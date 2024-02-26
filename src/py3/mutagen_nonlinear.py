# -*- coding: utf-8 -*-
"""in silico mutagenesis using MENTER
"""
import argparse
import math
import pyfasta
import torch
from torch.utils.serialization import load_lua
import numpy as np
import pandas as pd
import xgboost as xgb
from six.moves import reduce
import sys
from sklearn.isotonic import IsotonicRegression
from sklearn.externals import joblib

# silent future warning
import warnings
warnings.simplefilter("ignore", category = FutureWarning)
import pandas as pd

# parser
parser = argparse.ArgumentParser(description='In silico mutagenesis using non-linear ML from Chromatin effects.')

# General Parameters
parser.add_argument('--peak_window_size', type=int, default=100000,
                    help="Sequence window size around peak (one side); Current best model is +/-100-kb (then default is 100000)")
parser.add_argument('--max_InDel_size', type=int, default=250,
                    help="Maximum length of InDel. [bp]")
parser.add_argument('--buffer', type=int, default=2000,
                    help="Buffer length for making left and righ seq [bp].")
parser.add_argument('--ref_fa', required=True,
                    help = "Reference fasta file")
parser.add_argument('--threads', type=int, default=5,
                    help="Number of threads (CPU).")
parser.add_argument('--input_cols', nargs="+", default=['CHR', 'POS', 'REF', 'ALT'],
                    help="Column index names for input file.")
parser.add_argument('--modelList', required=True,
                    help="A list of paths of binary xgboost model files")
parser.add_argument('--calib_modelList', required=True,
                    help="A list of paths of calibration model files")
parser.add_argument('--modelList_prefix', required=True,
                    help="Prefix for the modelList and calib_modelList files (i.e. where MENTR is installed.)")
parser.add_argument('--verbose', action="store_true")

# input/output setting
parser.add_argument('--multifiles', action="store_true",
                    help = 'If True, MENTR use use --inFiles; else use --inputfile, --geneFile, and --out_cmn.')
parser.add_argument('--inFiles',
                    help = 'TSV file. col1: input file; col2: closestgene file; col3: out_cmn')
parser.add_argument('--inputfile',
                    help='Input file; see "--vcf-cols"')
parser.add_argument('--geneFile',
                    help = 'closestgene file')
parser.add_argument('--peakinfo', required=True,
                    help='Location of resources/cage_peak.txt.gz')
parser.add_argument('--out_cmn',
                    help='Output suffix')

# For DeepSEA
parser.add_argument('--deepsea', required=True,
                    help="Location of deepsea.beluga.2002.cpu (please download from https://github.com/FunctionLab/ExPecto)")
parser.add_argument('--nfeatures', type=int, default=2002,
                    help="This depends on trainded-DeepSEA window size.; default = 2002 and this would not be changed.")
parser.add_argument('--cuda', action='store_true')
parser.add_argument('--inputsize', type=int, default=2000,
                    help="The input sequence window size for neural network. This depends on trainded-DeepSEA window size.; default = 2kb and this would not be changed.")
parser.add_argument('--batchsize', dest="batchsize", type=int, default=50,
                    help="Batch size for running deepsea.")

args = parser.parse_args()
inputsize = args.inputsize
nfeatures = args.nfeatures
peak_w_Size = args.peak_window_size
batchSize=args.batchsize

# OpenMP threads
torch.set_num_threads(args.threads)

# Number of bin (each side from peak); default 500
n_bin_a_side = int(peak_w_Size / 200.0)

# ref_fa
genome = pyfasta.Fasta(args.ref_fa)

# deepsea model
model = load_lua(args.deepsea)
model.evaluate()
if args.cuda:
    model.cuda()

# Length of sequence to use trancate function
# default peak_window_size(100kb) + buffer; each side from TSS
#  buffer: 900bp (window size used in DeepSEA) +
#          2000bp (> max_indel_length (250bp)), which is enough for short InDel analysis
bffr_DeepSEA = 900
use_seq_len = peak_w_Size + bffr_DeepSEA + args.buffer

# init
CHRS = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
        'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
        'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX','chrY']

def encodeSeq(seq):
    """Convert sequences to 0-1 encoding.
    The output concatenates the forward and reverse complement sequence encodings.

    Args:
        seq: sequence

    Returns:
        numpy array of dimension: 2 x 4 x len(seq)

    2 means the concatenation of forward and reverse complement sequences.
    4 means base sequence corresponding to 'ATGC'
    """
    seqsnp = np.zeros((1, 4, len(seq)), np.bool_)
    mydict = {'A': np.asarray([1, 0, 0, 0]), 'G': np.asarray([0, 1, 0, 0]),
              'C': np.asarray([0, 0, 1, 0]), 'T': np.asarray([0, 0, 0, 1]),
              'N': np.asarray([0, 0, 0, 0]), 'H': np.asarray([0, 0, 0, 0]),
              'U': np.asarray([0, 0, 0, 0]), 'R': np.asarray([0, 0, 0, 0]),
              'Y': np.asarray([0, 0, 0, 0]), 'K': np.asarray([0, 0, 0, 0]),
              'M': np.asarray([0, 0, 0, 0]), 'S': np.asarray([0, 0, 0, 0]),
              'W': np.asarray([0, 0, 0, 0]), 'B': np.asarray([0, 0, 0, 0]),
              'D': np.asarray([0, 0, 0, 0]), 'V': np.asarray([0, 0, 0, 0]),
              'a': np.asarray([1, 0, 0, 0]), 'g': np.asarray([0, 1, 0, 0]),
              'c': np.asarray([0, 0, 1, 0]), 't': np.asarray([0, 0, 0, 1]),
              'n': np.asarray([0, 0, 0, 0]), '-': np.asarray([0, 0, 0, 0])}
    for i, c in enumerate(seq):
        seqsnp[0, :, i] = mydict[c]
    # get the complementary sequences (ex. ATGC -> GCAT (complement and flip)); default in DeepSEA/ExPecto
    dataflip = seqsnp[:, ::-1, ::-1]
    seqsnp = np.concatenate([seqsnp, dataflip], axis=0)
    return seqsnp.astype(np.float32)

def truncate(seq, l_or_r, tr_size):
    '''Truncate sequence

    Args:
        seq: sequnce including buffer regions for Indel
        l_or_r: 'left' or 'right'
        tr_size: size
                
    Returns:
        seq
        
    '''
    if l_or_r == 'right':
        # Left side is exact.
        st = 0
        en = tr_size
    elif l_or_r == 'left':
        # Right side is exact.
        st = len(seq) - tr_size
        en = len(seq)
    else:
        raise ValueError("l_or_r must be either right/left")
    assert st >= 0 and en >= 0
    out = seq[st:en]

    assert len(out) == tr_size

    return out

def fetchSeqs(chr, pos, ref, alt, tss_pos, bffr_tss_window = use_seq_len, bffr_DeepSEA = 900):
    """Fetches sequences from the genome.

    Retrieves sequences centered at the representative tss position with the given bffr_tss_window (100-kb + plenty of buffer sequence).
    Returns both reference and alternative allele sequences.
    
    Args:
        chr: the chromosome name that must matches one of the names in CHRS.
        pos: chromosome coordinate (1-based).
        ref: the reference allele.
        alt: the alternative allele.
        tss_pos: the representative tss position
        bffr_tss_window: the targeted sequence length including buffer region (one-side) [bp]
        bffr_DeepSEA: window size used in DeepSEA [bp]
        
    Returns:
        reference sequnce (left side from TSS)
        alternative sequence (left side from TSS)
        reference sequnce (right side from TSS)
        alternative sequence (right side from TSS)
        booling for checking reference allele (left sequence)
        booling for checking reference allele (right sequence)

    """
    matched_bool_left = False
    matched_bool_right = False

    # left of the TSS
    # first bin: {tss - 100k + 1, tss - 100k + 200} with plenty of buffer window
    # last bin: {tss - 1 - 200, tss} with right side-buffer region -> {tss + 1, tss + bffr_DeepSEA}; Must be exact position.
    left_st = tss_pos - bffr_tss_window + 1
    left_en = tss_pos + bffr_DeepSEA
    
    # right of the TSS
    # first bin: {tss, tss + 199} with left side-buffer region -> {tss - bffr_DeepSEA , tss - 1}; Must be exact position.
    # last bin: {tss + 100k - 199 - 1, tss + 100k - 1} with plenty of buffer window
    right_st = tss_pos - bffr_DeepSEA
    right_en = tss_pos + bffr_tss_window - 1
    
    # get sequence (reference); 1-based count
    seq_left = genome.sequence({'chr': chr, 'start': left_st, 'stop': left_en})
    seq_right = genome.sequence({'chr': chr, 'start': right_st, 'stop': right_en})
    
    # mutation spot (relative posision to *_st): use for slicing
    mutpos_left = pos - left_st
    mutpos_right = pos - right_st

    # mutagenesis (left)
    if mutpos_left >= 0 and pos <= left_en:
        # Reference check
        refseq_left = seq_left[:mutpos_left] + ref + seq_left[(mutpos_left + len(ref)):]
        altseq_left = seq_left[:mutpos_left] + alt + seq_left[(mutpos_left + len(ref)):] # len(ref) is not a mistake.
        matched_bool_left = (seq_left[mutpos_left:(mutpos_left + len(ref))].upper() == ref.upper())
    else:
        refseq_left = seq_left
        altseq_left = seq_left
        matched_bool_left = True

    # mutagenesis (right)
    if mutpos_right >= 0 and pos <= right_en:
        # Reference check
        refseq_right = seq_right[:mutpos_right] + ref + seq_right[(mutpos_right + len(ref)):]
        altseq_right = seq_right[:mutpos_right] + alt + seq_right[(mutpos_right + len(ref)):] # len(ref) is not a mistake.
        matched_bool_right = (seq_right[mutpos_right:(mutpos_right + len(ref))].upper() == ref.upper())
    else:
        refseq_right = seq_right
        altseq_right = seq_right
        matched_bool_right = True

    # truncate
    tran_len = peak_w_Size + (bffr_DeepSEA * 2)
    refseq_left = truncate(refseq_left, 'left', tran_len)
    altseq_left = truncate(altseq_left, 'left', tran_len)
    refseq_right = truncate(refseq_right, 'right', tran_len)
    altseq_right = truncate(altseq_right, 'right', tran_len)
    assert len(refseq_left) == tran_len
    assert len(altseq_left) == tran_len
    assert len(refseq_right) == tran_len
    assert len(altseq_right) == tran_len
    return refseq_left, altseq_left, refseq_right, altseq_right, matched_bool_left, matched_bool_right

def make_bins(left_seq, right_seq, inputSize = inputsize, n_bin_a_side = n_bin_a_side):
    '''Make bins from left and right sequence

    Args:
        left_seq: sequnce (left side from TSS which MUST BE already truncated)
        right_seq: sequnce (right side from TSS which MUST BE already truncated)
        inputsize: total sequence length for DeepSEA inference (default: 2-kb)
        n_bin_a_side: how many bins are used for a side. (default: 500)
       
    Returns:
        list of bins; left and right sequences are trimed and connected.

    '''
    # Make sequence bins with the inputSize
    bins = []

    # Size Check (to check whether the sequence was trancated)
    # assume: 100kb + 1.8kb = 200bp * 500 + 900bp-buffer * 2
    assert len(left_seq) == (200.0 * n_bin_a_side + (bffr_DeepSEA * 2))
    assert len(right_seq) == (200.0 * n_bin_a_side + (bffr_DeepSEA * 2))

    # left seq; the margin exists in the right of the TSS
    st = 0
    for _ in range(n_bin_a_side):
        _en = st + inputSize
        _bin = left_seq[st:_en]
        bins.append(_bin.upper())
        st += 200

    # right seq; the margin exists in the left of the TSS
    st = 0
    for _ in range(n_bin_a_side):
        _en = st + inputSize
        _bin = right_seq[st:_en]
        bins.append(_bin.upper())
        st += 200
    return bins

# for bin calculation
# n_bin_a_side = 500;
# 4: (left and right) * (f and r)
nums = [i for i in range(n_bin_a_side * 4)] 
evens = nums[0::2]
odds = nums[1::2]

def calc_deepsea(bins, batchSize = batchSize, evens = evens, odds = odds):
    '''Calculate epigenetic effects from DeepSEA model

    Args:
        bins: list of sequences of each bin (default: total 500 * 2 bins)
        batchSize: batch size for DeepSEA
        evens: [0, 2, ..., 1998] 
        odds: [1, 3, ..., 1999]
                
    Returns:
        Epigenetic effects, in which effects of forward and reverse&complemental sequences are averaged.
        
    '''
    # bins (list of sequence, left to right) -> encode ->
    #  chrom preds (list of sequence effect, #bin * [2, 2002 features])

    # make complementary seq and encode them
    bin_encodeds = []
    for i in range(len(bins)):
        bin_encodeds.append(encodeSeq(bins[i]))  
    bin_encodeds = np.vstack(bin_encodeds)

    # predict chromatin effets
    chrom_preds = []
    for i in range(int(1 + bin_encodeds.shape[0] / batchSize)):
        if bin_encodeds.shape[0] == int(i*batchSize):
            continue
        input = torch.from_numpy(bin_encodeds[int(i*batchSize):int((i+1)*batchSize),:,:]).unsqueeze(3)
        if args.cuda:
            input = input.cuda()
        chrom_preds.append(model.forward(input).cpu().numpy().copy())
    chrom_preds = np.vstack(chrom_preds)

    # merge forward and reverse effects; #bin * [1, 2002 features]
    seqEffects = (chrom_preds[evens, :] + chrom_preds[odds, :]) / 2.0

    return seqEffects

# load xgboost models
modelList = pd.read_csv(args.modelList, sep='\t', header=0)
models = []
for file in modelList['ModelName']:
    bst = xgb.Booster({'nthread': args.threads})
    bst.load_model(args.modelList_prefix + "/" + file.strip())
    models.append(bst)

# calibratoin models
if args.calib_modelList is None:
    raise ValueError("Set --calib_modelList for Probability Mode.")
calib_modelList = pd.read_csv(args.calib_modelList, sep='\t', header=0)
calib_models = []
for Tissue in modelList['Tissue']:
    file = calib_modelList[calib_modelList['Tissue'] == Tissue]['ModelName'].values[0]
    ir = joblib.load(args.modelList_prefix + "/" + file.strip())
    calib_models.append(ir)


def compute_effects(seqEffects, all_models, nfeatures = nfeatures, peak_w_Size = peak_w_Size):
    """Compute expression effects (log fold-change).

    Args:
        seqEffects: list of chromatin effect numpy arrays; default, 500 * 2 lists (bins; distance from TSS) of [1, 2002 features]
        all_models: list of ExPecto model files.
        nfeatures: number of chromatin/epigenomic features (default: 2002)
        peak_w_Size: My ML's window size, one-side (default: 100-kb)

    Returns:
        numpy array of size num_variants x num_models. Each value represents
        predicted log fold-change
    """
    # 500 * 2 list of distance from peak, each [1, 5]; strand is not considered due to treating enhancer
    max_bin = peak_w_Size - 100
    min_bin = -1 * max_bin
    Xreducedall = [np.vstack([
        np.exp(-0.01 * np.floor(np.abs(dist / 200.0))),
        np.exp(-0.02 * np.floor(np.abs(dist / 200.0))),
        np.exp(-0.05 * np.floor(np.abs(dist / 200.0))),
        np.exp(-0.1  * np.floor(np.abs(dist / 200.0))),
        np.exp(-0.2  * np.floor(np.abs(dist / 200.0)))
    ]).T for dist in list(range(min_bin, (max_bin + 200), 200))]  # the same order with bins(seqEffects)
    
    assert len(Xreducedall) == n_bin_a_side * 2
    assert Xreducedall[0].shape == (1, 5)
    assert len(seqEffects) == n_bin_a_side * 2
    assert seqEffects.shape == (n_bin_a_side * 2, nfeatures)

    # j: index of bin, used for list
    #
    # for each bin:
    #   np.tile(seqEffects[j, :], 5) * np.repeat(Xreducedall[j], nfeatures(=2002), axis=1)
    #    -> (1, 10010)                     -> (1, 10010)
    # 
    # compute gene expression level with models over the distances to peak, [1, 10010 (=2002 * 5)] from 200 bins of [1, 2002]
    # According to the distance to the peak, seqEffects[j] [1, 2002] in the bins were weighted by Xreducedall[j] [0, 0:5]
    # 10010 order is: [seqEffects[j][0, 0:2002] * Xreducedall[j][0, 0], seqEffects[j][0, 0:2002] * Xreducedall[j][0, 1], ...]
    X = reduce(lambda x, y: x + y, [np.tile(seqEffects[j, :], 5) * 
                                    np.repeat(Xreducedall[j], nfeatures, axis=1) for j in range(len(Xreducedall))])
    Xtest = xgb.DMatrix(X)
    effect = np.zeros((1, len(all_models)))
    for j in range(len(all_models)):
        #effect[:, j] = all_models[j].predict(Xtest, ntree_limit = all_models[j].best_ntree_limit)
        # see:
        ## https://github.com/dmlc/xgboost/issues/805
        ## https://stackoverflow.com/questions/43534219/xgboost-what-is-the-difference-among-bst-best-score-bst-best-iteration-and-bst
        best_itr = int(all_models[j].attributes()["best_iteration"])
        # num_parallel_tree was not set (default:1)
        num_parallel_tree = 1
        effect[:, j] = all_models[j].predict(Xtest, ntree_limit = (best_itr + 1) * num_parallel_tree)
    return effect

# main
if __name__ == '__main__':

    if args.multifiles:
        inFiles = pd.read_csv(args.inFiles, sep='\t',header=None)
        assert inFiles.shape[1] == 3
        inputfiles = inFiles.iloc[:, 0].values
        geneFiles = inFiles.iloc[:, 1].values
        out_cmns = inFiles.iloc[:, 2].values
        max_k = inFiles.shape[0]
    else:
        inputfiles = [args.inputfile]
        geneFiles = [args.geneFile]
        out_cmns = [args.out_cmn]
        max_k = 1

    for k in range(max_k):

        # logging
        print("=== Running File {}/{}".format(k, max_k))
        sys.stderr.write("=== Running File {}/{}\n".format(k, max_k))

        #== single input vcf file
        vcf = pd.read_csv(inputfiles[k], sep='\t', header=None, comment='#')
        vcf.columns = args.input_cols

        # standardize
        vcf.iloc[:, 0] = 'chr' + vcf.iloc[:, 0].map(str).str.replace('chr', '')
        vcf = vcf[vcf.iloc[:, 0].isin(CHRS)]

        # chr
        chr_uniq = vcf.iloc[:, 0].unique()
        if len(chr_uniq) == 1:
            chr = chr_uniq.reshape(-1).tolist()[0]
        else:
            raise NotImplementedError("Currently input file with multiple chromosome is not supported.")

        #== gene and distance to TSS information
        gene = pd.read_csv(geneFiles[k],sep='\t',header=None)
        geneinds = pd.match(vcf.iloc[:,0].map(str).str.replace('chr','')+' '+vcf.iloc[:,1].map(str),
                    gene.iloc[:,0].map(str).str.replace('chr','')+' '+gene.iloc[:,2].map(str))
        if np.any(geneinds==-1):
            raise ValueError("Genefile does not match the input file.")

        dist = - np.asarray(gene.iloc[geneinds,-1]) # -1 is for output file
        genename = np.asarray(gene.iloc[geneinds,-2])
        strand= np.asarray(gene.iloc[geneinds,-3])
        if len(np.unique(genename)) != 1:
            print("Unique genename was:")
            print(np.unique(genename))
            raise NotImplementedError("Currently geneFile with multiple genes is not supported.")

        gene_id = np.unique(genename)[0]

        # get gene info gene (paek_data)
        peakinfo = pd.read_csv(args.peakinfo, sep='\t',header=None)
        tss = peakinfo[peakinfo[0] == gene_id].iloc[:, -1].values[0]
        # verbose check
        assert all(gene_id == genename)

        #- main
        diffs = []
        refs = []
        alts = []
        ref_matched_bools_left = []
        ref_matched_bools_right = []

        # each variant
        for i in range(vcf.shape[0]):
            if args.verbose:
                print("Running: {}/{}".format(i, vcf.shape[0]))
            elif i % 100 == 0:
                print("Running: {}/{}".format(i, vcf.shape[0]))
            refseq_left, altseq_left, refseq_right, altseq_right, ref_matched_bool_left, ref_matched_bool_right = fetchSeqs(vcf['CHR'][i], vcf['POS'][i], vcf['REF'][i], vcf['ALT'][i], tss)
            ref_matched_bools_left.append(ref_matched_bool_left)
            ref_matched_bools_right.append(ref_matched_bool_right)
            # Make sequence bins with the inputSize
            bins_ref = make_bins(refseq_left, refseq_right)
            bins_alt = make_bins(altseq_left, altseq_right)
            # run deepsea (including averaging reverse&complement sequence)
            seqEffects_ref = calc_deepsea(bins_ref)
            seqEffects_alt = calc_deepsea(bins_alt)
            # comptue expression effects
            seqExpEffects_ref = compute_effects(seqEffects_ref, models).reshape(-1)
            seqExpEffects_alt = compute_effects(seqEffects_alt, models).reshape(-1)
            # calibration
            for j in range(len(calib_models)):
                seqExpEffects_ref[j] = calib_models[j].transform(seqExpEffects_ref[j:(j+1)])
                seqExpEffects_alt[j] = calib_models[j].transform(seqExpEffects_alt[j:(j+1)])
            # summary
            diff = seqExpEffects_alt - seqExpEffects_ref
            refs.append(seqExpEffects_ref)
            alts.append(seqExpEffects_alt)
            diffs.append(diff)
        refs = np.asarray(refs)
        alts = np.asarray(alts)
        diffs = np.asarray(diffs)
        ref_matched_bools_left = np.asarray(ref_matched_bools_left)
        ref_matched_bools_right = np.asarray(ref_matched_bools_right)

        # check
        if np.sum(ref_matched_bools_left) != len(ref_matched_bools_left):
            print("(Left seq): Number of variants with reference allele matched with reference genome:")
            print(np.sum(ref_matched_bools_left), "(k = ", k," , i = ", i, ")")
            print("(Left seq): Number of input variants (left side):")
            print(len(ref_matched_bools_left), "(k = ", k," , i = ", i, ")")
        if np.sum(ref_matched_bools_right) != len(ref_matched_bools_right):
            print("(Right seq): Number of variants with reference allele matched with reference genome:")
            print(np.sum(ref_matched_bools_right), "(k = ", k," , i = ", i, ")")
            print("(Right seq): Number of input variants (right side):")
            print(len(ref_matched_bools_right), "(k = ", k," , i = ", i, ")")

        # write output
        def write_res(res, out_name):
            snpExpEffects_df = vcf
            snpExpEffects_df['dist'] = dist
            snpExpEffects_df['gene'] = genename
            snpExpEffects_df['strand'] = strand
            snpExpEffects_df['ref_matched'] = ref_matched_bools_left * ref_matched_bools_right
            snpExpEffects_df=pd.concat([snpExpEffects_df.reset_index(), pd.DataFrame(res, columns = modelList.iloc[:,1])], axis=1, ignore_index =False)
            snpExpEffects_df.to_csv(out_name, header = True, index = False, sep='\t', compression = 'gzip')

        write_res(diffs, out_cmns[k] + '_diff.txt.gz')
        write_res(refs, out_cmns[k] + '_refs.txt.gz')
        write_res(alts, out_cmns[k] + '_alts.txt.gz')
