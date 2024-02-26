# -*- coding: utf-8 -*-

"""Predict epigenetic features from reference sequence.
"""

import argparse
import torch
from torch.utils.serialization import load_lua
import numpy as np
from six.moves import reduce
from subprocess import check_output as system

# silent future warning
import warnings
warnings.simplefilter("ignore", category = FutureWarning)
#import h5py
import pandas as pd

parser = argparse.ArgumentParser(
    description='Predict variant epigenetic features from reference sequence using DeepSEA Beluga')
# General Parameters
parser.add_argument('--inputsize', type=int, default=2000,
                    help="Sequence window size for Deep Learning")
parser.add_argument('--peak_window_size', type=int, default=100000,
                    help="Sequence window size around peak (one side)")
parser.add_argument('--threads', type=int, default=1,
                    help="Number of threads.")
parser.add_argument('--output', required=True,
                    help="Prefix for output files")
parser.add_argument('--ref', type=str, required=True,
                    help='File location for reference sequence')
parser.add_argument('--peak_file', type=str, required=True,
                    help = 'Peak file [cage_peak.txt.gz]')

# For DeepSEA
parser.add_argument('--model', type=str, required=True,
                    help="Path for DeepSea model [deepsea.beluga.2002.cpu]")
parser.add_argument('--cuda', action='store_true',
                    help="If you use GPU, please activate this. ** Strongly Recommended **")
parser.add_argument('--nfeatures', type=int, default=2002)
parser.add_argument('--batchsize', dest="batchsize", type=int, default=50,
                    help="Batch size for neural network predictions.")

args = parser.parse_args()
inputSize = args.inputsize
peak_w_Size = args.peak_window_size
nfeatures=args.nfeatures
batchSize=args.batchsize

# Number of bin (each side from peak)
n_bin_a_side = int(peak_w_Size / 200.0)
# Length of sequence for using DL
use_seq_len = peak_w_Size + inputSize - 200

# OpenMP threads
torch.set_num_threads(args.threads)

# path for reference sequence
ref = args.ref

# load peak info file
peak_data = pd.read_table(args.peak_file, header = None, sep = '\t', names = ['ID', 'peak'], dtype = {'ID': str, 'peak': int})

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
        raise ValueError
    if st < 0 or en < 0:
        return None

    out = seq[st:en]
    if len(out) != tr_size:
        return None
    else:
        return out

def get_seq(st, en, l_or_r, chr, tr_size = use_seq_len, ref = ref):
    '''Get sequence and trancate it.

    Args:
        st: start position
        en: end position
        l_or_r: 'left' or 'right'
        chr: chromosome (without 'chr')

    Returns:
        sequence (truncated)

    '''
    if st < 0:
        return None
    else:
        pipe1 = "samtools faidx {0} chr{1}:{2}-{3}".format(ref, chr, st, en)
        pipe2 = "grep -v \"^>\""
        seq = system(pipe1 + " | " + pipe2, shell = True).decode('utf-8').split("\n")
        seq = ''.join(seq[0:(len(seq) - 1)])
        # truncate
        seq = truncate(seq, l_or_r, tr_size)
        return seq

def make_bins(left_seq, right_seq, inputSize = inputSize, n_bin_a_side = n_bin_a_side):
    '''Make bins from left and right sequence

    Args:
        left_seq: sequnce (left side from peak; already truncated)
        right_seq: sequnce (right side from peak; already truncated)
        
    Returns:
        list of bins; left and right sequences are trimed and connected.

    '''
    # Make sequence bins with the inputSize
    bins = []
    # left seq; the margin exists in the right of the peak
    st = 0
    for _ in range(n_bin_a_side):
        _en = st + inputSize
        _bin = left_seq[st:_en]
        bins.append(_bin.upper())
        st += 200
    # right seq; the margin exists in the left of the peak
    st = 0
    for _ in range(n_bin_a_side):
        _en = st + inputSize
        _bin = right_seq[st:_en]
        bins.append(_bin.upper())
        st += 200
    return bins

def get_left_right_pos(peak_pos, peak_w_Size = peak_w_Size):
    """Fetches sequences from the genome.

    Retrieves sequences centered at the representative peak position with the given bffr_peak_window (including plenty of buffer sequence).
    Returns both reference and alternative allele sequences . 

    Args:
        peak_pos: the representative peak position

    """
    bffr_peak_window = peak_w_Size + 900
    bffr_expecto = 900
    # left of the peak
    # start: bin {peak - 100k + 1, peak - 100k + 200} with plenty of buffer window
    # end: bin {peak - 1 - 200, peak} with right side-buffer region {peak + 1, peak + 900}; Must be exact position.
    left_st = peak_pos - bffr_peak_window + 1
    left_en = peak_pos + bffr_expecto
    # right of the peak
    # start: bin {peak, peak + 199} with left side-buffer region {peak - 900 , peak - 1}; Must be exact position.
    # end: bin {peak + 100k - 199 - 1, peak + 100k - 1} with plenty of buffer window
    right_st = peak_pos - bffr_expecto
    right_en = peak_pos + bffr_peak_window - 1
    return left_st, left_en, right_st, right_en

def encodeSeq(seq):
    """Convert sequences to 0-1 encoding.
    The output concatenates the forward and reverse complement sequence encodings.

    Args:
        seq: sequence

    Returns:
        numpy array of dimension: 2 x 4 x len(seq)

          2: the concatenation of forward and reverse complement sequences.
          4: base sequence corresponding to 'ATGC'
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

# for bin calculation
nums = [i for i in range(n_bin_a_side * 4)] 
evens = nums[0::2]
odds = nums[1::2]

# Load Model
model = load_lua(args.model)
model.evaluate()
if args.cuda:
    model.cuda()

def calc_deepsea(bins, batchSize = batchSize, evens = evens, odds = odds):
    '''Calculate epigenetic effects from DeepSEA model

    Args:
        bins: list of sequences of each bin (total 200 bins)
        batchSize
        evens: [0, 2, ..., 198] 
        odds: [1, 3, ..., 199]
                
    Returns:
        Epigenetic effects, in which effects of forward and reverse&complemental sequences are averaged.
        
    '''
    # bins (list of sequence, left to right) -> encode -> chrom preds (list of sequence effect, #bin * [2, 2002 features])
    bin_encodeds = []
    for i in range(len(bins)):
        bin_encodeds.append(encodeSeq(bins[i]))  
    bin_encodeds = np.vstack(bin_encodeds)
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

def compute_effects(seqEffects, nfeatures = nfeatures, peak_w_Size = peak_w_Size):
    """Compute summarized chromatin effects
    
    Args:
        seqEffects: list of chromatin effect numpy arrays; default, 200 lists (bins; distance from peak) of [1, 2002 features]
        nfeatures: number of chromatin/epigenomic features.

    Returns:

    """
    # lists of distance from peak, each [1, 5]; strand is not considered due to treating enhancer
    max_bin = peak_w_Size - 100
    min_bin = -1 * max_bin
    Xreducedall = [np.vstack([
        np.exp(-0.01 * np.floor(np.abs(dist / 200.0))),
        np.exp(-0.02 * np.floor(np.abs(dist / 200.0))),
        np.exp(-0.05 * np.floor(np.abs(dist / 200.0))),
        np.exp(-0.1  * np.floor(np.abs(dist / 200.0))),
        np.exp(-0.2  * np.floor(np.abs(dist / 200.0)))
    ]).T for dist in list(range(min_bin, (max_bin + 200), 200))]  # the same order with bins(seqEffects)
    X = reduce(lambda x, y: x + y, [np.tile(seqEffects[j, :], 5) * 
                                    np.repeat(Xreducedall[j], nfeatures, axis=1) for j in range(len(Xreducedall))])
    return X

def calc_eff(seqChromEffects_list, left_st, left_en, right_st, right_en, chr, nfeatures = nfeatures):
    left_seq = get_seq(left_st, left_en, 'left', chr)
    right_seq = get_seq(right_st, right_en, 'right', chr)
    if left_seq is None or right_seq is None:
        seqEffects = []
        for _i in range(n_bin_a_side * 2):
            seqEffects.append(np.full((1, nfeatures), -9999999999.0))
        seqChromEffects = np.full((1, nfeatures*5), -9999999999.0)
    else:
        # Make sequence bins with the inputSize
        bins = make_bins(left_seq, right_seq)
        # run deepsea (including averaging reverse&complement sequence)
        seqEffects = calc_deepsea(bins)
        assert all(sum(seqEffects < 0) == 0)
        # comptue chrom effects
        seqChromEffects = compute_effects(seqEffects)
    #seqEffects_list.append(seqEffects)
    seqChromEffects_list.append(seqChromEffects)
    return seqChromEffects_list #, seqEffects_list

# main
if __name__ == '__main__':
    #seqEffects_list = []
    seqChromEffects_list = []
    for i in range(len(peak_data)):
        if i % 100 == 0:
            print("Running: {}/{}".format(i, len(peak_data)))
        chr = peak_data.iloc[i,0].split(':')[0].replace('chr', '')
        peak = peak_data.iloc[i,1]
        left_st, left_en, right_st, right_en = get_left_right_pos(peak)
        seqChromEffects_list = calc_eff(seqChromEffects_list, left_st, left_en, right_st, right_en, chr)

    # write output (reduced chrom effects, txt.gz)
    np.savetxt(args.output + '.peak_window_' + str(peak_w_Size) + '.reduced_wo_strand.txt.gz', np.array(seqChromEffects_list).reshape([len(peak_data), nfeatures*5]), delimiter = '\t')

    print("===END===")
