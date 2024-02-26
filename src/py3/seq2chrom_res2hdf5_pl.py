# -*- coding: utf-8 -*-

"""Make HDF5 format file for ML

Args for inputs:
    --chromFile
    --peakFile
    --expFile
    --clusterFile

Args for outputs:
    --output

Please see --help.

All of the merging processes are based on 'inner join' by the cluster ID key.

"""

import argparse
import numpy as np
import os
import pandas as pd

import warnings
warnings.simplefilter("ignore", category = FutureWarning)
import h5py

def mkdir_p(path):
    if not os.path.isdir(path):
        os.makedirs(path)

parser = argparse.ArgumentParser(description='Make HDF5 format file for ML')
parser.add_argument('--output', required=True,
                    help = 'Output directory. All of the output files are used in the following ML.')
parser.add_argument('--chromFile', required=True,
                    help = 'Reduced X file from seq2chrom_ref.py.')
parser.add_argument('--peakFile', required=True,
                    help = 'Information of --chromFile. The order must be the same with the chromFile')
parser.add_argument('--expFile', required=True,
                    help = 'CAGE transcriptome data (TSV); 1st column is cluster ID (the header must be clusterID) and 2nd- columns are normalized but antilogarithm expression levels')
parser.add_argument('--clusterFile', required=True,
                    help = 'CAGE cluster information. [F5.cage_cluster.hg19.info.tsv.gz]')
parser.add_argument('--lowmem', action='store_true',
                    help="Prepare HDF5 file in low memory mode. This mode is useful for large data. Default: False.")
parser.add_argument('--chunksize', type=int, default=10000,
                    help="Chunk size for low memory mode. Default: 10000")
parser.add_argument('--polars', action='store_true',
                    help="Use Polars instead of Pandas. Default: False")
args = parser.parse_args()

# make output dir
mkdir_p(args.output)

# read resources

if args.polars:
    raise NotImplementedError("Mode: Polars is not implemented yet.")

    # TODO:

    ## expression data
    Expr_head = pl.read_csv(args.expFile, delimiter = '\t', header=None, n_rows = 1)
    _dtype = {'prmtrID': str}
    for i in range(1, Expr_head.shape[1]):
        _dtype[Expr_head.iloc[0, i]] = 'float32'
    
    Expr = pl.read_csv(args.expFile, delimiter = '\t', dtype = _dtype)
    Expr.rename(columns = {'prmtrID': 'clusterID'}, inplace = True)

    ## Peak file
    Peak = pl.read_csv(args.peakFile, delimiter = '\t', names = ['clusterID', 'peak'], dtype = {'clusterID': str, 'peak': int}, header=None)
    ## Cluster file
    _dtype = {'clusterID': str,
            'clusterName': str,
            'type': str,
            'mask': str,
            'F5_tag_count': int,
            'geneNum': int,
            'trnscptIDStr': str,
            'geneIDStr': str,
            'geneNameStr': str,
            'geneClassStr': str,
            'F5_anno': str}
    Cluster = pl.read_csv(args.clusterFile, delimiter = '\t', dtype = _dtype)
    assert Expr.shape[0] == Cluster.shape[0]

    # Inner join between Cluster and Expr
    Cluster_cols = Cluster.columns
    Expr_cols = Expr.columns
    _Expr = Cluster.join(Expr, on = 'clusterID', how = 'inner')
    assert _Expr.shape[0] == Cluster.shape[0]
    assert _Expr.shape[1] == (Cluster.shape[1] + Expr.shape[1] - 1)    

else:
    ## expression data
    Expr_head = pd.read_csv(args.expFile, delimiter = '\t', header=None, nrows = 1)
    _dtype = {'prmtrID': str}
    for i in range(1, Expr_head.shape[1]):
        _dtype[Expr_head.iloc[0, i]] = 'float32'

    Expr = pd.read_csv(args.expFile, delimiter = '\t', dtype = _dtype)
    Expr.rename(columns = {'prmtrID': 'clusterID'}, inplace = True)

    ## Peak file
    Peak = pd.read_csv(args.peakFile, delimiter = '\t', names = ['clusterID', 'peak'], dtype = {'clusterID': str, 'peak': int}, header=None)
    ## Cluster file
    _dtype = {'clusterID': str,
            'clusterName': str,
            'type': str,
            'mask': str,
            'F5_tag_count': int,
            'geneNum': int,
            'trnscptIDStr': str,
            'geneIDStr': str,
            'geneNameStr': str,
            'geneClassStr': str,
            'F5_anno': str}
    Cluster = pd.read_csv(args.clusterFile, delimiter = '\t', dtype = _dtype)
    assert Expr.shape[0] == Cluster.shape[0]

    # Inner join between Cluster and Expr
    Cluster_cols = Cluster.columns
    Expr_cols = Expr.columns
    _Expr = pd.merge(Cluster, Expr, on = 'clusterID')
    assert _Expr.shape[0] == Cluster.shape[0]
    assert _Expr.shape[1] == (Cluster.shape[1] + Expr.shape[1] - 1)

# Find NA rows in chromatin effects
def find_na(x):
    na_val = -9999999999.0
    return all(x == na_val)

def get_strand(x):
    try:
        strand = x.split(',')[1]
    except IndexError:
        strand = None
    return strand
    
## chromatin effect (float32; torch default)
if args.lowmem:

    # Use chunksize for low memory mode
    readers = pd.read_csv(args.chromFile, dtype = 'float32', delimiter = '\t', header=None, chunksize=args.chunksize)
    # Index for extracting the same rows from Peak
    idx_start = 0
    for Xreducedall_chunk in readers:

        chunked_size = Xreducedall_chunk.shape[0]
        print('Running chunk: {}-{}'.format(idx_start, idx_start + chunked_size))

        # Prepare index for extracting the same rows from Peak
        idx_end = idx_start + chunked_size
        Peak_chunk = Peak.iloc[idx_start:idx_end, :]

        _names = []
        for i in range(Xreducedall_chunk.shape[1]):
            _names.append('chrom{}'.format(i))

        Xreducedall_chunk.columns = _names

        # Check data size before removing NA in chromatin effect
        assert Xreducedall_chunk.shape[0] == Peak_chunk.shape[0]

        na_rows_chunk = Xreducedall_chunk.apply(find_na, axis = 1)
        na_rows = na_rows_chunk if 'na_rows' not in locals() else pd.concat([na_rows, na_rows_chunk])

        # Remove NA rows from Xreducedall_chunk and Peak_chunk
        Xreducedall_chunk = Xreducedall_chunk.loc[na_rows_chunk == False, :]
        Peak_chunk = Peak_chunk.loc[na_rows_chunk == False, :]

        # Join Peak_chunk and Xreducedall_chunk for inner join with Expr and Cluster
        Xreducedall_chunk_cols = Xreducedall_chunk.columns
        Peak_chunk_cols = Peak_chunk.columns
        _Xreducedall_chunk = pd.concat([Peak_chunk, Xreducedall_chunk], axis = 1)
        assert _Xreducedall_chunk.shape[0] == Peak_chunk.shape[0]

        # Inner join between _Xreducedall_chunk and _Expr
        all_chunk = pd.merge(_Xreducedall_chunk, _Expr, on = 'clusterID')
        all_chunk_shape = all_chunk.shape

        # add chr, strand column from clusterID
        _chr_chunk = [x.split(':')[0] for x in all_chunk['clusterID'].values]
        _strand_chunk = [get_strand(x) for x in all_chunk['clusterID'].values]
        all_chunk['chr'] = _chr_chunk
        all_chunk['strand'] = _strand_chunk

        # Split data into train, test, others
        trainind_chunk = np.asarray(~all_chunk.loc[:, 'chr'].isin(['chrX', 'chrY', 'chrM', 'chr8']))
        testind_chunk = np.asarray(all_chunk.loc[:, 'chr'] == 'chr8')
        otherind_chunk = np.asarray(all_chunk.loc[:, 'chr'].isin(['chrX', 'chrY', 'chrM']))
        X_train_chunk = all_chunk[Xreducedall_chunk_cols].loc[trainind_chunk]
        Y_train_chunk = all_chunk[Expr_cols[1:]].loc[trainind_chunk]
        info_train_chunk = all_chunk.drop(columns=Xreducedall_chunk_cols).drop(columns=Expr_cols[1:]).loc[trainind_chunk]
        X_test_chunk = all_chunk[Xreducedall_chunk_cols].loc[testind_chunk]
        Y_test_chunk = all_chunk[Expr_cols[1:]].loc[testind_chunk]
        info_test_chunk = all_chunk.drop(columns=Xreducedall_chunk_cols).drop(columns=Expr_cols[1:]).loc[testind_chunk]
        X_other_chunk = all_chunk[Xreducedall_chunk_cols].loc[otherind_chunk]
        Y_other_chunk = all_chunk[Expr_cols[1:]].loc[otherind_chunk]
        info_other_chunk = all_chunk.drop(columns=Xreducedall_chunk_cols).drop(columns=Expr_cols[1:]).loc[otherind_chunk]

        # Save X and Y as hdf5 format (save object is not torch object)
        if idx_start == 0:
            with h5py.File(args.output + 'X_Y.h5', 'w') as f:
                # train
                f.create_group('train')
                f['train'].create_dataset('X', data = X_train_chunk, compression = "gzip", compression_opts = 9, maxshape=(None, X_train_chunk.shape[1]))
                f['train'].create_dataset('y', data = Y_train_chunk, compression = "gzip", compression_opts = 9, maxshape=(None, Y_train_chunk.shape[1]))
                # test
                f.create_group('test')
                f['test'].create_dataset('X', data = X_test_chunk, compression = "gzip", compression_opts = 9, maxshape=(None, X_test_chunk.shape[1]))
                f['test'].create_dataset('y', data = Y_test_chunk, compression = "gzip", compression_opts = 9, maxshape=(None, Y_test_chunk.shape[1]))
                # others
                f.create_group('other')
                f['other'].create_dataset('X', data = X_other_chunk, compression = "gzip", compression_opts = 9, maxshape=(None, X_other_chunk.shape[1]))
                f['other'].create_dataset('y', data = Y_other_chunk, compression = "gzip", compression_opts = 9, maxshape=(None, Y_other_chunk.shape[1]))
        else:
            with h5py.File(args.output + 'X_Y.h5', 'a') as f:
                # train
                ## X
                f['train']['X'].resize((f['train']['X'].shape[0] + X_train_chunk.shape[0]), axis = 0)
                f['train']['X'][-X_train_chunk.shape[0]:] = X_train_chunk
                ## y
                f['train']['y'].resize((f['train']['y'].shape[0] + Y_train_chunk.shape[0]), axis = 0)
                f['train']['y'][-Y_train_chunk.shape[0]:] = Y_train_chunk
                # test
                ## X
                f['test']['X'].resize((f['test']['X'].shape[0] + X_test_chunk.shape[0]), axis = 0)
                f['test']['X'][-X_test_chunk.shape[0]:] = X_test_chunk
                ## y
                f['test']['y'].resize((f['test']['y'].shape[0] + Y_test_chunk.shape[0]), axis = 0)
                f['test']['y'][-Y_test_chunk.shape[0]:] = Y_test_chunk
                # others
                ## X
                f['other']['X'].resize((f['other']['X'].shape[0] + X_other_chunk.shape[0]), axis = 0)
                f['other']['X'][-X_other_chunk.shape[0]:] = X_other_chunk
                ## y
                f['other']['y'].resize((f['other']['y'].shape[0] + Y_other_chunk.shape[0]), axis = 0)
                f['other']['y'][-Y_other_chunk.shape[0]:] = Y_other_chunk
        
        # Save info, X column, Y column, common rows as txt.gz file
        if idx_start == 0:
            pd.DataFrame(X_train_chunk.columns).to_csv(args.output + 'X_train_col.txt.gz', header = False, compression = 'gzip', sep = "\t", index=False)
            pd.DataFrame(X_test_chunk.columns).to_csv(args.output + 'X_test_col.txt.gz', header = False, compression = 'gzip', sep = "\t", index=False)
            pd.DataFrame(X_other_chunk.columns).to_csv(args.output + 'X_other_col.txt.gz', header = False, compression = 'gzip', sep = "\t", index=False)
            pd.DataFrame(Y_train_chunk.columns).to_csv(args.output + 'Y_train_col.txt.gz', header = False, compression = 'gzip', sep = "\t", index=False)
            pd.DataFrame(Y_test_chunk.columns).to_csv(args.output + 'Y_test_col.txt.gz', header = False, compression = 'gzip', sep = "\t", index=False)
            pd.DataFrame(Y_other_chunk.columns).to_csv(args.output + 'Y_other_col.txt.gz', header = False, compression = 'gzip', sep = "\t", index=False)

        info_train = info_train_chunk if 'info_train' not in locals() else pd.concat([info_train, info_train_chunk])
        info_test = info_test_chunk if 'info_test' not in locals() else pd.concat([info_test, info_test_chunk])
        info_other = info_other_chunk if 'info_other' not in locals() else pd.concat([info_other, info_other_chunk])

        # data metrics
        metrics_chunk = pd.DataFrame([[len(X_train_chunk), len(X_test_chunk), len(X_other_chunk)]])
        metrics_chunk.columns = ['N_train', 'N_test', 'N_other']
        if 'metrics' not in locals():
            metrics = metrics_chunk
        else:
            metrics = metrics + metrics_chunk
        
        # Update index
        idx_start += chunked_size

elif args.polars:
    raise NotImplementedError("Mode: Polars is not implemented yet.")

    # TODO:

else:
    Xreducedall = pd.read_csv(args.chromFile, dtype = 'float32', delimiter = '\t', header=None)
    _names = []
    for i in range(Xreducedall.shape[1]):
        _names.append('chrom{}'.format(i))

    Xreducedall.columns = _names

    # Check data size before removing NA in chromatin effect
    assert Xreducedall.shape[0] == Peak.shape[0]
    
    na_rows = Xreducedall.apply(find_na, axis = 1)

    # Remove NA rows from Xreducedall and Peak
    Xreducedall = Xreducedall.loc[na_rows == False, :]
    Peak = Peak.loc[na_rows == False, :]

    # Join Peak and Xreducedall for inner join with Expr and Cluster
    Xreducedall_cols = Xreducedall.columns
    Peak_cols = Peak.columns
    _Xreducedall = pd.concat([Peak, Xreducedall], axis = 1)
    assert _Xreducedall.shape[0] == Peak.shape[0]

    # Inner join between _Xreducedall and _Expr 
    all = pd.merge(_Xreducedall, _Expr, on = 'clusterID')
    all_shape = all.shape

    # add chr, strand column from clusterID
    _chr = [x.split(':')[0] for x in all['clusterID'].values]
    _strand = [get_strand(x) for x in all['clusterID'].values]
    all['chr'] = _chr
    all['strand'] = _strand

    # Split data into train, test, others
    trainind = np.asarray(~all.loc[:, 'chr'].isin(['chrX', 'chrY', 'chrM', 'chr8']))
    testind = np.asarray(all.loc[:, 'chr'] == 'chr8')
    otherind = np.asarray(all.loc[:, 'chr'].isin(['chrX', 'chrY', 'chrM']))
    X_train = all[Xreducedall_cols].loc[trainind]
    Y_train = all[Expr_cols[1:]].loc[trainind]
    info_train = all.drop(columns=Xreducedall_cols).drop(columns=Expr_cols[1:]).loc[trainind]
    X_test = all[Xreducedall_cols].loc[testind]
    Y_test = all[Expr_cols[1:]].loc[testind]
    info_test = all.drop(columns=Xreducedall_cols).drop(columns=Expr_cols[1:]).loc[testind]
    X_other = all[Xreducedall_cols].loc[otherind]
    Y_other = all[Expr_cols[1:]].loc[otherind]
    info_other = all.drop(columns=Xreducedall_cols).drop(columns=Expr_cols[1:]).loc[otherind]

    # Save X and Y as hdf5 format (save object is not torch object)
    with h5py.File(args.output + 'X_Y.h5', 'w') as f:
        # train
        f.create_group('train')
        f['train'].create_dataset('X', data = X_train, compression = "gzip", compression_opts = 9)
        f['train'].create_dataset('y', data = Y_train, compression = "gzip", compression_opts = 9)
        # test
        f.create_group('test')
        f['test'].create_dataset('X', data = X_test, compression = "gzip", compression_opts = 9)
        f['test'].create_dataset('y', data = Y_test, compression = "gzip", compression_opts = 9)
        # others
        f.create_group('other')
        f['other'].create_dataset('X', data = X_other, compression = "gzip", compression_opts = 9)
        f['other'].create_dataset('y', data = Y_other, compression = "gzip", compression_opts = 9)

    # Save X column, Y column, common rows as txt.gz file
    pd.DataFrame(X_train.columns).to_csv(args.output + 'X_train_col.txt.gz', header = False, compression = 'gzip', sep = "\t", index=False)
    pd.DataFrame(X_test.columns).to_csv(args.output + 'X_test_col.txt.gz', header = False, compression = 'gzip', sep = "\t", index=False)
    pd.DataFrame(X_other.columns).to_csv(args.output + 'X_other_col.txt.gz', header = False, compression = 'gzip', sep = "\t", index=False)
    pd.DataFrame(Y_train.columns).to_csv(args.output + 'Y_train_col.txt.gz', header = False, compression = 'gzip', sep = "\t", index=False)
    pd.DataFrame(Y_test.columns).to_csv(args.output + 'Y_test_col.txt.gz', header = False, compression = 'gzip', sep = "\t", index=False)
    pd.DataFrame(Y_other.columns).to_csv(args.output + 'Y_other_col.txt.gz', header = False, compression = 'gzip', sep = "\t", index=False)

    # data metrics
    metrics = pd.DataFrame([[len(X_train), len(X_test), len(X_other)]])
    metrics.columns = ['N_train', 'N_test', 'N_other']

print("#na_rows: ", sum(na_rows))

# Write info
info_train.to_csv(args.output + 'info_train.txt.gz', header = True, compression = 'gzip', sep = "\t", index=False)
info_test.to_csv(args.output + 'info_test.txt.gz', header = True, compression = 'gzip', sep = "\t", index=False)
info_other.to_csv(args.output + 'info_other.txt.gz', header = True, compression = 'gzip', sep = "\t", index=False)


# Write common rows for train, test, others
pd.DataFrame(info_train['clusterID']).to_csv(args.output + 'common_train_rows.txt.gz', header = False, compression = 'gzip', sep = "\t", index=False)
pd.DataFrame(info_test['clusterID']).to_csv(args.output + 'common_test_rows.txt.gz', header = False, compression = 'gzip', sep = "\t", index=False)
pd.DataFrame(info_other['clusterID']).to_csv(args.output + 'common_other_rows.txt.gz', header = False, compression = 'gzip', sep = "\t", index=False)

# Write metrics
metrics.to_csv(args.output + 'metrics.txt', header = True, sep = "\t", index=False)

print('=== Finish ===')