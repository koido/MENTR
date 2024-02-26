# -*- coding: utf-8 -*-

"""Make HDF5 format file for ML

Args for inputs:
    --chromFile
    --peakFile
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
                    help = 'Reduced X file from ./cd4t/src/py3/seq2chrom_ref.py.')
parser.add_argument('--peakFile', required=True,
                    help = 'Information of --chromFile. The order must be the same with the chromFile. No header.')
parser.add_argument('--clusterFile', required=True,
                    help = 'CAGE cluster information. Header is required. In this version, two columns, clusterID and type (promoter or enhancer), must be required. No further columns are not accepted.')

args = parser.parse_args()

# make output dir
mkdir_p(args.output)

# read resources

## chromatin effect (float32; torch default)
Xreducedall = pd.read_csv(args.chromFile, dtype = 'float32', delimiter = '\t', header=None)
_names = []
for i in range(Xreducedall.shape[1]):
    _names.append('chrom{}'.format(i))

Xreducedall.columns = _names

## Peak file
Peak = pd.read_csv(args.peakFile, delimiter = '\t', names = ['clusterID', 'peak'], dtype = {'clusterID': str, 'peak': int}, header=None)

## Cluster file
### NOTE: In this version, no header.
_dtype = {'clusterID': str,
          'type': str}
Cluster = pd.read_csv(args.clusterFile, delimiter = '\t', dtype = _dtype, header=None, names = ['clusterID', 'type'])

# Check data size before removing NA in chromatin effect
assert Xreducedall.shape[0] == Peak.shape[0]

# Find NA rows in chromatin effects
def find_na(x):
    na_val = -9999999999.0
    return all(x == na_val)

na_rows = Xreducedall.apply(find_na, axis = 1)
print("#na_rows: ", sum(na_rows))

# Remove NA rows from Xreducedall and Peak
Xreducedall = Xreducedall.loc[na_rows == False, :]
Peak = Peak.loc[na_rows == False, :]

# Concatenate Peak and Xreducedall for inner join with Cluster
Xreducedall_cols = Xreducedall.columns
Peak_cols = Peak.columns
_Xreducedall = pd.concat([Peak, Xreducedall], axis = 1)
assert _Xreducedall.shape[0] == Peak.shape[0]

# Inner join between _Xreducedall and Cluster 
all = pd.merge(_Xreducedall, Cluster)
all_shape = all.shape

# add chr, strand column from clusterID
def get_strand(x):
    try:
        strand = x.split(',')[1]
    except IndexError:
        strand = None
    return strand

_chr = [x.split(':')[0] for x in all['clusterID'].values]
_strand = [get_strand(x) for x in all['clusterID'].values]
all['chr'] = _chr
all['strand'] = _strand

# Split data into train, test, others
trainind = np.asarray(~all.loc[:, 'chr'].isin(['chrX', 'chrY', 'chrM', 'chr8']))
testind = np.asarray(all.loc[:, 'chr'] == 'chr8')
otherind = np.asarray(all.loc[:, 'chr'].isin(['chrX', 'chrY', 'chrM']))
X_train = all[Xreducedall_cols].loc[trainind]
info_train = all.drop(columns=Xreducedall_cols).loc[trainind]
X_test = all[Xreducedall_cols].loc[testind]
info_test = all.drop(columns=Xreducedall_cols).loc[testind]
X_other = all[Xreducedall_cols].loc[otherind]
info_other = all.drop(columns=Xreducedall_cols).loc[otherind]

# NOTE: This train samples should be further splited into train and valid in a training script

# Save X and Y as hdf5 format for memory saving mode (save object is not torch object)
with h5py.File(args.output + 'X_Y.h5', 'w') as f:
    # train
    f.create_group('train')
    f['train'].create_dataset('X', data = X_train, compression = "gzip", compression_opts = 9)
    # test
    f.create_group('test')
    f['test'].create_dataset('X', data = X_test, compression = "gzip", compression_opts = 9)
    # others
    f.create_group('other')
    f['other'].create_dataset('X', data = X_other, compression = "gzip", compression_opts = 9)

# Save info, X column, Y column, common rows as txt.gz file
info_train.to_csv(args.output + 'info_train.txt.gz', header = True, compression = 'gzip', sep = "\t", index=False)
info_test.to_csv(args.output + 'info_test.txt.gz', header = True, compression = 'gzip', sep = "\t", index=False)
info_other.to_csv(args.output + 'info_other.txt.gz', header = True, compression = 'gzip', sep = "\t", index=False)
pd.DataFrame(X_train.columns).to_csv(args.output + 'X_train_col.txt.gz', header = False, compression = 'gzip', sep = "\t", index=False)
pd.DataFrame(X_test.columns).to_csv(args.output + 'X_test_col.txt.gz', header = False, compression = 'gzip', sep = "\t", index=False)
pd.DataFrame(X_other.columns).to_csv(args.output + 'X_other_col.txt.gz', header = False, compression = 'gzip', sep = "\t", index=False)
pd.DataFrame(info_train['clusterID']).to_csv(args.output + 'common_train_rows.txt.gz', header = False, compression = 'gzip', sep = "\t", index=False)
pd.DataFrame(info_test['clusterID']).to_csv(args.output + 'common_test_rows.txt.gz', header = False, compression = 'gzip', sep = "\t", index=False)
pd.DataFrame(info_other['clusterID']).to_csv(args.output + 'common_other_rows.txt.gz', header = False, compression = 'gzip', sep = "\t", index=False)

# data metrics
metrics = pd.DataFrame([[len(X_train), len(X_test), len(X_other)]])
metrics.columns = ['N_train', 'N_test', 'N_other']
metrics.to_csv(args.output + 'metrics.txt', header = True, sep = "\t", index=False)

print('=== Finish ===')
