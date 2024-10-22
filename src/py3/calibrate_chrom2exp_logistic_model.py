# -*- coding: utf-8 -*-
"""Calibrate MENTR
"""

import argparse
import numpy as np
import os
import sys

# sklearn
from sklearn.isotonic import IsotonicRegression
from sklearn.externals import joblib

# silent future warning
import warnings
warnings.simplefilter("ignore", category = FutureWarning)
import pandas as pd

parser = argparse.ArgumentParser(description='Calibrate MENTR models.')

# input setting
parser.add_argument('--xgb_res_path', required=True,
                    default='Specify MENTR result files. Pask train/test/other by "{}". e.g., "..._evalauc_type.{}.expr_pred.txt.gz"')

# output setting
parser.add_argument('--out_dir', required=True,
                    help='Output directory')
parser.add_argument('--out_suffix',
                    help='Suffix for output files.')

# training parameters
parser.add_argument('--filterStr', type=str, default = 'all',
                    help='all or promoter or enhancer. This must be the same for training MENTR model.')
parser.add_argument('--seed', type=int, default=12345,
                    help='This must be the same seed for training MENTR model.')
parser.add_argument('--verbose', action = 'store_true')

args = parser.parse_args()
if args.out_suffix is None:
    args.out_suffix = ""
out_cmn=args.out_dir + '/' + args.out_suffix
xgb_res_path = args.xgb_res_path

# numpy setting
np.random.seed(args.seed)

# load xgboost-predicted data
## define dtype
train_data_header = pd.read_csv(xgb_res_path.format('train'), dtype = 'str', delimiter = '\t', nrows=1, header=None).values.reshape(-1).tolist()
dtypes = {}
for x in train_data_header:
    if x in ['expr', 'pred']:
        dtypes[x] = 'float64'
    else:
        dtypes[x] = 'str'
# load
train_data = pd.read_csv(xgb_res_path.format('train'), dtype = dtypes, delimiter = '\t')
y_train = train_data['expr']
p_train = train_data['pred']

# cluster type-based selection
if args.filterStr == 'promoter':
    filt_train = np.asarray(train_data.loc[:, 'type'].isin(['promoter']))
elif args.filterStr == 'enhancer':
    filt_train = np.asarray(train_data.loc[:, 'type'].isin(['enhancer']))
elif args.filterStr == 'all':
    filt_train = np.asarray(train_data.loc[:, 'type'].isin(['promoter', 'enhancer']))
    assert filt_train.shape[0] == train_data.shape[0]
else:
    raise ValueError('filterStr has to be one of all, promoter, enhancer')

# gene QC
filt_train = filt_train * \
    np.isfinite(y_train).reshape(-1)

# update X_train and y_train
if args.verbose:
    print("Updating y_train and p_train.")
    print("Shape of y_train (before filtering):", y_train.shape)
    print("Shape of p_train (before filtering):", p_train.shape)
y_train = y_train[filt_train]
p_train = p_train[filt_train]

if args.verbose:
    print("Shape of y_train (after filtering):", y_train.shape)
    print("Shape of p_train (after filtering):", p_train.shape)

# split train data into new train and validation data for tuning [80%, 20%]
indices = np.random.permutation(y_train.shape[0])
train_idx, valid_idx = indices[:int(len(indices) * 0.8)], indices[int(len(indices) * 0.8):]
y_train_valid = y_train[valid_idx]
p_train_valid = p_train[valid_idx]

# binarize, y
y_train_valid = (y_train_valid > 0).astype(int)

# Isotonic Regression
ir = IsotonicRegression(y_min = 0, y_max = 1, out_of_bounds = 'clip')
ir.fit(p_train_valid, y_train_valid)

# save joblib
joblib.dump(ir, out_cmn+'_IsotonicRegression.joblib')

## Save calibrated values
for type in ['train', 'test', 'other']:

    if args.verbose:
        print("Calibrating {} data".format(type))

    ## define dtype
    train_data_header = pd.read_csv(xgb_res_path.format(type), dtype = 'str', delimiter = '\t', nrows=1, header=None).values.reshape(-1).tolist()
    dtypes = {}
    for x in train_data_header:
        if x in ['expr', 'pred']:
            dtypes[x] = 'float64'
        else:
            dtypes[x] = 'str'

    # load
    data = pd.read_csv(xgb_res_path.format(type), dtype = dtypes, delimiter = '\t')
    y = data['expr']
    p = data['pred']

    # calibration & order check
    p_c = ir.transform(p)

    # binarize, y
    y_binary = (y > 0).astype(int)

    # update data
    data['pred'] = p_c

    # write
    out_name = out_cmn + '_type.' + type + '_expr_pred_c.txt.gz'
    data.to_csv(out_name, header=True, index=False, sep='\t', compression = 'gzip')

print('=== Finish ===')
