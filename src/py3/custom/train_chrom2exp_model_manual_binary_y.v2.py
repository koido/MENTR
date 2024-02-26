# -*- coding: utf-8 -*-
"""Train MENTR (customized for CD4T PJ)
"""
import argparse
import xgboost as xgb
import numpy as np
import os
import sys

# silent future warning
import warnings
warnings.simplefilter("ignore", category = FutureWarning)
import h5py
import pandas as pd

def mkdir_p(path):
    if not os.path.isdir(path):
        os.makedirs(path)

parser = argparse.ArgumentParser(description='Train MENTR models.')

# settings for y
parser.add_argument('--sample_name', required=True,
                    help='Name of training sample. e.g., C0_Naive.')

# output setting
parser.add_argument('--out_dir', required=True,
                    help='Output directory')
parser.add_argument('--out_suffix',
                    help='Suffix for output files.')

# setting for X and y
parser.add_argument('--infile_dir', required=True,
                    help='Input directory containing X_Y.h5, and info_*.txt.gz')
parser.add_argument('--y', required=True,
                    help='Case and Control ID list (col1: cluster ID, col2: case or control (int)); no header.')

# training parameters
parser.add_argument('--filterStr', default='all', choices=['all', 'promoter', 'enhancer'],
                    help='Filtering for training; all or promoter or enhancer')
parser.add_argument('--num_round', type=int, default=500,
                    help='Max round for training ML models.')
parser.add_argument('--early_stopping_rounds', type=int, default=10,
                    help='Validation score needs to improve at least every "early_stopping_rounds" round(s) to continue training.')
parser.add_argument('--l1', type=int, default=0,
                    help='Alpha (L1 for regularization)')
parser.add_argument('--l2', type=int, default=50,
                    help='Lambda (L2 for regularization)')
parser.add_argument('--eta', type=float, default=0.01,
                    help='Eta (learning rate)')
parser.add_argument('--max_depth', type=int, default=4,
                    help='Max depth for GBT. Used if --gbt is set.')
parser.add_argument('--base_score', type=float, default=0.5,
                    help='Initial pred_y value')
parser.add_argument('--evalauc', action='store_true',
                    help='AUROC is used for evaluation for validation data (it work with --logistic); MENTR use this.')
parser.add_argument('--gbt', action='store_true',
                    help='Gradient Boosted Trees mode; MENTR use this.')
parser.add_argument('--threads', type=int, default=1,
                    help='Number of CPUs for xgboost.')
parser.add_argument('--seed', type=int, default=12345)
parser.add_argument('--verbose', action='store_true')

args = parser.parse_args()

if args.verbose:
    print(args)

if args.out_suffix is None:
    args.out_suffix = ""

# numpy setting
np.random.seed(args.seed)

# Input datasets (h5 format)
hdf5_f = h5py.File(args.infile_dir + '/X_Y.h5', 'r')

objective = 'binary:logistic'
objective_name = 'logistic'

# check and setting
if args.gbt:
    if args.verbose:
        print("Gradient boosted trees mode.")
    booster = 'gbtree'
else:
    if args.verbose:
        print("Linear effects only.")
    booster = 'gblinear'

# read resources for training
if args.verbose:
    print("Loading X train...")
X_train = hdf5_f['train']['X'].value

# row names (cluster ID)
X_train_names = pd.read_csv(args.infile_dir + '/common_train_rows.txt.gz', delimiter = '\t', names = ['clusterID'], dtype = {'clusterID': str}, header=None)
assert X_train.shape[0] == X_train_names.shape[0]

if args.verbose:
    print("Loading y train...")

# y
## col1: cluster ID, col2: case (1) or control (0) (int)
y = pd.read_csv(args.y, delimiter = '\t', names = ['clusterID', 'label'], dtype = {'clusterID': str, 'label': int}, header=None)

if args.verbose:
    print("head of y after loading:")
    print(y.head())

info_train = pd.read_csv(args.infile_dir + '/info_train.txt.gz', dtype = 'str', delimiter = '\t')

assert X_train.shape[0] == info_train.shape[0]

# cluster type-based selection
if args.filterStr == 'promoter':
    filt_train = np.asarray(info_train.loc[:, 'type'].isin(['promoter']))
elif args.filterStr == 'enhancer':
    filt_train = np.asarray(info_train.loc[:, 'type'].isin(['enhancer']))
elif args.filterStr == 'all':
    filt_train = np.asarray(info_train.loc[:, 'type'].isin(['promoter', 'enhancer']))
    assert filt_train.shape[0] == info_train.shape[0]
else:
    raise ValueError('filterStr has to be one of all, promoter, enhancer')

# Add case/control status to X_train
# col1: cluster ID, col2: case (1) or control (0) (int)
# colurmn order is the same as X_train
y_train = pd.merge(X_train_names, y, on = 'clusterID', how = 'left')
assert y_train.shape[0] == X_train.shape[0]
# Check: the colurmn order of y_train is the same as X_train
assert np.all(y_train.loc[:, 'clusterID'] == X_train_names.loc[:, 'clusterID'])

# get clusterID values from y_train with label is not nan
filt_non_nan = ~np.isnan(y_train.loc[:, 'label'])
filt_train = filt_train * filt_non_nan.reshape(-1)

# Filtering X_train and y_train by the above QC
if args.verbose:
    print("Updating X_train and y_train.")
    print("Shape of X_train (before filtering):", X_train.shape)
    print("Shape of y_train (before filtering):", y_train.shape)
X_train = X_train[filt_train, :]
y_train = y_train[filt_train]

if args.verbose:
    print("Shape of X_train (after filtering):", X_train.shape)
    print("Shape of y_train (after filtering):", y_train.shape)

# Convert y_train into numpy array
# but only use column 2 (case/control status)
y_train = y_train.loc[:, 'label'].values

# split train data into new train and validation data for tuning [80%, 20%]
indices = np.random.permutation(X_train.shape[0])
train_idx, valid_idx = indices[:int(len(indices) * 0.8)], indices[int(len(indices) * 0.8):]
X_train_train = X_train[train_idx]
y_train_train = y_train[train_idx]
X_train_valid = X_train[valid_idx]
y_train_valid = y_train[valid_idx]

if args.verbose:
    print("Shape of X_train_train:", X_train_train.shape)
    print("Shape of y_train_train:", y_train_train.shape)
    print("Shape of X_train_valid:", X_train_valid.shape)
    print("Shape of y_train_valid:", y_train_valid.shape)
    print("Head of y_train_train:")
    print(y_train_train[:10])


# make DMatrix
if args.verbose:
    print("Converting DMatrix...")
dtrain = xgb.DMatrix(X_train_train)
dtest = xgb.DMatrix(X_train_valid)
if args.verbose:
    print("Setting labels...")
dtrain.set_label(y_train_train)
dtest.set_label(y_train_valid)

# setting
evallist = [(dtrain, 'train'), (dtest, 'eval')]

# training xgboost
def training_xgb(l1, l2, eta, max_depth = None, args = args, dtrain = dtrain, booster = booster, evalauc = args.evalauc):
    param = {'booster': booster,
             'nthread': args.threads,
             'seed': args.seed,
             'base_score': args.base_score,
             'alpha': l1,
             'lambda': l2,
             'eta': eta}

    param['objective'] = 'binary:logistic'


    if booster == 'gbtree':
        param['max_depth'] = max_depth
    if evalauc:
        param['eval_metric'] = 'auc'

    if args.verbose:
        print(param)

    # run
    evals_result = {}
    if evalauc == False:
        bst = xgb.train(param, dtrain, args.num_round, evallist, evals_result = evals_result, early_stopping_rounds = args.early_stopping_rounds)
    else:
        bst = xgb.train(param, dtrain, args.num_round, evallist, evals_result = evals_result, early_stopping_rounds = args.early_stopping_rounds, maximize=True)
    return bst, evals_result

# No tuning
if args.verbose:
    print("No tuning.")

# default parameters
default_param = {
    'l1': args.l1, # alpha
    'l2': args.l2, # lambda
    'eta': args.eta, # learning rate
}

# save ...
cmn_header = args.out_dir + '/' + args.out_suffix + \
    booster + '_' + objective_name
cmn_header = cmn_header + \
    '_sample_name.' + args.sample_name
cmn_header = cmn_header + \
    '_filterStr.' + args.filterStr + \
    '_num_round.' + str(int(args.num_round)) + \
    '_early_stopping_rounds.' + str(int(args.early_stopping_rounds)) + \
    '_base_score.' + str(args.base_score) + \
    '_l1.' + str(args.l1) + \
    '_l2.' + str(args.l2) + \
    '_eta.' + str(args.eta) + \
    '_seed.' + str(int(args.seed))

# gbtree mode
if booster == 'gbtree':
    default_param['max_depth'] = int(args.max_depth)
    cmn_header = cmn_header + '_max_depth.' + str(int(args.max_depth))
else:
    default_param['max_depth'] = None

if args.evalauc:
    cmn_header = cmn_header + '_evalauc'

# train using the default param
bst, evals_result = training_xgb(l1 = default_param['l1'], l2 = default_param['l2'],
                                 eta = default_param['eta'], max_depth = default_param['max_depth'])

# save model
outfile_path = cmn_header + '_model'
bst.save_model(outfile_path + '.save')
bst.dump_model(outfile_path + '.dump')

# save prediction accuracies
if args.evalauc:
    out_auc = pd.DataFrame({'train': evals_result['train']['auc'],
                            'test': evals_result['eval']['auc']})
    out_auc.to_csv(cmn_header + '.auc.txt.gz', header=True, index=False, sep='\t', compression = 'gzip')
else:
    out_error = pd.DataFrame({'train': evals_result['train']['error'],
                            'test': evals_result['eval']['error']})
    out_error.to_csv(cmn_header + '.error.txt.gz', header=True, index=False, sep='\t', compression = 'gzip')

if args.verbose:
    print("Training is end.")

# save expr and its prediction for all data in a file
for type in ['train', 'test', 'other']:

    if args.verbose:
        print("Evaluating {} data".format(type))

    # load X
    X = hdf5_f[type]['X'].value
    X_names = pd.read_csv(args.infile_dir + '/common_' + type + '_rows.txt.gz', delimiter = '\t', names = ['clusterID'], dtype = {'clusterID': str}, header=None)
    assert X.shape[0] == X_names.shape[0]

    # load info
    info = pd.read_csv(args.infile_dir + '/info_' + type + '.txt.gz', dtype = 'str', delimiter = '\t')
    # preprocess
    y_label = pd.merge(X_names, y, on = 'clusterID', how = 'left')
    assert y_label.shape[0] == X.shape[0]
    # Check: the colurmn order of y_train is the same as X_train
    assert np.all(y_label.loc[:, 'clusterID'] == X_names.loc[:, 'clusterID'])

    dall = xgb.DMatrix(X)
    # infer
    ypred = bst.predict(dall, ntree_limit = bst.best_ntree_limit)
    info['expr'] = y_label['label']
    info['pred'] = ypred
    # write
    out_name = cmn_header + '_type.' + type + '.expr_pred.txt.gz'
    info.to_csv(out_name, header=True, index=False, sep='\t', compression = 'gzip')

print('=== Finish ===')
