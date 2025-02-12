# -*- coding: utf-8 -*-
"""Train MENTR
"""
import argparse
import xgboost as xgb
import numpy as np
import sys

# silent future warning
import warnings
warnings.simplefilter("ignore", category = FutureWarning)
import h5py
import pandas as pd

parser = argparse.ArgumentParser(description='Train MENTR models.')

# output setting
parser.add_argument('--out_dir', required=True,
                    help='Output directory')
parser.add_argument('--out_suffix',
                    help='Suffix for output files.')

# setting for X and y
parser.add_argument('--infile_dir', required=True,
                    help='Input directory containing X_Y.h5, Y_*_col.txt.gz, and info_*.txt.gz')

# training parameters
parser.add_argument('--pseudocount', type=float, default=0.0001,
                    help='Pseudo-count for y; Used if --logistic is not specified.')
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
parser.add_argument('--logistic', action='store_true',
                    help='logistic regression mode (or classification mode); MENTR use this.')
parser.add_argument('--evalauc', action='store_true',
                    help='AUROC is used for evaluation for validation data (it work with --logistic); MENTR use this.')
parser.add_argument('--gbt', action='store_true',
                    help='Gradient Boosted Trees mode; MENTR use this.')
parser.add_argument('--threads', type=int, default=1,
                    help='Number of CPUs for xgboost.')
parser.add_argument('--seed', type=int, default=12345)
parser.add_argument('--shorter_name', action='store_true',
                    help='Shorter name for output files.')

# args for custom training
parser.add_argument('--sample_ID',
                    help='Sample ID in train_clusterID_Y')
# resources/GSE92743/mk/GSE92743_Broad_GTEx_L1000_Level3_Q2NORM_n3176x12320_lincs2020annot_unique_fantom5_annot.tsv.gz
parser.add_argument('--train_clusterID_Y', required=True,
                    help='ClusterID + Y TSV file to use for training(+validation); Need header.; col1 should be clusterID.; col2- should be Y.')
parser.add_argument('--landmark', action='store_true',
                    help='Use landmark only for training.')
parser.add_argument('--highest_promoter', action='store_true',
                    help='Use highest_promoter only for training.')
parser.add_argument('--no_logarithm', action='store_true',
                    help='Do not use logarithm for training.')
parser.add_argument('--verbose', action='store_true')

args = parser.parse_args()
if args.out_suffix is None:
    args.out_suffix = ""

# numpy setting
np.random.seed(args.seed)

# Input datasets (h5 format)
hdf5_f = h5py.File(args.infile_dir + '/X_Y.h5', 'r')

if args.logistic == False and args.evalauc == True:
    raise ValueError("Please use --evalauc with --logistic")

if args.logistic:
    if args.verbose:
        print("Logistic regression mode.")
    assert args.base_score >= 0 and args.base_score <= 1
    objective = 'binary:logistic'
    objective_name = 'logistic'
else:
    if args.verbose:
        print("Linear regression mode.")
    objective = 'reg:linear'
    objective_name = 'linear'

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
X_train = hdf5_f['train']['X'][()]
info_train_FANTOM5 = pd.read_csv(args.infile_dir + '/info_train.txt.gz', dtype = 'str', delimiter = '\t')
filt_train_df = pd.DataFrame({'clusterID': info_train_FANTOM5['clusterID']})

if args.verbose:
    print("Loading y train...")

# read train_clusterID_Y
Y_train = pd.read_csv(args.train_clusterID_Y, delimiter='\t', header=0)
# rename the first column name
Y_train.columns.name = 'clusterID'
# Check if the sample ID is in the Y_train
assert args.sample_ID in Y_train.columns
# Check if clusterID is primary key
assert Y_train['clusterID'].is_unique, "clusterID is not unique in train_clusterID_Y"

if args.landmark:
    Y_train = Y_train[Y_train['feature_space'] == 'landmark']

if args.highest_promoter:
    Y_train = Y_train[Y_train['highest_promoter'] == 1]

# Extract "clusterID" and args.sample_ID from Y_train in a DataFrame
y_train_df = Y_train[['clusterID', args.sample_ID]]
if args.verbose:
    print("Shape of y_train_df:", y_train_df.shape)
    # basic statistics of "args.sample_ID" column
    print("Basic statistics of {} column:".format(args.sample_ID))
    print(y_train_df[args.sample_ID].describe())
y_train_filter = pd.DataFrame({'clusterID': y_train_df['clusterID'], 'filter': True})

# update filt_train_df (but keeping the order of filt_train_df)
filt_train_df = filt_train_df.merge(y_train_filter, on='clusterID', how='left')
filt_train_df['filter'] = filt_train_df['filter'].fillna(False)
filt_train = np.asarray(filt_train_df['filter'])

# Filtering X_train and y_train by the above QC
if args.verbose:
    print("Updating X_train.")
    print("Shape of X_train (before filtering):", X_train.shape)
    print("Shape of info_train_FANTOM5 (before filtering):", info_train_FANTOM5.shape)
X_train = X_train[filt_train, :]
info_train_FANTOM5 = info_train_FANTOM5.loc[filt_train]

if args.verbose:
    print("Shape of X_train (after filtering):", X_train.shape)
    print("Shape of info_train_FANTOM5 (after filtering):", info_train_FANTOM5.shape)

# reorder of y_train_df (the same order of clusterID in the filtered info_train_FANTOM5)
y_train = info_train_FANTOM5.loc[:, ["clusterID"]].merge(y_train_df, on='clusterID', how='left')
# remove the clusterID column
y_train = y_train.drop(columns=['clusterID'])
# Then nunber of column should be 1
assert y_train.shape[1] == 1, "Number of column should be 1"
# as numpy array (1D)
y_train = y_train.values.reshape(-1)

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

if args.logistic:
    # binarize
    if args.verbose:
        print("Binarizing...")
    y_train_train = (y_train_train > 0).astype(int)
    y_train_valid = (y_train_valid > 0).astype(int)

if args.verbose:
    print("Converting DMatrix...")

# make DMatrix
dtrain = xgb.DMatrix(X_train_train)
dvalid = xgb.DMatrix(X_train_valid)

if args.logistic:
    # make y
    dtrain.set_label(y_train_train)
    dvalid.set_label(y_train_valid)
else:
    if args.no_logarithm:
        if args.verbose:
            print("No logarithm.")
        dtrain.set_label(y_train_train + args.pseudocount)
        dvalid.set_label(y_train_valid + args.pseudocount)
    else:
        # make y
        dtrain.set_label(np.log(y_train_train + args.pseudocount))
        dvalid.set_label(np.log(y_train_valid + args.pseudocount))

# setting
evallist = [(dtrain, 'train'), (dvalid, 'eval')]

# training xgboost
def training_xgb(l1, l2, eta, max_depth = None, args = args, dtrain = dtrain, logistic = args.logistic, booster = booster, evalauc = args.evalauc):
    param = {'booster': booster,
             'nthread': args.threads,
             'seed': args.seed,
             'base_score': args.base_score,
             'alpha': l1,
             'lambda': l2,
             'eta': eta}
    if logistic:
        param['objective'] = 'binary:logistic'
    else:
        param['objective'] = 'reg:linear'
    if booster == 'gbtree':
        param['max_depth'] = max_depth
    if evalauc:
        param['eval_metric'] = 'auc'

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
    booster + '_' + objective_name + \
    '_sample_ID.' + args.sample_ID

if not args.shorter_name:
    cmn_header = cmn_header + \
        '_num_round.' + str(int(args.num_round)) + \
        '_early_stopping_rounds.' + str(int(args.early_stopping_rounds)) + \
        '_base_score.' + str(args.base_score) + \
        '_l1.' + str(args.l1) + \
        '_l2.' + str(args.l2) + \
        '_eta.' + str(args.eta) + \
        '_seed.' + str(int(args.seed))

if args.logistic == False:
    cmn_header = cmn_header + \
        '_pseudocount.' + str(args.pseudocount)

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
elif args.logistic:
    out_error = pd.DataFrame({'train': evals_result['train']['error'],
                            'test': evals_result['eval']['error']})
    out_error.to_csv(cmn_header + '.error.txt.gz', header=True, index=False, sep='\t', compression = 'gzip')
else:
    out_rmse = pd.DataFrame({'train': evals_result['train']['rmse'],
                            'test': evals_result['eval']['rmse']})
    out_rmse.to_csv(cmn_header + '.rmse.txt.gz', header=True, index=False, sep='\t', compression = 'gzip')

if args.verbose:
    print("Training is end.")

# save expr and its prediction for all data in a file
for type in ['train', 'test', 'other']:

    if args.verbose:
        print("Evaluating {} data".format(type))

    # load
    X = hdf5_f[type]['X'][()]
    info = pd.read_csv(args.infile_dir + '/info_' + type + '.txt.gz', dtype = 'str', delimiter = '\t')

    dall = xgb.DMatrix(X)
    # infer
    ypred = bst.predict(dall)
    info['pred'] = ypred
    # write
    out_name = cmn_header + '_type.' + type + '.expr_pred.txt.gz'
    info.to_csv(out_name, header=True, index=False, sep='\t', compression = 'gzip')

print('=== Finish ===')
