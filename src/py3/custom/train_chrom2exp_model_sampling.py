# -*- coding: utf-8 -*-
"""Train MENTR
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
parser.add_argument('--sample_ontology_ID',
                    help='e.g., CL:0000034')
parser.add_argument('--sample_ontology_File', required=True,
                    help='Path for supp_table_10.sample_ontology_information.tsv.gz')
parser.add_argument('--CAGE_lib_ID',
                    help='e.g., CNhs10608; If this is given, sample ontology ID is NOT used.')

# output setting
parser.add_argument('--out_dir', required=True,
                    help='Output directory')
parser.add_argument('--out_suffix',
                    help='Suffix for output files.')

# setting for X and y
parser.add_argument('--infile_dir', required=True,
                    help='Input directory containing X_Y.h5, Y_*_col.txt.gz, and info_*.txt.gz')

# training parameters
parser.add_argument('--filterStr', default='all', choices=['all', 'promoter', 'enhancer'],
                    help='Filtering for training; all or promoter or enhancer')
parser.add_argument('--geneClassStr', default='all',
                    help='Filtering for training')
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

# args for sampling in train
parser.add_argument('--seed_sampling', type=int, default=12345)
parser.add_argument('--n_sampling', type=int, default=0,
                    help='Number of samples to use for training.; 0 means all samples are used.')
parser.add_argument('--verbose', action='store_true')

args = parser.parse_args()
if args.out_suffix is None:
    args.out_suffix = ""

# numpy setting
np.random.seed(args.seed)

# Input datasets (h5 format)
hdf5_f = h5py.File(args.infile_dir + '/X_Y.h5', 'r')

# check and setting
if args.CAGE_lib_ID is None and args.sample_ontology_ID is None:
    raise ValueError("Please specify either CAGE_lib_ID or sample_ontology_ID")

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

if args.CAGE_lib_ID is None:

    assert args.sample_ontology_ID is not None

    if args.verbose:
        print("Sample Ontology ID is used.")

    # CAGE sample information
    sample_ontology = pd.read_csv(args.sample_ontology_File, dtype = 'str', delimiter = '\t')
    sample_ontology_term = sample_ontology[sample_ontology['sample_ontology_ID'] == args.sample_ontology_ID].loc[:, 'sample_ontology_term'].values.tolist()
    assert len(sample_ontology_term) == 1
    sample_ontology_term = sample_ontology_term[0]
    sample_ontology_type = sample_ontology[sample_ontology['sample_ontology_ID'] == args.sample_ontology_ID].loc[:, 'sample_ontology_type'].values.tolist()
    assert len(sample_ontology_type) == 1
    sample_ontology_type = sample_ontology_type[0]
    CAGE_lib_ID = sample_ontology[sample_ontology['sample_ontology_ID'] == args.sample_ontology_ID].loc[:, 'CAGE_lib_ID'].values.tolist()
    CAGE_lib_ID = CAGE_lib_ID[0].split(',')

    # select samples to aggregate
    Y_train_col = pd.read_csv(args.infile_dir + '/Y_train_col.txt.gz', dtype = 'str', delimiter = '\t', header=None).values.reshape(-1).tolist()
    assert len(set(CAGE_lib_ID) & set(Y_train_col)) == len(CAGE_lib_ID)
    use_Y_train_col = [i for i, x in enumerate(Y_train_col) if x in CAGE_lib_ID]
    assert len(set([Y_train_col[i] for i in use_Y_train_col]) & set(CAGE_lib_ID)) == len(CAGE_lib_ID)

else:

    if args.verbose:
        print("CAGE library ID is used.")

    # select samples to aggregate
    Y_train_col = pd.read_csv(args.infile_dir + '/Y_train_col.txt.gz', dtype = 'str', delimiter = '\t', header=None).values.reshape(-1).tolist()
    use_Y_train_col = [i for i, x in enumerate(Y_train_col) if x == args.CAGE_lib_ID]

# read resources for training
if args.verbose:
    print("Loading X train...")
X_train = hdf5_f['train']['X'].value

if args.verbose:
    print("Loading y train...")

if args.CAGE_lib_ID is None:
    y_train = hdf5_f['train']['y'][:, use_Y_train_col].mean(axis = 1)
else:
    y_train = hdf5_f['train']['y'][:, use_Y_train_col]
info_train = pd.read_csv(args.infile_dir + '/info_train.txt.gz', dtype = 'str', delimiter = '\t')

assert X_train.shape[0] == y_train.shape[0]
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

# geneClassStr filter
if args.geneClassStr != 'all':
    filt_train = filt_train * \
        np.asarray(info_train.loc[:, 'geneClassStr'].isin([args.geneClassStr]))

if args.logistic:
    # gene QC
    filt_train = filt_train * \
        np.isfinite(y_train).reshape(-1)
else:
    # gene QC
    filt_train = filt_train * \
        np.isfinite(np.log(y_train + args.pseudocount)).reshape(-1) * \
        np.isfinite(y_train).reshape(-1)

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

# split train data into new train and validation data for tuning [80%, 20%]
if args.n_sampling == 0:
    indices = np.random.permutation(X_train.shape[0])
    train_idx, valid_idx = indices[:int(len(indices) * 0.8)], indices[int(len(indices) * 0.8):]
else:
    # local seed
    #rng = np.random.default_rng(args.seed_sampling)
    rng = np.random.RandomState(args.seed_sampling)
    indices = rng.choice(X_train.shape[0], size = args.n_sampling, replace = False)
    # From sampled indices, split train (80%) and validation (20%) data
    train_idx = indices[:int(len(indices) * 0.8)]
    valid_idx = indices[int(len(indices) * 0.8):]
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
dtest = xgb.DMatrix(X_train_valid)

if args.logistic:
    # make y
    dtrain.set_label(y_train_train)
    dtest.set_label(y_train_valid)
else:
    # make y
    dtrain.set_label(np.log(y_train_train + args.pseudocount))
    dtest.set_label(np.log(y_train_valid + args.pseudocount))

# setting
evallist = [(dtrain, 'train'), (dtest, 'eval')]

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
    booster + '_' + objective_name
if args.CAGE_lib_ID is None:
    cmn_header = cmn_header + \
        '_sample_ontology_ID.' + args.sample_ontology_ID
else:
    cmn_header = cmn_header + \
        '_CAGE_lib_ID.' + args.CAGE_lib_ID

if not args.shorter_name:
    cmn_header = cmn_header + \
        '_filterStr.' + args.filterStr + \
        '_num_round.' + str(int(args.num_round)) + \
        '_early_stopping_rounds.' + str(int(args.early_stopping_rounds)) + \
        '_base_score.' + str(args.base_score) + \
        '_l1.' + str(args.l1) + \
        '_l2.' + str(args.l2) + \
        '_eta.' + str(args.eta) + \
        '_seed.' + str(int(args.seed))

if args.n_sampling > 0:
    cmn_header = cmn_header + \
        '_n_sampling.' + str(int(args.n_sampling)) + \
        '_seed_sampling.' + str(int(args.seed_sampling))
if args.logistic == False:
    cmn_header = cmn_header + \
        '_pseudocount.' + str(args.pseudocount)
if args.geneClassStr != 'all':
    cmn_header = cmn_header + \
        '_geneClassStr.' + args.geneClassStr

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

    Y_col = pd.read_csv(args.infile_dir + '/Y_' + type + '_col.txt.gz', dtype = 'str', delimiter = '\t', header=None).values.reshape(-1).tolist()
    if args.CAGE_lib_ID is None:      
        assert len(set(CAGE_lib_ID) & set(Y_col)) == len(CAGE_lib_ID)
        use_Y_col = [i for i, x in enumerate(Y_col) if x in CAGE_lib_ID]
        assert len(set([Y_col[i] for i in use_Y_train_col]) & set(CAGE_lib_ID)) == len(CAGE_lib_ID)
    else:
        use_Y_col = [i for i, x in enumerate(Y_col) if x == args.CAGE_lib_ID]

    # load
    X = hdf5_f[type]['X'].value
    y_all = hdf5_f[type]['y'][:, use_Y_col]
    info = pd.read_csv(args.infile_dir + '/info_' + type + '.txt.gz', dtype = 'str', delimiter = '\t')
    # preprocess
    if args.logistic:
        if args.CAGE_lib_ID is None:
            y = y_all.mean(axis = 1)
        else:
            y = y_all
    else:
        if args.CAGE_lib_ID is None:
            y = np.log(y_all.mean(axis = 1) + args.pseudocount)
        else:
            y = np.log(y_all + args.pseudocount)
    dall = xgb.DMatrix(X)
    # infer
    ypred = bst.predict(dall)
    if args.logistic:
        info['expr'] = y
    else:
        info['log_expr'] = y
    info['pred'] = ypred
    # write
    out_name = cmn_header + '_type.' + type + '.expr_pred.txt.gz'
    info.to_csv(out_name, header=True, index=False, sep='\t', compression = 'gzip')

print('=== Finish ===')
