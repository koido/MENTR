# -*- coding: utf-8 -*-
"""Predict MENTR using REF
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
parser.add_argument('--use_best_itr', action="store_true",
                    help="Use best iteration for xgboost model prediction. default: False")
parser.add_argument('--threads', type=int, default=1,
                    help='Number of CPUs for xgboost.')
parser.add_argument('--seed', type=int, default=12345)
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

info_train = pd.read_csv(args.infile_dir + '/info_train.txt.gz', dtype = 'str', delimiter = '\t')

# reconstract cmn_header
cmn_header = args.out_dir + '/' + args.out_suffix + \
    booster + '_' + objective_name
if args.CAGE_lib_ID is None:
    cmn_header = cmn_header + \
        '_sample_ontology_ID.' + args.sample_ontology_ID
else:
    cmn_header = cmn_header + \
        '_CAGE_lib_ID.' + args.CAGE_lib_ID
cmn_header = cmn_header + \
    '_filterStr.' + args.filterStr + \
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
    cmn_header = cmn_header + '_max_depth.' + str(int(args.max_depth))

if args.evalauc:
    cmn_header = cmn_header + '_evalauc'

# load model
bst = xgb.Booster({'nthread': args.threads})
bst.load_model(cmn_header + '_model' + '.save')

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
    if args.use_best_itr:
        ## https://github.com/dmlc/xgboost/issues/805
        ## https://stackoverflow.com/questions/43534219/xgboost-what-is-the-difference-among-bst-best-score-bst-best-iteration-and-bst
        best_itr = int(bst.attributes()["best_iteration"])
        # num_parallel_tree was not set (default:1)
        num_parallel_tree = 1
        ypred = bst.predict(dall, ntree_limit = (best_itr + 1) * num_parallel_tree)
        if type == 'train':
            print("Best iteration (1-based): {}".format((best_itr + 1) * num_parallel_tree))
        out_name = cmn_header + '_type.' + type + '.best_itr.expr_pred.txt.gz'
    else:
        ypred = bst.predict(dall)
        out_name = cmn_header + '_type.' + type + '.expr_pred.txt.gz'
    if args.logistic:
        info['expr'] = y
    else:
        info['log_expr'] = y
    info['pred'] = ypred
    # write
   
    info.to_csv(out_name, header=True, index=False, sep='\t', compression = 'gzip')

print('=== Finish ===')
