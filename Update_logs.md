# Update logs

## MENTR v1.0.0

 - Initial release.

## MENTR v1.0.1

 - Citation information was updated (bioRxiv preprint → Koido et al. Nat. Biomed. Eng. 2022).
 - The others are the same as the initial release (v.1.0.0).

## MENTR v1.0.2

 - Update options and fix typos in shell scripts
 - Update dependency for conda for Mac Silicon users (CPU mode)
 - Add Python script for probability calibration
 - Pre-trained models are the same as before.
 - Then, these changes do not affect the results of in silico mutagenesis.

## MENTR on branch dev (<font color="red">not incorporated into the main branch yet</font>)

 - Update shell scripts to increase freedom of input and output files
 - Set `ntree_limit` for XGBoost models to use the best iteration, which may change the results of in silico mutagenesis (but maybe better than before)
 - Rename: `seq2chrom_hg19_ref.py` -> `seq2chrom_ref.py`
 - Add:
    - `src/py3/custom/seq2chrom_res2hdf5_nocount.py`
        - Compared to the normal script (`src/py3/seq2chrom_res2hdf5.py`), this script accounts for the situation when users do not use expression levels
        - When we use this to make HDF5 database, the following two files are required.
    - `src/py3/custom/train_chrom2exp_model_manual_binary_y.v2.py
        - Train xgboost model without providing expression levels
        - Users have to provide on/off transcript id information (see help)
    - `src/py3/custom/calibrate_chrom2exp_logistic_model_manual_binary_y.py`
        - Calibrate xgboost model without providing expression levels

### TODO:
 - Run `src/sh/debug/res2hdf5_lowmem.sh`
 - Update wiki





 

