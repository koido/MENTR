# MENTR (mutation effect prediction on ncRNA transcription)

Please cite the following papers when you use MENTR:

 - [MENTR paper](https://www.nature.com/articles/s41551-022-00961-8): M Koido et al. Prediction of the cell-type-specific transcription of non-coding RNAs from genome sequences via machine learning. Nat. Biomed. Eng. 2022. [SharedIt URL](https://rdcu.be/c26Uj), [News&View in NBME](https://rdcu.be/c26TZ)
 - FANTOM5 papers ([FANTOM CAT](http://dx.doi.org/10.1038/nature21374), [Promoter atlas](https://doi.org/10.1038/nature13182), [Enhancer atlas](https://doi.org/10.1038/nature12787)) and [LCL CAGE QTL paper](https://doi.org/10.1038/s41467-017-01467-7), whose resources we used for training MENTR ML models   
 - [ExPecto paper](https://doi.org/10.1038/s41588-018-0160-6) (We used the pre-trained DeepSEA Beluga model, which this repository does NOT contain, and script for utilizing the model)

This README briefly describes how to run MENTR (in silico mutagenesis).
The details of output files, how to customize parameters and settings, how to train MENTR ML models, are written in [our wiki](https://github.com/koido/MENTR/wiki).

# 1. How to perform in silico mutagenesis

## 1.1. Dependencies

### Software and libraries

 - python3 (see mentr_env.txt for dependent python libraries)

You can create `mentr` environment in conda:

```bash
conda create --name mentr -c conda-forge -c bioconda --file mentr_env.txt python=3.6.4
```

NOTE: Docker or Singularity image including the depending software is available via Docker Hub (https://hub.docker.com/repository/docker/mkoido/mentr).
You can find the same structure of this MENTR repository under `/MENTR`.
In addition to the image, you need to get the external resources (see [here](https://github.com/koido/MENTR#external-resources)).

### External resources

 - hg19.fa
 - deepsea.beluga.2002.cpu (from https://github.com/FunctionLab/ExPecto, especially [here]( https://github.com/FunctionLab/ExPecto/blob/447737793d8d21e50e82379feca44fb9465fdc79/download_resources.sh); `wget http://deepsea.princeton.edu/media/code/expecto/resources.tar.gz`)

These are required even if you use Docker or Singularity image.

## 1.2. Usage

### Step 0) Set environment

```bash
# Please specify
cmn_dir="where you want to write results"
deepsea="where deepsea.beluga.2002.cpu is"
hg19_fa="where hg19.fa is"
```

### Step 1) Prepare input file

TSV file with:
col1: Chromosome number
col2: Position (hg19)
col3: Reference allele (hg19)
col4: Alternate allele (hg19)

The default setting assumes no header line.
see example file (./resources/example/input.txt.gz)

If you want to analyze multiallelic sites, please split them into multiple rows.
For col3 and col4, SNVs and short InDels are acceptable.

```bash
# input file
in_f=./resources/example/input.txt.gz
```

### Step 2) Run MENTR

Here is the single command, quick run method.
 - Using default parameters and settings of MENTR
 - Need 1-GPU and 1-CPU.

```bash
bash bin/quick.mutgen.sh -i ${in_f} -o ${cmn_dir} -m ${deepsea} -f ${hg19_fa}
```

You can get all_qcd_res.txt.gz.

See [our wiki](https://github.com/koido/MENTR/wiki/in-silico-mutagenesis) to know the output files and alternative ways to run MENTR with other parameters/settings or with much faster.

# 2. License

This project is licensed under GNU GPL v3.

