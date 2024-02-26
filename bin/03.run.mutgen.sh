#!/bin/bash
#
# 03.run.mutgen.sh
#   Run in silico mutagenesis
#
#   Usage:
#     03.run.mutgen.sh [-o cmn_dir] ...
#       -o TEXT  output dir [required]
#       -t INT   Number of CPUs for xgboost. [default:1]
#       -j INT   Number of jobs per 1-GPU. If >1 value is specified, ${j} jobs will use single GPU, simultaneously.
#                **Must be the same with 02.collect.inputs.sh** [default:1]
#       -m TEXT  File of DeepSEA Beluga model (deepsea.beluga.2002.cpu) [required]
#       -M TEXT  Directory for xgboost models
#       -f TEXT  File of reference fasta file [required]
#       -l TEXT  xgboost logisit regression model file [required]
#       -q TEXT  calibrated model list file [required]
#       -P TEXT  peak info file [required]
#       -d       Dry Run mode
#       -c       Only CPU mode
#       -h       Show help (this message)

# Move to MENTR/bin dir
cd `dirname $0`

# == utils ==
source ../src/sh/utils.sh

# == arg parser ==

# init
cpu=0
DryRun=0
n_cpu=1
job_1gpu=1
ml_dir=../

# get
while getopts ":o:t:j:m:M:f:l:q:P:cdh" optKey; do
  case "${optKey}" in
    o)
      cmn_dir="${OPTARG}";;
    t)
      n_cpu="${OPTARG}";;
    j)
      job_1gpu="${OPTARG}";;
    m)
      deepsea="${OPTARG}";;
    M)
      ml_dir="${OPTARG}";;
    f)
      ref_fa="${OPTARG}";;
    l)
      model_list_log="${OPTARG}";;
    q)
      calib_modelList="${OPTARG}";;
    P)
      peak_info="${OPTARG}";;
    c)
      cpu=1;;
    d)
      DryRun=1;;
    h)
      help `basename $0`;;
    :)
      echo -e "Error: Undefined options were observed.\n"
      help `basename $0`;;
    \?)
      echo -e "Error: Undefined options were specified.\n"
      help `basename $0`;;
  esac
done

# required args
req_arg ${cmn_dir} "o"
req_arg ${deepsea} "m"
req_arg ${ref_fa} "f"
req_arg ${model_list_log} "l"
req_arg ${calib_modelList} "q"
req_arg ${peak_info} "P"

# INT check
int_chk ${n_cpu} "-t"
int_chk ${job_1gpu} "-j"

# file check
err_miss ${ref_fa}
err_miss ${deepsea}
err_miss ${model_list_log}
err_miss ${calib_modelList}
err_miss ${peak_info}

# main scripts
if [ ${cpu} -eq 1 ]; then
  sh_script=../src/sh/insilico_mutgen_run_cpu.sh
else
  sh_script=../src/sh/insilico_mutgen_run.sh
fi
py_script=../src/py3/mutgen_nonlinear.py

# Show args
echo -e "# == Args for 03.run.mutgen.sh =="
echo -e "# cmn_dir: ${cmn_dir}"
echo -e "# ml_dir: ${ml_dir}"
echo -e "# deepsea: ${deepsea}"
echo -e "# ref_fa: ${ref_fa}"
echo -e "# job_1gpu: ${job_1gpu}"
echo -e "# n_cpu: ${n_cpu}"
echo -e "# model_list_log: ${model_list_log}"
echo -e "# calib_modelList: ${calib_modelList}"
echo -e "# peak_info: ${peak_info}"
echo -e "# DryRun: ${DryRun}"
echo -e "# cpu: ${cpu}"
echo -e "# py_script: ${py_script}"
echo -e "# sh_script: ${sh_script}; This was automatically set."
echo -e ""

# == Dry RUN ==
set -eu
if [ ${DryRun} -eq 1 ]; then
  echo -e "# ** Dry Run **"
  echo -e "bash ${sh_script} ${cmn_dir} ${ml_dir} ${ref_fa} ${deepsea} ${job_1gpu} ${n_cpu} ${model_list_log} ${calib_modelList} ${peak_info} ${py_script}"
  exit 1
fi

# == RUN ==
bash ${sh_script} ${cmn_dir} ${ml_dir} ${ref_fa} ${deepsea} ${job_1gpu} ${n_cpu} ${model_list_log} ${calib_modelList} ${peak_info} ${py_script}
