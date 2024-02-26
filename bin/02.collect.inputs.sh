#!/bin/bash
#
# 02.collect.inputs.sh
#   Preprocessing: collect information of input files
#
#   Usage:
#     002.collect.inputs.sh [-o cmn_dir] ...
#       -o TEXT  output dir [required]
#       -g INT   Number of GPUs which you have in the node. If >1 value is specified, ${g} GPUs on the node will be used. [default:1]
#       -j INT   Number of jobs per 1-GPU. If >1 value is specified, ${j} jobs will use single GPU, simultaneously. [default:1]
#       -d       Dry Run mode
#       -h       Show help (this message)

# Move to MENTR/bin dir
cd `dirname $0`

# == utils ==
source ../src/sh/utils.sh

# == arg parser ==

# init
DryRun=0
n_gpu=1
job_1gpu=1

# get
while getopts ":o:g:j:dh" optKey; do
  case "${optKey}" in
    o)
      cmn_dir="${OPTARG}";;
    g)
      n_gpu="${OPTARG}";;
    j)
      job_1gpu="${OPTARG}";;
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

# INT check
int_chk ${n_gpu} "-g"
int_chk ${job_1gpu} "-j"

# Show args
echo -e "# == Args for 02.collect.inputs.sh =="
echo -e "# cmn_dir: ${cmn_dir}"
echo -e "# n_gpu: ${n_gpu}"
echo -e "# job_1gpu: ${job_1gpu}"
echo -e "# DryRun: ${DryRun}"
echo -e ""

# == Dry RUN ==
set -eu
if [ ${DryRun} -eq 1 ]; then
  echo -e "# ** Dry Run **"
  echo -e "bash ../src/sh/insilico_mutgen_prep_env.sh ${cmn_dir} ${n_gpu} ${job_1gpu}"
  exit 1
fi

# == RUN ==
bash ../src/sh/insilico_mutgen_prep_env.sh ${cmn_dir} ${n_gpu} ${job_1gpu}
