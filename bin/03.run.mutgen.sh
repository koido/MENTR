#!/bin/bash
#
# 03.run.mutgen.sh
#   Run in silico mutagenesis
#
#   Usage:
#     03.run.mutgen.sh [-o cmn_dir] ...
#       -o TEXT  output dir [required]
#       -c INT   Number of CPUs for xgboost. [default:1]
#       -j INT   Number of jobs per 1-GPU. If >1 value is specified, ${j} jobs will use single GPU, simultaneously. **Must be the same with 02.collect.inputs.sh** [default:1]
#       -m TEXT  File of DeepSEA Beluga model (deepsea.beluga.2002.cpu)
#       -f TEXT  File of hg19.fa
#       -d       Dry Run mode
#       -h       Show help (this message)

# Move to MENTR/bin dir
cd `dirname $0`

# == utils ==
source ../src/sh/utils.sh

# == arg parser ==

# init
DryRun=0
n_cpu=1
job_1gpu=1

# get
while getopts ":o:c:j:m:f:dh" optKey; do
  case "${optKey}" in
    o)
      cmn_dir="${OPTARG}";;
    c)
      n_cpu="${OPTARG}";;
    j)
      job_1gpu="${OPTARG}";;
    m)
      deepsea="${OPTARG}";;
    f)
      hg19_fa="${OPTARG}";;
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
if [ -z "${cmn_dir}" ]; then
  echo "Error: -o is required"
  help `basename $0`
fi
if [ -z "${deepsea}" ]; then
  echo "Error: -m is required"
  help `basename $0`
fi
if [ -z "${hg19_fa}" ]; then
  echo "Error: -f is required"
  help `basename $0`
fi

# INT check
int_chk ${n_cpu} "-c"
int_chk ${job_1gpu} "-j"

# == common args ==
ml_dir=../

# file check
if [ ! -e ${hg19_fa} ]; then
  echo -e "Error: ${hg19_fa} is missing. Please see the Resources in README."
  exit 1
fi
if [ ! -e ${deepsea} ]; then
  echo -e "Error: ${deepsea} is missing. Please see the Resources in README."
  exit 1
fi

# == Dry RUN ==
set -eu
if [ ${DryRun} -eq 1 ]; then
  echo -e "# ** Dry Run **"
  echo -e "# cmn_dir: ${cmn_dir}"
  echo -e "# deepsea: ${deepsea}"
  echo -e "# hg19_fa: ${hg19_fa}"
  echo -e "# job_1gpu: ${job_1gpu}"
  echo -e "# n_cpu: ${n_cpu}"
  echo -e "bash ../src/sh/insilico_mutgen_run.sh ${cmn_dir} ${ml_dir} ${hg19_fa} ${deepsea} ${job_1gpu} ${n_cpu}"
  exit 1
fi

# == RUN ==
bash ../src/sh/insilico_mutgen_run.sh ${cmn_dir} ${ml_dir} ${hg19_fa} ${deepsea} ${job_1gpu} ${n_cpu}

