#!/bin/bash
#
# quick.mutgen.sh
#   Run in silico mutagenesis with default parameters
#
#   Usage:
#     quick.mutgen.sh [-i input_file] [-o cmn_dir] ...
#       -i TEXT  input variant file [required]
#       -o TEXT  output dir [required]
#       -m TEXT  File of DeepSEA Beluga model (deepsea.beluga.2002.cpu) [required]
#       -M TEXT  Directory for xgboost models
#       -f TEXT  File of reference fasta [required]
#       -t INT   Number of threads. [default:1]
#       -w INT   window size to find variant-CAGE_ID pairs [bp] [default:100000]
#       -b TEXT  input bed file [default: ALL FANTOM5 CAGE IDs]
#                 col1: chr
#                 col2: start (TSS) or (midpoint - 1)
#                 col3: end (TSS or midpoint)
#                 col4: strand [+/-/N]
#                 col5: CAGE_peak_ID
#       -p TEXT  CAGE_ID list to find variant-CAGE_ID pairs. If many, please write IDs separated by semicolon [default: ALL; -> All FANTOM5 CAGE IDs]
#       -l TEXT  xgboost logisic regression model file
#       -q TEXT  calibrated model list file
#       -P TEXT  peak info file
#       -c       CPU mode
#       -h       Show help (this message)

# Move to MENTR/bin dir
cd `dirname $0`

# == utils ==
source ../src/sh/utils.sh

# == arg parser ==

# init
# Only CPU mode
cpu=0
# Number of CPU
n_cpu=1
# window size [bp]; default is (Promoter TSS / Enhancer midpoint) +/- 100-kb) but you can change this parameter
window_n=100000
# F5 promoter TSS and enhancer mid position (maximum set for ML models; that is, this file do NO include transcripts in edge region o
in_bed=../resources/mutgen_cmn/F5.cage_cluster.hg19.win_100kb.sort.bed
# CAGE ID list for evaluation
cage_list=ALL
# model parameters
ml_dir=../
model_list_log=${ml_dir}/resources/mutgen_cmn/modellist_xgb.txt
calib_modelList=${ml_dir}/resources/mutgen_cmn/modellist_calib.txt
peak_info=${ml_dir}/resources/mutgen_cmn/cage_peak.txt.gz

# get
while getopts ":i:o:m:M:f:t:w:b:p:l:q:P:ch" optKey; do
  case "${optKey}" in
    i)
      in_f="${OPTARG}";;
    o)
      cmn_dir="${OPTARG}";;
    m)
      deepsea="${OPTARG}";;
    M)
      ml_dir="${OPTARG}";;
    f)
      ref_fa="${OPTARG}";;
    t)
      n_cpu="${OPTARG}";;
    w)
      window_n="${OPTARG}";;
    b)
      in_bed="${OPTARG}";;
    p)
      cage_list="${OPTARG}";;
    l)
      model_list_log="${OPTARG}";;
    q)
      calib_modelList="${OPTARG}";;
    P)
      peak_info="${OPTARG}";;
    c)
      cpu=1;;
    h)
      help;;
    :)
      echo -e "Error: Undefined options were observed.\n"
      help;;
    \?)
      echo -e "Error: Undefined options were specified.\n"
      help;;
  esac
done

# required args
req_arg ${in_f} "i"
req_arg ${cmn_dir} "o"
req_arg ${deepsea} "m"
req_arg ${ref_fa} "f"

# INT check
int_chk ${n_cpu} "-t"
int_chk ${window_n} "-w"

if [ ${cpu} -eq 1 ]; then
  echo -e "Running in silico mutagenesis using only CPU(s)...\n\n"
  bash ./03.run.mutgen.sh -o ${cmn_dir} -m ${deepsea} -M ${ml_dir} -f ${ref_fa} -t ${n_cpu} -l ${model_list_log} -q ${calib_modelList} -P ${peak_info} -c
else
  echo -e "Running in silico mutagenesis...\n\n"
  bash ./03.run.mutgen.sh -o ${cmn_dir} -m ${deepsea} -M ${ml_dir} -f ${ref_fa} -t ${n_cpu} -l ${model_list_log} -q ${calib_modelList} -P ${peak_info}
fi

echo -e "Postprocessing...\n\n"
bash ./04.post.mutgen.sh -o ${cmn_dir}

cd - 1> /dev/null 2>&1
