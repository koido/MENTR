#!/bin/bash
#
# 01.get.pairs.sh
#   Preprocessing: get variant-promoter/enhancer pairs 
#
#   Usage:
#     01.get.pairs.sh [-i input_file] [-o cmn_dir] ...
#       -i TEXT  input variant file [required]
#       -o TEXT  output dir [required]
#       -w INT   window size [bp] [default: 100000]
#       -b TEXT  input bed file to find variant-CAGEID pairs [default: ALL FANTOM5 CAGE IDs ../resources/mutgen_cmn/F5.cage_cluster.hg19.win_100kb.sort.bed]
#                 col1: chr
#                 col2: start (TSS - 1) or (midpoint - 1)
#                 col3: end (TSS or midpoint)
#                 col4: strand [+/-/N] N is for enhancer
#                 col5: CAGE_peak_ID
#       -t INT   Number of threads. [default:1]
#       -d       Dry Run mode
#       -h       Show help (this message)

# Move to MENTR/bin dir
cd `dirname $0`

# == utils ==
source ../src/sh/utils.sh

# == arg parser ==

# Default parameters
DryRun=0
# window size [bp]; default is (Promoter TSS / Enhancer midpoint) +/- 100-kb)
window_n=100000
# Number of CPU
n_cpu=1
# F5 promoter TSS and enhancer mid position (maximum set for ML models (hg19)
in_bed=../resources/mutgen_cmn/F5.cage_cluster.hg19.win_100kb.sort.bed

# get
while getopts ":i:o:w:b:t:dh" optKey; do
  case "${optKey}" in
    i)
      in_f="${OPTARG}";;
    o)
      cmn_dir="${OPTARG}";;
    w)
      window_n="${OPTARG}";;
    b)
      in_bed="${OPTARG}";;
    t)
      n_cpu="${OPTARG}";;
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
req_arg ${in_f} "i"
req_arg ${cmn_dir} "o"

# INT check
int_chk ${window_n} "-w"
int_chk ${n_cpu} "-t"

# == common args ==

# whether input file contains header or not [T/F]
header=F

# Show args
echo -e "# == Args for 01.get.pairs.sh =="
echo -e "# in_f: ${in_f}"
echo -e "# cmn_dir: ${cmn_dir}"
echo -e "# window_n: ${window_n}"
echo -e "# in_bed: ${in_bed}"
echo -e "# DryRun: ${DryRun}"
echo -e "# n_cpu: ${n_cpu}"
#echo -e "# header: ${header}"
echo -e ""

if [ ${n_cpu} -gt 1 ]; then
  parallel=T

  # check availability of GNU parallel
  if ! type parallel > /dev/null 2>&1; then
    echo -e "Error: GNU parallel is not available."
    exit 1
  fi

else
  parallel=F
fi

# == Dry RUN ==
if [ ${DryRun} -eq 1 ]; then
  echo -e "# ** Dry Run **"
  echo -e "bash ../src/sh/insilico_mutgen_make_closest_gene_prospective.sh ${cmn_dir} ${in_bed} ${in_f} ${window_n} ${header} ${parallel} ${n_cpu}"
  exit 1
fi

# == RUN ==
bash ../src/sh/insilico_mutgen_make_closest_gene_prospective.sh \
  ${cmn_dir} \
  ${in_bed} \
  ${in_f} \
  ${window_n} \
  ${header} \
  ${parallel} ${n_cpu} && \
  bash ../src/sh/insilico_mutgen_targets_count.sh ${cmn_dir}
