#!/bin/bash
#
# 01.get.pairs.sh
#   Preprocessing: get variant-promoter/enhancer pairs 
#
#   Usage:
#     01.get.pairs.sh [-i input_file] [-o cmn_dir] ...
#       -i TEXT  input file [required]
#       -o TEXT  output dir [required]
#       -w INT   window size [bp] [default:100000]
#       -d       Dry Run mode
#       -h       Show help (this message)

# Move to MENTR/bin dir
cd `dirname $0`

# == utils ==
source ../src/sh/utils.sh

# == arg parser ==

# init
DryRun=0
# window size [bp]; default is (Promoter TSS / Enhancer midpoint) +/- 100-kb) but you can change this parameter
window_n=100000

# get
while getopts ":i:o:w:dh" optKey; do
  case "${optKey}" in
    i)
      in_f="${OPTARG}";;
    o)
      cmn_dir="${OPTARG}";;
    w)
      window_n="${OPTARG}";;
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
if [ -z "${in_f}" ]; then
  echo "Error: -i is required"
  help `basename $0`
fi
if [ -z "${cmn_dir}" ]; then
  echo "Error: -o is required"
  help `basename $0`
fi

# INT check
int_chk ${window_n} "-w"

# == common args ==

# F5 promoter TSS and enhancer mid position (maximum set for ML models; that is, this file do NO include transcripts in edge region o
in_bed=../resources/mutgen_cmn/F5.cage_cluster.hg19.win_100kb.sort.bed

# whether input file contains header or not [T/F]
header=F

# == Dry RUN ==
if [ ${DryRun} -eq 1 ]; then
  echo -e "# ** Dry Run **"
  echo -e "# in_f: ${in_f}"
  echo -e "# cmn_dir: ${cmn_dir}"
  echo -e "bash ../src/sh/insilico_mutgen_make_closest_gene_prospective.sh ${cmn_dir} ${in_bed} ${in_f} ${window_n} ${header}"
  exit 1
fi

# == RUN ==
bash ../src/sh/insilico_mutgen_make_closest_gene_prospective.sh \
  ${cmn_dir} \
  ${in_bed} \
  ${in_f} \
  ${window_n} \
  ${header} && \
  bash ../src/sh/insilico_mutgen_targets_count.sh ${cmn_dir}
