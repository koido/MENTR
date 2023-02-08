#!/bin/bash
#
# quick.mutgen.sh
#   Run in silico mutagenesis with default parameters
#
#   Usage:
#     quick.mutgen.sh [-i input_file] [-o cmn_dir] ...
#       -i TEXT  input file [required]
#       -o TEXT  output dir [required]
#       -m TEXT  File of DeepSEA Beluga model (deepsea.beluga.2002.cpu) [required]
#       -f TEXT  File of hg19.fa [required]
#       -t INT   Number of threads. [default:1]
#       -w INT   window size to find variant-CAGE_ID pairs [bp] [default:100000]
#       -b TEXT  input bed file to find variant-CAGE_ID pairs [default: ALL FANTOM5 CAGE IDs]
#       -p TEXT  CAGE_ID list to find variant-CAGE_ID pairs. If many, please write IDs separated by semicolon [default: ALL FANTOM5 CAGE IDs]
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

# get
while getopts ":i:o:m:f:t:w:b:p:ch" optKey; do
  case "${optKey}" in
    i)
      in_f="${OPTARG}";;
    o)
      cmn_dir="${OPTARG}";;
    m)
      deepsea="${OPTARG}";;
    f)
      hg19_fa="${OPTARG}";;
    t)
      n_cpu="${OPTARG}";;
    w)
      window_n="${OPTARG}";;
    b)
      in_bed="${OPTARG}";;
    p)
      cage_list="${cage_list}";;
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
if [ -z "${in_f}" ]; then
  echo "Error: -i is required"
  help
fi
if [ -z "${cmn_dir}" ]; then
  echo "Error: -o is required"
  help
fi
if [ -z "${deepsea}" ]; then
  echo "Error: -m is required"
  help
fi
if [ -z "${hg19_fa}" ]; then
  echo "Error: -f is required"
  help
fi

# INT check
int_chk ${n_cpu} "-t"
int_chk ${window_n} "-w"

if [ ${cage_list} != "ALL" ];
  join -t$'\t' -1 5 \
    <(sort -t$'\t' -k5,5 ${in_bed}) \
    <(echo ${cage_list} | sed -e "s/;/\n/g" | sort) | \
    awk -F"\t" 'BEGIN{OFS="\t"}{print $2,$3,$4,$5,$1}' > ${cmn_dir}/custom.F5.cage_cluster.hg19.win_100kb.sort.bed
  in_bed=${cmn_dir}/custom.F5.cage_cluster.hg19.win_100kb.sort.bed
  echo -e "\nUse custom bed file (${in_bed}) for finding variant-promoter/enhancer pairs.\n"
fi

echo -e "Preprocessing for getting variant-promoter/enhancer pairs...\n\n"
bash ./01.get.pairs.sh -i ${in_f} -o ${cmn_dir} -w ${window_n} -b ${in_bed}

echo -e "Preprocessing for collecting information of input files...\n\n"
bash ./02.collect.inputs.sh -o ${cmn_dir}

if [ ${cpu} -eq 1 ]; then
  echo -e "Running in silico mutagenesis using only CPUs...\n\n"
  bash ./03.run.mutgen.sh -o ${cmn_dir} -m ${deepsea} -f ${hg19_fa} -t ${n_cpu} -c
else
  echo -e "Running in silico mutagenesis...\n\n"
  bash ./03.run.mutgen.sh -o ${cmn_dir} -m ${deepsea} -f ${hg19_fa} -t ${n_cpu}
fi

echo -e "Postprocessing...\n\n"
bash ./04.post.mutgen.sh -o ${cmn_dir}

cd - 1> /dev/null 2>&1
