#!/bin/bash
#
# quick.mutgen.sh
#   Run in silico mutagenesis with default parameters
#
#   Usage:
#     quick.mutgen.sh [-i input_file] [-o cmn_dir] ...
#       -i TEXT  input file [required]
#       -o TEXT  output dir [required]
#       -m TEXT  File of DeepSEA Beluga model (deepsea.beluga.2002.cpu)
#       -f TEXT  File of hg19.fa
#       -h       Show help (this message)

# Move to MENTR/bin dir
cd `dirname $0`

# == utils ==
source ../src/sh/utils.sh

# == arg parser ==

# init

# get
while getopts ":i:o:m:f:h" optKey; do
  case "${optKey}" in
    i)
      in_f="${OPTARG}";;
    o)
      cmn_dir="${OPTARG}";;
    m)
      deepsea="${OPTARG}";;
    f)
      hg19_fa="${OPTARG}";;
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

echo -e "Preprocessing for getting variant-promoter/enhancer pairs...\n\n"
bash ./01.get.pairs.sh -i ${in_f} -o ${cmn_dir}

echo -e "Preprocessing for collecting information of input files...\n\n"
bash ./02.collect.inputs.sh -o ${cmn_dir}

echo -e "Running in silico mutagenesis...\n\n"
bash ./03.run.mutgen.sh -o ${cmn_dir} -m ${deepsea} -f ${hg19_fa}

echo -e "Postprocessing...\n\n"
bash ./04.post.mutgen.sh -o ${cmn_dir}

cd - 1> /dev/null 2>&1
