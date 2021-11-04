#!/bin/bash
#
# 04.post.mutgen.sh
#   Check output files and collect all results in a file
#
#   Usage:
#     04.post.mutgen.sh [-o cmn_dir] ...
#       -o TEXT  output dir [required]
#       -j INT   Number of jobs per 1-GPU. If >1 value is specified, ${j} jobs will use single GPU, simultaneously. **Must be the same with 02.collect.inputs.sh** [default:1]
#       -h       Show help (this message)

# Move to MENTR/bin dir
cd `dirname $0`

# == utils ==
source ../src/sh/utils.sh

# == arg parser ==

# init
job_1gpu=1

# get
while getopts ":o:j:dh" optKey; do
  case "${optKey}" in
    o)
      cmn_dir="${OPTARG}";;
    j)
      job_1gpu="${OPTARG}";;
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

# INT check
int_chk ${job_1gpu} "-j"

# == RUN ==
set -eu

# check outputs
function chk_output(){

    cmn_dir=${1}
    index=${2}

    type=logistic

    for i in `cat ${cmn_dir}/input/multifile_${index}_${type}.txt | cut -f3`
    do
      for j in diff refs alts
      do
        if [ ! -e ${i}_${j}.txt.gz ]; then
            echo -e "Check ${i}"
        elif ! gzip -t ${i}_${j}.txt.gz; then
          echo -e "Check ${i}"
        fi
      done
    done
}
export -f chk_output

if [ `expr ${job_1gpu} - 1` -eq 0 ]; then
  max_j=0
else
  max_j=`expr ${job_1gpu} - 1`
fi

for j in `seq 0 ${max_j}`
do
    chk_output ${cmn_dir} ${j}
done

# Collect all results in a file.
bash ../src/sh/insilico_mutgen_postproess.sh ${cmn_dir}
