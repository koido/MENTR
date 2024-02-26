#!/usr/bin/bash
#SBATCH --job-name mutgen
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1

set -u

out_dir_cmn=$1
ml_dir=$2
ref_fa=$3
deepsea=$4
# Number of jobs per 1 GPU
job_1gpu=$5
# number of thread (CPU)
n_core=$6
# mutgen parameters
model_list_log=$7
calib_modelList=$8
# peak file
peak_info=$9
# Main Python script
py_script=${10}

# Header name of input file which you prepared; comma separated
in_cols="CHR,POS,REF,ALT"

# ML type [logistic/regression] (Github's version only support logistic)
type=logistic

# function for run probability models
function run_mutgen(){

    set -u

    cmn_dir=${1}
    type=${2}
    py_script=${3}
    n_core=${4}
    model_list_log=${5}
    index=${6}
    ref_fa=${7}
    peak_info=${8}
    deepsea=${9}
    in_cols=${10}
    ml_dir=${11}

    log_cmn=${cmn_dir}/output/${type}_${index}

    python -u ${py_script} \
      --ref_fa ${ref_fa} \
      --peakinfo ${peak_info} \
      --deepsea ${deepsea} \
      --threads ${n_core} \
      --modelList ${model_list_log} \
      --multifiles \
      --input_cols `echo ${in_cols} | sed -e "s/,/ /g"` \
      --inFiles ${cmn_dir}/input/multifile_${index}_${type}.txt \
      --calib_modelList ${calib_modelList} \
      --modelList_prefix ${ml_dir} 1> ${log_cmn}.std.log 2> ${log_cmn}.err.log

}
export -f run_mutgen

# init
export CUDA_DEVICE_ORDER=PCI_BUS_ID
max_j=`expr ${job_1gpu} - 1`

# run
for j in `seq 0 ${max_j}`
do
    run_mutgen ${out_dir_cmn} ${type} ${py_script} ${n_core} ${model_list_log} ${j} ${ref_fa} ${peak_info} ${deepsea} "${in_cols}" ${ml_dir} &
done

wait;
