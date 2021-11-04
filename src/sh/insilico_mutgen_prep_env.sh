#!/usr/bin/sh

#-------------------------------------------
#
# Prepare in silico mutagenesis
#
# Input arguments:
#   arg1. working dir
#   arg2. how many GPUs do you have? -> 0 - ($3 - 1) of CUDA_VISIBLE_DEVICES are assigned.
#   arg3. Number of jobs per 1 GPU
#
#-------------------------------------------

set -u

# args
cmn_dir=$1
max_gpu=$2
job_1gpu=$3

# depreciated variable
type=logistic

# cmn function for cleaning before first running
function init_run(){

    set -u

    type=${1}
    cmn_dir=${2}
    index=${3}

    log_cmn=${cmn_dir}/output/${type}_${index}
    job_name=${type}_${index}
    rm -f ${log_cmn}.std.log
    rm -f ${log_cmn}.err.log
    rm -f ${cmn_dir}/input/multifile_${index}_${type}.txt
}
export -f init_run

# preparation of input files
function prep_run(){

    set -u
    
    clstr_id=${1}
    cmn_dir=${2}
    index=${3}
    type=${4}

    chr=`echo ${clstr_id} | awk -F":" '{print $1}' | sed -e "s/^chr//g"`

    closest_gene_dir=${cmn_dir}/input/closest_gene/chr${chr}
    out_dir=${cmn_dir}/output/chr${chr}
    mkdir -p ${out_dir}
    mkdir -p ${out_dir}/${type}

    # (for Re-RUN)
    if [ -e ${closest_gene_dir}/${clstr_id}.closestgene.gz ]; then
        wc_closest=`zcat ${closest_gene_dir}/${clstr_id}.closestgene.gz | head -n1 | wc -w`
    else
        wc_closest=99
    fi

    if [ ! -e ${out_dir}/${type}/${clstr_id}_diff.txt.gz ] && [ ${wc_closest} -ne 5 ]; then
        if [ -e ${closest_gene_dir}/${clstr_id}.txt.gz ]; then
            n_records=`zcat ${closest_gene_dir}/${clstr_id}.txt.gz | wc -l`
            echo -e "${closest_gene_dir}/${clstr_id}.txt.gz\t${closest_gene_dir}/${clstr_id}.closestgene.gz\t${out_dir}/${type}/${clstr_id}" >> ${cmn_dir}/input/multifile_${index}_${type}.txt
        fi
    fi
}
export -f prep_run

max_index=`expr ${max_gpu} - 1`
max_ij=`expr ${max_gpu} \* ${job_1gpu}`

# init
ij=0
for index in `seq 0 ${max_index}`
do
    for j in `seq 1 ${job_1gpu}`
    do
        rm -f ${cmn_dir}/output/multifile_${ij}_${type}.txt
        ((ij++))
    done
done

# make lists of multiple input files
ij=0
for index in `seq 0 ${max_index}`
do
    for j in `seq 1 ${job_1gpu}`
    do
        arg1s=`ls -1d ${cmn_dir}/input/closest_gene/chr*/*.txt.gz | \
            xargs -I % bash -c "basename %" | \
            sed -e "s/.txt.gz$//" | \
            sort | uniq | \
            awk -v ind=${ij} -v n_para=${max_ij} '{if((NR % n_para)==ind){print}}'`
        for arg1 in ${arg1s[@]}
        do
            prep_run ${arg1} ${cmn_dir} ${ij} ${type}
        done
        ((ij++))
    done
done
