#!/usr/bin/sh

#-------------------------------------------
#
# Make closest gene files (for unknown corresponding CAGE peak ID mode)
#
# Input arguments:
#   arg1. working dir
#   arg2. bed file (chr#NUM, st_peak, end_peak, strand (+/-/N), CAGE_peak_ID)
#   arg3. input file name (CHR, POS, REF, ALT, ...; tsv file format)
#   arg4. buffer region for ML (I set 100kb + 250 for safety; but later, I use only variants whose distance <=100kb)
#   arg5. does the input file contain header? [T/F]
#
#-------------------------------------------

set -eu

# args
cmn_dir=$1
in_bed=$2
INFILE=$3
buffer=$4
header=$5

# init
out_cmn=${cmn_dir}/input/closest_gene; mkdir -p ${out_cmn}

# judge gzip or not for the input file 
suffix_in=${INFILE##*.}
if [ "${suffix_in}" == "gz" ]; then
    cat_cmd=zcat
else
    cat_cmd=cat
fi

# get all target chromosomes
if [ "${header}" == "T" ]; then
    chrs=(`${cat_cmd} ${INFILE} | tail -n +2 | cut -f1 | sort | uniq`)
elif [ "${header}" == "F" ]; then
    chrs=(`${cat_cmd} ${INFILE} | cut -f1 | sort | uniq`)
else
    echo "Please specify arg6 by using [T/F]"
    exit
fi

# make closest files for each promoter/enhancer ID
function make_inputs(){

    set -u
    clstr_id=${1}
    out_cmn=${2}
    in_bed=${3}
    buffer=${4}
    INFILE=${5}
    chr=${6}
    cat_cmd=${7}
    tmp_INFILE=${8}
    header=${9}

    closest_file=${out_cmn}/chr${chr}/${clstr_id}.closestgene

    # closest file (minimum set)
    cat ${in_bed} | awk -v gen=${clstr_id} -v chr=chr${chr} '{if($5 == gen && $1 == chr){print $0}}' | sed -e "s/^chr//" > ${out_cmn}/chr${chr}/${clstr_id}.bed
    closest-features --delim '\t' \
        --closest \
        --dist \
        ${tmp_INFILE}.chr${chr} ${out_cmn}/chr${chr}/${clstr_id}.bed | \
            awk -F"\t" -v buffer=${buffer} '{if(abs($9)<=buffer){print $0}} function abs(x){return (x>0)? x:-x}' > ${closest_file}

    if [ `cat ${closest_file} | wc -l` -eq 0 ]; then
        rm -f ${closest_file} ${out_cmn}/chr${chr}/${clstr_id}.bed
    else

        # make input file for selected clstr_id
        cat ${closest_file} | awk '{print $3}' | sort -k1,1 -g > ${out_cmn}/chr${chr}/${clstr_id}.tmp.sort
        min_inf=`cat ${out_cmn}/chr${chr}/${clstr_id}.tmp.sort | head -n 1`
        max_inf=`cat ${out_cmn}/chr${chr}/${clstr_id}.tmp.sort | tail -n 1`
        if [ "${header}" == "T" ]; then
            ${cat_cmd} ${INFILE} | \
                tail -n +2 | \
                awk -v chr=${chr} -v min_inf=${min_inf} -v max_inf=${max_inf} '
                {
                    if($1 == chr && $2 >= min_inf && $2 <= max_inf){
                        print $0
                    }
                }' | \
                gzip -c > ${out_cmn}/chr${chr}/${clstr_id}.txt.gz
        elif [ "${header}" == "F" ]; then
            ${cat_cmd} ${INFILE} | \
                awk -v chr=${chr} -v min_inf=${min_inf} -v max_inf=${max_inf} '
                {
                    if($1 == chr && $2 >= min_inf && $2 <= max_inf){
                        print $0
                    }
                }' | \
                gzip -c > ${out_cmn}/chr${chr}/${clstr_id}.txt.gz
        fi
        rm -f ${out_cmn}/chr${chr}/${clstr_id}.bed

        # clstr_ids list for deepsea & expecto
        if [ `zcat ${out_cmn}/chr${chr}/${clstr_id}.txt.gz | wc -l` -gt 0 ]; then
            # gzip, closest_file
            cat ${closest_file} | gzip -c > ${closest_file}.gz
        else
            # clean
            rm -f ${out_cmn}/chr${chr}/${clstr_id}.txt.gz
        fi

        # clean
        rm -f ${closest_file} ${out_cmn}/chr${chr}/${clstr_id}.tmp.sort

    fi
}
export -f make_inputs

# make common temp file
tmp_INFILE=${INFILE}_`date +"%Y%m%d%I%M%S"`
if [ "${header}" == "T" ]; then
    ${cat_cmd} ${INFILE} | tail -n +2 | awk '{printf $1"\t"$2-1"\t"$2"\n"}' | sort-bed - > ${tmp_INFILE}
elif [ "${header}" == "F" ]; then
    ${cat_cmd} ${INFILE} | awk '{printf $1"\t"$2-1"\t"$2"\n"}' | sort-bed - > ${tmp_INFILE}
fi


# RUN
for chr in ${chrs[@]}
do
    if [ "${chr}" != "Y" ]; then
        mkdir -p ${out_cmn}/chr${chr}
        echo -e "Running variants on chr${chr} ..."

        cat ${tmp_INFILE} | awk -F"\t" -v chr=${chr} 'BEGIN{OFS="\t"} {if($1 == chr){print $0}}' | sort-bed - > ${tmp_INFILE}.chr${chr}

        arg1s=`cat ${in_bed} | awk -F"\t" -v chr=chr${chr} '{if($1 == chr){print $5}}' | sort | uniq`
        for arg1 in ${arg1s[@]}
        do
            make_inputs ${arg1} ${out_cmn} ${in_bed} ${buffer} ${INFILE} ${chr} ${cat_cmd} ${tmp_INFILE} ${header}
        done
    fi
    # clean
    rm -f ${tmp_INFILE}.chr${chr}
done

# clean
rm -f ${tmp_INFILE}

echo -e "End.\n"
