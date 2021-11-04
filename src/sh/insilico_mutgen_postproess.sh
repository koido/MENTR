#!/bin/bash

set -eu

cmn_dir=$1

cd ${cmn_dir}/output
out_f=./all_qcd_res.txt

# collect all results
for chr in `seq 1 22` X
do
    if [ -e ./chr${chr} ]; then
        for i in `ls ./chr${chr}/logistic/*_diff.txt.gz`
        do
            if [ ! -e ${out_f} ]; then
                zcat ${i} | cut -f2- > ${out_f}
            else
                zcat ${i} | tail -n +2 | cut -f2- >> ${out_f}
            fi
        done
    fi
done

# exclude refmatch False records
#  -> If exist, the records' predictions are not exact due to too complex conditions (e.g., InDel near TSS)
cat ${out_f} | \
    awk -F"\t" '
        BEGIN{
            OFS="\t"
        }
        {
            if(NR==1){
                print
            }else if($8 == "True"){
                print
            }
        }' | gzip -c > ${out_f}.gz
rm -f ${out_f}

cd -
