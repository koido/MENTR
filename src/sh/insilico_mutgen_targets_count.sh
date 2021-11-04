#!/usr/bin/sh

set -eu

# args
cmn_dir=$1

# init
out_cmn=${cmn_dir}/input/closest_gene

# check
vt_pairs=0
targets=0
for tmp in `ls -1d ${out_cmn}/chr*`
do
    vt_pair=`zcat ${tmp}/*closestgene.gz | wc -l`
    target=`ls ${tmp}/*closestgene.gz | wc -l`
    vt_pairs=`expr ${vt_pairs} + ${vt_pair}`
    targets=`expr ${targets} + ${target}`
done

echo -e "Total number of variant-transcripts pairs was: ${vt_pairs}\n"
echo -e "Total number of target transcripts was: ${targets}\n"
