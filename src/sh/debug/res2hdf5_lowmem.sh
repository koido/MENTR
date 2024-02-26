#!/bin/bash

# Preparation of toy examples for debugging.
rm -Rf debug
mkdir -p debug/{input,output}
mkdir -p debug/output/{standard,lowmem}

# chromatin effect (100000 rows 100 columns)
## Use numpy
## Output: debug/input/toy.chrom.txt.gz
python -c "import numpy as np; np.random.seed(0); np.savetxt('debug/input/toy.chrom.txt.gz', np.random.rand(100000, 100), delimiter='\t')"
## Very small toy example with missing values
#echo -e "0.1\t0.2\t0.3\t0.4\t0.5\t0.6\t0.7\t0.8\t0.9\t0.10" > debug/input/toy.chrom.txt
#for i in {1..9989}; do echo -e "0.1\t0.2\t0.3\t0.4\t0.5\t0.6\t0.7\t0.8\t0.9\t0.10" >> debug/input/toy.chrom.txt; done
#for i in {1..10}; do echo -e "-9999999999.0\t-9999999999.0\t-9999999999.0\t-9999999999.0\t-9999999999.0\t-9999999999.0\t-9999999999.0\t-9999999999.0\t-9999999999.0\t-9999999999.0" >> debug/input/toy.chrom.txt; done
#for i in {1..9989}; do echo -e "0.1\t0.2\t0.3\t0.4\t0.5\t0.6\t0.7\t0.8\t0.9\t0.10" >> debug/input/toy.chrom.txt; done
#for i in {1..10}; do echo -e "-9999999999.0\t-9999999999.0\t-9999999999.0\t-9999999999.0\t-9999999999.0\t-9999999999.0\t-9999999999.0\t-9999999999.0\t-9999999999.0\t-9999999999.0" >> debug/input/toy.chrom.txt; done
#for i in {1..80001}; do echo -e "0.1\t0.2\t0.3\t0.4\t0.5\t0.6\t0.7\t0.8\t0.9\t0.10" >> debug/input/toy.chrom.txt; done
#gzip debug/input/toy.chrom.txt

# peak file
#echo -e "clusterID peak" > debug/input/toy.peak.txt
touch debug/input/toy.peak.txt
for i in {1..100000}; do echo -e "cluster${i}\t$i" >> debug/input/toy.peak.txt; done
gzip debug/input/toy.peak.txt

# expression file
echo -e "clusterID\texpr1\texpr2\texpr3\texpr4\texpr5" > debug/input/toy.exp.txt
for i in {1..100000}; do echo -e "cluster${i}\t1\t2\t3\t4\t5" >> debug/input/toy.exp.txt; done
gzip debug/input/toy.exp.txt

# cluster file
echo -e "clusterID\tclusterName\ttype\tmask\tF5_tag_count\tgeneNum\ttrnscptIDStr\tgeneIDStr\tgeneNameStr\tgeneClassStr\tF5_anno" > debug/input/toy.cluster.txt
for i in {1..100000}; do echo -e "cluster${i}\tcluster${i}\ttype\tmask\t1000\t1000\ttrnscptIDStr\tgeneIDStr\tgeneNameStr\tgeneClassStr\tF5_anno" >> debug/input/toy.cluster.txt; done
gzip debug/input/toy.cluster.txt

# Run python script with measuring time for each job
/usr/bin/time -v python src/py3/seq2chrom_res2hdf5.py --chromFile debug/input/toy.chrom.txt.gz --peakFile debug/input/toy.peak.txt.gz --expFile debug/input/toy.exp.txt.gz --clusterFile debug/input/toy.cluster.txt.gz --output debug/output/standard/ 1> debug/output/standard/log.txt 2> debug/output/standard/err.txt
/usr/bin/time -v python src/py3/seq2chrom_res2hdf5.py --chromFile debug/input/toy.chrom.txt.gz --peakFile debug/input/toy.peak.txt.gz --expFile debug/input/toy.exp.txt.gz --clusterFile debug/input/toy.cluster.txt.gz --output debug/output/lowmem/ --lowmem --chunksize 100 1> debug/output/lowmem/log.txt 2> debug/output/lowmem/err.txt

# Check the consistency of the output files
for i in info_train.txt.gz info_test.txt.gz info_other.txt.gz X_train_col.txt.gz X_test_col.txt.gz X_other_col.txt.gz Y_train_col.txt.gz Y_test_col.txt.gz Y_other_col.txt.gz common_train_rows.txt.gz common_test_rows.txt.gz common_other_rows.txt.gz metrics.txt
do
  diff_rows=`zdiff debug/output/standard/$i debug/output/lowmem/$i | wc -l`
  if [ $diff_rows -ne 0 ]; then
    echo "Error: $i"
  fi
done

# Check X_Y.h5
echo "=== h5diff ==="
h5diff -c debug/output/standard/X_Y.h5 debug/output/lowmem/X_Y.h5

# Show time and memory usage
echo "=== time (s) ==="
echo "standard: `cat debug/output/standard/err.txt | grep "User time" | cut -d":" -f2`"
echo "lowmem: `cat debug/output/lowmem/err.txt | grep "User time" | cut -d":" -f2`"
echo "=== memory (MB) ==="
echo "standard: `cat debug/output/standard/err.txt | grep "Maximum resident set size" | cut -d":" -f2 | awk '{print int($1/1024)}'`"
echo "lowmem: `cat debug/output/lowmem/err.txt | grep "Maximum resident set size" | cut -d":" -f2 | awk '{print int($1/1024)}'`"
# rm -Rf debug