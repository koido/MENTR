# Preparation of toy examples for debugging.
mkdir -p debug/input

# chromatin effect
echo "0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.10" > debug/input/toy.chrom.txt
for i in {1..9999}; do echo "0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.10" >> debug/input/toy.chrom.txt; done
gzip debug/input/toy.chrom.txt

# peak file
echo "clusterID peak" > debug/input/toy.peak.txt
for i in {1..10000}; do echo "cluster${i} $i" >> debug/input/toy.peak.txt; done
gzip debug/input/toy.peak.txt

# expression file
echo "clusterID expr1 expr2 expr3 expr4 expr5" > debug/input/toy.exp.txt
for i in {1..10000}; do echo "cluster${i} 1 2 3 4 5" >> debug/input/toy.exp.txt; done
gzip debug/input/toy.exp.txt

# cluster file
echo "clusterID clusterName type mask F5_tag_count geneNum trnscptIDStr geneIDStr geneNameStr geneClassStr F5_anno" > debug/input/toy.cluster.txt
for i in {1..10000}; do echo "cluster${i} cluster${i} type mask 1000 1000 trnscptIDStr geneIDStr geneNameStr geneClassStr F5_anno" >> debug/input/toy.cluster.txt; done
gzip debug/input/toy.cluster.txt

# Run python script with measuring time for each job
time python src/py3/seq2chrom_res2hdf5.py --chromFile debug/input/toy.chrom.txt.gz --peakFile debug/input/toy.peak.txt.gz --expFile debug/input/toy.exp.txt.gz --clusterFile debug/input/toy.cluster.txt.gz --output debug/output/standard 1> debug/output/log.txt 2> debug/output/err.txt
time python src/py3/seq2chrom_res2hdf5.py --chromFile debug/input/toy.chrom.txt.gz --peakFile debug/input/toy.peak.txt.gz --expFile debug/input/toy.exp.txt.gz --clusterFile debug/input/toy.cluster.txt.gz --output debug/output/lowmem --lowmem 1> debug/output/log.lowmem.txt 2> debug/output/err.lowmem.txt

# Check the consistency of the output files
for i in X_Y.h5 info_train.txt.gz info_test.txt.gz info_other.txt.gz X_train_col.txt.gz X_test_col.txt.gz X_other_col.txt.gz Y_train_col.txt.gz Y_test_col.txt.gz Y_other_col.txt.gz common_train_rows.txt.gz common_test_rows.txt.gz common_other_rows.txt.gz metrics.txt
do
  diff_rows=`diff debug/output/standard/$i debug/output/lowmem/$i | wc -l`
  if [ $diff_rows -ne 0 ]; then
    echo "Error: $i"
  fi
done
# rm -Rf debug