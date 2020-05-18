#!/bin/sh


rm -r /mnt/tmp
mkdir /mnt/tmp
export TMPDIR=/mnt/tmp


# # # # # # # # # # # # # # # # # # # # # # # # #


mkdir ~/seqs

qiime tools import \
--type SampleData[PairedEndSequencesWithQuality] \
--input-path ~/manifest.txt \
--input-format PairedEndFastqManifestPhred33V2 \
--output-path ~/seqs/combined.qza


# # # # # # # # # # # # # # # # # # # # # # # # #


qiime cutadapt trim-paired \
--i-demultiplexed-sequences ~/seqs/combined.qza \
--p-front-f AGGATTAGATACCCTGGTA \
--p-front-r CRRCACGAGCTGACGAC \
--p-error-rate 0.20 \
--output-dir ~/trim \
--verbose


# # # # # # # # # # # # # # # # # # # # # # # # #


mkdir ~/reports

qiime demux summarize \
--i-data  ~/trim/trimmed_sequences.qza \
--o-visualization  ~/reports/trimmed_seqs.qzv


# # # # # # # # # # # # # # # # # # # # # # # # #


qiime dada2 denoise-paired \
--i-demultiplexed-seqs ~/trim/trimmed_sequences.qza \
--p-trunc-len-f 210 \
--p-trunc-len-r 166 \
--p-n-threads 0 \
--output-dir ~/dada2out \
--verbose


# # # # # # # # # # # # # # # # # # # # # # # # #


qiime feature-table summarize \
--i-table ~/dada2out/table.qza \
--m-sample-metadata-file ~/metadata.txt \
--o-visualization ~/reports/table.qzv \
--verbose

qiime feature-table tabulate-seqs \
--i-data ~/dada2out/representative_sequences.qza \
--o-visualization ~/reports/rep_seqs.qzv \
--verbose

qiime metadata tabulate \
--m-input-file ~/dada2out/denoising_stats.qza \
--o-visualization ~/reports/denoising_stats.qzv \
--verbose


# # # # # # # # # # # # # # # # # # # # # # # # #


qiime feature-classifier classify-sklearn \
--i-classifier ~/silva132_classifier.qza \
--i-reads ~/dada2out/representative_sequences.qza \
--p-n-jobs -8 \
--output-dir ~/taxonomy \
--verbose


# # # # # # # # # # # # # # # # # # # # # # # # #


qiime metadata tabulate \
--m-input-file ~/taxonomy/classification.qza \
--o-visualization ~/reports/taxonomy.qzv \
--verbose


# # # # # # # # # # # # # # # # # # # # # # # # #


qiime taxa filter-seqs \
--i-sequences ~/dada2out/representative_sequences.qza \
--i-taxonomy ~/taxonomy/classification.qza \
--p-exclude fbf10cb48bb23a060d4b0168b8fa5a9c,0ae23b19ed2134869ab1dbfdb6270834 \
--p-mode 'exact' \
--o-filtered-sequences ~/dada2out/rep_seqs_filtered.qza \
--verbose


# # # # # # # # # # # # # # # # # # # # # # # # #


qiime taxa filter-table \
--i-table ~/dada2out/table.qza \
--i-taxonomy ~/taxonomy/classification.qza  \
--p-exclude Mitochondria,Chloroplast \
--p-include D_1__ \
--o-filtered-table ~/dada2out/table_filtered.qza \
--verbose


# # # # # # # # # # # # # # # # # # # # # # # # #


qiime feature-table summarize \
--i-table ~/dada2out/table_filtered.qza \
--m-sample-metadata-file ~/metadata.txt \
--o-visualization ~/reports/table_filtered.qzv \
--verbose


# # # # # # # # # # # # # # # # # # # # # # # # #


mkdir ~/aligned

qiime alignment mafft \
--i-sequences ~/dada2out/rep_seqs_filtered.qza \
--p-n-threads 15 \
--o-alignment ~/aligned/aligned_rep_seqs.qza \
--verbose


# # # # # # # # # # # # # # # # # # # # # # # # #


qiime alignment mask \
--i-alignment ~/aligned/aligned_rep_seqs.qza \
--o-masked-alignment ~/aligned/masked_aligned_rep_seqs.qza \
--verbose


# # # # # # # # # # # # # # # # # # # # # # # # #


mkdir ~/tree

qiime phylogeny fasttree \
--i-alignment ~/aligned/masked_aligned_rep_seqs.qza \
--p-n-threads 1 \
--o-tree ~/tree/unrooted_tree.qza \
--verbose


# # # # # # # # # # # # # # # # # # # # # # # # #


qiime phylogeny midpoint-root \
--i-tree ~/tree/unrooted_tree.qza \
--o-rooted-tree ~/tree/rooted_tree.qza \
--verbose


# # # # # # # # # # # # # # # # # # # # # # # # #


qiime taxa barplot \
--i-table ~/dada2out/table_filtered.qza \
--i-taxonomy ~/taxonomy/classification.qza \
--m-metadata-file ~/metadata.txt \
--o-visualization ~/reports/barchart.qzv \
--verbose


# # # # # # # # # # # # # # # # # # # # # # # # #


mkdir ~/emperor/

qiime diversity beta \
--i-table ~/dada2out/table_filtered.qza \
--p-metric 'braycurtis' \
--p-n-jobs -8 \
--o-distance-matrix ~/emperor/bray.qza \
--verbose

qiime diversity pcoa \
--i-distance-matrix ~/emperor/bray.qza \
--o-pcoa ~/emperor/pcoa.qza \
--verbose

qiime emperor plot \
--i-pcoa ~/emperor/pcoa.qza \
--m-metadata-file ~/metadata.txt \
--o-visualization ~/reports/emperor.qzv \
--verbose


# # # # # # # # # # # # # # # # # # # # # # # # #


mkdir ~/output

qiime tools export \
--input-path ~/tree/unrooted_tree.qza \
--output-path ~/output/

qiime tools export \
--input-path ~/taxonomy/classification.qza \
--output-path ~/output/

qiime tools export \
--input-path ~/dada2out/table_filtered.qza \
--output-path ~/output/

biom convert  -i ~/output/feature-table.biom \
--to-tsv \
-o ~/output/table.tsv


# # # # # # # # # # # # # # # # # # # # # # # # #


# Use UNIX commands to tidy up the output files

# Remove biom file
rm ~/output/feature-table.biom

# Remove header from otu table
sed -i '1d' ~/output/table.tsv

# Remove header from taxonomy table
sed -i '1,2d' ~/output/taxonomy.tsv

# Remove confidence values from taxonomy table
cut -f 1,2 ~/output/taxonomy.tsv > ~/output/tax.tsv

# Remove old taxonomy file
rm ~/output/taxonomy.tsv

# Replace semi-colons in taxonomy table with tabs
sed -i 's/;/\t/g' ~/output/tax.tsv

# Fix weird name
sed -i 's/D_6__Corynebacterium doosanense CAU 212 = DSM 45436/D_6__Corynebacterium doosanense/g' ~/output/tax.tsv

# # # # # # # # # # # # # # # # # # # # # # # # #
