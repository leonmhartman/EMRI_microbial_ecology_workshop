# load R packages used in the EMRI workshop

library(ape)
library(phyloseq)
library(vegan)
library(decontam)
library(ggplot2)
library(microbiome)
library(gridExtra)
library(car)
library(mvabund)

#--------------------------------

# import data

emri_otu <- read.table(
  "table.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1)

emri_tax <- read.table(
  "tax.tsv",
  sep = "\t",
  fill = TRUE,
  row.names = 1)

emri_met <- read.table(
  "meta.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1)

emri_tre <- read_tree("tree.nwk")

#--------------------------------

# convert to phyloseq

# add levels of taxonomy to taxonomy table
colnames(emri_tax) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# convert data frames to matrices for compatibility with phyloseq
emri_otu_mat <- as.matrix(emri_otu)
emri_tax_mat <- as.matrix(emri_tax)

# combine OTU, taxonomy, and metadata files into a phyloseq object
phy <- phyloseq(otu_table(emri_otu_mat, taxa_are_rows = T),
                tax_table(emri_tax_mat),
                sample_data(emri_met))

# add tree data
phy <- merge_phyloseq(phy, emri_tre)

# remove zero sum ASVs
phy <- prune_taxa((taxa_sums(phy) > 0), phy)
# 542 ASVs

rm(emri_otu_mat)
rm(emri_tax_mat)
