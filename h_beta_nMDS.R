# 2D ordination to assess beta diversity

# subset anemone samples
anemones <- subset_samples(phy, sample_type == "anemone")
# 530 ASVs

anemones <- prune_taxa((taxa_sums(anemones) > 0), anemones)
# 434 ASVs

# convert samplingDay varaible from numeric to factor
anemones@sam_data$samplingDay <- factor(anemones@sam_data$samplingDay)

anemones.nmds <- ordinate(anemones, method = "NMDS", k = 2, distance = "bray")
# stress 0.1129381

p <- plot_ordination(
  anemones,
  ordination = anemones.nmds,
  type = "samples",
  axes = c(1, 2),
  color = "samplingDay",
  shape = "genotype"
  ) +
  theme_bw() +
  geom_point(size = 8) +
  stat_ellipse()

print(p)

rm(anemones)
rm(anemones.nmds)
rm(p)
