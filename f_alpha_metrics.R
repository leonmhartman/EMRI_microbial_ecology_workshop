# diversity metrics

# subset anemone samples
anemones <- subset_samples(phy, sample_type == "anemone")
# 530 ASVs

anemones <- prune_taxa((taxa_sums(anemones) > 0), anemones)
# 434 ASVs

# rarefy at a level equal to the sample with the fewest reads 
anemones <- rarefy_even_depth(anemones, min(sample_sums(anemones)), rngseed = 1)
# 439 ASVs

# - - - - - - - - - - - - - - - -

# create diversity metric plot object of "average observed ASVs" per genotype/day
ob <- plot_richness(anemones, "grouping", measures = c("Observed"))

# specify box-plot
ob <- ob + geom_boxplot(data=ob$data, aes(x=grouping, y=value), show.legend = FALSE) +
  theme_bw() +
  theme(axis.title = element_blank())

# reorder data along x-axis
ob$data$grouping <- factor(ob$data$grouping, levels = c("a1_1", "a4_1", "a1_2", "a4_2"))

# - - - - - - - - - - - - - - - -

# create diversity metric plot object of "Simpson Index" evenness per genotype/day
sm <- plot_richness(anemones, "grouping", measures = c("Simpson"))

# specify box-plot
sm <- sm + geom_boxplot(data=sm$data, aes(x=grouping, y=value), show.legend = FALSE) +
  theme_bw() +
  theme(axis.title = element_blank())

# reorder data along x-axis
sm$data$grouping <- factor(sm$data$grouping, levels = c("a1_1", "a4_1", "a1_2", "a4_2"))

# - - - - - - - - - - - - - - - -

# create diversity metric plot object of "Shannon Index" alpha diversity per genotype/day
sh <- plot_richness(anemones, "grouping", measures = c("Shannon"))

# specify box-plot
sh <- sh + geom_boxplot(data=sh$data, aes(x=grouping, y=value), show.legend = FALSE) +
  theme_bw() +
  theme(axis.title = element_blank())

# reorder data along x-axis
sh$data$grouping <- factor(sh$data$grouping, levels = c("a1_1", "a4_1", "a1_2", "a4_2"))

# - - - - - - - - - - - - - - - -

# print plots on one page
grid.arrange(ob, sm, sh, nrow = 1)

rm(ob)
rm(sm)
rm(sh)
rm(anemones)
