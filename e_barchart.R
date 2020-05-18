# Subset anemone samples
anemones <- subset_samples(phy, sample_type == "anemone")

# Transform counts to percentages
anemonesPer <- transform_sample_counts(anemones, function(x) 100 * x/sum(x))

# Collapse taxa to specified level of taxonomy
anemonesClass <- tax_glom(anemonesPer, "Class", NArm = FALSE)

# Melt phyloseq data
q <- psmelt(anemonesClass)

# change stacking so most abundant (overall) are at top
q$Class <- reorder(q$Class, q$Abundance)
q$Class <- factor(q$Class, levels = rev(levels(q$Class)))

# change order of samples along x axis
q$Sample <- factor(q$Sample, levels = c("a1_1_01","a1_1_02","a1_1_03","a1_1_04","a1_1_05",
                                        "a4_1_01","a4_1_02","a4_1_03","a4_1_04","a4_1_05",
                                        "a1_2_01","a1_2_02","a1_2_03","a1_2_04","a1_2_05",
                                        "a4_2_01","a4_2_02","a4_2_03","a4_2_04","a4_2_05"))

# Generate ggplot2 object
p <- ggplot(q, aes_string(x = "Sample", y = "Abundance", fill = "Class"))
# Customise and output ggplot2 object
p +  
  geom_bar(stat = 'identity', position = "stack") +
  ylab("Reads assigned to Class") +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = barColours3) +
  guides(fill = guide_legend(ncol = 4)) +
  theme(axis.text.x = element_text(angle = 90))

rm(anemones)
rm(anemonesClass)
rm(anemonesPer)
rm(q)
rm(p)
rm(barColours1)
rm(barColours2)
rm(barColours3)
