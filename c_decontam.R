# identify the contaminants
consList <- isContaminant(seqtab = phy, neg = "neg", method = "prevalence")

# get the names of contaminants
cons <- rownames(consList)[consList$contaminant=="TRUE"]
# there are 12

# - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# what percentage were the contaminants in the anemone samples?
anemones <- subset_samples(phy, sample_type == "anemone")

# merge the anemone samples
anemones <- merge_samples(anemones, "sample_type", fun = sum)

# transform counts to percentages
anemonesPer <- transform_sample_counts(anemones, function(x) 100 * x/sum(x))

# subset the contaminants
conSub <- prune_taxa(x = anemonesPer, taxa = cons)

# write otu table to dataframe
conData <- data.frame(t(conSub@otu_table))

# add tax data to dataframe
conData <- cbind(conData, conSub@tax_table[,2:7])

# # write data to file for reference
# write.csv(file = "contaminants.csv", x = conData, row.names = FALSE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# remove the contaminants from the main phyloseq file
phy <- remove_taxa(phy, taxa = cons)

rm(consList)
rm(cons)
rm(anemones)
rm(anemonesPer)
rm(conSub)
rm(conData)
