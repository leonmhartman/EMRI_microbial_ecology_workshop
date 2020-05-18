# Here, we test for sgnificant differences between the bacterial community compositions of groups of samples.
# mvabund fits data to a generalised linear model to overcome issues with count data (zeroes, heteroscedasticity, skew).
# Other approaches (e.g. ANOSIM, permANOVA) deal with this by sub-sampling the data, but this results in a loss of power.

# subset anemone samples
anemones <- subset_samples(phy, sample_type == "anemone")
# 530 ASVs

anemones <- prune_taxa((taxa_sums(anemones) > 0), anemones)
# 434 ASVs

# collapse the data to Genus to reduce noise and accelerate computation for this demo
anemonesGen <- tax_glom(anemones, 'Genus', NArm = FALSE)
# 220 ASVs

# get the OTU table & transpose it
otuTab <- t(anemonesGen@otu_table)

# subset the metadata
mvaMeta <- emri_met[1:20,2:3]

# combine the metadata & OTU table
combo <- cbind(mvaMeta, otuTab)

# create mvabund object from count data & check it
comboMva <- mvabund(combo[,3:220]) # I am using data from all rows (i.e. all samples), but only columns 3-220 (i.e. only the ASV counts)
is.mvabund(comboMva)

# check mean vs variance relationship
meanvar.plot(comboMva)
# variance increases with the mean = heteroscedasticity
# no worries, mvabund can cope with this :)

# create model of data
mod <- manyglm(comboMva ~ combo$genotype * combo$samplingDay, family="poisson")

# check that the family gives us a well-fitted model
plot(mod)
# nope, we get a horrible pattern. Let's try again

mod <- manyglm(comboMva ~ combo$genotype * combo$samplingDay, family="negative_binomial")

# check that the family gives us a well-fitted model
plot(mod)
# that's better

# OK, let's run an anova on the model
anova(mod)

#                                      Res.Df   Df.diff   Dev       Pr(>Dev)    
# (Intercept)                          19                           
# combo$genotype                       18       1         439.3     0.010 **  
# combo$samplingDay                    17       1        1033.1     0.001 ***
# combo$genotype:combo$samplingDay     16       1         207.0     0.031 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# The data differ significantly based on genotype & sampling-day.
# However, there is also a significant interaction.
# I'm going to look at the differences between the genotypes on each sampling-day

comboDay1 <- subset(combo, samplingDay == "1")
comboMva1 <- mvabund(comboDay1[,3:220])
mod1 <- manyglm(comboMva1 ~ comboDay1$genotype, family="negative_binomial")
plot(mod1)
anova(mod1)
#                         Res.Df  Df.diff   Dev       Pr(>Dev)  
# (Intercept)             9                         
# comboDay1$genotype      8       1         528.8     0.019 *
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# On Day 1, the genotypes were sig different!
# What about Day 2?

comboDay2 <- subset(combo, samplingDay == "2")
comboMva2 <- mvabund(comboDay2[,3:220])
mod2 <- manyglm(comboMva2 ~ comboDay2$genotype, family="negative_binomial")
plot(mod2)
anova(mod2)
#                         Res.Df  Df.diff   Dev     Pr(>Dev)
# (Intercept)             9                       
# comboDay2$genotype      8       1         180.2     0.233

# By Day 2, the genotypes were not sig different any more!

rm(anemones)
rm(anemonesGen)
rm(otuTab)
rm(mvaMeta)
rm(combo)
rm(comboMva)
rm(mod)
rm(comboDay1)
rm(comboMva1)
rm(mod1)
rm(comboDay2)
rm(comboMva2)
rm(mod2)

