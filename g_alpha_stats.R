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

# create subset of original metadata so we have somewhere to put the anemone diversity metric values
diversityMetrics <- emri_met[1:20,]

# get values for diversity indices and add them to the new metadata file
diversityMetrics <- cbind(diversityMetrics, estimate_richness(anemones, measures = c("Observed","Simpson","Shannon")))

# - - - - - - - - - - - - - - - -

# all values should be checked for assumptions of normality
# here's a test on the a1_1 samples only:
shapiro.test((diversityMetrics$Observed)[1:5])
# W = 0.95661, p-value = 0.7842
# Look at the p-value - the 'Observed' data for the a1_1 samples does not differ significantly from normality :)

# all values should also be checked for assumptions of homogeneity of variance
# here's a test on all the 'Observed' values:
leveneTest(Observed ~ grouping, data = diversityMetrics)
#         Df    F value   Pr(>F)
# group    3    0.4288    0.7351
#         16
# look at the p-value - the variance of the 'Observed' data does not differ significantly from homogeneity :)

# if your data fail these tests, linear models (e.g. ANOVA) may not be appropriate
# instead, you may need a non-linear test e.g. GLS

# - - - - - - - - - - - - - - - -

# tests of OVERALL differences between the diversity metric values

summary(aov(Observed ~ genotype * samplingDay, data = diversityMetrics))
#                       Df    Sum Sq  Mean Sq   F value   Pr(>F)    
# genotype              1     110     110       0.561     0.465    
# samplingDay           1     14742   14742     74.844    1.98e-07 ***
# genotype:samplingDay  1     61      61        0.311     0.585    
# Residuals             16    3152    197                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# No sig dif based on 'Genotype'
# Sig dif based on 'Day'
# No sig interaction


summary(aov(Simpson ~ genotype * samplingDay, data = diversityMetrics))
#                       Df  Sum Sq    Mean Sq     F value   Pr(>F)  
# genotype              1   0.00121   0.001205    0.554     0.4675  
# samplingDay           1   0.01459   0.014587    6.706     0.0198 *
# genotype:samplingDay  1   0.00031   0.000307    0.141     0.7120  
# Residuals             16  0.03480   0.002175                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


summary(aov(Shannon ~ genotype * samplingDay, data = diversityMetrics))
#                       Df  Sum Sq  Mean Sq   F value   Pr(>F)    
# genotype              1   0.074   0.074     0.522     0.480568    
# samplingDay           1   3.636   3.636     25.549    0.000117 ***
# genotype:samplingDay  1   0.167   0.167     1.176     0.294253    
# Residuals             16  2.277   0.142                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# - - - - - - - - - - - - - - - -

# which pairs of samples were significantly different?
# Because the data meet the normality & homogeneity of variance assumptions, we can simply do pair-wise t-tests

pairwise.t.test(x = diversityMetrics$Observed, g = diversityMetrics$grouping, p.adjust.method = "holm")
#         a1_1        a1_2      a4_1   
# a1_2    3.6e-05     -         -      
# a4_1    0.73863     0.00013   -      
# a4_2    3.4e-05     0.89415   0.00013

pairwise.t.test(x = diversityMetrics$Simpson, g = diversityMetrics$grouping, p.adjust.method = "holm")
#         a1_1    a1_2    a4_1
# a1_2    0.55    -       -   
# a4_1    0.88    0.19    -   
# a4_2    0.63    0.88    0.26

pairwise.t.test(x = diversityMetrics$Shannon, g = diversityMetrics$grouping, p.adjust.method = "holm")
#         a1_1      a1_2      a4_1  
# a1_2    0.0379    -         -     
# a4_1    0.4393    0.0043    -     
# a4_2    0.0297    0.8011    0.0030

# Do these results make seem resonable when we look at the box-plots?

rm(anemones)
rm(diversityMetrics)

