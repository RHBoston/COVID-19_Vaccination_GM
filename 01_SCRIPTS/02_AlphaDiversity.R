# set a working directory
setwd("~/Desktop/COVID-19_Vaccination_GM")

# source R settings
source("01_SCRIPTS/Settings.R")

# read pre-processed files back in
meta_new <- readRDS("02_RESULTS/Modified_metadata.rds")
Tax_table_new <- readRDS("02_RESULTS/Modified_tax_table.rds")
Otu_table <- readRDS("02_RESULTS/Modified_otu_table.rds")
ps <- readRDS("02_RESULTS/Phyloseq.rds")

# transform takes abundances (not = to 1 within sample) and makes them relative abundances (=1 within sample)
ps_main_composition <- microbiome::transform(ps, "compositional")


## Assess the alpha diversity of the samples
# Make df long for eventual plotting in facets
meta_new_long <- pivot_longer(meta_new, cols = c("chao1", "diversity_shannon"),
                              names_to = "Diversity_Measure",
                              values_to = "Value")

#The above makes values = characters, so move them back to numeric
meta_new_long$Value <- as.numeric(meta_new_long$Value)

#  ggplot to produce results + create figures
Fig1B <- ggplot(meta_new_long, 
                aes(x = Cohort, 
                    y = Value,
                    colour = TimepointF,
                    na.rm = T)) +
  geom_boxplot() +
  geom_point(position=position_jitterdodge(jitter.width = 0), size = 4) +
  facet_wrap(~Diversity_Measure, scales = 'free_y') +
  labs(x= "Cohort", y= "Diversity Measure") +
  scale_colour_manual(values = timepoint_colours, name = "Timepoint") +
  theme(axis.text = element_text(size=16), 
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size=16),
        strip.text.x = element_text(size = 20)) #+ 
# geom_pwc(aes(group = Cohort),method = "wilcox_test", p.adjust.method = "fdr", tip.length = 0) 
# adds space for p-val, padj wrong, adjust accordingly

Supp_Fig1A <- ggplot(meta_new_long, 
                     aes(x = Cohort, 
                         y = Value,
                         colour = Cohort,
                         na.rm = T)) +
  geom_boxplot() +
  geom_point(position=position_jitterdodge(jitter.width = .1), size = 4) +
  facet_wrap(~Diversity_Measure, scales = 'free_y') +
  labs(x= "Cohort", y= "Diversity Measure") +
  scale_colour_manual(values = cohort_colours, name = "Cohort") +
  theme(axis.text = element_text(size=16), 
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size=16),
        strip.text.x = element_text(size = 20)) #+ 
# geom_pwc(aes(group = Cohort),method = "wilcox_test", p.adjust.method = "fdr", tip.length = 0) 
# adds space for p-val, padj wrong, adjust accordingly

# save pdfs
ggsave("03_FIGURES/Figure1B.pdf", Fig1B)
ggsave("03_FIGURES/Supp_Fig1A.pdf", Supp_Fig1A)


## loop to test normality of the alpha div measures 
# create empty df for normality testing
alpha_diversity_normality <- c()

# subset out the alpha div measures from the metadata
subset_alpha <- meta_new[, c("chao1", "diversity_shannon")]

# loop through each column of the data frame and test for normality
for (col in names(subset_alpha)) {
  p_value <- shapiro.test(subset_alpha[[col]])$p.value
  
  if (p_value < 0.05) {
    is_normal <- FALSE
  } else {
    is_normal <- TRUE
  }
  alpha_diversity_normality <- rbind(alpha_diversity_normality, c(col, is_normal))
}

# view the results
alpha_diversity_normality <- as.data.frame(alpha_diversity_normality)
alpha_diversity_normality 

## Statistical testing of alpha div values for graphing 
# format values to be numeric
meta_new$chao1 <- as.numeric(meta_new$chao1)
meta_new$diversity_shannon <- as.numeric(meta_new$diversity_shannon)

#take the unique values within the cohort/timepoint combinations 
grp <- unique(meta_new$Cohort_Timepoint)

# make a new file for the results of the stats test to go into
chao1_shannon_statistics <- c()
passed_tests <- character(0)

# use a wilcox test to see significance of alpha div between cohort timepoints the cohort (change cohort_timepoint to cohort after)
for (i in grp) {
  this_i <- meta_new[meta_new$Cohort_Timepoint == i,]
  for (j in grp) {
    this_j <- meta_new[meta_new$Cohort_Timepoint == j,]
    this_p_shannon <- wilcox.test(this_i$diversity_shannon, this_j$diversity_shannon, exact = FALSE)$p.value
    this_p_chao1 <- wilcox.test(this_i$chao1, this_j$chao1, exact = FALSE)$p.value
    if (i != j & length(grep(paste(j, i, sep = "__"), passed_tests)) == 0) {
      this_output_shannon <- c(i, j, this_p_shannon)
      chao1_shannon_statistics <- rbind(chao1_shannon_statistics, this_output_shannon)
      this_output_chao1 <- c(i, j, this_p_chao1)
      chao1_shannon_statistics <- rbind(chao1_shannon_statistics, this_output_chao1)
      
      passed_tests[length(passed_tests) + 1] <- paste(i, j, sep = "__")
    }
  }
}

# edit loop output to be able to adjust for fdr
chao1_shannon_statistics <- as.data.frame(chao1_shannon_statistics)
colnames(chao1_shannon_statistics) <- c("Group_1", "Group_2", "P_value")
chao1_shannon_statistics$FDR <- p.adjust(chao1_shannon_statistics$P_value, method = 'fdr') 

# save the resulting statistics document
saveRDS(chao1_shannon_statistics, "02_RESULTS/chao1_shannon_statistics.rds")

# write out data 
# Fig1B_SuppFig1A_data <- cbind(SampleID = rownames(meta_new), meta_new[, c(1, 2, 3, 4, 34, 35)])
# saveRDS(Fig1B_SuppFig1A_data, "00_DATA/Fig1B_SuppFig1A_data.rds")
