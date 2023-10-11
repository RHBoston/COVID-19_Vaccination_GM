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


## Beta-diversity 
## compute the PCA 
res <- as.data.frame(as.matrix(ps_main_composition@otu_table))  %>% t()
res <- res[, colSums(res) != 0]
res.pca <- prcomp(res, scale = TRUE)

# take the factors out of the metadata which will be used for shape/colour in the PCA
meta <- as.data.frame(as.matrix(ps_main_composition@sam_data))
timepoints <- as.factor(meta$TimepointF)
cohorts <- as.factor(meta$Cohort)

# factor the timepoints = defined order 
meta$TimepointF<-factor(meta$Timepoint, levels=c("Pre-Dose",
                                                 "Acute",
                                                 "Late"))

# Plot the PCA, not in fviz (alpha = 0), but using geom_point() to utilise factors pulled out above
Fig1c <- fviz_pca_biplot(res.pca, label = "none",
                     alpha = 0,
                     invisible = "quali",
                     addEllipses = F, 
                     pointshape = 1,
                     pointsize = 1,
                     legend.title = "Cohort/Timepoint",
                     ggtheme = theme_minimal(),
                     title = " ") +
  labs(x = "PC1 (2.9%)", y = "PC2 (2.41%)") + #make sure to run with and without labels to get right values
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=16),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size=16)) +
  geom_point(aes(colour = factor(meta$TimepointF), shape = factor(meta$Cohort)), size = 4, alpha = 0.65) +
  scale_colour_manual(values = timepoint_colours) +
  scale_shape_manual(values = c(15, 16, 17)) 

# save the pca_plot
ggsave("03_FIGURES/Figure1C.pdf", Fig1c)

# write out data used to make pca
# Fig1C_meta <- cbind(SampleID = rownames(meta_new), meta_new[, c(1, 2, 3, 4)])
# Fig1C_data <- merge(Fi1C_meta, res, by = "row.names")
# saveRDS(Fig1C_data, "00_DATA/Fig1C_data.rds")


## Determine the influence of timepoints on PCA plot
# take the metadata we have and create a new column for our mixed effects linear modelling 
pca_df <- meta_new[,c("Cohort", "Patient", "Vaccine", "Timepoint")]
pca_df$tp_ord <- NA
pca_df[pca_df$Timepoint == "Pre-Dose", "tp_ord"] <- 0
pca_df[pca_df$Timepoint == "Acute", "tp_ord"] <- 1
pca_df[pca_df$Timepoint == "Late", "tp_ord"] <- 2

# write in the PC values for PC1-5 into the dataframe
pca_df$PC1 <- res.pca$x[,1]
pca_df$PC2 <- res.pca$x[,2]
pca_df$PC3 <- res.pca$x[,3]
pca_df$PC4 <- res.pca$x[,4]
pca_df$PC5 <- res.pca$x[,5]

# give the list of cohorts for the loop
cohorts <- c("HC", "ICP", "PID")  

# create df for eventual p values to store into after the loop
p_value_df <- data.frame(Cohort = character(0), PC = character(0), P_Value = numeric(0))

# for each cohort, rotate between PC1-5, fitting the baseline model, then the mixed effect and compare
for (cht in cohorts) {
  for (pc_num in 1:5) {  
    baseline <- lmer(paste0("PC", pc_num, " ~ 1 + (1 | Patient)"), data = subset(pca_df, Cohort == cht))
    rand.intercept.slope <- lmer(paste0("PC", pc_num, " ~ 1 + tp_ord + (tp_ord | Patient)"), data = subset(pca_df, Cohort == cht))
    comparison_result <- model.comparison(baseline, rand.intercept.slope)
    
    # then extract the p-values from the model summary
    p_values <- as.numeric(comparison_result$statistics$p)
    
    # produce the resulting dataframe describing if baseline model is improved on by adding the effect of timepoints
    for (p_value in p_values) {
      if (!is.na(p_value)) {
        p_value_df <- rbind(p_value_df, data.frame(Cohort = cht, PC = paste0("PC", pc_num), P_Value = p_value))
      }
    }
  }
}

# reformat the final data frame with p-values 
model_pvalues <- pivot_wider(p_value_df, names_from = PC, values_from = P_Value)

# print final table as a figure
Supp_Fig1C <- tableGrob(model_pvalues)
grid.arrange(Supp_Fig1C)

# save the figure
ggsave("03_FIGURES/Supp_Fig1C.pdf", grid.arrange(Supp_Fig1C))

# write out the data 
# saveRDS(pca_df, "00_DATA/Supp_Fig1C_data.rds")
