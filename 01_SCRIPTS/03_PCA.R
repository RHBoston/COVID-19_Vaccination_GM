# set a working directory 
setwd("~/Desktop/COVID-19_Vaccination_GM")

# source R settings
source("01_SCRIPTS/00_Settings.R")

# read pre-processed files back in
meta_new <- readRDS("02_RESULTS/Modified_metadata.rds")
Tax_table_new <- readRDS("02_RESULTS/Modified_tax_table.rds")
Otu_table <- readRDS("02_RESULTS/Modified_otu_table.rds")
ps <- readRDS("02_RESULTS/Phyloseq.rds")

# transform takes abundances (not = to 1 within sample) and makes them relative abundances (=1 within sample)
ps_main_composition <- microbiome::transform(ps, "compositional")


## Beta-diversity 
## compute the PCA 
res <- as.data.frame(as.matrix(ps_main_composition@otu_table)) %>% t()
res <- res[, colSums(res) != 0]
res.pca <- prcomp(res, scale = TRUE)


# take the factors out of the metadata which will be used for shape/colour in the PCA
meta <- as.data.frame(as.matrix(ps_main_composition@sam_data))
meta$TimepointF<-factor(meta$Timepoint, levels=c("Pre-Dose",
                                                 "Acute",
                                                 "Late"))

# Plot the PCA, not in fviz (alpha = 0), but using geom_point() to utilise factors pulled out above
Fig1c <- fviz_pca_ind(res.pca, label = "none",
                        # select.var = list(contrib = 107),
                     alpha = 0,
                     invisible = "quali",
                     addEllipses = F, 
                     pointshape = 1,
                     pointsize = 1,
                     ggtheme = theme_minimal(),
                     title = " ") +
  labs(x = "PC1 (2.9%)", y = "PC2 (2.4%)", shape="Cohort", colour="Timepoint") + #make sure to run with and without labels to get computed values
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
cohorts <- unique(pca_df$Cohort)

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
Supp_Fig2C <- tableGrob(model_pvalues)
grid.arrange(Supp_Fig2C)

# save the figure
ggsave("03_FIGURES/Supp_Fig2C.pdf", grid.arrange(Supp_Fig2C))


## PERMANOVA
otu <- as.data.frame(t(as.matrix(otu_table(ps_main_composition))))
meta <- as.data.frame(as.matrix(sample_data(ps_main_composition)))

dist <- as.matrix(vegdist(otu, method = "euclidean"))

#run the test
permanova <- adonis2(formula = dist ~ Patient*Vaccine*Timepoint + Sex + Age + Therapy + Disease,
                    data = meta, permutations=999)

# print final table as a figure
perma_plot <- tableGrob(permanova)
title <- textGrob("adonis2(formula = dist ~ Patient*Vaccine*Timepoint + Sex + Age + Therapy + Disease, 
                  data = meta, permutations = 999, method = euclidean)", gp=gpar(fontsize=14))

grid.arrange(title, perma_plot, ncol = 1, heights = c(15,15))

# print final table as a figure
perma_table <- tableGrob(permanova)
grid.arrange(perma_table)

ggsave("03_FIGURES/Supp_Table2.pdf", grid.arrange(perma_table))
saveRDS(permanova, "02_RESULTS/permanova.rds")



## NMDS b-diversity ---- EXTRA ---- 
# run the ordination 
ord.nmds.bray <- ordinate(ps_main_composition, method="NMDS", distance="bray")

# extract NMDS1 and NMDS2 values
nmds_data <- data.frame(ord.nmds.bray$points)
# merge with sample data from the phyloseq object
nmds_df <- merge(nmds_data, meta, by = 0, all = TRUE)

meta$TimepointF<-factor(meta$TimepointF, 
                                   levels=c("Pre-Dose", "Acute", "Late"))

# create the NMDS plot
nmds_ct <-  ggplot(nmds_df,
                   aes(x = MDS1, y = MDS2)) +
  geom_point(aes(colour = factor(meta$TimepointF), shape = factor(meta$Cohort)), size = 4, alpha = 0.65) +
  scale_colour_manual(values = timepoint_colours) +
  scale_shape_manual(values = c(15, 16, 17)) +
  labs(x = "NMDS1", y = "NMDS2", shape="Cohort", colour="Timepoint") +
  theme(axis.text = element_text(size=16), 
        axis.title = element_text(size=16),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size=16))

nmds_coh <- ggplot(nmds_df,
                   aes(x = MDS1, y = MDS2)) +
  geom_point(aes(colour = factor(meta$Cohort), shape = factor(meta$Cohort)), size = 4, alpha = 0.65) +
  scale_colour_manual(values = cohort_colours) +
  scale_shape_manual(values = c(15, 16, 17)) +
  labs(x = "NMDS1", y = "NMDS2", shape="Cohort", colour="Cohort") +
  theme(axis.text = element_text(size=16), 
        axis.title = element_text(size=16),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size=16))

ggsave("03_FIGURES/nmds_cohort_timepoint.pdf", nmds_ct)
ggsave("03_FIGURES/nmds_cohort.pdf", nmds_coh)


## Determine the influence of timepoints on NMDS plot
pca_df$NMDS1 <- nmds_data[,1]
pca_df$NMDS2 <- nmds_data[,2]

# create df for eventual p values to store into after the loop
p_value_df2 <- data.frame(Cohort = character(0), NMDS = character(0), P_Value = numeric(0))

# for each cohort, rotate between NMDS1 + NMDS2, fitting the baseline model, then the mixed effect and compare
for (cht in cohorts) {
  for (nmds_num in 1:2) {  
    baseline <- lmer(paste0("NMDS", nmds_num, " ~ 1 + (1 | Patient)"), data = subset(pca_df, Cohort == cht))
    rand.intercept.slope <- lmer(paste0("NMDS", nmds_num, " ~ 1 + tp_ord + (tp_ord | Patient)"), data = subset(pca_df, Cohort == cht))
    comparison_result <- model.comparison(baseline, rand.intercept.slope)
    
    # then extract the p-values from the model summary
    p_values <- as.numeric(comparison_result$statistics$p)
    
    # produce the resulting dataframe describing if baseline model is improved on by adding the effect of timepoints
    for (p_value in p_values) {
      if (!is.na(p_value)) {
        p_value_df2 <- rbind(p_value_df2, data.frame(Cohort = cht, NMDS = paste0("NMDS", nmds_num), P_Value = p_value))
      }
    }
  }
}

# reformat the final data frame with p-values 
model_pvalues2 <- pivot_wider(p_value_df2, names_from = NMDS, values_from = P_Value)

# print final table as a figure
nmds_plot <- tableGrob(model_pvalues2)
grid.arrange(nmds_plot)
ggsave("03_FIGURES/nmds_stats.pdf", grid.arrange(nmds_plot))


# NMDS like permanova
dist <- as.matrix(vegdist(otu, method = "bray"))
#run the test
permanova_bray <- adonis2(formula = dist ~ Patient*Vaccine*Timepoint + Sex + Age + Therapy + Disease,
                     data = meta, permutations=999)

# print final table as a figure
perma_plot_bray <- tableGrob(permanova_bray)
grid.arrange(perma_plot_bray)

ggsave("03_FIGURES/nmds_perma_stats.pdf", grid.arrange(perma_plot_bray))
