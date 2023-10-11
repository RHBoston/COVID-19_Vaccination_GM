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

## Pre-Processing to get appropriate paired analysis dataframes
# add variable to subset data by 
tmp <- phyloseq::sample_data(ps_main_composition)
tmp$Patient_Vaccine <- paste0(tmp$Patient, "_", tmp$Vaccine)
sample_data(ps_main_composition) <- phyloseq::sample_data(tmp)

# remove wrong timepoint samples
ps_main_composition <- subset_samples(ps_main_composition, !(rownames(sample_data(ps_main_composition)) %in% "ICP22_426"))
ps_main_composition <- subset_samples(ps_main_composition, !(rownames(sample_data(ps_main_composition)) %in% "ICP5_791"))

# subset dataframes to be paired 
aPD_A <- subset_samples(ps_main_composition, (Timepoint == "Pre-Dose" | Timepoint == "Acute"))
aA_L <- subset_samples(ps_main_composition, (Timepoint == "Acute" | Timepoint == "Late"))
aPD_L <- subset_samples(ps_main_composition, (Timepoint == "Pre-Dose" | Timepoint == "Late"))

aPD_A_sample_counts <- table(sample_data(aPD_A)$Patient_Vaccine)
aA_L_sample_counts <- table(sample_data(aA_L)$Patient_Vaccine)
aPD_L_sample_counts <- table(sample_data(aPD_L)$Patient_Vaccine)

aPD_A_final <- subset_samples(aPD_A, Patient_Vaccine %in% names(aPD_A_sample_counts[aPD_A_sample_counts == 2]))
aA_L_final <- subset_samples(aA_L, Patient_Vaccine %in% names(aA_L_sample_counts[aA_L_sample_counts == 2]))
aPD_L_final <- subset_samples(aPD_L, Patient_Vaccine %in% names(aPD_L_sample_counts[aPD_L_sample_counts == 2]))

# remove what you don't need
rm(aPD_A_sample_counts, aA_L_sample_counts, aPD_L_sample_counts)
rm(aPD_A, aA_L, aPD_L, tmp)

# merge final phyloseqs together
patient_linked <- merge_phyloseq(aPD_A_final, aA_L_final, aPD_L_final)
patient_linked <- microbiome::transform(patient_linked, "compositional")


## Alpha diversity
patient_linked_alpha <- sample_data(patient_linked) %>% as.matrix() %>% as.data.frame()
patient_linked_alpha_long <- pivot_longer(patient_linked_alpha, cols = c("chao1", "diversity_shannon"),
                             names_to = "Diversity_Measure",
                             values_to = "Value")

patient_linked_alpha_long$Value <- as.numeric(patient_linked_alpha_long$Value)
patient_linked_alpha_long$TimepointF <- factor(patient_linked_alpha_long$Timepoint, levels=c("Pre-Dose", "Acute","Late"))

Supp_Fig1B <- ggplot(data = subset(patient_linked_alpha_long),
                   aes(x=TimepointF, 
                       y=Value,
                       group = Patient_Vaccine,
                       colour = Cohort)) + 
              geom_point(size = 4) + geom_line() +
              facet_wrap(~Diversity_Measure, scales = 'free_y', ncol =2) +
              scale_colour_manual(values = cohort_colours) +
              theme(axis.text = element_text(size=16),
                    axis.title = element_blank(),
                    legend.position = "bottom",
                    legend.text = element_text(size = 16),
                    legend.title = element_text(size=16),
                    strip.text.x = element_text(size = 16))

ggsave("03_FIGURES/Supp_Fig1B.pdf", Supp_Fig1B)


# statistics
patient_linked_alpha_stats <- sample_data(patient_linked) %>% as.matrix() %>% as.data.frame() %>% rownames_to_column(var = "Sample")
patient_linked_alpha_stats$chao1 <- as.numeric(patient_linked_alpha_stats$chao1)
patient_linked_alpha_stats$diversity_shannon <- as.numeric(patient_linked_alpha_stats$diversity_shannon)

cohorts <- unique(patient_linked_alpha_stats$Cohort)
all_pvalues <- c()
t1 <- "Pre-Dose"
diversity_metrics <- c("chao1", "diversity_shannon")

for (k in cohorts) {
  sub <- subset(patient_linked_alpha_stats, Cohort == k)
  
  for (t2 in c("Acute", "Late")) {
    print(paste(t1, t2))
    
    for (metric in diversity_metrics) {
      this_wlc <- sub[sub$Timepoint %in% c(t1, t2), ]
      this_count <- this_wlc %>% group_by(Patient_Vaccine) %>% 
        summarise(n = length(Sample))
      
      this_count <- this_count[this_count$n > 1, ]
      this_wlc <- this_wlc[this_wlc$Patient_Vaccine %in% this_count$Patient_Vaccine, ]
      this_stat <- t.test(this_wlc[this_wlc$Timepoint == t1, ][[metric]],
                          this_wlc[this_wlc$Timepoint == t2, ][[metric]], 
                          alternative = "greater", paired = T)
      
      all_pvalues <- rbind(all_pvalues, c(k, t1, t2, metric, this_stat$p.value))
    }
  }
}

colnames(all_pvalues) <- c("Cohort", "Timepoint_1", "Timepoint_2", "Diversity_Metric", "P_value")
all_pvalues <- as.data.frame(all_pvalues)
all_pvalues$P_adj <- p.adjust(all_pvalues$P_value, method = "fdr")

# write out data
# SuppFig1B_data <- patient_linked_alpha[, c(1:7,34:35)]
# saveRDS(SuppFig1B_data, "00_DATA/Supp_Fig1B_data.rds")



## Composition
# merge and melt on phylum level
patient_linked_phy <- tax_glom(patient_linked, taxrank = "Phylum")
patient_linked_phy_melt <- psmelt((patient_linked_phy))

Supp_Fig1D <-  ggplot(subset(patient_linked_phy_melt, Phylum %in%c("p__Firmicutes", 
                                                     "p__Bacteroidetes", 
                                                     "p__Proteobacteria", 
                                                     "p__Actinobacteria", 
                                                     "p__Verrucomicrobia",
                                                     "p__Euryarchaeota")),
                     aes(x=TimepointF, 
                         y=Abundance,
                         group = Patient_Vaccine,
                         colour = Cohort)) + 
                geom_point(size =4) + geom_line() +
                facet_wrap(~Phylum, scales = 'free_y') +
                scale_colour_manual(values = cohort_colours) +
                labs(y = "Relative Abundance") + 
                theme(axis.text = element_text(size=16),
                      axis.title.x = element_blank(),
                      axis.title.y = element_text(size=20),
                      legend.position = "bottom",
                      legend.text = element_text(size = 16),
                      legend.title = element_text(size=16),
                      strip.text.x = element_text(size = 16)) 

ggsave("03_FIGURES/Supp_Fig1B.pdf", Supp_Fig1D)

# write out data
# SuppFig1E_data <- patient_linked_phy_melt[, c(1:7,40:41)]
# saveRDS(SuppFig1E_data, "00_DATA/Supp_Fig1E_data.rds")

# statistics 
patient_linked_phy_stats <- unique(patient_linked_phy_melt[, c("Sample", "Patient_Vaccine", "Timepoint", "Cohort", "Phylum", "Abundance")])

cmp_phyla <- unique(patient_linked_phy_stats$Phylum)
cohorts <- unique(patient_linked_phy_stats$Cohort)
all_phyla <- c()
t1 <- "Pre-Dose"

for (k in cohorts){
  sub <- subset(patient_linked_phy_stats, Cohort == k)
  
  for (t2 in c("Acute", "Late")) {
    print(paste(t1, t2))
    this_wlc <- sub[sub$Timepoint %in% c(t1, t2), ]
    
    for (p in cmp_phyla) {
      this_wlc_p <- this_wlc[this_wlc$Phylum == p, ]
      this_count <- this_wlc_p %>% group_by(Patient_Vaccine) %>% 
        summarise(n = length(Sample))
      
      this_count <- this_count[this_count$n > 1, ]
      this_wlc_p <- this_wlc_p[this_wlc_p$Patient_Vaccine %in% this_count$Patient_Vaccine, ]
      this_stat <- wilcox.test(this_wlc_p[this_wlc_p$Timepoint == t1, ]$Abundance,
                          this_wlc_p[this_wlc_p$Timepoint == t2, ]$Abundance,
                          alternative = "greater", paired = T)
      all_phyla <- rbind(all_phyla, c(k, t1, t2, p, this_stat$p.value))
    }
  }
}
colnames(all_phyla) <- c("Cohort","Timepoint_1", "Timepoint_2", "Phylum", "P_value")
all_phyla <- as.data.frame(all_phyla)
all_phyla$P_adj <- p.adjust(all_phyla$P_value, method = "fdr")

saveRDS(all_pvalues, "02_RESULTS/patient_linked_phylum_statistics.rds")


# DESeq 

sample_data(patient_linked)$Timepoint <- as.factor(sample_data(patient_linked)$Timepoint)
patient_linked_DESeq <- patient_linked

tax_table(patient_linked_DESeq) <- tax_table(patient_linked)[,c(1, 2 ,7)]
ps.D <- tax_glom(patient_linked_DESeq, taxrank = 'Species', NArm = F)

# pairwise comparison between two things of choice = chose a cohort and leave out one of the timepoints
ps.D.sub <- subset_samples(ps.D, Cohort %in% "PID" & Timepoint != "Late")

# filter sparse features, with > 95% zeros
ps.D.sub <- prune_taxa(rowSums(otu_table(ps.D.sub) == 0) < ncol(otu_table(ps.D.sub)) * 0.95, ps.D.sub)

#convert otu  to integers 
otu_table(ps.D.sub) <- otu_table(ps.D.sub)*100000

# convert to DESeq
ps_ds = phyloseq_to_deseq2(ps.D.sub, ~ Patient_Vaccine + Timepoint)

# use alternative estimator on a condition of "every gene contains a sample with a zero"
ds <- estimateSizeFactors(ps_ds, type="poscounts")
ds = DESeq(ds, test="Wald", fitType="local")


alpha = 0.05 
res = results(ds, alpha=alpha, cooksCutoff = FALSE)
res = res[order(res$pvalue, na.last=NA), ]
taxa_sig = rownames(res[1:35, ])
taxa_sig_check = res[(res$padj < alpha), ]

