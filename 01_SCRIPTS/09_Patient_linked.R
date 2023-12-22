# processing -------
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


saveRDS(patient_linked, "02_RESULTS/Patient_linked_ps.rds")

## Alpha diversity ----

# dataframe formatting
patient_linked_alpha <- sample_data(patient_linked) %>% as.matrix() %>% as.data.frame() %>% rownames_to_column(var = "Sample")
patient_linked_alpha_long <- pivot_longer(patient_linked_alpha, cols = c("chao1", "diversity_shannon"),
                             names_to = "Diversity_Measure",
                             values_to = "Value")

patient_linked_alpha_long$Value <- as.numeric(patient_linked_alpha_long$Value)
patient_linked_alpha_long$TimepointF <- factor(patient_linked_alpha_long$Timepoint, 
                                               levels=c("Pre-Dose", "Acute","Late"))
patient_linked_alpha_long$Vaccine_TimepointF<-factor(patient_linked_alpha_long$Vaccine_Timepoint, 
                                        levels=c("1_Pre-Dose", "1_Acute", "1_Late", "2_Pre-Dose", "2_Acute", "2_Late", "3_Pre-Dose", "3_Acute", "3_Late"))


Supp_Fig2B <- ggplot(data = patient_linked_alpha_long,
                   aes(x=Vaccine_TimepointF, 
                       y=Value,
                       group = Patient_Vaccine,
                       colour = Cohort)) + 
              geom_point(size = 4) + geom_line() +
              scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
              facet_wrap(~Diversity_Measure, scales = 'free_y', ncol =1) +
              scale_colour_manual(values = cohort_colours) +
              theme(axis.text = element_text(size=16),
                    axis.title = element_blank(),
                    legend.position = "bottom",
                    legend.text = element_text(size = 16),
                    legend.title = element_text(size=16),
                    strip.text.x = element_text(size = 16))  

ggsave("03_FIGURES/Supp_Fig2B.pdf", Supp_Fig2B, height = 8.57, width = 10.8, units = "in")


# statistics
patient_linked_alpha_stats <- patient_linked_alpha 
patient_linked_alpha_stats$chao1 <- as.numeric(patient_linked_alpha_stats$chao1)
patient_linked_alpha_stats$diversity_shannon <- as.numeric(patient_linked_alpha_stats$diversity_shannon)
patient_linked_alpha_stats$Vaccine <- as.character(patient_linked_alpha_stats$Vaccine)

cohorts <- unique(patient_linked_alpha_stats$Cohort)
vaccine <- unique(patient_linked_alpha_stats$Vaccine)
timepoint <- unique(patient_linked_alpha_stats$Timepoint)

all_pvalues <- c()
diversity_metrics <- c("chao1", "diversity_shannon")

for (v in vaccine) {
  v_sub <- subset(patient_linked_alpha_stats, Vaccine == v)
  
  for (k in cohorts) {
    sub <- subset(v_sub, Cohort == k)
    
    for (t1 in c("Pre-Dose", "Acute")) {
      for (t2 in c("Acute", "Late")) {
        print(paste(t1, t2))
        
        for (metric in diversity_metrics) {
          this_wlc <- sub[sub$Timepoint %in% c(t1, t2), ]
          this_count <- this_wlc %>% group_by(Patient_Vaccine) %>% 
            summarise(n = length(Sample))
          
          this_count <- this_count[this_count$n > 1, ]
          this_wlc <- this_wlc[this_wlc$Patient_Vaccine %in% this_count$Patient_Vaccine, ]
          if(nrow(this_wlc) < 6) next
          
          n1 <- length(unique(this_wlc[this_wlc$Timepoint == t1, ]$Sample))
          
          this_stat <- wilcox.test(this_wlc[this_wlc$Timepoint == t1, ][[metric]],
                              this_wlc[this_wlc$Timepoint == t2, ][[metric]], exact = F,
                              paired = T)
          
          all_pvalues <- rbind(all_pvalues, c(metric, v, k, t1, t2, n1, this_stat$p.value))
        }
      }
    }
  }
}

colnames(all_pvalues) <- c("Diversity_Measure", "Vaccine", "Cohort", "Timepoint_1", "Timepoint_2", "Samples", "P_value")
all_pvalues <- as.data.frame(all_pvalues)
final_alpha_stats <- c()

for (v in vaccine) {
  v_stats <- subset(all_pvalues, Vaccine == v)
  
  for (c in cohorts) {
    cohort_rows <- subset(v_stats, Cohort == c)
    
    chao1_rows <- cohort_rows$Diversity_Measure == "chao1"
    cohort_rows$fdr[chao1_rows] <- p.adjust(cohort_rows$P_value[chao1_rows], method = "bonferroni")
    
    diversity_shannon_rows <- cohort_rows$Diversity_Measure == "diversity_shannon"
    cohort_rows$fdr[diversity_shannon_rows] <- p.adjust(cohort_rows$P_value[diversity_shannon_rows], method = "bonferroni")
    
    final_alpha_stats <- rbind(final_alpha_stats, cohort_rows)
  }
}

final_alpha_stats <- final_alpha_stats[order(final_alpha_stats$Vaccine, final_alpha_stats$Diversity_Measure),]


# print final table as a figure
Supp_Tab1 <- tableGrob(final_alpha_stats)
grid.arrange(Supp_Tab1)

ggsave("03_FIGURES/Supp_Table1.pdf", grid.arrange(Supp_Tab1))
saveRDS(final_alpha_stats, "02_RESULTS/patient_linked_alpha_statistics.rds")


## Composition -----
# merge and melt on phylum level
patient_linked_phy <- tax_glom(patient_linked, taxrank = "Phylum")
patient_linked_phy_melt <- psmelt((patient_linked_phy))

patient_linked_phy_melt$Vaccine_TimepointF<-factor(patient_linked_phy_melt$Vaccine_Timepoint, 
                                                     levels=c("1_Pre-Dose", "1_Acute", "1_Late", "2_Pre-Dose", "2_Acute", "2_Late", "3_Pre-Dose", "3_Acute", "3_Late"))


Supp_Fig2E <-  ggplot(subset(patient_linked_phy_melt, Phylum %in%c("p__Lentisphaerae",  
                                                     "p__Bacteroidetes", 
                                                     "p__Proteobacteria", 
                                                     "p__Actinobacteria", 
                                                     "p__Verrucomicrobia",
                                                     "p__Euryarchaeota")),
                     aes(x=Vaccine_TimepointF, 
                         y=Abundance,
                         group = Patient_Vaccine,
                         colour = Cohort)) + 
                geom_point(size =4) + geom_line() +
                facet_wrap(~Phylum, scales = 'free_y') +
                scale_colour_manual(values = cohort_colours) +
                labs(y = "Relative Abundance") + 
                theme(axis.text.y = element_text(size=16),
                      axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size = 16),
                      axis.title.x = element_blank(),
                      axis.title.y = element_text(size=20),
                      legend.position = "bottom",
                      legend.text = element_text(size = 16),
                      legend.title = element_text(size=16),
                      strip.text.x = element_text(size = 16)) 

ggsave("03_FIGURES/Supp_Fig2E.pdf", Supp_Fig2E, height = 8, width = 12.5, units = "in")


# statistics 
patient_linked_phy_stats <- unique(patient_linked_phy_melt[, c("Sample", "Vaccine", "Patient", "Patient_Vaccine", "Timepoint", "Cohort", "Phylum", "Abundance")])
vaccine <- unique(patient_linked_phy_stats$Vaccine)
cohorts <- unique(patient_linked_phy_stats$Cohort)
cmp_phyla <- unique(patient_linked_phy_stats$Phylum)

all_phyla <- c()

for (v in vaccine) {
  v_sub <- subset(patient_linked_phy_stats, Vaccine == v)
  
for (k in cohorts){
  sub <- subset(v_sub, Cohort == k)
  
  for (t1 in c("Pre-Dose", "Acute")) {
  for (t2 in c("Acute", "Late")) {
    print(paste(t1, t2))
    this_wlc <- sub[sub$Timepoint %in% c(t1, t2), ]
    
    for (p in cmp_phyla) {
      this_wlc_p <- this_wlc[this_wlc$Phylum == p, ]
      this_count <- this_wlc_p %>% group_by(Patient_Vaccine) %>% 
        summarise(n = length(Sample))
      
      this_count <- this_count[this_count$n > 1, ]
      this_wlc_p <- this_wlc_p[this_wlc_p$Patient_Vaccine %in% this_count$Patient_Vaccine, ]
      if(length(unique(this_wlc_p$Timepoint)) < 2) next
      
      a1 <- this_wlc_p[this_wlc_p$Timepoint == t1, ]$Abundance
      a2 <- this_wlc_p[this_wlc_p$Timepoint == t2, ]$Abundance
      
      a3 <- median((a1/a2), na.rm = TRUE)
      
      n1 <- length(unique(this_wlc_p[this_wlc_p$Timepoint == t1, ]$Sample))
      n2 <- length(unique(this_wlc_p[this_wlc_p$Timepoint == t2, ]$Sample))
      if (n1 != n2) stop("not paired")
      
      this_stat <- wilcox.test(this_wlc_p[this_wlc_p$Timepoint == t1, ]$Abundance,
                          this_wlc_p[this_wlc_p$Timepoint == t2, ]$Abundance, exact = F, 
                          paired = T)
      
      length(unique(this_wlc_p[this_wlc_p$Timepoint == t1, ]$Sample))
      
      all_phyla <- rbind(all_phyla, c(v, k, t1, t2, p, a3, n1, this_stat$p.value))
        }
      }
    }
  }
}

colnames(all_phyla) <- c("Vaccine", "Cohort","Timepoint_1", "Timepoint_2", "Phylum", "foldchange", "comparisons", "P_value")
all_phyla <- as.data.frame(all_phyla)


final_phylum_paired <- all_phyla[0,]
final_phylum_paired[1,] <- NA
final_phylum_paired$padj <- NA
final_phylum_paired <- final_phylum_paired[0,]

for (v in vaccine) {
  v_stats <- subset(all_phyla, (Vaccine == v &
                                  !is.na(foldchange) &
                                  abs(as.numeric(foldchange)) != Inf &
                                  as.numeric(foldchange) != 0))
  
  for (c in cohorts) {
    cohort_rows <- subset(v_stats, Cohort == c)
    cohort_rows$padj <- p.adjust(cohort_rows$P_value, method = "bonferroni")
    final_phylum_paired <- rbind(final_phylum_paired, cohort_rows)
  }
}

final_phylum_paired$log2FC <- log2(as.numeric(final_phylum_paired$foldchange))
final_phylum_paired$full_id <- paste(final_phylum_paired$Vaccine, final_phylum_paired$Cohort, final_phylum_paired$Timepoint_1, final_phylum_paired$Timepoint_2, final_phylum_paired$Phylum, sep = "_")
Supp_Fig2F <- EnhancedVolcano(final_phylum_paired, 
                                 lab = final_phylum_paired$full_id, 
                                 x = "log2FC", y = "padj", 
                                 pCutoff = 0.05,
                                 title = NULL,
                                 subtitle = NULL,
                                 caption = NULL,
                                 pointSize = 4)
ggsave("03_FIGURES/Supp_Fig2F.pdf", Supp_Fig2F, height = 5, width = 5, units = "in")



saveRDS(final_phylum_paired, "02_RESULTS/patient_linked_phylum_statistics.rds")
presented_final_phylum_paired <- subset(final_phylum_paired, padj <1)
presented_final_phylum_paired <- presented_final_phylum_paired[order(presented_final_phylum_paired$padj),]

saveRDS(presented_final_phylum_paired, "02_RESULTS/presented_patient_linked_phylum_statistics.rds")

# print final table as a figure
Supp_Tab3 <- tableGrob(presented_final_phylum_paired)
grid.arrange(Supp_Tab3)

ggsave("03_FIGURES/Supp_Table3.pdf", grid.arrange(Supp_Tab3))


# DESeq  ----

sample_data(patient_linked)$Timepoint <- as.factor(sample_data(patient_linked)$Timepoint)
patient_linked_DESeq <- patient_linked

tax_table(patient_linked_DESeq) <- tax_table(patient_linked)[,c(1, 2 ,7)]
ps.D <- tax_glom(patient_linked_DESeq, taxrank = 'Species', NArm = F)

# pairwise comparison between two things of choice = chose a cohort and leave out one of the timepoints
ps.D.sub <- subset_samples(ps.D, Cohort %in% "PID" & Timepoint != "Acute")

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
res$padj <- ifelse(is.na(res$padj), 1, res$padj)
taxa_sig = rownames(res[1:35, ])
taxa_sig_check = res[(res$padj < alpha), ]
taxa_sig_check

sigtab = cbind(as(taxa_sig_check, "data.frame"), as(tax_table(ps.D.sub)[rownames(taxa_sig_check), ], "matrix"))



patient_linked_spec <- tax_glom(patient_linked, taxrank = "Species")
patient_linked_spec_melt <- psmelt((patient_linked_spec))

patient_linked_spec_melt$Vaccine_TimepointF<-factor(patient_linked_spec_melt$Vaccine_Timepoint, 
                                                   levels=c("1_Pre-Dose", "1_Acute", "1_Late", "2_Pre-Dose", "2_Acute", "2_Late", "3_Pre-Dose", "3_Acute", "3_Late"))

Supp_Fig3E <- ggplot(data = subset(patient_linked_spec_melt, Species %in% "s__ Enterobacter sp."), 
                 aes(x=Cohort, 
                     y=log(Abundance),
                     colour = TimepointF)) + 
  geom_boxplot() +
  geom_point(position=position_jitterdodge(jitter.width = .1), size = 4) +
  labs(y = "log(Relative Abundance)") +
  scale_colour_manual(values = timepoint_colours, name = "Timepoint") +
  theme(axis.text = element_text(size=16), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size=16),
        strip.text.x = element_text(size = 20)) 

ggsave("03_FIGURES/Supp_Fig3E.pdf", Supp_Fig3E, height = 8, width = 12.5, units = "in")

