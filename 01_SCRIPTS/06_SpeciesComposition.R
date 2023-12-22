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

# merge on species for analysis
ps_main_composition_spec = tax_glom(ps_main_composition, taxrank="Species") 


## T15 species in each cohort - see effect of timepoints
# take OTU out of phyloseq object 
cohort_sub <- subset_samples(ps_main_composition_spec, Cohort %in% "ICP")
sub_otu_table <- as.data.frame(as.matrix(cohort_sub@otu_table))

#makes relative abundance between samples
df_rabs <- data.frame(mOTU_ID = row.names(sub_otu_table),
                      RA = (rowSums(sub_otu_table))/(ncol(sub_otu_table)))

## takes present Otu, gets rid of what's not there, then ranks them based on RA
df_rabs_trim <- df_rabs %>% filter(RA > 0)
df_rabs_trim$Rank <- rank(-(df_rabs_trim$RA))

# Take data of average abundance of each taxa within all samples, ranks, filters the top 15
t15rabs <- df_rabs_trim %>% arrange(Rank) %>% filter(Rank %in% c(1:15))

# add species name from Tax_table_new
tax_table_temp <- as.data.frame(Tax_table_new)

add_name_df <- data.frame(mOTU_ID = rownames(tax_table_temp),
                          Species = tax_table_temp$Species)
add_name_df_sub <- add_name_df[add_name_df$mOTU_ID %in% rownames(t15rabs), ]
t15rabs <- merge(t15rabs, add_name_df)
t15rabs <- t15rabs %>% arrange(Rank) 
otu_lst15 <- t15rabs$mOTU_ID

# gather meta data to merge
sub_meta <- as.data.frame(as.matrix(cohort_sub@sam_data))
sub_meta$Sample_ID <- rownames(sub_meta)

top15rabs <- merge(t15rabs, sub_meta)

# creates plot of top 15 species abundance, black line
p1 <- ggplot(data = top15rabs,
             aes(x=Species, 
                 y=(log(RA+1e-5)),
                 group = TimepointF)) + geom_point() + geom_line() +
  scale_x_discrete(limits = t15rabs$Species) + 
  scale_y_continuous(limits = c(-12,0)) +
  labs(x = "Species", y = "log(Relative Abundance)") + 
  theme_minimal_hgrid(color = "gray88", line_size = 0) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 8))

# new Otu_table for top20
top15_otu <- sub_otu_table[rownames(sub_otu_table) %in% t15rabs$mOTU_ID, ] %>% as.data.frame()
top15_otu$mOTU_ID <- rownames(top15_otu)

# merge documents to get species names added
top15_otu_long <-  pivot_longer(top15_otu, cols = !mOTU_ID, names_to = "Sample_ID", values_to = "RA")
final <- merge(top15_otu_long, sub_meta)
final <- merge(final, add_name_df)

# factor Timepoints
final$TimepointF <- factor(final$Timepoint, levels=c("Pre-Dose",
                                                     "Acute", 
                                                     "Late"))
# makes scatterplot for each species at each timepoint
p2 <- ggplot(final,
             aes(x = Species, 
                 y = log(RA+1e-5),
                 colour = TimepointF,
                 na.rm = T)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = .02), size = 0.8) +
  scale_colour_manual(values = timepoint_colours, name = "Timepoint") +
  scale_x_discrete(limits = t15rabs$Species) + 
  scale_y_continuous(limits = c(-12, 0),position = "right") +
  theme_minimal_hgrid(color = "gray88", line_size = 0.3) + 
  theme(axis.text.x = element_blank(),
        axis.title = element_blank())

aligned_plots <- align_plots(p1, p2, align="hv", axis="tblr")
Fig2E <- ggdraw(aligned_plots[[2]]) + draw_plot(aligned_plots[[1]])

ggsave("03_FIGURES/Figure2E.pdf", Fig2E)

## Abundances of species level in scatter plots
# melt the phyloseq = pivot longer on species
ps_main_full <-psmelt(ps_main_composition) 

# Look at specific species 
Fig2C <- ggplot(data = subset(ps_main_full, Cohort %in% "ICP" & Species %in% "s__Klebsiella pneumoniae"),
                aes(x=Cohort, 
                    y=log(Abundance),
                    colour = TimepointF)) + 
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

# Look at specific species 
Fig2D <- ggplot(data = subset(ps_main_full, Cohort %in% "ICP" & Species %in% "s__Butyrivibrio crossotus"),
                aes(x=Cohort, 
                    y=log(Abundance),
                    colour = TimepointF)) + 
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


Fig2F <-  ggplot(data = subset(ps_main_full, Species %in% "s__Faecalibacterium prausnitzii"), 
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
        strip.text.x = element_text(size = 20)) +
  geom_pwc(aes(group = Timepoint), method = "wilcox_test", p.adjust.method = "bonferroni", tip.length = 0,
           p.adjust.by = c("group"), 
           label = "p.adj.format") 

Fig2G <-  ggplot(data = subset(ps_main_full, Species %in% "s__Akkermansia muciniphila"), 
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
        strip.text.x = element_text(size = 20)) +
  geom_pwc(aes(group = Timepoint), method = "wilcox_test", p.adjust.method = "bonferroni", tip.length = 0,
           p.adjust.by = c("group"), 
           label = "p.adj.format") 
  

Fig2H <-  ggplot(data = subset(ps_main_full, Species %in% "s__Escherichia coli"), 
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
        strip.text.x = element_text(size = 20)) +
  geom_pwc(aes(group = Timepoint), method = "wilcox_test", p.adjust.method = "bonferroni", tip.length = 0,
           p.adjust.by = c("group"), 
           label = "p.adj.format") 


# save figures
ggsave("03_FIGURES/Figure2C.pdf", Fig2C)
ggsave("03_FIGURES/Figure2D.pdf", Fig2D)
ggsave("03_FIGURES/Figure2E.pdf", Fig2E)
ggsave("03_FIGURES/Figure2F.pdf", Fig2F)
ggsave("03_FIGURES/Figure2G.pdf", Fig2G)
ggsave("03_FIGURES/Figure2H.pdf", Fig2H)


# statistics 
stats_sub <- ps_main_full
p_values_spec <- data.frame(Species = character(),
                       Cohort = character(),
                       Timepoint1 = character(),
                       Timepoint2 = character(),
                       p_value = numeric())

# Get unique cohorts, timepoints and species
cohorts <- unique(stats_sub$Cohort)
timepoints <- unique(stats_sub$Timepoint)
spec_list <- unique(stats_sub$Species)

# Loop through each species, cohort, and unique timepoint combinations within each cohort
for (species in spec_list) {
  for (cohort in cohorts) {
    unique_timepoints <- unique(stats_sub$Timepoint[stats_sub$Cohort == cohort])
    for (i in 1:(length(unique_timepoints) - 1)) {
      for (j in (i + 1):length(unique_timepoints)) {
        timepoint1 <- unique_timepoints[i]
        timepoint2 <- unique_timepoints[j]
        
        # Create subsets for the two timepoints within the same cohort
        group1 <- subset(stats_sub, Species == species & Cohort == cohort & Timepoint == timepoint1)
        group2 <- subset(stats_sub, Species == species & Cohort == cohort & Timepoint == timepoint2)
        
        # Perform Wilcoxon signed-rank test
        p_value <- wilcox.test(group1$Abundance, group2$Abundance)$p.value
        
        # Store the results in the data frame
        result <- data.frame(Species = species,
                             Cohort = cohort,
                             Timepoint1 = timepoint1,
                             Timepoint2 = timepoint2,
                             p_value = p_value)
        
        p_values_spec <- rbind(p_values_spec, result)
      }
    }
  }
}

p_values_spec$padj <- p.adjust(p_values_spec$p_value, method = "bonferroni") 

saveRDS(p_values_spec, "02_RESULTS/p_values_species.rds")

