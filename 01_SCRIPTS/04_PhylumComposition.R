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

# Phylum composition analysis 
# merge the OTUs - the phyloseq function tax_glom merges the OTUs with the same taxonomy, summing the abundances:
ps_main_composition_phylum = tax_glom(ps_main_composition, taxrank="Phylum") 

# Composition barplot 
# the relative abundance of each phylum in each sample
Fig1D <- plot_bar(ps_main_composition_phylum, fill="Phylum") + 
            facet_grid(~Cohort_TimepointF,  scales="free", space = "free", labeller = label_wrap_gen()) +  
            scale_fill_manual(values = getOI(18)) +
            labs(y = "Relative Abundance") + 
            theme(axis.text.x = element_blank(),
                  legend.position = "bottom") 

# save the final figure - plots stitched together in Adobe Illustrator
ggsave("03_FIGURES/Figure1D.pdf", Fig1D, width = 75, height = 24, units = "cm")


# Scatterplots 
# melt the phyloseq so that you can plot each phylum in separate facets
ps_main_full_phylum <-psmelt(ps_main_composition_phylum) 

# plot first for samples from each cohort 
Supp_Fig1D <- ggplot(data = subset(ps_main_full_phylum, Phylum %in%c("p__Firmicutes", 
                                                                     "p__Bacteroidetes", 
                                                                     "p__Proteobacteria", 
                                                                     "p__Actinobacteria", 
                                                                     "p__Verrucomicrobia",
                                                                     "p__Euryarchaeota")),
                     aes(x=Cohort, 
                         y=Abundance,
                         colour = Cohort)) +
  geom_boxplot() +
  geom_point(position=position_jitterdodge(jitter.width = .1), size = 4) +
  facet_wrap(~Phylum, ncol = 3, scales = "free", labeller = label_wrap_gen(2)) +
  theme(axis.text = element_text(size=16), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size=16),
        strip.text.x = element_text(size = 20)) +
  labs(y = "Relative Abundance") +
  scale_colour_manual(values = cohort_colours, name = "Cohort") #+
# geom_pwc(aes(group = Cohort),method = "wilcox_test", p.adjust.method = "fdr", tip.length = 0)
# adds space for p-val, padj wrong, adjust accordingly

# plot then for timepoints within each cohort
Fig1E <- ggplot(data = subset(ps_main_full_phylum, Phylum %in%c("p__Firmicutes", 
                                                       "p__Bacteroidetes", 
                                                       "p__Proteobacteria", 
                                                       "p__Actinobacteria", 
                                                       "p__Verrucomicrobia",
                                                       "p__Euryarchaeota")),
       aes(x=Cohort, 
           y=Abundance,
           colour = TimepointF)) +
  geom_boxplot() +
  geom_point(position=position_jitterdodge(jitter.width = .1), size = 4) +
  facet_wrap(~Phylum, ncol = 4, scales = "free", labeller = label_wrap_gen(2)) +
  theme(axis.text = element_text(size=16), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size=16),
        strip.text.x = element_text(size = 20)) +
  labs(y = "Relative Abundance") +
  scale_colour_manual(values = timepoint_colours, name = "Timepoint") #+
  # geom_pwc(aes(group = Timepoint),method = "wilcox_test", p.adjust.method = "fdr", tip.length = 0)
  # adds space for p-val, padj wrong, adjust accordingly

# save the final figures
ggsave("03_FIGURES/Supp_Fig1D.pdf", Supp_Fig1D)
ggsave("03_FIGURES/Figure1E.pdf", Fig1E)


# Statistical testing
# Cohort_Timepoint is first assessed, change to just Cohort after 
stats_sub <- ps_main_full_phylum
grp <- unique(stats_sub$Cohort_Timepoint)
phy_list <- unique(stats_sub$Phylum)

# make new file for the results of the stats test to go into
p_values <- c()
passed_tests <- c()

# loop the stats test 
for (k in phy_list){
  sub <- subset(stats_sub, Phylum == k)
  for (i in grp) {
    this_i = sub[sub$Cohort_Timepoint == i,]
    for (j in grp) {
      if (i != j & length(grep(paste(j, i, sep="__"), passed_tests)) == 0){
        this_j = sub[sub$Cohort_Timepoint == j,]
        this_p <- wilcox.test(this_i$Abundance, this_j$Abundance, exact = F)$p.value
        
        this_output <- c(k, i, j, this_p)
        p_values <- rbind(p_values, this_output)
        passed_tests[length(passed_tests)+1] <- paste(i, j, sep = "__")
      }
    }
  }
}
p_values <- as.data.frame(p_values)
colnames(p_values) <- c("Function", "Group_1", "Group_2", "P_value")
p_values$FDR <- p.adjust(p_values$P_value, method = "fdr") 
p_values$P_value <- as.numeric(p_values$P_value)

# final table of statistics = p_values + save
saveRDS(p_values, "02_RESULTS/p_values_phylum.rds")

# write out data
# Fig1DE_SuppFig1D_data <- ps_main_full_phylum[, c(1:7, 39,40)]
# saveRDS(Fig1DE_SuppFig1D_data, "00_DATA/Fig1DE_SuppFig1D_data.rds")