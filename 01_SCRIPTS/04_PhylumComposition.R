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
Supp_Fig2D <- ggplot(data = subset(ps_main_full_phylum, Phylum %in%c("p__Firmicutes", 
                                                                     "p__Bacteroidetes", 
                                                                     "p__Proteobacteria", 
                                                                     "p__Actinobacteria", 
                                                                     "p__Verrucomicrobia",
                                                                     "p__Euryarchaeota")),
                     aes(x=Cohort, 
                         y=Abundance,
                         colour = Cohort)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = .1), size = 4) +
  facet_wrap(~Phylum, ncol = 3, scales = "free", labeller = label_wrap_gen(2)) +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 20)) + 
  labs(y = "Relative Abundance") +
  scale_colour_manual(values = cohort_colours, name = "Cohort") +
  geom_pwc(aes(group = Cohort), method = "wilcox_test", p.adjust.method = "bonferroni", tip.length = 0,
           p.adjust.by = c("group"), 
           label = "p.adj.format") 
# adds space for p-val, padj wrong, adjust accordingly

ggsave("03_FIGURES/Supp_Fig2D.pdf", Supp_Fig2D, height = 7.5, width = 12.5, units = "in")

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
  geom_boxplot(outlier.shape = NA) +
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
  scale_colour_manual(values = timepoint_colours, name = "Timepoint") +
  geom_pwc(aes(group = Timepoint), method = "wilcox_test", p.adjust.method = "bonferroni", tip.length = 0,
           p.adjust.by = c("group"), 
           label = "p.adj.format") 

# save the final figures
ggsave("03_FIGURES/Supp_Fig2D.pdf", Supp_Fig2D)
ggsave("03_FIGURES/Figure1E.pdf", Fig1E)

