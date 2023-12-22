# set a working directory
setwd("~/Desktop/COVID-19_Vaccination_GM")

# source R settings
source("01_SCRIPTS/00_Settings.R")

# read pre-processed files back in
meta_new <- readRDS("02_RESULTS/Modified_metadata.rds")

## Assess the alpha diversity of the samples ##
# Make df long for eventual plotting in facets
meta_new_long <- pivot_longer(meta_new, cols = c("chao1", "diversity_shannon"),
                              names_to = "Diversity_Measure",
                              values_to = "Value")

#The above makes values = characters, so move them back to numeric
meta_new_long$Value <- as.numeric(meta_new_long$Value)

#  ggplot to produce results + create figures
Supp_Fig2A <- ggplot(subset(meta_new_long),
                     aes(x = Cohort, 
                         y = Value,
                         colour = Cohort,
                         na.rm = T)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = .1), size = 4) +
  facet_wrap(~Diversity_Measure, scales = 'free_y', ncol = 2) +
  labs(x= "Cohort", y= "Diversity Measure") +
  scale_colour_manual(values = cohort_colours, name = "Cohort") +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 20)) + 
  geom_pwc(aes(group = Cohort), method = "wilcox_test", p.adjust.method = "bonferroni", tip.length = 0,
           p.adjust.by = c("group"), 
           label = "p.adj.format") 


Fig1B <- ggplot(meta_new_long, 
                aes(x = Cohort, 
                    y = Value,
                    colour = TimepointF,
                    na.rm = T)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0), size = 4) +
  facet_wrap(~Diversity_Measure, scales = 'free_y') +
  labs(x= "Cohort", y= "Diversity Measure") +
  scale_colour_manual(values = timepoint_colours, name = "Timepoint") +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 20)) + 
  geom_pwc(aes(group = Timepoint), method = "wilcox_test", p.adjust.method = "bonferroni", tip.length = 0,
           p.adjust.by = c("group"), 
           label = "p.adj.format") 


# save pdfs
ggsave("03_FIGURES/Supp_Fig2A.pdf", Supp_Fig2A, height = 4.5, width = 8, units = "in")
ggsave("03_FIGURES/Figure1B.pdf", Fig1B, height = 4.5, width = 8, units = "in")
