# set a working directory
setwd("~/Desktop/COVID-19_Vaccination_GM")

# source R settings
source("01_SCRIPTS/00_Settings.R")

# read pre-processed files back in - contains alpha diversity measures 
meta_new <- readRDS("02_RESULTS/Modified_metadata.rds")

# How vaccine efficacy (neutralisation of live virus (NT) & spike-specific antibodies) correlates with diversity
### NT @ 2nd vs Shannon
Fig3A <- ggplot(subset(meta_new, Vaccine == "2"),
                    aes(x = log(v2D21),
                        y = diversity_shannon,
                        colour = Cohort,
                        na.rm = TRUE)) +
              geom_point(size =4) +
              labs(x = "log(NT) at v2D21", y = "Shannon Diversity") +
              facet_wrap(~ TimepointF) +
              scale_x_continuous() +
              theme(axis.text = element_text(size=16), 
                    axis.title = element_text(size=20),
                    legend.position = "bottom",
                    legend.text = element_text(size = 16),
                    legend.title = element_text(size=16),
                    strip.text.x = element_text(size = 16)) +
              stat_cor(method = "spearman", cor.coef.name = "rho") +
              geom_smooth(method = "lm", se = FALSE)


### NT @ 3rd vs Shannon
Fig3B <- ggplot(subset(meta_new, Vaccine == "3"),
               aes(x = log(v3D28),
                   y = diversity_shannon,
                   colour = Cohort,
                   na.rm = TRUE)) +
          geom_point(size =4) +
          labs(x = "log(NT) at v3D28", y = "Shannon Diversity") +
          facet_wrap(~ TimepointF) +
          scale_x_continuous() +
          theme(axis.text = element_text(size=16), 
                axis.title = element_text(size=20),
                legend.position = "bottom",
                legend.text = element_text(size = 16),
                legend.title = element_text(size=16),
                strip.text.x = element_text(size = 16)) +
        stat_cor(method = "spearman", cor.coef.name = "rho") +
          geom_smooth(method = "lm", se = FALSE)


# anti-spike IgG antibodies @ 2nd vs Shannon
Supp_Fig4 <- ggplot(subset(meta_new, Vaccine == "2"),
                   aes(x = log(Spike),
                       y = diversity_shannon,
                       colour = Cohort,
                       na.rm = TRUE)) +
              geom_point(size =4) +
              labs(x = "log(anti-spike IgG antibodies at v2D21)", y = "Shannon Diversity") +
              facet_wrap(~ TimepointF) +
              scale_x_continuous() +
              theme(axis.text = element_text(size=16), 
                    axis.title = element_text(size=20),
                    legend.position = "bottom",
                    legend.text = element_text(size = 16),
                    legend.title = element_text(size=16),
                    strip.text.x = element_text(size = 16)) +
              stat_cor(method = "spearman", cor.coef.name = "rho") +
              geom_smooth(method = "lm", se = FALSE)


# save pdfs
ggsave("03_FIGURES/Figure3A.pdf", Fig3A)
ggsave("03_FIGURES/Figure3B.pdf", Fig3B)
ggsave("03_FIGURES/Supp_Fig4.pdf", Supp_Fig4)
