# set a working directory
setwd("~/Desktop/COVID-19_Vaccination_GM")

# source R settings
source("01_SCRIPTS/Settings.R")

# read pre-processed files back in - contains alpha diversity measures 
meta_new <- readRDS("02_RESULTS/Modified_metadata.rds")

## Functional analysis
# Read in breakdowns of what each COG represents
EGGNOG_Groups <- read_excel("00_DATA/EGGNOG_functions.xlsx")

# read in the data analysis of the main functional groups
functional_maingroups <- read.delim("00_DATA/EGGNOG_MainFunctionalGroups.txt") 

# Starting with main_groups
functional_maingroups <- read.delim("00_DATA/EGGNOG_MainFunctionalGroups.txt") %>%
  remove_rownames() %>%
  column_to_rownames(var = "Datasets") %>%
  t() %>%
  as.data.frame() 

# dataframe manipulation for plotting
functional_maingroups <- functional_maingroups / rowSums(functional_maingroups)
functional_maingroups[is.na(functional_maingroups)] <- 0
functional_maingroups_meta <- merge(meta_new, functional_maingroups, by = "row.names")
functional_maingroups_meta <- functional_maingroups_meta %>% mutate(Sample = Row.names) %>% remove_rownames %>% column_to_rownames(var="Row.names")
functional_maingroups_meta$Sample <- as.character(functional_maingroups_meta$Sample)
functional_maingroups_meta_long <- functional_maingroups_meta %>% pivot_longer(cols="Information Storage and Processing":"Metabolism", 
                                                                      names_to='Functional_Main', 
                                                                      values_to='Functional_Abundance')

Fig4A <-  ggplot(functional_maingroups_meta_long,
               aes(x = Cohort, 
                   y = Functional_Abundance,
                   colour = TimepointF)) +
          geom_boxplot()+
          geom_point(position=position_jitterdodge(jitter.width = .1), size = 4) +
          facet_wrap(~Functional_Main, labeller = label_wrap_gen(), scales = 'free_y') +
          labs(y="Relative Abundance") + 
          scale_colour_manual(values = timepoint_colours, name = "Timepoint") +
          theme(axis.text = element_text(size=16), 
                axis.title.y = element_text(size=20),
                axis.title.x = element_blank(),
                legend.position = "bottom",
                legend.text = element_text(size = 16),
                legend.title = element_text(size=16),
                strip.text.x = element_text(size = 14)) #+
        #  geom_pwc(aes(group = Cohort),method = "wilcox_test", p.adjust.method = "fdr", tip.length = 0)
        #  adds space for p-val, padj wrong, adjust accordingly

Supp_Fig4A <-  ggplot(functional_maingroups_meta_long,
                 aes(x = Cohort, 
                     y = Functional_Abundance,
                     colour = Cohort)) +
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width = .1), size = 4) +
  facet_wrap(~Functional_Main, labeller = label_wrap_gen(), scales = 'free_y') +
  labs(y="Relative Abundance") + 
  scale_colour_manual(values = cohort_colours, name = "Cohort") +
  theme(axis.text = element_text(size=16), 
        axis.title.y = element_text(size=20),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size=16),
        strip.text.x = element_text(size = 14)) #+
#  geom_pwc(aes(group = Cohort),method = "wilcox_test", p.adjust.method = "fdr", tip.length = 0)
#  adds space for p-val, padj wrong, adjust accordingly

# save pdfs
ggsave("03_FIGURES/Figure4A.pdf", Fig4A)
ggsave("03_FIGURES/Supp_Fig4A.pdf", Supp_Fig4A)


# statistics 
mainfunc_stats_sub <- functional_maingroups_meta_long
p_values_mainfunc <- data.frame(Main_Func = character(),
                            Cohort = character(),
                            Timepoint1 = character(),
                            Timepoint2 = character(),
                            p_value = numeric())

# Get unique cohorts, timepoints and species
cohorts <- unique(mainfunc_stats_sub$Cohort)
timepoints <- unique(mainfunc_stats_sub$Timepoint)
mainfunc_list <- unique(mainfunc_stats_sub$Functional_Main)

# Loop through each species, cohort, and unique timepoint combinations within each cohort
for (func in mainfunc_list) {
  for (cohort in cohorts) {
    unique_timepoints <- unique(mainfunc_stats_sub$Timepoint[mainfunc_stats_sub$Cohort == cohort])
    for (i in 1:(length(unique_timepoints) - 1)) {
      for (j in (i + 1):length(unique_timepoints)) {
        timepoint1 <- unique_timepoints[i]
        timepoint2 <- unique_timepoints[j]
        
        # Create subsets for the two timepoints within the same cohort
        group1 <- subset(mainfunc_stats_sub, Functional_Main == func & Cohort == cohort & Timepoint == timepoint1)
        group2 <- subset(mainfunc_stats_sub, Functional_Main == func & Cohort == cohort & Timepoint == timepoint2)
        
        # Perform Wilcoxon signed-rank test
        p_value <- wilcox.test(group1$Functional_Abundance, group2$Functional_Abundance)$p.value
        
        # Store the results in the data frame
        result <- data.frame(Main_Func = func,
                             Cohort = cohort,
                             Timepoint1 = timepoint1,
                             Timepoint2 = timepoint2,
                             p_value = p_value)
        
        p_values_mainfunc <- rbind(p_values_mainfunc, result)
      }
    }
  }
}

p_values_mainfunc$FDR <- p.adjust(p_values_mainfunc$p_value, method = "fdr") 
saveRDS(p_values_mainfunc, "02_RESULTS/p_values_main-func.rds")


# read in the data analysis of the next level of functional grouping + manipulate for graphing
ENOG <- read.delim("00_DATA/EGGNOG_FunctionalGroups.txt")
ENOG <- ENOG %>% remove_rownames %>% column_to_rownames(var="Datasets_Datasets")
ENOG <-  ENOG %>% t() %>% as.data.frame()
ENOG <- ENOG/rowSums(ENOG)
ENOG[is.na(ENOG)] <- 0 

# merge metadata + functional data
ENOG_meta <- merge(meta_new, ENOG, by = "row.names")
colnames(ENOG_meta)[colnames(ENOG_meta) == "Row.names"] <- "Sample"
ENOG_meta$Sample <- as.character(ENOG_meta$Sample)

# pivot longer for graphing
ENOG_meta_long <- ENOG_meta %>% pivot_longer(cols="[A] RNA processing and modification":"[Q] Secondary metabolites biosynthesis, transport and catabolism", 
                                              names_to='Functional_Group', 
                                              values_to='Functional_Abundance')
ENOG_meta_long$Functional_Abundance <- as.numeric(ENOG_meta_long$Functional_Abundance)

## functional composition barplot
Fig4B <- ggplot(ENOG_meta_long,
                 aes(x = Sample, 
                     y = Functional_Abundance,
                     fill = Functional_Group)) +
            geom_bar(position="stack", stat="identity", width =1) +
            scale_fill_manual(values = getOI(22)) +
            facet_grid(~Cohort_TimepointF, scales='free_x', space = "free") + 
            labs(x = "Sample", y="Relative Abundance") + 
            theme(axis.text.x = element_blank(),
                  axis.ticks.x=element_blank(),
                  legend.position = "bottom")

ggsave("03_FIGURES/Figure4B.pdf", Fig4B, width = 75, height = 24, units = "cm")


# dot plot of the most abundant functions in each main category
Fig4C <- ggplot(subset(ENOG_meta_long, Functional_Group %in% c("[M] Cell wall/membrane/envelope biogenesis","[E] Amino acid transport and metabolism", "[L] Replication, recombination and repair")),
                 aes(x = Cohort, 
                     y = Functional_Abundance,
                     colour = TimepointF)) +
            geom_boxplot() +
            geom_point(position=position_jitterdodge(jitter.width = .01), size = 4) +
            facet_wrap(~Functional_Group, labeller = label_wrap_gen(), scales='free_y') + 
            labs(y="Relative Abundance") + 
            scale_colour_manual(values = timepoint_colours, name = "Timepoint") +
            theme(axis.text = element_text(size=16), 
                  axis.title.y = element_text(size=20),
                  axis.title.x = element_blank(),
                  legend.position = "bottom",
                  legend.text = element_text(size = 16),
                  legend.title = element_text(size=16),
                  strip.text.x = element_text(size = 14)) 

# save it
ggsave("03_FIGURES/Figure4C.pdf", Fig4C)


# design a facet to depict the three rows of main functional groups
design <- "DKLOPQRS
           CEFGHMN#
           ABIJ####"

# remaining functional groups (n=19) 
Supp_Fig4B <- ggplot(subset(ENOG_meta_long, !Functional_Group %in% c("[M] Cell wall/membrane/envelope biogenesis","[E] Amino acid transport and metabolism", "[L] Replication, recombination and repair")),
                           aes(x = Cohort, 
                               y = Functional_Abundance,
                               colour = TimepointF)) +
                geom_boxplot() +
                geom_point(position=position_jitterdodge(jitter.width = .01)) +
                facet_manual(~Functional_Group, design = design, labeller = label_wrap_gen(), scales='free_y') +
                scale_colour_manual(values = timepoint_colours, name = "Timepoint") +
                labs(y="Abundance") + 
                theme(strip.background = element_rect(),
                      axis.title.x = element_blank(),
                      legend.position = "bottom")
# save it
ggsave("03_FIGURES/Supp_Fig4B.pdf", Supp_Fig4B)

# statistics 
func_stats_sub <- ENOG_meta_long
p_values_func <- data.frame(Function = character(),
                            Cohort = character(),
                            Timepoint1 = character(),
                            Timepoint2 = character(),
                            p_value = numeric())

# Get unique cohorts, timepoints and species
cohorts <- unique(func_stats_sub$Cohort)
timepoints <- unique(func_stats_sub$Timepoint)
func_list <- unique(func_stats_sub$Functional_Group)

# Loop through each species, cohort, and unique timepoint combinations within each cohort
for (func in func_list) {
  for (cohort in cohorts) {
    unique_timepoints <- unique(func_stats_sub$Timepoint[func_stats_sub$Cohort == cohort])
    for (i in 1:(length(unique_timepoints) - 1)) {
      for (j in (i + 1):length(unique_timepoints)) {
        timepoint1 <- unique_timepoints[i]
        timepoint2 <- unique_timepoints[j]
        
        # Create subsets for the two timepoints within the same cohort
        group1 <- subset(func_stats_sub, Functional_Group == func & Cohort == cohort & Timepoint == timepoint1)
        group2 <- subset(func_stats_sub, Functional_Group == func & Cohort == cohort & Timepoint == timepoint2)
        
        # Perform Wilcoxon signed-rank test
        p_value <- wilcox.test(group1$Functional_Abundance, group2$Functional_Abundance)$p.value
        
        # Store the results in the data frame
        result <- data.frame(Func = func,
                             Cohort = cohort,
                             Timepoint1 = timepoint1,
                             Timepoint2 = timepoint2,
                             p_value = p_value)
        
        p_values_func <- rbind(p_values_func, result)
      }
    }
  }
}

p_values_func$FDR <- p.adjust(p_values_func$p_value, method = "fdr") 
saveRDS(p_values_func, "02_RESULTS/p_values_func.rds")


# lowest functional annotation (COG) level 
functional <- read.delim("00_DATA/EGGNOG_FunctionalCOGs.txt")
functional <- functional %>% remove_rownames %>% column_to_rownames(var="X.Datasets")
functional <-  functional %>% t() %>% as.data.frame()
functional <- functional/rowSums(functional)
functional[is.na(functional)] <- 0 

# merge with metadata
functional_meta <- merge(meta_new, functional, by = "row.names")
functional_meta <- functional_meta %>% mutate(Sample = Row.names) %>% remove_rownames %>% column_to_rownames(var="Row.names")

# pivot longer for statistics
COG <- functional_meta %>% pivot_longer("COG1949 3'-to-5' exoribonuclease specific for small oligoribonucleotides (By similarity)":"COG4910 Dehydratase, small subunit", 
                                      names_to='Functions', 
                                      values_to='Functional_Abundance')

# add meta of which preceding functional groups to each COG
COG <- merge(COG, EGGNOG_Groups)
COG$Functional_Abundance <- as.numeric(COG$Functional_Abundance)

# statistics 
COG_stats_sub <- COG
p_values_COG <- data.frame(COG = character(),
                            Cohort = character(),
                            Timepoint1 = character(),
                            Timepoint2 = character(),
                            p_value = numeric())

# Get unique cohorts, timepoints and species
cohorts <- unique(COG_stats_sub$Cohort)
timepoints <- unique(COG_stats_sub$Timepoint)
cog_list <- unique(COG_stats_sub$Functions)

# Loop through each species, cohort, and unique timepoint combinations within each cohort
for (cog in cog_list) {
  for (cohort in cohorts) {
    unique_timepoints <- unique(COG_stats_sub$Timepoint[COG_stats_sub$Cohort == cohort])
    for (i in 1:(length(unique_timepoints) - 1)) {
      for (j in (i + 1):length(unique_timepoints)) {
        timepoint1 <- unique_timepoints[i]
        timepoint2 <- unique_timepoints[j]
        
        # Create subsets for the two timepoints within the same cohort
        group1 <- subset(COG_stats_sub, Functions == cog & Cohort == cohort & Timepoint == timepoint1)
        group2 <- subset(COG_stats_sub, Functions == cog & Cohort == cohort & Timepoint == timepoint2)
        
        # Perform Wilcoxon signed-rank test
        p_value <- wilcox.test(group1$Functional_Abundance, group2$Functional_Abundance)$p.value
        
        # Store the results in the data frame
        result <- data.frame(Species = cog,
                             Cohort = cohort,
                             Timepoint1 = timepoint1,
                             Timepoint2 = timepoint2,
                             p_value = p_value)
        
        p_values_COG <- rbind(p_values_COG, result)
      }
    }
  }
}

p_values_COG$FDR <- p.adjust(p_values_COG$p_value, method = "fdr") 

saveRDS(p_values_spec, "02_RESULTS/p_values_species.rds")