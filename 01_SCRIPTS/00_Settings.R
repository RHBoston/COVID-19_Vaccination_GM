# load in required packages 
library("readxl")
library("writexl")
library("tidyverse")
library("ggplot2")
library("ggpubr") 
library("phyloseq") 
library("microbiome") 
library("vegan")
library("factoextra") 
library("lme4") 
library("flexplot") 
library("RColorBrewer") 
library("gridExtra")
library("DESeq2")
library("ComplexHeatmap") 
library("cowplot") 
library("ggh4x") 
library("rstatix") 
library("EnhancedVolcano")

# define a colour palette
palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                               "#0072B2", "#D55E00", "#CC79A7", "#999999")
# create a function to define number of colours requested using the defined palette as a template
getOI = colorRampPalette(palette_OkabeIto)

# define colours for specific aspects of data analysis                                
cohort_colours <- c("#F8766D",  "#00BA38", "#619CFF")
timepoint_colours <- c("#E69F00", "#CC79A7", "#009E73")

# set theme               
theme_set(theme_bw())
                               