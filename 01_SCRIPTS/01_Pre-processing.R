# set a working directory 
setwd("~/Desktop/COVID-19_Vaccination_GM")

# source R settings
source("01_SCRIPTS/00_Settings.R")

# read in data from hpc
merged <- read.delim("00_DATA/merged_final.txt")

# extract tax table from the merged document
Tax_table <- unique(merged[,1:8])
rownames(Tax_table) <- Tax_table$mOTU_ID
Tax_table <- Tax_table[, -1]

# make species column more readable
Tax_table_full <- Tax_table %>% as.data.frame() %>% 
  separate(Species, c("Species", "Species_Other"), sep= "\\s+\\[" , 
           extra = "drop", fill = "right", remove=TRUE)
Tax_table_new <- as.matrix(Tax_table_full[1:7])
rm(Tax_table_full, Tax_table)
   
# extract OTU table from the merged document
Otu_table.temp <- merged[, c("mOTU_ID", "Sample", "RA")]

# pivot the table
Otu_table<- Otu_table.temp %>% pivot_wider(names_from ='Sample', 
                                           values_from ='RA') 
# remove temp file
rm(Otu_table.temp)

# remove the column of OTU things and make it rownames
Otu_table <- Otu_table %>% remove_rownames %>% column_to_rownames(var="mOTU_ID")
Otu_table <- as.data.frame(as.matrix(Otu_table))

# change the NAs in the OTU table to 0
Otu_table[is.na(Otu_table)] <- 0
Otu_table <- apply(Otu_table, 2, function(x) x/sum(x))

# read in metadata 
Metadata <- read_excel("00_DATA/metadata_final.xlsx")
Metadata <- Metadata %>% remove_rownames %>% column_to_rownames(var="Sample_ID")

# add factors to order graphics later
Metadata$TimepointF<-factor(Metadata$Timepoint, 
                            levels=c("Pre-Dose", "Acute", "Late"))


# concat cohort and timepoint into a single variable 
Metadata$Cohort_Timepoint<-paste0(Metadata$Cohort, "_", Metadata$Timepoint)
Metadata$Cohort_Vaccine<-paste0(Metadata$Cohort, "_", Metadata$Vaccine)
Metadata$Vaccine_Timepoint<-paste0(Metadata$Vaccine, "_", Metadata$Timepoint)
Metadata$Vaccine_Cohort<-paste0(Metadata$Vaccine, "_", Metadata$Cohort)

# Factor based on a new variable
Metadata$Cohort_TimepointF<-factor(Metadata$Cohort_Timepoint, 
                                   levels=c("HC_Pre-Dose", "HC_Acute", 
                                            "HC_Late", "ICP_Pre-Dose", 
                                            "ICP_Acute", "ICP_Late", 
                                            "PID_Pre-Dose", "PID_Acute",  
                                            "PID_Late"))

Metadata <- Metadata %>% mutate(Health = case_when(Cohort %in% c("ICP","PID") 
                                                   ~ "Immunocompromised",
                                                   Cohort %in% c("HC") 
                                                   ~ "Healthy Controls"))

Metadata <- Metadata %>% mutate(NTresponse = case_when(v2D21 > 2481 ~ "Response",
                                                       v2D21 <= 1480 ~ "Low response",
                                                       v2D21 <= 5 ~ "No response",
                                                  is.na(v2D21) ~ "No data"))       

# Create Phyloseq object
ps <- phyloseq(otu_table(Otu_table, taxa_are_rows = TRUE), 
               sample_data(Metadata), 
               tax_table(Tax_table_new))

# convert the class of columns 
sample_data(ps)$Vaccine <-  as.character(sample_data(ps)$Vaccine)
sample_data(ps)$v3D28 <- as.numeric(sample_data(ps)$v3D28) 
sample_data(ps)$v2D21 <- as.numeric(sample_data(ps)$v2D21) 
sample_data(ps)$Spike <- as.numeric(sample_data(ps)$Spike) 
sample_data(ps)$Age <- as.numeric(sample_data(ps)$Age) 


# Alpha Diversity
alpha_div <- alpha(ps, index = c('shannon', 'chao1'))
new_Metadata <-  merge(sample_data(ps), alpha_div, by=0, all.x=TRUE)
new_Metadata <- new_Metadata %>% remove_rownames() %>% column_to_rownames(var="Row.names")

# add new Metadata back into ps
sample_data(ps) = new_Metadata
rm(Metadata) # doesn't have alpha diversity measures

# save modified files for processing
saveRDS(new_Metadata, "02_RESULTS/Modified_metadata.rds")
saveRDS(Otu_table, "02_RESULTS/Modified_otu_table.rds")
saveRDS(Tax_table_new, "02_RESULTS/Modified_tax_table.rds")
saveRDS(ps, "02_RESULTS/Phyloseq.rds")


