# set a working directory
setwd("~/Desktop/COVID-19_Vaccination_GM")

# source R settings
source("01_SCRIPTS/00_Settings.R")

# read in the patient_linked df
patient_linked <- readRDS("02_RESULTS/Patient_linked_ps.rds")
patient_linked_df <- sample_data(patient_linked) %>% as.matrix() %>% as.data.frame() %>% rownames_to_column(var = "Sample")

# formatting for loops to acquire sample numbers for power calculations
vaccine <-unique(patient_linked_df$Vaccine)
cohorts <- unique(patient_linked_df$Cohort)
timepoint_pairs <- unique(patient_linked_df$Timepoint)

df_list <- data.frame()
final_power_calc <- data.frame()

for (v in vaccine) {
  linked_sub <- subset(patient_linked_df, Vaccine == v)
  
  for (tp in timepoint_pairs) {
  timepoint_sub <- subset(linked_sub, Timepoint != tp)
  
    for (c in cohorts) {
    sub_df <-  subset(timepoint_sub, Cohort == c)
    sample_counts <- table(sub_df$Patient_Vaccine)
    final_def <- subset(sub_df, Patient_Vaccine %in% names(sample_counts[sample_counts == 2]))
    
    pairs <- length(unique((final_def$Patient_Vaccine)))
    timepoint_comp <- timepoint_pairs[timepoint_pairs != tp]
    
    df_list <- rbind(df_list, data.frame(Vaccine = v, Timepoints = paste0(timepoint_comp[1], sep = '_', timepoint_comp[2]), Cohort = c, Unique_pairs =pairs))
    }
  }
}

## while final_power_calc is a subsetted df of paired sample counts by which each cohort comparision of vaccine timepoints in total

final_power_calc <- df_list %>%
  group_by(Timepoints) %>%
  summarize(Unique_pairs = sum(Unique_pairs))

