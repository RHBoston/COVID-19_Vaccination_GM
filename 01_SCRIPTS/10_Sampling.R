# load in package to plot
library(ggplot2)

# set a working directory 
setwd("~/Desktop/COVID-19_Vaccination_GM")

# read in the excel file which details how each participant sampled
df <- read_excel("00_DATA/actual_sampling.xlsx")

# create a empty new dataframe that will be utilized 
df2 <- data.frame(PatientID = character(), Cohort = character(), Timepoint = character(), Value = numeric())

#take out variables for orders of the x (sample timepoints) and y (patient order) axis in eventual plot 
col_order <- df$PatientID
timepoint_order <- colnames(df[,3:11])

# loop through the read in excel file dataframe (df) to create a new long dataframe (df2) detailing only the timepoints in which each patient sampled
for (i in 1:length(rownames(df))){
  for (j in 4:length(colnames(df))){
    if (df[i,j] == 1){
      df2 <- rbind(df2,
                   data.frame(PatientID = df[i,2], Cohort = df[i,1], Timepoint = colnames(df)[j], Value = 1))
    }
  }
}

# set the patient order as a factor for graphing
df2$PatientID <- factor(df2$PatientID, levels = rev(col_order))

# plot 
sampling <- ggplot(df2, aes(Timepoint, 
                            PatientID, 
                            group = PatientID, 
                            colour = Cohort)) + 
              geom_line() + geom_point(size = 4) + 
              scale_x_discrete(limits = timepoint_order) 

# results in a abacus plot detailing each sample from each participant in the study, colours depict cohorts.
ggsave("03_FIGURES/Sampling.pdf", sampling)
