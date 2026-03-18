# This Script -------------------------------------------------------------
  # creates dataframes from Columbus Instrument + AutoTrack open-field locomotion sessions and project metadata, then pre-processes, plots, and analysis the data
  # used here for 15mo old Dat-cre Neuromelanin R01 mice
# Load Libraries
# Create Raw Dataframe
# Pre-Processing & Cleaning
  # drop files with recording errors (z-axis not working in some chambers or any other issue)
# Summary Dataframes
# all summary types are ran for virus/timepoint and again for virus/sex/timepoint
  # drop missed injections or other animal specific issues to remove them before analysis
  # Time In Zone (TIZ) summary
  # Locomotion summaries (resting, distance, rearing, etc)
# Analysis 
# Load Libraries ---------------------------------------------------------------
pacman::p_load('dplyr', 'tidyverse', 'ggplot2', 'here', 'purrr', 'readr', 'viridis', 'scatterplot3d', 'plotly','lme4','emmeans', 'writexl') # auto installs required packages not yet installed
source(here('Loco_functions.R'))

here::i_am("AutoTrack_openField.R")


# to collect pack citations, uncomment and run line below:
# cite_packages()

# to collect pack versions, uncomment and run line below:
# report_packages()

############################### Create Raw Dataframe ##########################################

# open files by timepoint - make sure these order animal in the same way the animal vector is set up
file_paths <- list(
  timepoint_0 = list.files(path = here("Locomotion", "wk_0", "Bins"), full.names= TRUE),
  timepoint_1 = list.files(path = here("Locomotion", "wk_1", "Bins"), full.names= TRUE),
  timepoint_2 = list.files(path = here("Locomotion", "wk_2", "Bins"), full.names= TRUE)
)
#print(file_paths)

# create metadata vectors
animal <- c('1CM', '2AM', '3AF', '4CF','5BM', '6BF') 
sex <- c( 'M', 'M', 'F', 'F', 'M', 'F')  
treatment <- c('C','A','A','C','B','B') 
timepoints <- c(0,1,2)


# Set an empty list to store dataframes
loco_list <-list()

expected_cols <- c(
  "Bin #", "Start", "Finish", "Records", "Zone", "Entrances",
  "Ambulatory (s)", "Stereotypic (s)", "Resting (s)", "Undetected (s)",
  "Rearing Events (Z-Axis)", "Rearing Events (V-Axis)",
  "Rearing Time Z (s)", "Rearing Time V (s)",
  "Distance (cm)", "Average Speed (cm/s)",
  "animal", "sex", "treatment", "timepoint")


# Loop through timepoints and files: outer loop iterates over timepoints
for (tp in seq_along(file_paths)) {
  #get file paths for current timepoint
  files <- file_paths[[tp]]
  
  # Loop through the files for this timepoint and skip first 21 lines: inner loop processes each file for the current timepoint
  for (i in seq_along(files)) {
    # read csv files, skipping the first 21 lines
    loco <- readr::read_csv(files[i], skip = 21, col_names = TRUE) # header = TRUE indicates line 22 is column names
    
    # keep only expected cols
    loco <- loco %>% select(any_of(expected_cols))
    
    # add metadata columns
    loco$animal <- animal[i]
    loco$sex <- sex[i]
    loco$treatment <- treatment[i]
    loco$timepoint <- timepoints[tp]
    
    # Append to the list
    loco_list <-append(loco_list, list(loco))
  }
}

# check columns before rbind
all_cols <- unique(unlist(lapply(loco_list, colnames)))
print(all_cols)

# combing all data to a dataframe
loco <- do.call(rbind, loco_list)

# Rearrange column names so metadata is first ###
loco <- select(loco,"animal", "sex", "treatment", "timepoint", everything())

# clean column names to lowercase and _
loco <- clean_names(loco)
#View(loco)




########################### Pre-processing & Cleaning #################################################################
## Use F-N baselineSums to collect timepoint == 0 summaries : used to counter-balance viral groups pre-injection
# pre_inj_totals <- baselineSums(
# df = loco, 
# group = "animal", 
# timepoint = "timepoint",
# bin_col = "bin_number", bin_range = 1:12,
# "distance_m"
#)
#write.csv(pre_inj_totals, file = here("Locomotion", "Dataframes", "preInjection_summary.csv"))

## Use F-N dropErrors to make corrections to recording errors. Un-comment below to run example. 
# loco <- dropErrors(loco, "animal", "1CM", "rearing_events_z_axis")  # add / adjust variables as needed


# convert distance in cm to m
loco <- loco %>%
  mutate(distance_cm = distance_cm / 100) %>% # convert
  rename(distance_m = distance_cm) # rename column

#View(loco)

write.csv(loco,file = here("Locomotion","Dataframes", "raw_locomotion.csv"))


##### Summary Dataframes ####
loco <- read.csv(file = here("Locomotion","Dataframes", "raw_locomotion.csv"))

# F-N animalSums will sum over full session (60min), while the added bin_range = 1:6 will sum first 30min, so either can be used in analysis and graphing
  # Similarly, F-N timeInZones will be used to calculate time spent in the center vs edge of the open field chamber

#1) Summarize time in zones (TIZ) (ambulatory, stereotypic, and resting) full session 
# by animal
TIZ_animal <- timeInZones(loco, bin_col = "bin_number", bin_range = 1:6, animal, treatment, timepoint, zone)
TIZ_animal <- TIZ_animal %>% filter(zone != 'None') # filter out 'none' observations
write.csv(TIZ_animal,file = here("Locomotion","Dataframes","timeInZone_animal.csv")) 

# by treatment
TIZ <- timeInZones(loco, bin_col= "bin_number", bin_range = 1:6, treatment, timepoint, zone)
TIZ <- TIZ %>% filter(zone != 'None') 
write.csv(TIZ,file = here("Locomotion","Dataframes","timeInZone_summary.csv")) 

# by sex and treatment
TIZ_2 <- timeInZones(loco, bin_col= "bin_number", bin_range = 1:6, treatment, sex, timepoint, zone)
TIZ_2 <- TIZ_2 %>% filter(zone != 'None')
write.csv(TIZ_2,file = here("Locomotion","Dataframes","timeInZone_sex_summary.csv"))


#2) Time in center (TIC) & Time in edge summaries (TIE)
TIC_animal <- TIZ_animal %>% filter(zone == "Center") #reduce to time in center
write.csv(TIC_animal,file = here("Locomotion","Dataframes","timeInCenter_byAnimal.csv"))

TIE_animal <- TIZ_animal %>% filter(zone == "Edge") #reduce to time in edge
write.csv(TIE_animal,file = here("Locomotion","Dataframes","timeInEdge_byAnimal.csv"))

#3) Get sum of all variable per animal (average speed handled outside of function)
# automate with variable vector
sum_variables <- c("ambulatory_s", "stereotypic_s", "resting_s", "rearing_events_z_axis", "rearing_time_z_s", "distance_m")

# Apply animalSums FN
summed <- animalSums(
  df = loco, 
  zone_col = "zone", 
  zones = c("Center", "Edge"), 
  bin_col = "bin_number", 
  bin_range = 1:6, 
  variables = sum_variables,
  animal, treatment, timepoint
)

# Handle speed separately
speed_summed <- loco %>%
  group_by(animal, treatment, timepoint, sex) %>%
  summarise(
    avg_speed_total = sum(`distance_m`[zone %in% c("Center", "Edge")], na.rm = TRUE) / 3600,
    avg_speed_bin = sum(`distance_m`[`bin_number` %in% 1:6 & zone %in% c("Center", "Edge")], na.rm = TRUE) / 1800)


# merge into df
animal_sums <- full_join(summed, speed_summed, by = c("animal", "treatment", "timepoint"))

# reorder columns to logical order
animal_sums <- relocate(animal_sums, 'avg_speed_total', .before = 'ambulatory_s_bin')
animal_sums <- relocate(animal_sums, 'sex', .before = 'treatment')


# replace full NA cols from recording error drops, which sum to zero. Un-comment below to run example.  
# recording_errors <- c('1CM')

#animal_sums <- animal_sums %>%
#  mutate(
#    `rearing_events_z_axis_total` = if_else((animal %in% recording_errors), NA_real_, `rearing_events_z_axis_total`),
#    `rearing_events_z_axis_bin` = if_else((animal %in% recording_errors),NA_real_, `rearing_events_z_axis_bin`))

#View(animal_sums)
#write.csv(animal_sums, file = here("Locomotion", "Dataframes", "animal_Sums.csv"))


#3) get summary stats w/ sex
sex_summary <- animal_sums %>%
  group_by(treatment, timepoint, sex) %>%
  summarise(
    across(
      c(ends_with("_total"), ends_with("_bin")),
      list(mean = ~ mean(.x, na.rm = TRUE),
           sem  = ~ sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))),
           n = ~ sum(!is.na(.x))
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )
write.csv(sex_summary, file = here("Locomotion", "Dataframes", "sex_summary.csv"))


# summary w/out sex
treatment_summary <- animal_sums %>%
  group_by(treatment, timepoint) %>%
  summarise(
    across(
      c(ends_with("_total"), ends_with("_bin")),
      list(mean = ~ mean(.x, na.rm = TRUE),
           sem = ~ sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))),
           n = ~ sum(!is.na(.x))
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

#View(treatment_summary)
write.csv(treatment_summary, file = here("Locomotion", "Dataframes", "treatment_summary.csv"))


################################ Analysis ################################
# Check Assumptions
# 1) check for extreme outliers
outliers <- identify_outliers(animal_sums, resting_s_total)  #replace w/ each variable to test: definitely outliers across groups due to viral effect (will remove if injections were confirmed misses)
View(outliers)

# 2) Normality check
normResult <- animal_sums %>%
  shapiro_test(resting_s_total)  
is_normal <- all(normResult$p > 0.05)
print(normResult)                         # not normal -> report p.adj
print(is_normal)


# Linear mixed-effects model (lmer from lme4)
  # mixed model: repeated measures
  # treatment and timepoint are fixed effect, animal is random effect


# locomotion variables
# set time as factor
animal_sums$timepoint <- factor(animal_sums$timepoint,
                                levels = c("0","1","2"))
# set reference terms explicitly 
animal_sums$timepoint <- relevel(factor(animal_sums$timepoint), ref = "0") # establishes baseline as with-in subject reference  
animal_sums$treatment     <- relevel(factor(animal_sums$treatment), ref = "A")  # establishes control treatment as between subject reference 

# run simple model: random intercept (animal) only
model <- lmer(distance_m_total ~ treatment * timepoint + (1 | animal), data = animal_sums) # adjust variable for stats on each measure of interest (resting mean, distance mean, etc.)
summary(model) # posthoc tests not necessary when running summary(model), multiple comparisons accounted for

# Time in Center
TIC_animal$timepoint <- relevel(factor(TIC_animal$timepoint), ref = "0") # establishes baseline as with-in subject reference  
TIC_animal$treatment     <- relevel(factor(TIC_animal$treatment), ref = "A")  # establishes control treatment as between subject reference 

model2 <- lmer(total_mean ~ treatment * timepoint + (1 | animal), data = TIC_animal)
summary(model2)
# help("pvalues",package="lme4") # pulls up R help page for model summary interpretation

# to extract traditional inference (p-values), use emmeans from emmean library 
emm <- emmeans(model2, ~ treatment * timepoint) # update model to run variables (model) vs TIC (model2)
comp <- pairs(emm, by = "timepoint", adjust = "holm") # Compare each treatment to reference treatment at each timepoint
comp_df <- as.data.frame(summary(comp, infer = c(TRUE, TRUE))) # gives CIs and p-values and stores in df
comp_df[] <- lapply(comp_df, function(x) if(is.numeric(x)) round(x, 4) else x) # round values to 4 decimal places. just note, this convers <0.001 to 0
write_xlsx(comp_df, "contrasts_results.xlsx") # save to an excel file to current wd 

# plot residuals if desired
#plot(model, type = c("p","smooth"))
#qqmath(model, id = 0.05)


######################################### Data Visualization #########################################
# Plot Open Field Measures with F-N plot_openField
  # Make plot adjustments in F-N script
  # can plot full session or subset by providing correct y_labels

# Plot Time in Center w/ viral summary TIZ df (TIC)
TIC <- TIZ%>% filter(zone == "Center") #reduce to time in center  

# plot by treatment
time_in_center <- plot_openField(
  data = TIC,
  y_var = "total_mean", 
  y_error = "total_sem",
  title = 'Time In Center',
  y_label = "Center Time (s)"
)
ggsave(here("Locomotion", "Graphs", "locomotor plots", "TimeInCenterPlot.pdf"),
       plot = time_in_center, width = 6, height = 4, dpi = 300)

# Plot Time in Edges (TIE)
TIE <- TIZ %>% filter(zone == "Edge") #reduce to time in edges
#print(TIE)


time_in_edge <- plot_openField(
  data = TIE,
  y_var = "total_mean", 
  y_error = "total_sem",
  title = "Time in Edges (60min)",
  y_label = "Time (s)"
)
ggsave(here("Locomotion", "Graphs", "locomotor plots", "TimeInEdgePlot.pdf"),
       plot = time_in_edge, width = 6, height = 4, dpi = 300)


# Any variable by treatment (Ex. ambulatory)
ambulatory <- plot_openField(
  data = aged_summary,
  y_var = "ambulatory_s_total_mean", 
  y_error = "ambulatory_s_total_sem",
  title = "Ambulatory Behavior",
  y_label = "Time (s)"
)
ggsave(here("Locomotion", "Graphs", "locomotor plots", "ambulatoryPlot.pdf"),
       plot = ambulatory, width = 6, height = 4, dpi = 300)

# Any variable by sex (Ex. ambulatory)
amb2 <- plot_openField2(
  data = aged_sex_summary,
  y_var = "ambulatory_s_total_mean", 
  y_error = "ambulatory_s_total_se",
  title = "Ambulatory Behavior",
  y_label = "Time (s)"
)
ggsave(here("Locomotion", "Graphs", "locomotor plots", "ambulatory_bySex.pdf"),
       plot = amb2, width = 6, height = 4, dpi = 300)




