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
pacman::p_load('dplyr', 'tidyverse', 'ggplot2', 'here', 'purrr', 'readr', 'viridis', 'scatterplot3d', 'plotly','lme4','lmerTest','emmeans') # auto installs required packages not yet installed
source(here('my_functionsKC.R'))

here::i_am("AutoTrack_openField.R")



############################### Create Raw Dataframe ##########################################

# open files by timepoint - make sure these order animal in the same way the animal vector is set up
file_paths <- list(
  timepoint_0 = list.files(path = here("Locomotion", "Pre-inj. Run 2", "Bins"), full.names= TRUE),
  timepoint_1 = list.files(path = here("Locomotion", "1wk post", "Bins"), full.names= TRUE),
  timepoint_2 = list.files(path = here("Locomotion", "2wk post", "Bins"), full.names= TRUE),
  timepoint_4 = list.files(path = here("Locomotion" ,"4wk post", "Bins"), full.names= TRUE),
  timepoint_10 = list.files(path = here("Locomotion" ,"10wk post", "Bins"), full.names= TRUE),
  timepoint_16 = list.files(path = here("Locomotion", "16wk post","Bins"), full.names = TRUE)
)
#print(file_paths)

# create metadata vectors
animal <- c('90030R','90030LL','102247L','102247R','102247LR', #cohort1
            '90046L', '90046R', '107544L', '107544R', '107544LR', #cohort2
            '94406R','94406LR','94406LL','94407L','94407R','94407LR','94407LL', #cohort3
            '107545LLR', '107545LR', '107549R', '107549LL', '90032R', '108807L', '108807R', #cohort4
            '108829LL', '108829LRR', '108829R', '112935L', '90042LR') #cohort 5

sex <- c( 'M', 'M', 'M', 'M', 'M', #cohort1
          'F', 'F', 'M', 'M', 'M', #cohort2
          'M','M','M','F','F','F','F', #cohort3
          'F','F','F','F','M','M','M', #cohort4
          'F','F','F','F','F') #cohort5 

virus <- c('AAV5-hTyr-Full','AAV5-hTyr-1:100','AAV5-hTyr-1:100','AAV5-hTyr-Full','AAV5-hTyr-Full', #cohort1 
           'AAV5-hTyr-1:100', 'AAV5-EYFP', 'AAV5-EYFP', 'AAV5-hTyr-1:100', 'AAV5-EYFP', #cohort2
           'AAV5-EYFP','AAV5-hTyr-Full','AAV5-hTyr-1:100','AAV5-hTyr-Full','AAV5-hTyr-1:100','AAV5-EYFP','AAV5-hTyr-Full', #cohort3
           'AAV5-hTyr-1:100','AAV5-hTyr-Full','AAV5-EYFP','AAV5-hTyr-Full','AAV5-hTyr-Full','AAV5-hTyr-1:100','AAV5-EYFP', #cohort4
           'AAV5-hTyr-1:100','AAV5-EYFP','AAV5-EYFP','AAV5-hTyr-Full','AAV5-hTyr-1:100') #cohort5

timepoints <- c(0,1,2,4,10,16)


# Set an empty list to store dataframes
loco_list <-list()

expected_cols <- c(
  "Bin #", "Start", "Finish", "Records", "Zone", "Entrances",
  "Ambulatory (s)", "Stereotypic (s)", "Resting (s)", "Undetected (s)",
  "Rearing Events (Z-Axis)", "Rearing Events (V-Axis)",
  "Rearing Time Z (s)", "Rearing Time V (s)",
  "Distance (cm)", "Average Speed (cm/s)",
  "animal", "sex", "virus", "timepoint")


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
    loco$virus <- virus[i]
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
loco <- dplyr::select(loco,"animal", "sex", "virus", "timepoint", everything())

# clean column names to lowercase and _
loco <- clean_names(loco)
#View(loco)




########################### Pre-processing & Cleaning #################################################################

## Use my FN dropErrors to make corrections to recording errors 

loco <- dropErrors(loco, "animal", "102247L", "rearing_events_z_axis","rearing_time_z_s") #ch 5 broken,
loco <- dropErrors(loco, 'animal', '90030LL', "rearing_events_z_axis","rearing_time_z_s") #ch 4 rearing time wrong
loco <- dropErrors(loco, 'animal', '90046R', "rearing_events_z_axis","rearing_time_z_s") #ch 4 rearing time wrong
loco <- dropErrors(loco, 'animal', '107549LL', "rearing_events_z_axis","rearing_time_z_s") #ch 4 rearing time wrong
loco <- dropErrors(loco, 'animal', '94407L', "rearing_events_z_axis","rearing_time_z_s") #ch 4 rearing time wrong

loco <- loco %>% mutate(across(
  c("rearing_events_z_axis","rearing_time_z_s"), 
  ~if_else(animal == '90030R' & timepoint == 16, NA_real_, .) # this mouse ran in CH 4 on week 16 which wasn't recording rearing time correctly, so only this timepoint rearing dropped
))

# drop 107545LLR data at 10 weeks - euthanized for hindlimb paralysis just after session
loco <- loco %>% filter(!(animal == "107545LLR" & timepoint == 10))

#View(loco)

write.csv(loco,file = here("Locomotion","Dataframes", "aged_raw_locomotion.csv"))

loco <- read.csv(file = here("Locomotion","Dataframes", "aged_raw_locomotion.csv"))

##### Summary Dataframes ####

# drop missed injection animals (confirmed by IHC / ephys)
loco <- loco %>%
  filter(animal != '94407LL')



# custom function: animalSums will sum over full session (60min), while the added bin_range = 1:6 will sum first 30min so either can be used in analysis and graphing
# similarly, timeInZones will be used to calculate time spent in the center vs edge of the open field chamber


#1) Summarize time in zones (TIZ) (ambulatory, stereotypic, and resting) full session 
# by virus
TIZ <- timeInZones(loco, bin_col= "bin_number", bin_range = 1:6, virus, timepoint, zone)

# filter out 'none' observations
TIZ <- TIZ %>% filter(zone != 'None')
#View(TIZ)

write.csv(TIZ,file = here("Locomotion","Dataframes","timeInZone_summary.csv"))

# by sex and virus
TIZ_2 <- timeInZones(loco, bin_col= "bin_number", bin_range = 1:6, virus, sex, timepoint, zone)
TIZ_2 <- TIZ_2 %>% filter(zone != 'None')
#View(TIZ_2)

write.csv(TIZ_2,file = here("Locomotion","Dataframes","timeInZone_sex_summary.csv"))



#2) Using custom FN: animalSums - get sum of all variable per animal (average speed handled outside of function)
# automate with variable vector
sum_variables <- c("ambulatory_s", "stereotypic_s", "resting_s", "rearing_events_z_axis", "rearing_time_z_s", "distance_cm")


# Apply animalSums FN
summed <- animalSums(
  df = loco, 
  zone_col = "zone", 
  zones = c("Center", "Edge"), 
  bin_col = "bin_number", 
  bin_range = 1:6, 
  variables = sum_variables,
  animal, virus, timepoint
)

# Handle speed separately
speed_summed <- loco %>%
  group_by(animal, virus, timepoint, sex) %>%
  summarise(
    avg_speed_total = sum(`distance_cm`[zone %in% c("Center", "Edge")], na.rm = TRUE) / 3600,
    avg_speed_bin = sum(`distance_cm`[`bin_number` %in% 1:6 & zone %in% c("Center", "Edge")], na.rm = TRUE) / 1800)


# merge into df
animal_sums <- full_join(summed, speed_summed, by = c("animal", "virus", "timepoint"))

# move cols to logical order
animal_sums <- relocate(animal_sums, 'avg_speed_total', .before = 'ambulatory_s_bin')
animal_sums <- relocate(animal_sums, 'sex', .before = 'virus')


# replace full NA cols from recording error drops, which sum to zero  
  # could not resolve this issue in FN
recording_errors <- c('102247L', '90030LL', '90046R', '107549LL', '94407L')

animal_sums <- animal_sums %>%
  mutate(
    `rearing_events_z_axis_total` = if_else((animal %in% recording_errors) | (animal == '90030R' & timepoint == 16),
                                            NA_real_, `rearing_events_z_axis_total`),
    `rearing_events_z_axis_bin` = if_else((animal %in% recording_errors) | (animal == '90030R' & timepoint == 16),
                                          NA_real_, `rearing_events_z_axis_bin`),
    `rearing_time_z_s_total` = if_else((animal %in% recording_errors) | (animal == '90030R' & timepoint == 16),
                                       NA_real_, `rearing_time_z_s_total`),
    `rearing_time_z_s_bin` = if_else((animal %in% recording_errors) | (animal == '90030R' & timepoint == 16),
                                     NA_real_, `rearing_time_z_s_bin`)
  )

#View(animal_sums)
write.csv(animal_sums, file = here("Locomotion", "Dataframes", "aged_animal_Sums.csv"))


#3) get summary stats w/ sex
aged_sex_summary <- animal_sums %>%
  group_by(virus, timepoint, sex) %>%
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


#View(aged_sex_summary)
write.csv(aged_sex_summary, file = here("Locomotion", "Dataframes", "aged_sex_summary.csv"))

# summary w/out sex
aged_summary <- animal_sums %>%
  group_by(virus, timepoint) %>%
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

#View(aged_summary)
write.csv(aged_summary, file = here("Locomotion", "Dataframes", "aged_virus_summary.csv"))








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

# 2-way linear fixed effect model (lmer from lme4)
  # mixed model: repeated measures
  # virus and timepoint are fixed effect, animal is random effect
# choosing lmer over 2-way anova here because missing data, repeated measures, and non-homogeneous variance b/w groups

# set time as factor
animal_sums$timepoint <- factor(animal_sums$timepoint,
                                levels = c("0","1","2","4","10","16"))

# run 2 models: model1 = random intercept (animal), model2 = randome intercept (animal) + slope (time trajectory / animal)
model1 <- lmer(resting_s_total ~ virus * timepoint + (1 | animal), data = animal_sums) # adjust variable for stats on each measure of interest (resting mean, distance mean, etc.)
# model2 <- lmer(resting_s_total ~ virus * timepoint + (timepoint | animal), data = animal_sums) # is singular, so likely an overfit model for the data

# ANOVA table for fixed effects
anova(model1)
# anova(model2)

# Test which model to use, though singularity tells us model2 is likely an overfit
# anova(model1,model2) # higher logLik = better, p-value will tell if model2 is sig. better fit, lower AIC = better fit
# AIC(model1,model2)  # delta AIC <2: models ~ equal, 4-7: model w/ lower AIC may be better, >10: model w/ lower AIC strongly better

# Post-hoc: Dunnett contrast: compare every virus against the control w/in each timepoint 
# emmean = Estimated Marginal Means, (Least-Squares Means)
emm <- emmeans(model1, ~ virus | timepoint)
dunnet_results <- contrast(emm, method = 'dunnett', ref = 'AAV5-EYFP', adjust = 'holm')
dunnet_results

# NOTE: these stats are different than reported in Prism where a 2 Way ANOVA was applied with Tukey multiple measures correction.
  # they're nearly the same outcome but I personally think LMER is more appropriate since it accounts for missing values & individual trajectories over time. (9.10.25)



######################################### Data Vizualization #########################################

## Plot Open Field Measures w/ my FN plot_openField
  # Make plot adjustments in FN script
  # can plot full session or subset by providing correct y_labels (all 60min plots currently 4.22.25)


# Plot Time in Center (TIC)
TIC <- TIZ %>% filter(zone == "Center") #reduce to time in center
# print(TIC)

# plot by virus
time_in_center <- plot_openField(
  data = TIC,
  y_var = "total_mean", 
  y_error = "total_sem",
  title = 'Time In Center',
  y_label = "Center Time (s)"
)
ggsave(here("Locomotion", "Graphs", "locomotor plots", "aged_TimeInCenterPlot.pdf"),
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
#ggsave(here("Locomotion", "Graphs", "locomotor plots", "aged_TimeInEdgePlot.pdf"),
#       plot = time_in_edge, width = 6, height = 4, dpi = 300)


# ambulatory by virus
ambulatory <- plot_openField(
  data = aged_summary,
  y_var = "ambulatory_s_total_mean", 
  y_error = "ambulatory_s_total_sem",
  title = "Ambulatory Behavior",
  y_label = "Time (s)"
)
#ggsave(here("Locomotion", "Graphs", "locomotor plots", "aged_ambulatoryPlot.pdf"),
#       plot = ambulatory, width = 6, height = 4, dpi = 300)

# ambulatory by sex
# amb2 <- plot_openField2(
#   data = aged_sex_summary,
#   y_var = "ambulatory_s_total_mean", 
#   y_error = "ambulatory_s_total_se",
#   title = "Ambulatory Behavior",
#   y_label = "Time (s)"
# )
# ggsave(here("Locomotion", "Graphs", "locomotor plots", "aged_ambulatory_bySex.pdf"),
#        plot = amb2, width = 6, height = 4, dpi = 300)



# stereotypic by virus
# stereotypic <- plot_openField(
#   data = aged_summary,
#   y_var = "stereotypic_s_total_mean", 
#   y_error = "stereotypic_s_total_se",
#   title = "Stereotypic Behavior",
#   y_label = "Time (s)"
# )
# ggsave(here("Locomotion", "Graphs", "locomotor plots", "aged_stereotypicPlot.pdf"),
#        plot = stereotypic, width = 6, height = 4, dpi = 300)

# stereotypic by sex
# stereotypic2 <- plot_openField2(
#   data = aged_sex_summary,
#   y_var = "stereotypic_s_total_mean", 
#   y_error = "stereotypic_s_total_se",
#   title = "Stereotypic Behavior",
#   y_label = "Time (s)"
# )
# ggsave(here("Locomotion", "Graphs", "locomotor plots", "aged_stereotypic_bySex.pdf"),
#        plot = stereotypic2, width = 6, height = 4, dpi = 300)



# resting by virus
resting <- plot_openField(
  data = aged_summary,
  y_var = 'resting_s_total_mean',
  y_error = 'resting_s_total_sem',
  title = 'Time Resting',
  y_label = 'Time (s)'
)
ggsave(here("Locomotion", "Graphs", "locomotor plots", "aged_restingPlot.pdf"),
       plot = resting, width = 6, height = 4, dpi = 300)

# resting by sex
resting2 <- plot_openField2(
  data = aged_sex_summary,
  y_var = 'resting_s_total_mean',
  y_error = 'resting_s_total_sem',
  title = 'Time Resting',
  y_label = 'Time (s)'
)
ggsave(here("Locomotion", "Graphs", "locomotor plots", "aged_resting_bySex.pdf"),
       plot = resting2, width = 6, height = 4, dpi = 300)



# rearing events by virus
rearing <- plot_openField(
  data = aged_summary,
  y_var = 'rearing_events_z_axis_total_mean',
  y_error = 'rearing_events_z_axis_total_sem',
  title = 'Rearing Events',
  y_label = 'Rearing Count'
)
ggsave(here("Locomotion", "Graphs", "locomotor plots", "aged_rearingEventsPlot.pdf"),
       plot = rearing, width = 6, height = 4, dpi = 300)

# rearing events by sex
rearing2 <- plot_openField2(
  data = aged_sex_summary,
  y_var = 'rearing_events_z_axis_total_mean',
  y_error = 'rearing_events_z_axis_total_sem',
  title = 'Rearing Events',
  y_label = 'Rearing Count'
)
ggsave(here("Locomotion", "Graphs", "locomotor plots", "aged_rearingEvents_bySex.pdf"),
       plot = rearing2, width = 6, height = 4, dpi = 300)




# rearing time by virus
time_rearing <- plot_openField(
  data = aged_summary,
  y_var = 'rearing_time_z_s_total_mean',
  y_error = 'rearing_time_z_s_total_sem',
  title = 'Time Rearing',
  y_label = 'Time (s)'
)
ggsave(here("Locomotion", "Graphs", "locomotor plots", "aged_rearingTimePlot.pdf"),
       plot = time_rearing, width = 6, height = 4, dpi = 300)

# rearing time by sex
time_rearing2 <- plot_openField2(
  data = aged_sex_summary,
  y_var = 'rearing_time_z_s_total_mean',
  y_error = 'rearing_time_z_s_total_sem',
  title = 'Time Rearing',
  y_label = 'Time (s)'
)
ggsave(here("Locomotion", "Graphs", "locomotor plots", "aged_rearingTime_bySex.pdf"),
       plot = time_rearing, width = 6, height = 4, dpi = 300)




# distance by virus
distance <- plot_openField(
  data = aged_summary,
  y_var = 'distance_cm_total_mean',
  y_error = 'distance_cm_total_sem',
  title = 'Distance Traveled',
  y_label = 'Distance (cm)'
)
ggsave(here("Locomotion", "Graphs", "locomotor plots", "aged_distancePlot.pdf"),
       plot = distance, width = 6, height = 4, dpi = 300)

# distance by sex
distance2 <- plot_openField2(
  data = aged_sex_summary,
  y_var = 'distance_cm_total_mean',
  y_error = 'distance_cm_total_sem',
  title = 'Distance Traveled',
  y_label = 'Distance (cm)'
)
ggsave(here("Locomotion", "Graphs", "locomotor plots", "aged_distance_bySex.pdf"),
       plot = distance2, width = 6, height = 4, dpi = 300)




# ave speed by virus
aveSpeed <- plot_openField(
  data = aged_summary,
  y_var = 'avg_speed_total_mean',
  y_error = 'avg_speed_total_sem',
  title = 'Average Speed',
  y_label = 'Speed (cm/s)'
)
ggsave(here("Locomotion", "Graphs", "locomotor plots", "aged_aveSpeedPlot.pdf"),
       plot = aveSpeed, width = 6, height = 4, dpi = 300)

# ave speed by sex
aveSpeed2 <- plot_openField2(
  data = aged_sex_summary,
  y_var = 'avg_speed_total_mean',
  y_error = 'avg_speed_total_sem',
  title = 'Average Speed',
  y_label = 'Speed (cm/s)'
)
ggsave(here("Locomotion", "Graphs", "locomotor plots", "aged_aveSpeed_bySex.pdf"),
       plot = aveSpeed2, width = 6, height = 4, dpi = 300)




# Scrap  ------------------------------------------------------------------

## Use my FN baselineSums to collect timepoint == 0 summaries : used to counter-balance viral groups pre-injection
# pre_inj_totals <- baselineSums(
# df = loco, 
# group = "animal", 
# timepoint = "timepoint",
# bin_col = "bin_number", bin_range = 1:12,
# "distance_cm"
#)
#View(pre_inj_totals)


#pre_injection_30min <- baselineSums(
# df = loco, 
# group = "animal", 
# timepoint = "timepoint",
# bin_col = "bin_number", bin_range = 1:6,
# "distance_cm"
#)

# rename the sum column for 30min to distinguish them 
#names(pre_injection_30min)[names(pre_injection_30min) == "distance_cm_sum"] <- "30min_Distance_sum"
#print(pre_injection_30min)

# bind and write to csv to counter-balance animals pre-injection (based on forward locomotion)
#pre_inj_sums <- merge(pre_inj_totals, pre_injection_30min, by = "animal")
#View(pre_inj_sums)

#write.csv(pre_inj_sums, file = here("Locomotion", "Dataframes", "preInjection_summary.csv"))




# Step 3) Run analysis (no longer need to do this 6.3.25)

# check observations for anova (missing observations will break ANOVA so must remove until dataset complete)
#table(animal_sums$sex, animal_sums$virus, animal_sums$timepoint) 


# A) Look at trends for lab meeting. Remove this when dataset complete and *attempt* my FN which runs everything built in w/out removing timepoints
#animal_sums3<- animal_sums2[animal_sums2$timepoint!=16, ] #remove timepoints w/out observations
#animal_sums3$timepoint=factor(animal_sums3$timepoint) #set time as factor
#animal_sums3=animal_sums3[,c(1:11)] #select columns of interest (full time 1:11)
#   30min bins would be 1:4 & 12:18)
# View(animal_sums3)

# Anova per each variable
# filter for rearing events 
# animal_sums3 <- animal_sums3 %>%
#   filter(rearing_events_z_axis_total != 'NA')

