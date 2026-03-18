# This Script -------------------------------------------------------------
  # Analyzes MouseSpecifics Digigait time series data across multiple speeds
  # The speeds dataframe is build from manually recorded data while the gait indices dataframe is build from the 'Indices' readout for each animal at each timepoint.
  # summaries & analyses are broken down by speed so comparisons are made within a single speed only. 
# Load Libraries
# Build Indices Dataframe
  # post-process digigait 'Indices' files
  # clean df: remove additional row coercies from df build
  # convert to proper data types 
# Build Speed Dataframe
  # create vectors of manual data and metadata
  # build df
  # drop treatment errors & save
# Pre-Processing
  # set datatypes
  # drop injection or recording errors from indices df
# Summary Dataframes
  # subset dfs by belt speed
  # subset by paw / belt speed
  # select variables of interest
  # summarize fore paws together and hind paws together per belt speed
# Visualization
  # delta max speed / group (Facet, trajectory, and lineplot)
  # delta indice variable / group (lineplot, box plot, violin plot)
# Analysis 
  # summary stats - variable / paw & belt speed
  # Linear mixed effects model
  # Assumptions checks
  # LMER from LME4, run summary on model
  # holm post-hoc multiple comparison to extract CIs & p-values


# Load Libraries ----------------------------------------------------------
pacman::p_load('janitor','dplyr', 'tidyverse', 'ggplot2', 'here', 'purrr', 'readr', 'viridis', 'scatterplot3d', 'plotly','mosaic') # auto installs required packages not yet installed
source(here('Gait_Functions.R'))
here::i_am("Digigait_timeSeries.R")

# to collect pack citations, uncomment and run line below:
# cite_packages()

# to collect pack versions, uncomment and run line below:
# report_packages()



# Build Indices Dataframe -------------------------------------------------------
# open files by timepoint - make sure these order animal in the same way the animal vector is set up
file_paths <- list(
  timepoint_0 = list.files(path = here("Digigait","wk_0"), full.names= TRUE),
  timepoint_1 = list.files(path = here("Digigait", "wk_1"), full.names= TRUE),
  timepoint_2 = list.files(path = here("Digigait", "wl_2"), full.names= TRUE))

# print filepaths to ensure all timepoints are present
print(file_paths)

# set metadata : define treatment mapping so treatmentes are mapped to animal IDs
treatment_mapping <- data.frame(
  animal = c('1CM', '2AM', '3AF', '4CF','5BM', '6BF'),
  sex = c( 'M', 'M', 'F', 'F', 'M', 'F'),
  treatment = c('C','A','A','C','B','B'))

  

# remove whitespaces for smooth df left_join on animal IDs
treatment_mapping <- treatment_mapping %>%
  mutate(animal = as.character(trimws(animal))) %>%
  distinct(animal, .keep_all = TRUE)

# set time metadata
timepoint <- c('0','1','2')


# Set an empty list to store dataframes
digi_list <-list()

# outer loop assigns timepoints to each file
for (tp in names(file_paths)) {
  #get file paths for current timepoint
  files <- file_paths[[tp]]
  timepoint_label <- gsub("timepoint_","", tp)
  
  #inner loop processes each file for the current timepoint 
  for (i in seq_along(files)) {
    
    # read the first two rows (headers)
    header <- read.csv(files[i], nrows = 2, header = FALSE) # read
    col_names <- paste(header[1, ], header[2, ]) # combine two headers into one, separated by space
    col_names <- make_clean_names(col_names) #clean names to avoid base R adding unnecessary X and punctuation
    
    # now read csv files skipping header rows
    digi <- read.csv(files[i], skip = 2, header = FALSE)
    colnames(digi) <- col_names
    
    # add metadata columns
    digi$timepoint <- timepoint_label
    
    # merge with treatment mapping
    digi <- digi %>%
      mutate(animal = trimws(as.character(animal))) %>% #ensure proper joining by trimming and converting animal IDs to characters. 
      left_join(treatment_mapping %>% mutate(animal = trimws(as.character(animal))), by = 'animal')
    
    # coerce column data types
    digi <- digi %>%
      mutate(across(c(stance_width_cm,step_angle_deg,sw_var_cm,step_angle_var_deg,stance_width_cv_cv_percent, step_angle_cv_cv_percent,hind_limb_shared_stance_time_s,
                      percent_shared_stance_percent, stance_factor_real_number, tau_propulsion_real_number, paw_drag_real_number,weight_grams), ~suppressWarnings(as.numeric(.))))
    
    # Append to the list
    digi_list <-append(digi_list, list(digi))
    
    # Re-apply data type coerce to weights (still reading one as char)
    digi_list <- lapply(digi_list, function(df) {
      df %>% mutate(weight_grams = suppressWarnings(as.numeric(weight_grams)))
    })
    
  }
}


# combing all data to a dataframe
digigait_df <- bind_rows(digi_list)

# reorder columns so metadata is first
digigait_df <- digigait_df[, c("animal", "timepoint","treatment","sex", setdiff(names(digigait_df), c("animal", "timepoint","treatment","sex")))]

# clean white spaces from limb type
digigait_df <- digigait_df %>%
  mutate(limb = trimws(limb))

# check it worked
unique(digigait_df$limb) # returns an empty extra string



View(digigait_df)


# save df
write.csv(digigait_df, file = here("Digigait", "Dataframes", "total_indices_df.csv"))


# Build Speed Dataframe ---------------------------------------------------
# min and max speed df
# 1 observation/animal at all timepoints, then move to next animal. Ex. Animal 1 time 0 = 35, time 1 = 25 , time 2-16 = 0, then Animal 2...

max_speed <- c(35,25,0,   35,25,35,   25,20,30,   35,15,0,   45,25,15,   30,15,30)
min_speed <- c(15,15,0,   15,20,20,   15,15,20,   15,15,0,   20,15,15,   15,15,20)

# create metadata vectors
animal <- c('1CM', '2AM', '3AF', '4CF','5BM', '6BF') 
sex <- c( 'M', 'M', 'F', 'F', 'M', 'F')  
treatment <- c('C','A','A','C','B','B') 
timepoints <- c(0,1,2)


# Expand grid so that each timepoint is repeated across all animals first, then trials
speed_df <- expand.grid(animal = animal, timepoint = timepoint)

# map sex and treatment to each animal based on their index
speed_df$sex <- sex[match(speed_df$animal, animal)]
speed_df$treatment <- treatment[match(speed_df$animal, animal)]

# Sort by timepoint, then animal, then trial to ensure correct order in df
speed_df <- speed_df[order(speed_df$animal, speed_df$timepoint), ]

# Add speed data
speed_df$max_speed <- max_speed[1:nrow(speed_df)]
speed_df$min_speed <- min_speed[1:nrow(speed_df)]

# Reorder columns 
speed_df <- speed_df[, c('animal', 'sex', 'treatment','timepoint', 'min_speed', 'max_speed')]

#View(speed_df)


# export to CSV file
write.csv(speed_df, here("Digigait", "Dataframes", "speed_df.csv"))






# Pre-Processing ----------------------------------------------------------
# Indices df 
digigait_df <- digigait_df %>%
  mutate(timepoint = as.numeric(as.character(timepoint))) # set timepoint datatype to chr

# set animal to factor and arrange df to match metadata vector
digi <- digi %>%
  mutate(animal = factor(animal, levels = treatment_mapping$animal)) %>% 
  arrange(animal)

#View(digi)


# Summary Dataframes ------------------------------------------------------

# 1) Individual paws by belt speed
# subset digigait df to speeds
belt_15 <- subset(digi, gait_speed_cm_s == 15) 
belt_20 <- subset(digi, gait_speed_cm_s == 20)


# write each to csv
write.csv(belt_15, here("Digigait", "Dataframes", "Indices_15cms", "15cm_s_indices.csv"))

# subset each into all 4 paws to get summary dfs : each belt speed captures a diff. level of performance that varies by animal
LH_15 <- subset(belt_15, limb == ' Left Hind ')
RH_15 <- subset(belt_15, limb == 'Right Hind')
LF_15 <- subset(belt_15, limb == 'Left Fore')
RF_15 <- subset(belt_15, limb == 'Right Fore')

LH_20 <- subset(belt_20, limb == 'Left Hind')
RH_20 <- subset(belt_20, limb == 'Right Hind')
LF_20 <- subset(belt_20, limb == 'Left Fore')
RF_20 <- subset(belt_20, limb == 'Right Fore')


# str(belt_15) to view variables of interest for summaries (can add or subtract here)
variables <- c("percent_swing_stride_percent", "percent_brake_stride_percent", "percent_propel_stride_percent", "percent_stance_stride_percent",
               "stride_length_cm", "stride_frequency_steps_s", "absolute_paw_angle_deg", "paw_angle_variability_deg", "stance_width_cm", "sl_var_cm", 
               "sw_var_cm", "step_angle_var_deg", "step_angle_var_deg", "hind_limb_shared_stance_time_s", "percent_shared_stance_percent")

# use gaitSummaries FN to get summaries by each paw: repeat for all speeds
LH_summary_15 <- gaitSummaries(df = LH_15, variables = variables, treatment, timepoint) #can add sex here
RH_summary_15 <- gaitSummaries(df = RH_15, variables = variables, treatment, timepoint)
LF_summary_15 <- gaitSummaries(df = LF_15, variables = variables, treatment, timepoint)
RF_summary_15 <- gaitSummaries(df = RF_15, variables = variables, treatment, timepoint)

# export each to CSV file if desired
write.csv(RF_summary_15, here("Digigait", "Dataframes", "Indices_15cms", "RF_15_treatmentSummary.csv"))


# 3) Hind paws together, then fore paws together per speed: change speed integer for different belt speeds in study
# 15 cm/s
hind_15 <- subset(belt_15, limb == 'Left Hind' | limb == 'Right Hind')
fore_15 <- subset(belt_15, limb == 'Left Fore' | limb == 'Right Fore')

hind_summary_15 <- gaitSummaries(df = hind_15, variables = variables, treatment, timepoint)
fore_summary_15 <- gaitSummaries(df = fore_15, variables = variables, treatment, timepoint)

# export to CSV file
write.csv(hind_summary_15, here("Digigait", "Dataframes", "Indices_15cms", "hindpaws_15_Summary.csv"))
write.csv(fore_summary_15, here("Digigait", "Dataframes", "Indices_15cms", "forepaws_15_Summary.csv"))





########################### Data Visualization #####################

## Change in max speed plots
speed_df$timepoint <- factor(speed_df$timepoint, levels=c('0','1','2'))

# 1) Facet Plot - treatment group variability 
facet_box <- ggplot(speed_df, aes(x=timepoint, y=max_speed, fill=treatment)) +
  geom_boxplot() +
  facet_wrap(~treatment) +    #facet plot by treatment
  labs(
    title = 'Max Speed Variability',
    x= 'Time from Treatment',
    y = 'Belt Speed (cm/s)',
    fill = 'treatment Type'
  ) +
  scale_fill_manual(
    values=c('#A3E4D7','#F8C471','#E67E22', 'black'),
    labels = c('A','B-1:100','C')
  ) +
  scale_x_discrete(breaks = c('0','1','2')    #x-axis labeled as timepoints 
  ) +
  theme_bw() 
facet_box
ggsave(here("Digigait", "Graphs", "beltSpeed_FacetPlot.pdf"),
       plot = facet_box, width = 6, height = 4, dpi = 300)



#2 Trajectory Plot - individual change from baseline / each animal 
# subset and mutate df
delta_speed <- speed_df %>% 
  filter(!is.na(max_speed)) %>% # filter NAs
  group_by(animal, treatment) %>%
  mutate(
    baseline =max_speed[timepoint == 0][1], #explicitly state baseline
    data_change = max_speed - baseline, #sub baseline at each timepoint
    timepoint = factor(timepoint)) #ensures timepoints print as exact in x-axis
#View(delta_speed)


# create plot        
traj <-  ggplot(delta_speed, aes(x =timepoint, y = data_change, color = treatment, group = animal)) +
  geom_line(linewidth = 1.0) +    #Line connecting viral points
  labs(title= 'Individual Change From Baseline : Max Speed',
       x = "Time from Treatment", 
       y = "Speed (cm/s)", 
       color = "treatment Type", 
       linetype = "treatment Type") +
  theme_minimal() +
  scale_color_manual(values=c('#A3E4D7','#F8C471','#E67E22', 'black')) +
  scale_x_discrete(labels = levels(delta_speed$timepoint)) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
traj
ggsave(here("Digigait", "Graphs", "delta_maxSpeed_trajectory.pdf"),
       plot = traj, width = 6, height = 4, dpi = 300)




########################### Analysis #####################

# summary stats w/ mosaic package
favstats(~treatment, data = LH_15) 

# max speed uses the df: speed_df
# all other variables: use subset df by speed of choice
speed_df <- read.csv(here("Digigait", "Dataframes", "speed_df.csv"))

# A) Check for Outliers
identify_outliers(speed_df, variable = "max_speed")  # Full will always through outliers, but they're true data

# 1) clean df for analysis    
# filter for any 'NA' in repeated measures 
filtered <- speed_df %>%
  filter(speed_df$max_speed != 'NA')
View(filtered)


# B) normality check
normResult <- filtered %>%
  shapiro_test(max_speed)
is_normal <- all(normResult$p > 0.05)
print(normResult) # super not normal -> report p.adj


# Linear mixed-effects model (lmer from lme4)
# mixed model: repeated measures
# treatment and timepoint are fixed effect, animal is random effect

# set reference terms explicitly 
speed_df$timepoint <- relevel(factor(speed_df$timepoint), ref = "0") # establishes baseline as with-in subject reference  
speed_df$treatment     <- relevel(factor(speed_df$treatment), ref = "A")  # establishes control treatment as between subject reference 

model <- lmer(max_speed~treatment * timepoint + (1 | animal), data = speed_df) # (1 | animal) accounts for repeated measures of animals
summary(model)

# to extract traditional inference (p-values), use emmeans from emmean library 
emm <- emmeans(model, ~ treatment * timepoint) # update model to run variables (model) vs TIC (model2)
comp <- pairs(emm, by = "timepoint", adjust = "holm") # Compare each treatment to reference treatment at each timepoint
comp_df <- as.data.frame(summary(comp, infer = c(TRUE, TRUE))) # gives CIs and p-values and stores in df
comp_df[] <- lapply(comp_df, function(x) if(is.numeric(x)) round(x, 4) else x) # round values to 4 decimal places. just note, this convers <0.001 to 0
write_xlsx(comp_df, "contrasts_results.xlsx") # save to an excel file to current wd 


# plot residuals if desired
#plot(model, type = c("p","smooth"))
#qqmath(model, id = 0.05)




