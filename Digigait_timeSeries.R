# This Script -------------------------------------------------------------
  # This script analysis MouseSpecifics Digigait timeseries data across multiple speeds
  # The speeds dataframe is build from manually recorded data while the gait indices dataframe is build from the 'Indices' readout for each animal at each timepoint.
  # summaries & analyses are broken down by speed so comparisons are made within a single speed only. 
# Load Libraries
# Build Indices Dataframe
  # post-process digigait 'analysis'Indices' files
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
  # delta speed by treatment group and by treatment group * sex
  # subset dfs by belt speed
  # subset by paw / belt speed
  # select variables of interest
  # summaries fore paws together and hind paws together per belt speed
# Visualization
  # delta max speed / group (Facet, trajectory, and lineplot)
  # delta indice variable / group (lineplot, box plot, violin plot)
# Analysis 
  # summary stats - variable / paw & belt speed
  # 3 Way Repeated Measures ANOVA
    # Assumptions checks
    # Remove animals not yet ran at speeds or removed from study 
    # Post hoc tests


# Load Libraries ----------------------------------------------------------
pacman::p_load('janitor','dplyr', 'tidyverse', 'ggplot2', 'here', 'purrr', 'readr', 'viridis', 'scatterplot3d', 'plotly','mosaic') # auto installs required packages not yet installed
source(here::here('Gait_Functions.R'))
here::i_am("Digigait_timeSeries.R")




# Build Indice Dataframe -------------------------------------------------------
# open files by timepoint - make sure these order animal in the same way the animal vector is set up
file_paths <- list(
  timepoint_0 = list.files(path = here::here("Digigait","pre-inj"), full.names= TRUE),
  timepoint_1 = list.files(path = here::here("Digigait", "1wk"), full.names= TRUE),
  timepoint_2 = list.files(path = here::here("Digigait", "2wk"), full.names= TRUE),
  timepoint_4 = list.files(path = here::here("Digigait","4wk"), full.names= TRUE),
  timepoint_10 = list.files(path = here::here("Digigait","9-10wk"), full.names= TRUE),
  timepoint_16 = list.files(path = here::here("Digigait", "16wk"), full.names= TRUE)
)

# print filepaths to ensure all timepoints are present
print(file_paths)

# set metadata : define virus mapping so viruses are mapped to animal IDs
virus_mapping <- data.frame(
  animal = c('90030R','90030LL','102247L','102247R','102247LR', #cohort1
             '90046L', '90046R', '107544L', '107544R', '107544LR', #cohort2
             '94406R','94406LR','94406LL','94407L','94407R','94407LR','94407LL', #cohort3 
             '107545LLR', '107545LR', '107549R', '107549LL', '90032R', '108807L', '108807R', #cohort4
             '108829LL', '108829LRR', '108829R', '112935L', '90042LR'), #cohort5
  
  virus = c('AAV5-hTyr-Full','AAV5-hTyr-1:100','AAV5-hTyr-1:100','AAV5-hTyr-Full','AAV5-hTyr-Full', #cohort1
            'AAV5-hTyr-1:100', 'AAV5-EYFP', 'AAV5-EYFP', 'AAV5-hTyr-1:100', 'AAV5-EYFP', #cohort2
            'AAV5-EYFP', 'AAV5-hTyr-Full', 'AAV5-hTyr-1:100', 'AAV5-hTyr-Full', 'AAV5-hTyr-1:100', 'AAV5-EYFP', 'AAV5-hTyr-Full', #cohort3
            'AAV5-hTyr-1:100','AAV5-hTyr-Full','AAV5-EYFP','AAV5-hTyr-Full','AAV5-hTyr-Full','AAV5-hTyr-1:100','AAV5-EYFP', #cohort4
            'AAV5-hTyr-1:100','AAV5-EYFP','AAV5-EYFP','AAV5-hTyr-Full','AAV5-hTyr-1:100') #cohort5
)

# remove whitespaces for smooth df left_join on animal IDs
virus_mapping <- virus_mapping %>%
  mutate(animal = as.character(trimws(animal))) %>%
  distinct(animal, .keep_all = TRUE)

# set time metadata
timepoint <- c('0','1','2','4','10','16')


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
    col_names <- janitor::make_clean_names(col_names) #clean names to avoid base R adding unnecessary X and punctuation

    # now read csv files skipping header rows
    digi <- read.csv(files[i], skip = 2, header = FALSE)
    colnames(digi) <- col_names
    
    # add metadata columns
    digi$timepoint <- timepoint_label
    
    # merge with virus mapping
    digi <- digi %>%
      mutate(animal = trimws(as.character(animal))) %>% #ensure proper joining by trimming and converting animal IDs to characters. 
      left_join(virus_mapping %>% mutate(animal = trimws(as.character(animal))), by = 'animal')
    
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
digigait_df <- digigait_df[, c("animal", "timepoint","virus","gender", setdiff(names(digigait_df), c("animal", "timepoint","virus","gender")))]

# clean white spaces 
digigait_df <- digigait_df %>%
  mutate(across(where(is.character), trimws))

# check it worked
unique(digigait_df$limb) # returns an empty extra string
sum(digigait_df$limb == "") # 7 in fact, hmm... reading blanks

# lets clean them up
digigait_df <- digigait_df %>% filter(limb != "")

#preview
View(digigait_df)

# Virus NAs coerced for 1 mouse in cohort 5 at time = 10 and for another mouse at time =16... cannot resolve, so fill manually
#digigait_df <- digigait_df$virus[digigait_df$animal == '108829LL' & is.na(digigait_df$virus)] <- 'AAV5-hTyr-1:100'

# save df
write.csv(digigait_df, file = here("Digigait", "Dataframes", "total_indices_df.csv"))


# Build Speed Dataframe ---------------------------------------------------
# min and max speed df
  # 1 observation/animal at all timepoints, then move to next animal. Ex. 90030R time 0 = 35, time 1 = 25 , time 2-16 = 0, then 90030LL...

max_speed <- c(35,25,0,0,0,0,      35,25,35,30,20,30,  25,20,30,20,20,25,  35,15,0,0,0,0,      45,25,15,15,0,0, #cohort1
               30,15,30,35,30,25,  20,40,45,45,45,45,  25,0,30,45,45,40,   20,20,30,40,35,40,  25,20,30,40,35,25, #cohort2 
               35,35,35,25,35,40,  15,20,15,0,0,0,     30,30,45,40,35,40,  40,20,15,0,0,0,     45,40,45,45,45,40,   45,40,35,25,25,25,  45,40,30,15,30,35, #cohort3
               45,35,20,0,0,0,     25,35,40,30,25,45,  20,45,45,40,45,25,  40,30,30,45,NA,NA,  35,45,30,15,15,0,    40,40,45,45,45,45,  40,30,20,15,0,0, #cohort4
               40,25,30,0,20,20,   30,15,NA,35,35,15,  40,15,NA,45,45,30,  35,15,35,20,0,0,   35,30,40,35,40,25) #cohort5

min_speed <- c(15,15,0,0,0,0,      15,20,20,20,20,20,  15,15,20,15,15,15,  15,15,0,0,0,0,      20,15,15,15,0,0, #cohort1
               15,15,20,20,15,15,  15,25,15,15,15,15,  20,0,20,15,15,15,   15,15,25,15,15,15,  15,15,15,15,15,15, #cohort2
               15,15,15,15,15,15,  15,15,15,0,0,0,     15,20,25,15,15,15,  15,15,15,0,0,0,     15,15,15,15,15,15,   15,15,15,15,15,15,   15,15,20,15,15,15, #cohort3
               15,15,15,0,0,0,     15,15,25,15,15,25,  15,20,15,15,15,15,  15,15,15,15,NA,NA,  15,15,15,15,15,0,    15,15,25,15,15,15,   15,15,15,15,0,0, #cohort4
               15,15,15,0,15,15,   15,15,NA,20,15,15,  15,15,NA,15,20,20,  15,15,20,15,0,0,   15,15,15,15,15,15) #cohort5

animal <- c('90030R', '90030LL', '102247L', '102247R', '102247LR', #cohort 1
            '90046L', '90046R', '107544L', '107544R', '107544LR', #cohort 2
            '94406R', '94406LR', '94406LL', '94407L', '94407R', '94407LR', '94407LL', #cohort 3
            '90032R', '108807L', '108807R', '107545LLR','107545LR','107549R', '107549LL', #cohort 4
            '108829LL', '108829LRR', '108829R', '112935L', '90042LR') #cohort 5

sex <- c( 'M', 'M', 'M', 'M', 'M', 
          'F', 'F', 'M', 'M', 'M', 
          'M', 'M', 'M', 'F', 'F', 'F', 'F',
          'M','M','M','F','F','F','F',
          'F','F','F','F','F') 

virus <- c('AAV5-hTyr-Full','AAV5-hTyr-1:100','AAV5-hTyr-1:100','AAV5-hTyr-Full','AAV5-hTyr-Full', #cohort1
           'AAV5-hTyr-1:100', 'AAV5-EYFP', 'AAV5-EYFP', 'AAV5-hTyr-1:100', 'AAV5-EYFP', #cohort2
           'AAV5-EYFP', 'AAV5-hTyr-Full', 'AAV5-hTyr-1:100', 'AAV5-hTyr-Full', 'AAV5-hTyr-1:100', 'AAV5-EYFP', 'AAV5-hTyr-Full', #cohort3
           'AAV5-hTyr-Full', 'AAV5-hTyr-1:100', "AAV5-EYFP", 'AAV5-hTyr-1:100','AAV5-hTyr-Full','AAV5-EYFP', 'AAV5-hTyr-Full', #cohort4
           'AAV5-hTyr-1:100','AAV5-EYFP','AAV5-EYFP','AAV5-hTyr-Full','AAV5-hTyr-1:100') #cohort5

timepoint <- c('0', '1', '2', '4','10','16')

# Expand grid so that each timepoint is repeated across all animals first, then trials
speed_df <- expand.grid(animal = animal, timepoint = timepoint)

# map sex and virus to each animal based on their index
speed_df$sex <- sex[match(speed_df$animal, animal)]
speed_df$virus <- virus[match(speed_df$animal, animal)]

# Sort by timepoint, then animal, then trial to ensure correct order in df
speed_df <- speed_df[order(speed_df$animal, speed_df$timepoint), ]

# Add speed data
speed_df$max_speed <- max_speed[1:nrow(speed_df)]
speed_df$min_speed <- min_speed[1:nrow(speed_df)]

# Reorder columns 
speed_df <- speed_df[, c('animal', 'sex', 'virus','timepoint', 'min_speed', 'max_speed')]

#View(speed_df)

# drop animals w/ missed injections (94407LL 8.26.25)
speed_df <- speed_df %>%
  filter(animal != '94407LL')

# export to CSV file
write.csv(speed_df, here("Digigait", "Dataframes", "speed_df.csv"))






# Pre-Processing ----------------------------------------------------------
# Indices df 
digigait_df <- digigait_df %>%
  mutate(timepoint = as.numeric(as.character(timepoint))) # set timepoint datatype to chr

# set animal to factor and arrange df to match metadata vector
digi <- digi %>%
  mutate(animal = factor(animal, levels = virus_mapping$animal)) %>% 
  arrange(animal)

#View(digi)

# drop rows with comments -> these indicate paws that were too messy to include in analysis
digi_clean <- digi %>%
  filter(is.na(comments) | comments == "")

# drop animals w/ missed injections (94407LL 8.26.25)
digi_clean <- digi_clean %>%
  filter(animal != '94407LL')

#View(digi_clean)


# Summary Dataframes ------------------------------------------------------

# 1 Speed
maxSpeed_summary <- groupSummary(speed_df, max_speed, "mean_se", virus, timepoint) # virus by time
write.csv(maxSpeed_summary, file = here("Digigait", "Dataframes", "maxSpeed_summary.csv"))

maxSpeed_sexSummary <- groupSummary(speed_df, max_speed, "mean_se", virus, sex, timepoint) # virus by time by sex
write.csv(maxSpeed_sexSummary, file = here("Digigait", "Dataframes", "maxSpeed_sexSummary.csv"))


# 2) Individual paws by belt speed
# subset digigait df to speeds
belt_15 <- subset(digi_clean, gait_speed_cm_s == 15) 
belt_20 <- subset(digi_clean, gait_speed_cm_s == 20)
belt_25 <- subset(digi_clean, gait_speed_cm_s == 25)

# write each to csv
write.csv(belt_15, here("Digigait", "Dataframes", "Indices_15cms", "15cm_s_indices.csv"))

# subset each into all 4 paws to get summary dfs : each belt speed captures a diff. level of performance that varies by animal
LH_15 <- subset(belt_15, limb == 'Left Hind')
RH_15 <- subset(belt_15, limb == 'Right Hind')
LF_15 <- subset(belt_15, limb == 'Left Fore')
RF_15 <- subset(belt_15, limb == 'Right Fore')

LH_20 <- subset(belt_20, limb == 'Left Hind')
RH_20 <- subset(belt_20, limb == 'Right Hind')
LF_20 <- subset(belt_20, limb == 'Left Fore')
RF_20 <- subset(belt_20, limb == 'Right Fore')

LH_25 <- subset(belt_25, limb == 'Left Hind')
RH_25 <- subset(belt_25, limb == 'Right Hind')
LF_25 <- subset(belt_25, limb == 'Left Fore')
RF_25 <- subset(belt_25, limb == 'Right Fore')

# str(belt_15) to view variables of interest for summaries (can add or subtract here)
variables <- c("percent_swing_stride_percent", "percent_brake_stride_percent", "percent_propel_stride_percent", "percent_stance_stride_percent",
               "stride_length_cm", "stride_frequency_steps_s", "absolute_paw_angle_deg", "paw_angle_variability_deg", "stance_width_cm", "sl_var_cm", 
               "sw_var_cm", "step_angle_var_deg", "step_angle_var_deg", "hind_limb_shared_stance_time_s", "percent_shared_stance_percent")

# use gaitSummaries FN to get summaries by each paw: repeat for all speeds
LH_summary_25 <- gaitSummaries(df = LH_25, variables = variables, virus, timepoint) #can add sex here
RH_summary_25 <- gaitSummaries(df = RH_25, variables = variables, virus, timepoint)
LF_summary_25 <- gaitSummaries(df = LF_25, variables = variables, virus, timepoint)
RF_summary_25 <- gaitSummaries(df = RF_25, variables = variables, virus, timepoint)

# export each to CSV file if desired
write.csv(RF_summary_25, here("Digigait", "Dataframes", "Indices_25cms", "RF_25_virusSummary.csv"))



# 3) Hind paws together, then fore paws together per speed
  # 15 cm/s
hind_15 <- subset(belt_15, limb == 'Left Hind' | limb == 'Right Hind')
fore_15 <- subset(belt_15, limb == 'Left Fore' | limb == 'Right Fore')

hind_summary_15 <- gaitSummaries(df = hind_15, variables = variables, virus, timepoint)
fore_summary_15 <- gaitSummaries(df = fore_15, variables = variables, virus, timepoint)

# export to CSV file
write.csv(hind_summary_15, here("Digigait", "Dataframes", "Indices_15cms", "hindpaws_15_Summary.csv"))
write.csv(fore_summary_15, here("Digigait", "Dataframes", "Indices_15cms", "forepaws_15_Summary.csv"))

  # 20cm/s
hind_20 <- subset(belt_20, limb == 'Left Hind' | limb == 'Right Hind')
fore_20 <- subset(belt_20, limb == 'Left Fore' | limb == 'Right Fore')

hind_summary_20 <- gaitSummaries(df = hind_20, variables = variables, virus, timepoint)
fore_summary_20 <- gaitSummaries(df = fore_20, variables = variables, virus, timepoint)

# export each to CSV file
write.csv(hind_summary_20, file = here("Digigait", "Dataframes", "Indices_20cms", "hindpaws_20_Summary.csv"))
write.csv(fore_summary_20, file = here("Digigait", "Dataframes", "Indices_20cms", "forepaws_20_Summary.csv"))


  #25cm/s
hind_25 <- subset(belt_25, limb == 'Left Hind' | limb == 'Right Hind')
fore_25 <- subset(belt_25, limb == 'Left Fore' | limb == 'Right Fore')

hind_summary_25 <- gaitSummaries(df = hind_25, variables = variables, virus, timepoint)
fore_summary_25 <- gaitSummaries(df = fore_25, variables = variables, virus, timepoint)

# export each to CSV file
write.csv(hind_summary_25, file = here("Digigait", "Dataframes","Indices_25cms", "hindpaws_25_Summary.csv"))
write.csv(fore_summary_25, file = here("Digigait", "Dataframes","Indices_25cms", "forepaws_25_Summary.csv"))




########################### Data Visualization #####################

## Change in max speed plots
speed_df$timepoint <- factor(speed_df$timepoint, levels=c('0','1','2','4','10','16'))

# 1) Facet Plot - treatment group variability 
facet_box <- ggplot(speed_df, aes(x=timepoint, y=max_speed, fill=virus)) +
  geom_boxplot() +
  facet_wrap(~virus) +    #facet plot by virus
  labs(
    title = 'Max Speed Variability',
    x= 'Time from Injection (weeks)',
    y = 'Belt Speed (cm/s)',
    fill = 'Virus Type'
  ) +
  scale_fill_manual(
    values=c('#A3E4D7','#F8C471','#E67E22', 'black'),
    labels = c('AAV5-EYFP','AAV5-hTyr-1:100','AAV5-hTy-Full', 'NA')
  ) +
  scale_x_discrete(breaks = c('0','1','2','4','10','16')      #x-axis labeled as timepoints 
  ) +
  theme_bw() 
facet_box
ggsave(here("Digigait", "Graphs", "beltSpeed_FacetPlot.pdf"),
       plot = facet_box, width = 6, height = 4, dpi = 300)



#2 Trajectory Plot - individual change from baseline / each animal 
# subset and mutate df
delta_speed <- speed_df %>% 
  filter(!is.na(max_speed)) %>% # filter NAs
  group_by(animal, virus) %>%
  mutate(
    baseline =max_speed[timepoint == 0][1], #explicitly state baseline
    data_change = max_speed - baseline, #sub baseline at each timepoint
    timepoint = factor(timepoint)) #ensures timepoints print as exact in x-axis
View(delta_speed)

# create plot        
traj <-  ggplot(delta_speed, aes(x =timepoint, y = data_change, color = virus, group = animal)) +
  geom_line(linewidth = 1.0) +    #Line connecting viral points
  labs(title= 'Individual Change From Baseline : Max Speed',
       x = "Time from Injection (weeks)", 
       y = "Speed (cm/s)", 
       color = "Virus Type", 
       linetype = "Virus Type") +
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


#3) Line Plots - group summary and sex summary data
# max speed by virus
p <- plotGait(data=maxSpeed_summary,
               y_var = 'mean',
               y_error = 'se',
               title = "Max Speed",
               y_label = 'Speed (cm/s)')
ggsave(here("Digigait", "Graphs", "maxSpeed.pdf"), 
       plot = p, width = 6, height = 4, dpi = 300)


# max speed by sex
p2 <- plotGait2(data=maxSpeed_sexSummary,
              y_var = 'mean',
              y_error = 'se',
              title = "Max Speed",
              y_label = 'Speed (cm/s)')
p2
ggsave(here("Digigait", "Graphs", "maxSpeed_bySex.pdf"), 
       plot = p2, width = 6, height = 4, dpi = 300)


#4) various gait variables per speed by paw df
p <- plotGait(
  data = LH_summary_20,
  y_var = "mean_sw_var_cm", 
  y_error = "se_sw_var_cm",
  title = "Stance Width Variability",
  y_label = "Width (cm)"
)

ggsave(here("Digigait", "Graphs", "stanceWidth_20_line.pdf"),
       plot =p , width = 6, height = 4, dpi = 300)


#5) try with violin plot to plot inter-group variability : use paw subset dfs, not summary dfs 

p2 <- violinGait(
  data = LH_20,
  y_var = "percent_shared_stance_percent", 
  title = "Shared Stance: 20cm/s",
  y_label = "% of Swing"
)
ggsave(here("Digigait", "Graphs", "HindShare_15_violin.pdf"),
       plot = p2, width = 6, height = 4, dpi = 300)


# 6) try with box plot (prefer this 5.5.25 KC)
p3_20 <- boxGait(
  data = LH_20,
  y_var = "sw_var_cm", 
  title = "Stance Width Variability",
  y_label = "Width (cm)"
)
ggsave(here("Digigait", "Graphs", "stanceWidth_20_box.pdf"),
       plot = p3_20, width = 6, height = 4, dpi = 300)




########################### Analysis #####################

# summary stats w/ mosaic package
favstats(percent_shared_stance_percent~virus, data = LH_20)


# 1) check ANOVA assumptions
  # max speed uses the df: speed_df
  # all other variables: use subset df by speed of choice

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

# 2-way linear fixed effect model (lmer from lme4)
# mixed model: repeated measures
# virus and timepoint are fixed effect, animal is random effect
# choosing lmer over 2-way anova here because missing data, repeated measures, and non-homogeneous variance b/w groups

# set time as factor
filtered$timepoint <- factor(filtered$timepoint,
                                levels = c("0","1","2","4","10","16"))

# run LMER model = random intercept (animal)
model <- lmer(max_speed ~ virus * timepoint + (1 | animal), data = filtered) # adjust variable for stats on each measure of interest (resting mean, distance mean, etc.)

# ANOVA table for fixed effects
anova(model)


# Post-hoc: Dunnett contrast: compare every virus against the control w/in each timepoint 
  # emmean = Estimated Marginal Means, (Least-Squares Means)
emm <- emmeans(model, ~ virus | timepoint)
dunnet_results <- contrast(emm, method = 'dunnett', ref = 'AAV5-EYFP', adjust = 'holm')
dunnet_results

# NOTE: these stats are different than reported in Prism where a 2 Way ANOVA was applied with Tukey multiple measures correction.
# they're nearly the same outcome but I personally think LMER is more appropriate since it accounts for missing values & individual trajectories over time. (9.10.25)


# Alternatively.... 
# 2) Run 2-Way Repeated-Measures ANOVA
anova_result <- anova_test(
  data= reduced,
  dv = max_speed,
  wid = animal,
  between = c(virus,sex),
  within = timepoint
)
print(anova_result)

# reduce to only significance
anova_results <- as.data.frame(anova_result$ANOVA) 
anova_sig <- anova_results %>% filter(p<0.05) # filter out ns 

# Add stars 
anova_sig$stars <- case_when(
  anova_sig$p < 0.001 ~ "***", 
  anova_sig$p < 0.01 ~ "**",  
  anova_sig$p < 0.05 ~ "*",   
  TRUE ~ ""  # No stars if p >= 0.05
)
anova_sig


# main effect : time
kw_result <- maxSpeed_summary %>%
  kruskal_test(mean~timepoint)
print(kw_result)

# interaction : time x virus
kw_interaction <- maxSpeed_summary %>%
  kruskal_test(mean ~ interaction(virus, timepoint))
print(kw_interaction)

# 3) Dunn posthoc on ungrouped df
# max speed : speed_df
dunn_main <- dunn_test(max_speed ~ timepoint, data = reduced, p.adjust.method = 'holm')
dunn_main

# interaction: create col for Dunn test
reduced$interaction_var <- interaction(reduced$virus, reduced$timepoint)
posthoc <- reduced %>%
  dunn_test(max_speed ~ interaction_var, p.adjust.method = 'holm')
test_type <- "Non-parametric (Dunn’s test with sidak correction)"
View(posthoc)




