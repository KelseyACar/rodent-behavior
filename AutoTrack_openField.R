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


################################ Summary Dataframes ###############################
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


#3) Using FN: animalSums, get sum of all variable per animal (average speed handled outside of function)
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


################################ Analysis Prep ################################


# 1) Define variables to test in a loop
  # animal and sex level = data_animal
  # zone = string defined, NULL if not running TIC or TIE


model_targets <- list(
  list(var = "resting_s_total", data = animal_sums, label = "Resting_Time"),
  list(var = "rearing_events_z_axis_total", data = animal_sums, label = "Rearing_Events" ), 
  list(var = "rearing_time_z_s_total", data = animal_sums, label = "Rearing_Time"), 
  list(var = "distance_m_total", data = animal_sums, label = "Distance"), 
  list(var = "avg_speed_total", data = animal_sums, label = "Avg_Speed"),
  list(var = "total_sum", data = TIZ_animal, label = "Time_in_Center", zone = "Center")
)


# 2) Define a helper fn to prep a df for a given model target
model_prep <- function(df, var, zone_filter = NULL) {
  
  # filter to zone if specified
  if (!is.null(zone_filter)) {
    df <- df %>% filter(zone == zone_filter)
  }
  
  # drop NAs for variable and model terms to prevent errors/ not running
  df <- df %>%
    filter(!is.na(.data[[var]]), !is.na(virus), !is.na(timepoint), !is.na(animal), !is.na(sex))
  
  # set factors and model reference levels
  df$timepoint <- relevel(factor(df$timepoint,
                                 levels = c("0", "1", "2", "4", "10", "16")), ref = "0")
  df$virus <- relevel(factor(df$virus), ref = "AAV5-EYFP")
  df$sex <- factor(df$sex)
  
  # convert to clean df
  df <- as.data.frame(df)
  
  return(df)
}



# Outlier check (raw data, pre-model) -------------------------------------------------------
for (target in model_targets) {
  
  var <- target$var
  df <- target$data
  label <- target$label

  # across all groups (outliers expected due to viral effect)
  outliers_all <- identify_outliers(df, !!sym(var))
  cat("Outliers across all groups:\n")
  print(outliers_all)
  
  # within each viral group separately (more meaningful check)
  for (v in unique(df$virus)) {
    df_v <- df %>% filter(virus == v)
    outliers_v <- identify_outliers(df_v, !!sym(var))
    if (nrow(outliers_v) > 0) {
      cat("\nOutliers in", v, ":\n")
      print(outliers_v)
    }
  }
}


# Analysis ----------------------------------------------------------------


# Fit Models  ------------------------------------------

# Fit a Linear Mixed Effects Model (lmer from lme4)
  # mixed model: repeated measures
  # treatment and timepoint are fixed effect, animal is random effect
  # choosing lmer over 2-way anova here because missing data, repeated measures
  # NOTE: independence is satisfied by design, the random effect on animal accounts for within-subject correlation across timepoints. No indep. assumption test needed.
  # help("pvalues",package="lme4") # pulls up R help page for model summary interpretation


# store fitted models for use in inference loop
fitted_models <- list()

for (target in model_targets) {
  var <- target$var
  label <- target$label
  
  cat("MODEL:", label, "\n")
  
  # prep dfs
  df_clean <- as.data.frame(model_prep(target$data, var, target$zone))
   
  
  # build model formulas
  f_treatment <- as.formula(paste(var, "~ treatment * timepoint + (1|animal)"))
  f_sex <- as.formula(paste(var, "~ treatment * timepoint * sex + (1|animal)"))
  
  # fit model 
  m_treatment <- lmer(f_treatment, data = df_clean)
  m_sex <- lmer(f_sex, data = df_clean)
  
  cat("treatment model fit:", label, "\n")
  cat("Sex model fit:", label, "\n")
  
  # store both models and metadata
  fitted_models[[label]] <- list(
    m_treatment = m_treatment,
    m_sex = m_sex,
    var = var,
    label = label
  )
  }


# Summaries & Fixed Effects -------------------------------------------------
for (entry in fitted_models) {
  cat("SUMMARY:", entry$label, "\n")
  
  cat("\n ", entry$label, "treatment Model Summary \n")
  print(summary(entry$m_treatment))
  cat("\n ", entry$label, "treatment Fixed Effects (F-table) \n")
  print(anova(entry$m_treatment, ddf = "Satterthwaite"))
  
  cat("\n ", entry$label, "Sex Model Summary \n")
  print(summary(entry$m_sex))
  cat("\n ", entry$label, "Sex Fixed Effects (F-table) \n")
  print(anova(entry$m_sex, ddf = "Satterthwaite"))
}



# Model Assumption Checks -------------------------------------------------
for (entry in fitted_models) {
  cat("ASSUMPTIONS:", entry$label, "\n")

  for (m in list(list(mod = entry$m_treatment, tag = "treatment"),
                 list(mod = entry$m_sex,  tag = "sex"))) {
    
    cat("\n---", entry$label, m$tag, "---\n")
    
    # Normality of residuals
    sw_resid <- shapiro.test(residuals(m$mod))
    p_resid  <- ifelse(sw_resid$p.value < 0.0001, "<0.0001", 
                       round(sw_resid$p.value, 4))
    cat("Shapiro-Wilk residuals: W =", round(sw_resid$statistic, 4),
        "p =", p_resid, "\n")
    
    # Normality of random effects
    ranef_vals <- ranef(m$mod)$animal[[1]]
    if (length(unique(ranef_vals)) == 1) {
      cat("Shapiro-Wilk random effects: SKIPPED — singular fit,",
          "all BLUPs identical (random effect variance = 0)\n")
    } else {
      sw_ranef <- shapiro.test(ranef_vals)
      p_ranef  <- ifelse(sw_ranef$p.value < 0.0001, "<0.0001",
                         round(sw_ranef$p.value, 4))
      cat("Shapiro-Wilk random effects: W =", round(sw_ranef$statistic, 4),
          "p =", p_ranef, "\n")
    }
    
    # Q-Q plots
    par(mfrow = c(1,2))
    qqnorm(residuals(m$mod), main = paste(entry$label, m$tag, "residuals Q-Q"))
    qqline(residuals(m$mod))
    qqnorm(ranef_vals, main = paste(entry$label, m$tag, "random effects Q-Q"))
    qqline(ranef_vals)
    par(mfrow = c(1,1))
    
    # Homoscedasticity
    print(plot(m$mod, type = c("p","smooth"),
               main = paste(entry$label, m$tag, "residuals vs fitted")))
  }
}



# Model Inference ---------------------------------------------------------
for (entry in fitted_models) {
  cat("Inference:", entry$label, "\n")
  
  # treatment model — compare each treatment to reference at each timepoint
  emm_treatment <- emmeans(entry$m_treatment, ~ treatment * timepoint)
  comp_treatment <- pairs(emm_treatment, by = "timepoint", adjust = "holm") # Holm's p adj
  comp_treatment_df <- as.data.frame(summary(comp_treatment, infer = c(TRUE, TRUE))) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
    mutate(p.value = ifelse(p.value < 0.0001, "<0.0001", as.character(p.value)))
  
  write_xlsx(comp_treatment_df,
             path = here("Locomotion", "Dataframes", "Linear Mixed Effects results",
                         paste0(entry$label, "_treatment_LME_results.xlsx")))
  cat("Saved:", entry$label, "treatment results\n")
  
  # Sex model — compare M vs F within each treatment at each timepoint
  emm_sex  <- emmeans(entry$m_sex, ~ sex * timepoint | treatment)
  comp_sex <- pairs(emm_sex, by = c("treatment", "timepoint"), adjust = "holm")
  comp_sex_df <- as.data.frame(summary(comp_sex, infer = c(TRUE, TRUE))) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
    mutate(p.value = ifelse(p.value < 0.0001, "<0.0001", as.character(p.value)))
  
  write_xlsx(comp_sex_df, path = here("Locomotion", "Dataframes", "Linear Mixed Effects results", paste0(entry$label, "_sex_LME_results.xlsx")))
  cat("Saved:", entry$label, "sex results\n")
}


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




