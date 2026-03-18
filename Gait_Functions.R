# Source File for Digigait and raw gait functions ###


# Libraries Required ------------------------------------------------------
pacman::p_load('pacman', 'ggplot2','mosaic', 'dplyr', 'ggridges', 'plotly', 'hrbrthemes','tidyverse','patchwork','ggstatsplot','rstatix','dunn.test','FSA','ggprism','tidyr', 'janitor')



# 1) groupSummary 
groupSummary <- function(df, dependVar, type = c("mean_sd", "mean_se", "median_iqr"), ...) {
  
  # Set grouping variables
  grouping_vars <- enquos(...)
  
  # Assign type to type
  type <- match.arg(type)
  
  # Get summary stats
  groupSummary_df <- df %>%
    group_by(!!!grouping_vars) %>%
    get_summary_stats( {{ dependVar }}, type = type) %>%
    mutate(across(where(is.numeric), ~round(.x, digits = 4)))
  
  print(groupSummary_df)
  return(groupSummary_df)  # Return the result instead of just printing
}

# To Call
# treatment_summary <- groupSummary(df, time_to_fall, "mean_se", treatment, timepoint)
# View(treatment_summary)

# animal_summary <- groupSummary(df, time_to_fall, "mean_se", animal, treatment, timepoint)
# View(animal_summary)


# 2) gait variable summaries
gaitSummaries <- function(df, variables, ...) {
  grouping_vars <- enquos(...)
  
  summaries <- df %>%
    group_by(!!!grouping_vars) %>%
    summarise(
      n = n_distinct(animal),  # count animals, not rows
      across(all_of(variables), 
             list(
               mean = ~ mean(.x, na.rm = TRUE),
               se   = ~ sd(.x, na.rm = TRUE) / sqrt(n),
               min  = ~ min(.x, na.rm = TRUE),
               max  = ~ max(.x, na.rm = TRUE)
             ), .names = "{.fn}_{.col}"),
      .groups = "drop")
  return(summaries)
}



# To Call
# variables <- c(vector of variables)
# summed_data <- gaitSummaries(
# df = gait, 
# variables = variables,  # Pass list of variables directly
# animal, treatment, timepoint
# )





######### Stance by time summary stats #########

# 1) animal summary stats for single and any combination of paws down (must create binary df first)
# df: dataframe with time + paws as columns containing binary (0|1) paw up / down data
# paws: character vector with paw column names (must match df)

paw_stanceStats <- function(df, paws = c("lf", "rf", "lr", "rr")) {
  df$sum_down <- rowSums(df[, paws] == 0)
  
  cv <- function(x) {
    if (mean(x) == 0) return(NA)
    sd(x) / mean(x)
  }
  
  combos <- list(
    both_fore = (df$lf==0) & (df$rf==0),
    both_hind = (df$lr==0) & (df$rr==0),
    both_left = (df$lf==0) & (df$lr==0),
    both_right = (df$rf==0) & (df$rr==0),
    contralateral_lf = (df$lf==0) & (df$rr==0),
    contralateral_rf = (df$rf==0) & (df$lr==0),
    exactly_three = (df$sum_down == 3)
  )
  
  # Compute stats for each combination
  results <- lapply(combos, function(x) {
    sum_x <- sum(x)
    mean_x <- mean(x)
    cv_x <- cv(as.numeric(x))
    duration <- sum_x * 0.01  # Convert sum to seconds
    c(sum = sum_x, mean = mean_x, CV = cv_x, duration_sec = duration)
  })
  
  # Combine results into a tidy dataframe
  results_df <- bind_rows(results, .id = "combo")
  
  return(results_df)
}

# To Call:
#stats_df <- paw_stats(wide_df)
#print(stats_df)



# ACF summary stats -------------------------------------------------------


## To process a long format ("processed") per paw step count and duration files and collect CV of durations for all paws, fore paws, and hind paws
process_file <- function(filepath, speed_label) {
  df <- read.csv(filepath)
  
  # ID metadata
  meta <- df[1, c("animal","treatment","timepoint")] %>% as.list()
  
  # calculate step duration
  df <- df %>%
    mutate(duration = end - start)
  
  # CV per paw
  paw_cv <- df %>%
    group_by(paw) %>%
    summarise(
      mean_duration = mean(duration, na.rm=TRUE),
      sd_duration   = sd(duration, na.rm=TRUE),
      cv = sd_duration / mean_duration,
      .groups="drop"
    )
  
  # combine to fore/hind/total
  fore_cv  <- mean(paw_cv$cv[paw_cv$paw %in% c("lf","rf")], na.rm=TRUE)
  hind_cv  <- mean(paw_cv$cv[paw_cv$paw %in% c("lr","rr")], na.rm=TRUE)
  total_cv <- mean(paw_cv$cv[paw_cv$paw %in% c("lf","rf","lr","rr")], na.rm=TRUE)
  
  tibble(
    id = meta$animal,
    treatment = meta$treatment,
    timepoint = meta$timepoint,
    speed = speed_label,
    fore_cv = fore_cv,
    hind_cv = hind_cv,
    total_cv = total_cv
  )
}

## To Call
# cv_df <- process_file(file_path$label[1], speed_label = "15cms") #must index file_path bc it is a list
# Then loop through all files (in rawGait.R script using this FN)

# 1) stats on ACF_values 
acf_perPawSummary <- function(acf_df){
  acf_df %>%
    mutate(
      acf_lag1 = map_dbl(acf_vals, ~ .x[2]),
      acf_decay = map_dbl(acf_vals, ~ which(.x < 0.5)[1] * 0.01),
      acf_auc = map_dbl(acf_vals, ~ sum(.x[-1])),
      acf_peak_lag = map_dbl(acf_vals, ~ {
        peak_idx <- which(diff(sign(diff(.x))) == -2)
        if (length(peak_idx) > 0) peak_idx[1] * 0.01 else NA
      })
    )
}
# To Call
# acf_perPawSummary(acf)  #on any acf file, then copy & paste to a csv



# 2) This function utilizes the data frame built from the above ACF per paw function (acf_animalSummary)
# must first copy & paste each animals acf_animalSummary
# Results are grouped metrics (auc, decay, lag1, peak lag): determine grouping by replacing ... w/ group factors

acf_groupSummary <- function(df, ...) {
  df %>%
    group_by(...) %>%
    summarise(
      # ---- ALL PAWS ----
      total_mean_auc = mean(acf_auc,na.rm = TRUE),
      total_se_auc = sd(acf_auc)/sqrt(sum(!is.na(acf_auc))),
      total_mean_lag1 = mean(acf_lag1, na.rm = TRUE),
      total_se_lag1 = sd(acf_lag1)/sqrt(sum(!is.na(acf_lag1))),
      total_mean_decay = mean(acf_decay, na.rm = TRUE),
      total_se_decay = sd(acf_decay)/sqrt(sum(!is.na(acf_decay))),
      total_mean_peaklag = mean(acf_peak_lag, na.rm = TRUE),
      total_se_peaklag = sd(acf_peak_lag)/sqrt(sum(!is.na(acf_peak_lag))),
      
      # ---- FORE PAWS (lf, rf) ----
      fore_mean_auc = mean(acf_auc[paw %in% c("lf","rf")], na.rm = TRUE),
      fore_se_auc = sd(acf_auc[paw %in% c("lf","rf")]) / sqrt(sum(!is.na(acf_auc[paw %in% c("lf","rf")]))),
      fore_mean_lag1 = mean(acf_lag1[paw %in% c("lf","rf")], na.rm = TRUE),
      fore_se_lag1 = sd(acf_lag1[paw %in% c("lf","rf")]) / sqrt(sum(!is.na(acf_lag1[paw %in% c("lf","rf")]))),
      fore_mean_decay = mean(acf_decay[paw %in% c("lf","rf")], na.rm = TRUE),
      fore_se_decay = sd(acf_decay[paw %in% c("lf","rf")]) / sqrt(sum(!is.na(acf_decay[paw %in% c("lf","rf")]))),
      fore_mean_peaklag = mean(acf_peak_lag[paw %in% c("lf","rf")], na.rm = TRUE),
      fore_se_peaklag = sd(acf_peak_lag[paw %in% c("lf","rf")]) / sqrt(sum(!is.na(acf_peak_lag[paw %in% c("lf","rf")]))),
      
      # ---- HIND PAWS (lr, rr) ----
      hind_mean_auc = mean(acf_auc[paw %in% c("lr","rr")], na.rm = TRUE),
      hind_se_auc = sd(acf_auc[paw %in% c("lr","rr")]) / sqrt(sum(!is.na(acf_auc[paw %in% c("lr","rr")]))),
      hind_mean_lag1 = mean(acf_lag1[paw %in% c("lr","rr")], na.rm = TRUE),
      hind_se_lag1 = sd(acf_lag1[paw %in% c("lr","rr")]) / sqrt(sum(!is.na(acf_lag1[paw %in% c("lr","rr")]))),
      hind_mean_decay = mean(acf_decay[paw %in% c("lr","rr")], na.rm = TRUE),
      hind_se_decay = sd(acf_decay[paw %in% c("lr","rr")]) / sqrt(sum(!is.na(acf_decay[paw %in% c("lr","rr")]))),
      hind_mean_peaklag = mean(acf_peak_lag[paw %in% c("lr","rr")], na.rm = TRUE),
      hind_se_peaklag = sd(acf_peak_lag[paw %in% c("lr","rr")]) / sqrt(sum(!is.na(acf_peak_lag[paw %in% c("lr","rr")]))),
      .groups = "drop"
    )
}

# To Call:
# acf_animalStats <- acf_groupSummary(acf_perPaw df, animal, timepoint)
# acf_treatmentStats <- acf_groupSummary(acf_perPaw df, animal, treatment, timepoint)
# acf_sexStats <- acf_groupSummary(acf_perPaw df, animal, sex, treatment, timepoint)




# Plots -------------------------------------------------------------------


## baselineChange

baselineChange <- function(df, time_var, baseline_timepoint, ...) {
  # Set grouping variables
  grouping_vars <- enquos(...)
  
  # Identify all columns that contain "mean_" (for mean values) 
  mean_cols <- grep("(_total_mean|_bin_mean)$", names(df), value = TRUE)
  
  # Calculate change from baseline
  results <- df %>%
    group_by(!!!grouping_vars) %>%
    mutate(across(
      .cols = all_of(mean_cols),
      .fns = ~ .[!!sym(time_var) == baseline_timepoint] - .,
      .names = "{.col}_delta"
    )) %>%
    ungroup()
}

# To Call
# new_df <- baselineChange(df, "timepoint", baseline_timepoint = 0, treatment)



# Line plot for gait variables --------------------------------------------

# 1) plot_gait FN : base structure of all gait plots ---> make edits here
# data - df of choice (summary df)
# y_var - y variable
# y_error - y variable SE
# y_label - y-axis label
# title - graph title

# line plots for treatment groups
plotGait <- function(data, y_var, y_error, title, y_label, color_palette = NULL) {
  
  #Default color palette if non provided
  if(is.null(color_palette)) {
    color_palette <- c('A' = '#A3E4D7', 'B' = '#F8C471', 'C' = '#E67E22')
  }
  
  # Base Plot
  p <- ggplot(data, aes(x= timepoint, y = .data[[y_var]], group = treatment, color = treatment)) +
    geom_line(linewidth=1.3) +
    geom_errorbar(
      aes(ymin = .data[[y_var]] - .data[[y_error]], 
          ymax = .data[[y_var]] + .data[[y_error]]), #plot error bars
      width = 0.2
    ) +
    scale_color_manual(values = color_palette) +
    scale_x_discrete(
      breaks = c(0,1,2)) + 
    labs(
      title = title,
      x = "Time Since Treatment",
      y = y_label,
      color = "treatment"
    ) +
    theme_minimal() +
    theme(
      legend.position = 'right',
      plot.title = element_text(hjust = 0.5),     #center title
      panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.8))
  
  print(p)
  return(invisible(p))
}
# To Call
# swing <- plotGait(  
#   data = left_summary,
#   y_var = "swing_cm_mean`", 
#   y_error = "swing_cm_se`",
#   title = "Gait Swing",
#   y_label = "Swing (cm)"
# )


# 2) line plots by sex
plotGait2 <- function(data, y_var, y_error, title, y_label, color_palette = NULL) {
  
  #Default color palette if non provided
  if(is.null(color_palette)) {
    color_palette <- c('A' = '#A3E4D7', 'B' = '#F8C471', 'C' = '#E67E22')
  }
  
  # Base Plot
  p <- ggplot(data, aes(x= timepoint, y = .data[[y_var]], linetype = sex, color = treatment, group = interaction(treatment, sex))) +
    geom_line(linewidth=1.3) +
    geom_errorbar(
      aes(ymin = .data[[y_var]] - .data[[y_error]], 
          ymax = .data[[y_var]] + .data[[y_error]]), #plot error bars
      width = 0.2
    ) +
    scale_color_manual(values = color_palette) +
    scale_x_discrete(
      breaks = c(0,1,2)) + 
    labs(
      title = title,
      x = "Time Since Treatment",
      y = y_label,
      color = "treatment"
    ) +
    theme_minimal() +
    theme(
      legend.position = 'right',
      plot.title = element_text(hjust = 0.5),     #center title
      panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.8))
  
  print(p)
  return(invisible(p))
}

# 3) violin gait plot
violinGait <- function(data, y_var, title, y_label, color_palette = NULL) {
  
  #Default color palette if non provided
  if(is.null(color_palette)) {
    color_palette <- c('A' = '#A3E4D7', 'B' = '#F8C471', 'C' = '#E67E22')
  }

  # Base Plot
  p <- ggplot(data, aes(x= as.factor(timepoint), y = .data[[y_var]], fill = treatment)) +
    geom_violin(position = 'dodge') +
    scale_fill_manual(values = color_palette) +
    scale_x_discrete(
      breaks = c(0,1,2)) + 
    labs(
      title = title,
      x = "Time Since Treatment",
      y = y_label,
      color = "treatment"
    ) +
    theme_minimal() +
    theme(
      legend.position = 'right',
      plot.title = element_text(hjust = 0.5),     #center title
      panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.8))
  
  print(p)
  return(invisible(p))
}
# To Call 
# swing <- violinGait(  
#   data = LH_df,         # not summary df
#   y_var = "swing_cm_mean`", 
#   title = "Gait Swing",
#   y_label = "Swing (cm)"
# )


# 4) box n whisker plots
boxGait <- function(data, y_var, title, y_label, color_palette = NULL) {
  
  #Default color palette if non provided
  if(is.null(color_palette)) {
    color_palette <- c('A' = '#A3E4D7', 'B' = '#F8C471', 'C' = '#E67E22')
  }
  
  # Base Plot
  p <- ggplot(data, aes(x= as.factor(timepoint), y = .data[[y_var]], fill = treatment)) +
    geom_boxplot() +
    scale_fill_manual(values = color_palette) +
    #geom_jitter(position = 'dodge') +
    scale_x_discrete(
      breaks = c(0,1,2)) + 
    labs(
      title = title,
      x = "Time Since Treatment",
      y = y_label,
      color = "treatment"
    ) +
    theme_minimal() +
    theme(
      legend.position = 'right',
      plot.title = element_text(hjust = 0.5),     #center title
      panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.8))
  
  print(p)
  return(invisible(p))
}
# To Call 
# swing <- boxGait(  
#   data = LH_df,       # not summary df
#   y_var = "swing_cm_mean`", 
#   title = "Gait Swing",
#   y_label = "Swing (cm)"
# )






