### Source File for my Autotrack functions ###

# LIBRARIES REQUIRED #
pacman::p_load('pacman', 'ggplot2','mosaic', 'dplyr', 'ggridges', 'plotly', 'hrbrthemes','tidyverse','patchwork','ggstatsplot','rstatix','dunn.test','FSA','ggprism','tidyr', 'janitor')


############################ Open Field Locomotion  FNs #############################################

##   1  dropErrors FN : drop any openfield recording erros
# df - your df
# id_column - the column housing a specific individual 
# - specific_id - the specific individual who needs measures dropped
# ... - list of measures (columns) to drop for specific individual

dropErrors <- function(df, id_column, specific_id, ...) { 
  
  # set col names to drop
  drop_measures <- c(...)
  
  # drop measurements for specific animal_id
  edited_df <- df %>%
    mutate(across(all_of(drop_measures), ~if_else(.data[[id_column]] == specific_id, NA_real_, .)))
  return(edited_df)
}

# To Call
# df <- dropErrors(df, 'animal','90030R', "measure1", "measure2")
# View(df) 



##   2   baselineSums FN
# df - df of choice
# group - grouping factor (individuals)
# timepoint - baseline == 0 
# bin_col - session bin column
# bin_range - session range length

baselineSums <- function(df, group, timepoint, bin_col = NULL, bin_range = NULL, ...){
  
  # set columns to sum
  sum_cols <- c(...)
  
  sum_df <- df %>%
    filter(!!sym(timepoint) == 0) %>% # sum only timepoint = 0 data
    
    # Apply bin filtering if given in FN call
    { if (!is.null(bin_col) & !is.null(bin_range)) {
          filter(., .data[[bin_col]] %in% bin_range)
          } else {
            .
          }  
      } %>%
    
    #dynamically reference grouping column
    group_by(.data[[group]]) %>%  
    
    # sum across selected columns
    summarise(across(all_of(sum_cols), ~sum(.x,na.rm = TRUE), .names = '{.col}_sum'), .groups = "drop")
  
  return(sum_df)
}

# To call 
# pre_injection_30min <- baselineSums(
  # df = loco, 
  # group = "animal", 
  # timepoint = "timepoint",
  # bin_col = "Bin #", bin_range = 1:6,
  # "Ambulatory (s)", "Resting (s)", "Rearing Events (Z-Axis)", 
  # "Rearing Time Z (s)", "Distance (cm)"
# )





##    3  timeInZones FN   : Calculates total time in each zone by summing appropriate cols
## df of choice
# ... means mutable group_by variables
# bin_col means session bin column
# bin_range means session range length
# uses dplyr

timeInZones <- function(df, bin_col, bin_range, ..., na.rm = TRUE) {
  # Add TIZ column
  df <- df %>%
    mutate(
      TIZ = ambulatory_s + stereotypic_s + resting_s,
      !!bin_col := as.numeric(!!sym(bin_col))  # Ensure bin_col is numeric
    )
  
  # Convert bin_col to symbol
  bin_col_sym <- sym(bin_col)

  # Check grouping
  grouped_df <- df %>% group_by(...)
  
  # Compute grouped summaries
  result <- grouped_df %>%
    summarize(
      total_sum = sum(TIZ, na.rm = na.rm),
      total_mean = mean(TIZ, na.rm = na.rm),
      total_sem = sd(TIZ, na.rm = na.rm) / sqrt(n()),
      
      bin_sum = sum(if_else(!!bin_col_sym %in% bin_range, TIZ, 0), na.rm = na.rm),
      bin_mean = mean(if_else(!!bin_col_sym %in% bin_range, TIZ, NA_real_), na.rm = na.rm),
      bin_sem = ifelse(sum(!!bin_col_sym %in% bin_range, na.rm = TRUE) > 0, 
                      sd(if_else(!!bin_col_sym %in% bin_range, TIZ, NA_real_), na.rm = na.rm) / 
                        sqrt(sum(!!bin_col_sym %in% bin_range, na.rm = TRUE)), 
                      NA),
      .groups = 'drop'
    )
  
  print(result)
  return(result)
}


# To Call
# timeInZones(df, bin_col= "Bin #", bin_range = 1:6,"virus", "timepoint", "Zone")





## 4 animalSums. : sums variable totals / animal by summing zones per variable
                      # does not include average speed as it needs to be handled separately
# df - df of choice
# zone_col - zone column by name
# zones - zones in df
# bin_col - bin column by name
# bin_range - bin range for subsetting session
# variable - variable summing
# ... - group_by factors

animalSums <- function(df, zone_col, zones, bin_col, bin_range, variables, ...){
  # Set grouping variables
  grouping_vars <- enquos(...)
  # Find sum of variables
  sum_result <- df %>%
    group_by(!!!grouping_vars) %>%
    summarise(across(all_of(variables), ~ sum(.x[df[[zone_col]] %in% zones], na.rm = TRUE), .names = "{.col}_total"),
              across(all_of(variables), ~ sum(.x[df[[bin_col]] %in% bin_range & df[[zone_col]] %in% zones], na.rm = TRUE), .names = "{.col}_bin"),
              .groups = "drop")
  return(sum_result)
}


# To Call
# variables <- c(vector of variables)
# summed_data <- animalSums(
  # df = loco, 
  # zone_col = "Zone", 
  # zones = c("A", "B"), 
  # bin_col = "Bin #", 
  # bin_range = 1:12, 
  # variables = sum_variables,  # Pass list of variables directly
  # animal, virus, timepoint
# )




## 5 baselineChange

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
# new_df <- baselineChange(df, "timepoint", baseline_timepoint = 0, virus)





## 6 plot_openField FN : base structure of all openfield locomotion plots ---> make edits here
# data - df of choice (summary df)
# y_var - y variable
# y_error - y variable SE
# y_label - y-axis label
# title - graph title
    ## To make more flexible, add x_var and group_var to FN ---> it is currently set for this project specifically 
library(ggplot2)

plot_openField <- function(data, y_var, y_error, title, y_label, color_palette = NULL) {
  
  #Default color palette if non provided
  if(is.null(color_palette)) {
    color_palette <- c('AAV5-EYFP' = '#A3E4D7', 'AAV5-hTyr-1:100' = '#F8C471', 'AAV5-hTyr-Full' = '#E67E22')
  }
  
  # Base Plot
  p <- ggplot(data, aes(x= timepoint, y = .data[[y_var]], group = virus, color = virus)) +
    geom_line(linewidth=1.3) +
    geom_errorbar(
      aes(ymin = .data[[y_var]] - .data[[y_error]], 
          ymax = .data[[y_var]] + .data[[y_error]]), #plot error bars
      width = 0.2
    ) +
    scale_color_manual(values = color_palette) +
    scale_x_continuous(
      breaks = c(0,1,2,3,4,10,16)) + 
    labs(
      title = title,
      x = "Time Since Injection (weeks)",
      y = y_label,
      color = "Virus"
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
# ambulatory <- plot_openField(
#   data = aged_summary,
#   y_var = "`Ambulatory (s)_total_mean`", 
#   y_error = "`Ambulatory (s)_total_se`",
#   title = "Ambulatory Behavior",
#   y_label = "Time (s)"
# )


## repeat to plot sex differences
plot_openField2 <- function(data, y_var, y_error, title, y_label, color_palette = NULL) {
  
  #Default color palette if non provided
  if(is.null(color_palette)) {
    color_palette <- c('AAV5-EYFP' = '#A3E4D7', 'AAV5-hTyr-1:100' = '#F8C471', 'AAV5-hTyr-Full' = '#E67E22')
  }
  
  # Base Plot
  p <- ggplot(data, aes(x= timepoint, y = .data[[y_var]], color = virus, linetype = sex, group = interaction(virus,sex))) +
    geom_line(linewidth=1.3) +
    geom_errorbar(
      aes(ymin = .data[[y_var]] - .data[[y_error]], 
          ymax = .data[[y_var]] + .data[[y_error]]), #plot error bars
      width = 0.2
    ) +
    scale_color_manual(values = color_palette) +
    scale_x_continuous(
      breaks = c(0,1,2,3,4,10,16)) + 
    labs(
      title = title,
      x = "Time Since Injection (weeks)",
      y = y_label,
      color = "Virus"
    ) +
    theme_minimal() +
    theme(
      legend.position = 'right',
      plot.title = element_text(hjust = 0.5),     #center title
      panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.8))
  
  print(p)
  return(invisible(p))
}



## 7 plot_change FN : without error bars ---> for baseline change, if desired ---> or fix function to calculate SE but can't figure that out yet 3.19.25
# data - df of choice (summary df)
# y_var - y variable
# y_label - y-axis label
# title - graph title
## To make more flexible, add x_var and group_var to FN ---> it is currently set for this project specifically 

library(ggplot2)

plot_change <- function(data, y_var, title, y_label, color_palette = NULL) {
  
  #Default color palette if non provided
  if(is.null(color_palette)) {
    color_palette <- c('AAV5-EYFP' = '#A3E4D7', 'AAV5-hTyr-1:100' = '#F8C471', 'AAV5-hTyr-Full' = '#E67E22')
  }
  
  # Base Plot
  p <- ggplot(data, aes(x= timepoint, y = .data[[y_var]], group = virus, color = virus)) +
    geom_line(linewidth=1.3) +
    scale_color_manual(values = color_palette) +
    scale_x_continuous(
      breaks = c(0,1,2,3,4,10,16)) + 
    labs(
      title = title,
      x = "Time Since Injection (weeks)",
      y = y_label,
      color = "Virus"
    ) +
    theme_minimal() +
    theme(
      legend.position = 'right',
      plot.title = element_text(hjust = 0.5),     #center title
      panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.8))
  
  print(p)
  return(invisible(p))
}


















  
