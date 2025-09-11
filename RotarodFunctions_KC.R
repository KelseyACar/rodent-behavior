############################ Rotarod FNs #############################################

##     groupSummary FN : for one variable 
library(dplyr)
library(rstatix)

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
# virus_summary <- groupSummary(df, time_to_fall, "mean_se", virus, timepoint)
# View(virus_summary)

# animal_summary <- groupSummary(df, time_to_fall, "mean_se", animal, virus, timepoint)
# View(animal_summary)





## rotarod_lineplot FN 
# data - df of choice (summary df)
# y_var - y variable
# y_error - y variable SE
# y_label - y-axis label
# title - graph title
## To make more flexible, add x_var and group_var to FN ---> it is currently set for this project specifically 
library(ggplot2)

line_plot <- function(data, y_var, y_error, title, y_label, color_palette = NULL) {
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
    scale_x_discrete(
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
# aged <- line_plot(
#   data = aged_virus_summary,
#   y_var = "mean", 
#   y_error = "se",
#   title = "Rotarod",
#   y_label = "Time to Fall (s)")




## lineplot grouped by sex

line_plot2 <- function(data, y_var, y_error, title, y_label, color_palette = NULL) {
  #Default color palette if non provided
  if(is.null(color_palette)) {
    color_palette <- c('AAV5-EYFP' = '#A3E4D7', 'AAV5-hTyr-1:100' = '#F8C471', 'AAV5-hTyr-Full' = '#E67E22')
  }
  
  # Base Plot
  p <- ggplot(data, aes(x= timepoint, y = .data[[y_var]], linetype = sex, color = virus, group = interaction(virus,sex))) +
    geom_line(linewidth = 1.3) +
    geom_errorbar(
      aes(ymin = .data[[y_var]] - .data[[y_error]], 
          ymax = .data[[y_var]] + .data[[y_error]]), #plot error bars
      width = 0.2
    ) +
    scale_color_manual(values = color_palette) +
    scale_x_discrete(
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



