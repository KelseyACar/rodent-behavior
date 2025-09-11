# This Rscripts analysis Digigait raw stance_by_time data


# This Script -------------------------------------------------------------
  # Load Libraries
  # Preprocess & Cleaning
    # read stance_by_time files : 1/animal/timepoint
    # remove steps with NA as starting or ending time (removes steps not captured in data)
    # add metadata and colors / paw
    # pivot wider to capture start and stop times of each paw's step
    # save as processed csv
  # Plot Stance Time
    # read in processed file
    # plot vistime
  # Autocorrelation Prep
    # get max session length from processed file & set time sequence range by 0.01s intervals
    # for each paw, create binary df : 0 = paw down, 1 = paw up
    # combine to df & save as binary csv
  # Plot Correlograms
    # read binary csv : 1/animal/timepoint & set paws
    # run ACF fn : auto and cross covariance and correlation
      # convert lag to time
    # plot with ggplot
  # Analysis
    # read in bini


# Load Libraries ----------------------------------------------------------
pacman::p_load('pacman', 'dplyr', 'tidyr', 'janitor', 'readxl', 'here', 'stringr', 'purrr', 'vistime')
source(here::here('Gait_Functions.R'))
here::i_am("rawGait.R")




# Pre-processing & Cleaning ------------------------------------------------

# make step plots, 1 / animal / timepoint: pull from folder, update animal ID and timepoint in file name to read in

df <- read_excel(here('Digigait','stance by time','16wk', '15cms', '108829LRR_15cms_0degUP_16wk_Stance_time_info.xlsm'))
df <- df[-1, ] # remove units row (first row)
df <- clean_names(df)
#colnames(df)
#str(df)
#View(df)

# convert data types to numeric and set decimals to 2 to match raw data xlsm files
df <- df %>%
  mutate(across(everything(), as.numeric)) # replaces '-' with NA as well, which we want. ignore the warning!

# change col names to each paw start and end, id'd by _ 
df <- df %>%
  rename_with(~ gsub("_stance_", "_", .x), .col = -step_number)
#colnames(df)
#View(df)

# pivot df to long format
df_long <- df %>%
  pivot_longer(cols = -`step_number`,
               names_to = c("paw", "phase"),
               names_sep = "_",
               values_to = 'time'
)
  
# clean df to exclude steps for paws that have an NA in either their start or stop 
df_clean <- df_long %>%
  group_by(step_number, paw) %>%
  filter(!any(is.na(time))) %>%
  ungroup()

# add metdata: Update details for each new plot
df_clean <- df_clean %>%
  mutate(
    animal = '108829LRR', 
    timepoint = '16',    
    virus = 'AAV5-eYFP'
  )

# assign paw colors
paw_colors <- c(
  lf = '#ccccff',
  rf = '#ccffff',
  rr = '#99ffcc',
  lr = '#ffccff'
)

# map paw color to df
df_clean <- df_clean %>%
  mutate(color = paw_colors[paw])

# pivot back to wide to set start and stop columns
df_wide <- df_clean %>%
  pivot_wider(
    names_from = phase,
    values_from = time,
  )

# reorder col names
df_wide <- df_wide %>%
  select(animal, virus, timepoint, step_number, paw, start, end, color)
#View(df_wide)

# save file to csv
write.csv(df_wide, file = here::here("Digigait", "stance by time", "processed","15cms", "108829LRR_16_processed.csv"), row.names = FALSE) # update  animal ID & timepoint


# Create dfs for ACF: up-down duration ------------------------------------

# Calculate stats for fore & hind paws from processed files (Run ACF on this to detect decay of correlation): all files must be processed first to loop through
file_path <- list(
  speed_15cms = list.files(path = here::here("Digigait", "stance by time", "processed", "15cms"), full.names= TRUE),
  speed_20cms = list.files(path = here::here("Digigait", "stance by time", "processed", "20cms"), full.names= TRUE))

# FN process_file processes 1 speed file: extracts metdata, calculates paw down duration for each step & paw, calculates CV / paw, combines all paw CV, fore CV, hind CV
single_file <- process_file(file_path$speed_15cms[1], speed_label = "15cms") # must index file_path list and process 15cms & 20cms separately
single_file2 <- process_file(file_path$speed_20cms[2], speed_label = "20cms")

# loop over all files in speed paths & bind to one df
CV_df <-bind_rows(
  map(file_path$speed_15cms, ~process_file(.x, "15cms")),
  map(file_path$speed_20cms, ~process_file(.x, "20cms"))
)

# save CV dfs
write.csv(CV_df %>% filter(speed =="15cms"), 
          here("Digigait", "Dataframes", "autocorrelation", "stepDuration_CV_15cms.csv"), row.names = FALSE)

write.csv(CV_df %>% filter(speed == "20cms"),
          here("Digigait", "Dataframes", "autocorrelation", "stepDuration_CV_20cms.csv"), row.names = FALSE)

# drop misses: (94407LL right now, still comfirming through IHC 8.21.25)
CV_df <- CV_df %>%
  filter(id != '94407LL')

# find group stats for CV
group_cv <- CV_df %>%
  group_by(speed, virus, timepoint) %>%
  summarise(
    n = n(),
    mean_fore_cv = mean(fore_cv, na.rm = TRUE),
    sd_fore_cv = sd(fore_cv, na.rm = TRUE),
    mean_hind_cv = mean(hind_cv, na.rm = TRUE),,
    sd_hind_cv = sd(hind_cv, na.rm = TRUE),,
    mean_total_cv = mean(total_cv, na.rm = TRUE),,
    sd_total_cv = sd(total_cv, na.rm = TRUE)
  )

write.csv(group_cv, here("Digigait", "Dataframes", "autocorrelation", "grouped_CV_stats.csv"), row.names = FALSE)


# Plot Stance Time --------------------------------------------------------
# raw step time plots (with processed files)
df <- read.csv(here::here("Digigait", "stance by time", "processed","15cms", "108829LRR_16_processed.csv")) #select single file by name

# set factor levels for paws so plot output is graphically sensical 
df <- df %>%
  mutate(paw = factor(paw, levels = c('lr', 'lf', 'rf', 'rr')))

# Raw time / step plot
vistime(df,
        col.event = 'step_number',
        col.start = 'start',
        col.end = 'end',
        col.group = 'paw',
        col.color = 'color',
        title = 'eYFP: 16 (108829LRR)',  # update per file
        optimize_y = TRUE,
        show_labels = FALSE
        ) 

# save from viewer window ----->


# Autocorrelation Prep ----------------------------------------------------

# convert to binary for autocorrelation
max(df$end) # get max time / file
time_seq <- seq(0, 4.03, by = 0.01)  # change second value to reflect max time of each df, set to 0.01s intervals

# Unique paws
paws <- unique(df$paw)

# Initialize empty list to store dataframes
df_list <- list()

for (p in paws) {
  df_paw <- df[df$paw == p, ]
  binary_vec <- rep(1, length(time_seq)) # initialize with 1s
  
  for (i in 1:nrow(df_paw)) {
    binary_vec[time_seq >= df_paw$start[i] & time_seq <= df_paw$end[i]] <- 0 # set to 0s if paw is down during start to end time sequence
  }
  
  # Create a dataframe for this paw
  temp_df <- data.frame(
    time = time_seq,
    paw = p,
    value = binary_vec
  )
  
  df_list[[as.character(p)]] <- temp_df
}

# combine all paws to 1 df
combined_df <- do.call(rbind, df_list)

# reset row names
rownames(combined_df) <- NULL

  
# pivot wider so each paw is a column
binary_df <- combined_df %>%
  select(time, paw, value) %>%
  pivot_wider(names_from = paw, values_from = value, values_fill = 0)


# add metdata: update metadata per file
binary_df <- binary_df %>%
  mutate(
    animal = '108829LRR',
    timepoint = '16',
    virus = 'AAV5-eYFP'
  )

# reorder col names
binary_df <- binary_df %>%
  select(animal, virus, timepoint, lf, rf, lr, rr)
#View(binary_df)

# save file to csv
write.csv(binary_df, here("Digigait", "stance by time", "binary","15cms", "108829LRR_16wk_binary.csv"), row.names = FALSE)







# Run ACF & Plot Correlograms -------------------------------------------------------

# Doing 1 file at a time and pasting to a csv because I dont have time to make batch process code at the moment - 6.2.25 (and it's hard) (for me)
df <- read.csv(here::here("Digigait", "stance by time", "binary","20cms", "90042LR_16wk_binary.csv")) #select single file by name

paws <- c('lf','rf','lr','rr')

# run ACF - Auto and Cross Covariance & Correlation FN 

# lag_max determines lags to compute, lag 1 = 0.01sec (my time intervals), 100 = 1 second later, i.e. how far it goes out
acf <- df %>%
  pivot_longer(cols = all_of(paws), names_to = "paw", values_to = "binary") %>%
  group_by(paw) %>%
  summarise(acf_vals = list(acf(binary, plot = FALSE, lag.max = 100)$acf),
            .groups = "drop")

# Un-nest ACF values for plotting
acf_result <- acf %>%
  mutate(lag = map(acf_vals, ~ seq_along(.x) - 1),
         acf = acf_vals) %>%
  select(paw, lag, acf) %>%
  unnest(c(lag, acf))


# Plot correlogram
data <- acf_result %>%
  mutate(lag_time = lag * 0.01) # convert lag to seconds


p <- ggplot(data, aes(x = lag_time, y = acf, color = paw)) +
  geom_line() +
  labs(title = "Full: 16wk (108829LL)", x = "Lag", y = "ACF")
p

# save from viewer window ----> 

# Build Stats & Summary Dataframes ----------------------------------------------------------------

### Autocorrelation Stats

# 1) ACF per paw animal stats w/ acf_perPawSummary -> Copy & paste to build csv (acf_pawStats_15cms.csv or acf_pawStats_20cms.csv)
# lag1, auc, peak lag, decay per paw/animal/timepoint
stats_acf <- acf_perPawSummary(acf)
View(stats_acf)


# 2) ACF per animal stats w/ acf_perAnimalSummary
# lag1, auc, peak lag, decay across all paws, fore paw, and hind paws only per animal, virus, or sex

# read in df made from step 1
acf_df <- read.csv(here::here("Digigait", "Dataframes","autocorrelation", "15cms", "acf_pawStats_15cms.csv")) 

# filter out animals with comments (94407LL = likely a miss)
acf_df <- acf_df %>%
  filter(comments == '')

# 2.1) per animal
acf_animal_stats <- acf_groupSummary(acf_df, animal, virus, timepoint) 
write.csv(acf_animal_stats, here("Digigait", "Dataframes","autocorrelation", "15cms", "acf_animalSummary_15cms.csv"))


# 2.2) per virus
acf_virus_stats <- acf_groupSummary(acf_df, virus, timepoint)
write.csv(acf_virus_stats, here("Digigait", "Dataframes","autocorrelation", "15cms", "acf_virusSummary_15cms.csv"))

# 2.3) per sex / virus
acf_sex_stats <- acf_groupSummary(acf_df, virus, sex, timepoint)
write.csv(acf_sex_stats, here("Digigait", "Dataframes","autocorrelation", "15cms", "acf_sexSummary_15cms.csv"))



# 3) Stance Time Per Paw stats w/ paw_stats function -> Copy to build csv
df2 <- read.csv(here("Digigait", "stance by time", "binary","20cms", "90042LR_16wk_binary.csv")) #select single file by name

paws <- paw_stanceStats(df2)
View(paws) # copy to csv


### Stance Stats

# 4) Stance time per paw grouped stats
# read in df created from step 3
paws_df <- read.csv(here("Digigait", "Dataframes", "stance_time", "15cms", "stancePerPaw_15cms.csv"))

# filter out animals with comments (94407LL = likely a miss)
paws_df <- paws_df %>%
  filter(comments == '')


# filter df by paw measure desired (couldn't get function to calculate mean & se across all paw_combos to work)
  # variables: both_fore, both_hind, both_left, both_right, contralateral_lf, contralateral_rf, exactly_three 
paw_set <- paws_df %>% 
  filter(paw_combo == "both_hind")
#paw_set

paws_summary <- groupSummary(paw_set, duration_sec, "mean_se", virus, timepoint) # adjust group by variables as needed

# save file to csv
write.csv(paws_summary, file = here::here("Digigait", "Dataframes","stance_time", "15cms","summaries", "hind_virusSummary_15cms.csv"), row.names = FALSE) # name by paw combo chosen







# plots -------------------------------------------------------------------
# using plot fns made for gait analysis 
color_palette <- c('AAV5-eYFP' = '#A3E4D7', 'AAV5-hTyr-1:100' = '#F8C471', 'AAV5-hTyr-Full' = '#E67E22')


# line plot w/ summary stats
fore <- ggplot(paws_group_stats, aes(x= timepoint, y = mean, group = virus, color = virus)) +
  geom_line(linewidth=1.3) +
  geom_errorbar(
    aes(ymin = mean - se, 
        ymax = mean + se), #plot error bars
    width = 0.2
  ) +
  scale_color_manual(values = color_palette) +
  scale_x_discrete(
    breaks = c(0,1,2,3,4,10,0)) + 
  labs(
    title = "Mean Fore Paws",
    x = "Time Since Injection (weeks)",
    y = "",
    color = "Virus"
  ) +
  theme_minimal() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),     #center title
    panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.8))
peak.line
ggsave(here::here("Digigait", "Graphs", "stance by time", "ACF_peakLagPlot.pdf"),
       plot = peak.line, width = 6, height = 4, dpi = 300)










# box plot (if we want them)
lag1 <- ggplot(acf_df, aes(x= as.factor(timepoint), y = acf_lag1, fill = virus)) +
  geom_boxplot() +
  scale_fill_manual(values = color_palette) +
  scale_x_discrete(
    breaks = c(0,1,2,3,4,10,0)) + 
  labs(
    title = "Lag 1",
    x = "Time Since Injection (weeks)",
    y = "ACF at Lag 1",
    color = "Virus"
  ) +
  theme_minimal() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),     #center title
    panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.8))
lag1
ggsave(here::here("Digigait", "Graphs", "stance by time", "ACF_lag1Plot.pdf"),
       plot = lag1, width = 6, height = 4, dpi = 300)


decay <- ggplot(acf_df, aes(x= as.factor(timepoint), y = acf_decay, fill = virus)) +
  geom_boxplot() +
  scale_fill_manual(values = color_palette) +
  scale_x_discrete(
    breaks = c(0,1,2,3,4,10,0)) + 
  labs(
    title = "ACF Decay",
    x = "Time Since Injection (weeks)",
    y = "ACF Decay",
    color = "Virus"
  ) +
  theme_minimal() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),     #center title
    panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.8))
decay
ggsave(here::here("Digigait", "Graphs", "stance by time", "ACF_decayPlot.pdf"),
       plot = decay, width = 6, height = 4, dpi = 300)


auc <- ggplot(acf_df, aes(x= as.factor(timepoint), y = acf_auc, fill = virus)) +
  geom_boxplot() +
  scale_fill_manual(values = color_palette) +
  scale_x_discrete(
    breaks = c(0,1,2,3,4,10,0)) + 
  labs(
    title = "ACF AUC",
    x = "Time Since Injection (weeks)",
    y = "AUC",
    color = "Virus"
  ) +
  theme_minimal() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),     #center title
    panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.8))
auc
ggsave(here::here("Digigait", "Graphs", "stance by time", "ACF_AUCPlot.pdf"),
       plot = auc, width = 6, height = 4, dpi = 300)


peak <- ggplot(acf_df, aes(x= as.factor(timepoint), y = acf_peak_lag, fill = virus)) +
  geom_boxplot() +
  scale_fill_manual(values = color_palette) +
  scale_x_discrete(
    breaks = c(0,1,2,3,4,10,0)) + 
  labs(
    title = "ACF Peak Lag",
    x = "Time Since Injection (weeks)",
    y = "Peak Lag",
    color = "Virus"
  ) +
  theme_minimal() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),     #center title
    panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.8))
peak
ggsave(here::here("Digigait", "Graphs", "stance by time", "ACF_peakLagPlot.pdf"),
       plot = peak, width = 6, height = 4, dpi = 300)
