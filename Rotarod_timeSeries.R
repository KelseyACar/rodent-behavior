# This Script -------------------------------------------------------------
  # Visualizes and analyzes time series data for Rotarod behavior (would apply to any single response variable) for 3 factors
  # The three factors are: time, sex, treatment
  ## DEPENDENT VARIABLE: Time-to-Fall (seconds), can convert to RPMS at fall with following equation
  # RPMs at Fall: ((ending speed of full session - starting speed) / full session length) x (seconds at fall) + starting speed
  # RPMS for our study: ((40-4) / 300) X (secs at fall) + 4
# Load Libraries 
# Create Data Frames
  # Create from raw, manually collected data (rotarod time to fall: 4 trials / animal across 6 timepoints)
# Summary Data Frames
  # Drop animals that had missed injections after collecting animal summaries but before group summaries
  # Using custom functions, collect animal summaries (mean time to fall / timepoint / animal), treatment group summaries, and sex by treatment group summaries
# Analysis
  # Check assumptions, reduce files to cols of interest
  # LMER from LME4, run summary on model w/ holms post-hoc correction for multiple comparisons (time, treatment)
  # pull CIs and p-values (emmeans)
# Plot w/ ggplot2 


# Load Libraries ---------------------------------------------------------------
pacman::p_load('pacman', 'dplyr','ggplot2', 'here', 'rstatix', 'tidyr','lme4','emmeans') # auto installs required packages not yet installed
source(here::here('RotarodFunctions.R')) # source in custom functions script 
here::i_am("Rotarod_timeSeries.R")

# to collect pack citations, uncomment and run line below:
# cite_packages()

# to collect pack versions, uncomment and run line below:
# report_packages()

# Create Data Frames ----------------------------------------------------

# define vectors (columns) of metadata
# create metadata vectors
animal <- c('1CM', '2AM', '3AF', '4CF','5BM', '6BF') 
sex <- c( 'M', 'M', 'F', 'F', 'M', 'F')  
treatment <- c('C','A','A','C','B','B') 
timepoint <- c(0,1,2)
trial <- c('T1', 'T2', 'T3', 'T4') 


# create a metadata data frame work that integrates the above vectors
# Expand grid so that each timepoint is repeated across all animals first, then trials
rotarod <- expand.grid(animal = animal, timepoint = timepoint, trial = trial)

# map sex and treatment to each animal based on their index
rotarod$sex <- sex[match(rotarod$animal, animal)]
rotarod$treatment <- treatment[match(rotarod$animal, animal)]

# Sort by timepoint, then animal, then trial to ensure correct order in df
rotarod <- rotarod[order(rotarod$timepoint, rotarod$animal, rotarod$trial), ]


# create data vector (all 4 trial values/animal/timepoint)
time_to_fall <- c(208,167,218,194,105,131,204,94,90,185,236,215,289,292,300,300,135,105,235,249,136,152,273,165,   #timepoint 0 
                  196,168,248,258,189,264,289,246,199,264,180,NA,283,263,265,224,178,130,186,219,300,300,147,45, #timepoint 1   
                  18,17,32,24,183,NA,131,53,214,174,250,232,34,161,155,128,23,37,24,42,155,193,136,185) #timepoint 2

# add data vector as a column to metdata data frame
rotarod$time_to_fall <- time_to_fall[1:nrow(rotarod)] # applies time_to_fall vector rotarod df from index 1 through n rows 

# Reorder columns as specified: animal, sex, treatment, timepoint, trial, time_to_fall
rotarod <- rotarod[, c('animal', 'sex', 'treatment','timepoint', 'trial', 'time_to_fall')]

#View(rotarod)

# export to CSV file
write.csv(rotarod, here("Rotarod", "Dataframes", "rotarod.csv"))



############################################# Summary Data Frames ############################################################################

# Use custom function: groupSummary - get summary stats per animal, then by treatment, lastly, by sex and treatment

# 1) Animal Summary -> get mean of 4 trials/timepoint
animal_summary <- groupSummary(rotarod, time_to_fall, "mean_se", animal, sex, treatment, timepoint)
write.csv(animal_summary, here("Rotarod", "Dataframes", "animal_summary.csv"))


# 2) treatment Summary -> using mean of 4trials/animal calculated in above df, get mean/treatment group/timepoint
treatment_summary <- groupSummary(animal_summary, mean, "mean_se", treatment, timepoint)
treatment_summary <- treatment_summary %>% mutate(variable = ifelse(variable == 'mean', 'time_to_fall', variable)) #rename variable from animal mean to 'time_to_fall'
write.csv(treatment_summary, file = here("Rotarod", "Dataframes", "treatment_summary.csv"))

# summary by sex and treatments (interaction) (*not, only 1 of each sex in sample data, so SE = NA)
sex_summary <-groupSummary(animal_summary, mean, "mean_se", sex, treatment, timepoint)
sex_summary <- sex_summary %>% mutate(variable = ifelse(variable == 'mean', 'time_to_fall', variable))
write.csv(sex_summary, file = here("Rotarod", "Dataframes", "sex_summary.csv"))



############################################# ANALYSIS ############################################################################

# Check  Assumptions

# 1) Grubbs Test : Check for Outliers (rstatix)
identify_outliers(data = animal_summary, variable = "mean")  # no outliers

# 2) Reduce df : select columns of interest (1-4 = metdata, 7 = animal trial mean)
reduced_rotarod=animal_summary[,c(1,2,3,4,7)]
#View(reduced_rotarod)

# 3) Normality check
normResult <- reduced_rotarod %>%
  shapiro_test(mean)  # Use y_var_sym instead of y_var
is_normal <- all(normResult$p > 0.05)
print(normResult)                         # not normal -> report p.adj


# linear mixed effect model (lmer from lme4)
# treatment and timepoint are fixed effect, animal is random effect

# set reference terms explicitly 
reduced_rotarod$timepoint <- relevel(factor(reduced_rotarod$timepoint), ref = "0") # establishes baseline as with-in subject reference  
reduced_rotarod$treatment <- relevel(factor(reduced_rotarod$treatment), ref = "A")  # establishes control treatment as between subject reference 

model <- lmer(mean~treatment * timepoint + (1 | animal), data = reduced_rotarod) # (1 | animal) accounts for repeated measures of animals
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


############################################# DATA VIZUALIZATION ############################################################################


# 1) Line Plot : Time to Fall by treatment
color_palette <- c('A' = '#A3E4D7', 'B' = '#F8C471', 'C' = '#E67E22')

treatment_summary %>% filter(!is.na(mean))   # filter NAs before plotting

p <- line_plot(data=treatment_summary,
                  y_var = 'mean',
                  y_error = 'se',
                  title = "Rotarod",
                  y_label = 'Time to Fall (s)')


p  
ggsave(here("Rotarod", "Graphs", "rotarod.pdf"), 
       plot = p, width = 6, height = 4, dpi = 300)


# repeat plot by sex
sex_summary %>% filter(!is.na(mean), !is.na(treatment))   # filter NAs before plotting

p2 <- line_plot2(data=sex_summary,
                 y_var = 'mean',
                 y_error = 'se',
                 title = 'Rotarod',
                 y_label = 'Time to Fall (s)')
p2
ggsave(here("Rotarod", "Graphs", "rotarod_bySex.pdf"), 
       plot = p2, width = 6, height = 4, dpi = 300)



#### FACETED BOXPLOT: show spread w/in groups over time

rotarod$timepoint <- factor(rotarod$timepoint, levels=c('0','1','2'))    # ensure timepoint is a factor w/ ordered levels

facet_box <- ggplot(animal_summary, aes(x=as.factor(timepoint), y=mean, fill=treatment)) +
  geom_boxplot() +
  facet_wrap(~treatment) +    #facet plot by treatment
  labs(
    title = 'Rotarod Animal Variability',
    x= 'Time from Treatment',
    y = 'Time to Fall (sec)',
    fill = 'treatment Type'
  ) +
  scale_fill_manual(
    values=c('#A3E4D7','#F8C471','#E67E22'),
    labels = c('A','B','C')
  ) +
  scale_x_discrete(breaks = c('0','1','2')      #x-axis labeled as timepoints 
  ) +
  theme_bw() 
facet_box
ggsave(here("Rotarod", "Graphs", "FacetPlot.pdf"),
       plot = facet_box, width = 6, height = 4, dpi = 300)



#### TRAJECTORY: indiv. change from baseline / each animal 
change <- animal_summary %>% 
  filter(!is.na(mean)) %>% # filter NAs
  group_by(animal, treatment) %>%
  mutate(
    baseline =mean[timepoint == 0][1], #explicitly state baseline
    data_change = mean - baseline, #sub baseline at each timepoint
    timepoint = factor(timepoint)) #ensures timepoints print as exact in x-axis


# Create plot        
traj <-  ggplot(change, aes(x =timepoint, y = data_change, color = treatment, group = animal)) +
  geom_line(linewidth = 1) +    #Line connecting viral points
  labs(title= 'Individual Change Trajectory',
       x = "Time from Treatment", 
       y = "Change from Baseline (Time to Fall)", 
       color = "Treatment Type", 
       linetype = "Treatment Type") +
  theme_minimal() +
  scale_color_manual(values=c('#A3E4D7','#F8C471','#E67E22')) +
  scale_x_discrete(labels = levels(change$timepoint)) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
traj
ggsave(here("Rotarod", "Graphs", "TrajectoryPlot.pdf"),
       plot = traj, width = 6, height = 4, dpi = 300)



