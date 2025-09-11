# This script for NM grip strength time series but applicable to any single response variable across 3 factors (here: virus, time, sex)


# libraries ---------------------------------------------------------------
pacman::p_load('pacman','dplyr','ggpubr','tidyverse','ggplot2','here','purrr','plotly','mosaic','ggridges','hrbrthemes','patchwork','rstatix', 'dunn.test') # auto installs required packages not yet installed

source(here::here('RotarodFunctions_KC.R'))
here::i_am("Grip_Strength.R")



# Create Raw Dataframe ----------------------------------------------------

# define vectors (columns)
animal <- c('90030R', '90030LL', '102247L', '102247R', '102247LR', #cohort 1
            '90046L', '90046R', '107544L', '107544R', '107544LR', #cohort 2
            '94406R', '94406LR', '94406LL', '94407L', '94407R', '94407LR', '94407LL', #cohort 3
            '90032R', '108807L', '108807R', '107545LLR', '107545LR', '107549R', '107549LL', #cohort 4
            '108829LL', '108829LRR', '108829R', '112935L', '90042LR') #cohort 5

sex <- c( 'M', 'M', 'M', 'M', 'M', 
          'F', 'F', 'M', 'M', 'M', 
          'M', 'M', 'M', 'F', 'F', 'F', 'F',
          'M','M','M','F','F','F','F',
          'F', 'F', 'F', 'F','F') 

virus <- c('AAV5-hTyr-Full','AAV5-hTyr-1:100','AAV5-hTyr-1:100','AAV5-hTyr-Full','AAV5-hTyr-Full', #cohort 1
           'AAV5-hTyr-1:100', 'AAV5-EYFP', 'AAV5-EYFP', 'AAV5-hTyr-1:100', 'AAV5-EYFP', #cohort2
           'AAV5-EYFP', 'AAV5-hTyr-Full', 'AAV5-hTyr-1:100', 'AAV5-hTyr-Full', 'AAV5-hTyr-1:100', 'AAV5-EYFP', 'AAV5-hTyr-Full', #cohort 3
           'AAV5-hTyr-Full', 'AAV5-hTyr-1:100', "AAV5-EYFP", 'AAV5-hTyr-1:100', 'AAV5-hTyr-Full', 'AAV5-EYFP', 'AAV5-hTyr-Full', #cohort 4
           'AAV5-hTyr-1:100','AAV5-EYFP','AAV5-EYFP','AAV5-hTyr-Full','AAV5-hTyr-1:100')

timepoint <- c('0', '1', '2', '4','10','16')

trial <- c('T1', 'T2', 'T3')

# Expand grid so that each timepoint is repeated across all animals first, then trials
aged_gripStrength <- expand.grid(animal = animal, timepoint = timepoint, trial = trial)

# map sex and virus to each animal based on their index
aged_gripStrength$sex <- sex[match(aged_gripStrength$animal, animal)]
aged_gripStrength$virus <- virus[match(aged_gripStrength$animal, animal)]

# Sort by timepoint, then animal, then trial to ensure correct order in df
aged_gripStrength <- aged_gripStrength[order(aged_gripStrength$timepoint, aged_gripStrength$animal, aged_gripStrength$trial), ]

# add all 3 trial values/animal/timepoint. ex. 90030R pre-inj t1-t3, repeats for all animals of all cohorts at time=0. Then add for time=1 starting back at 90030R. linebreaks indicate timepoint breaks.
  # NOTE: this test was started at cohort 1 16wks-post injection / cohort 3 4wks post injection. NAs for timepoints not caught in certain cohorts.

grip_strength <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,            NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,              NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,        68,115,94,120,135,88,109,82,167,131,149,160,161,117,131,147,147,139,84,102,102,  75,29,46,100,42,50,60,60,84,56,42,28,61,19,44,    #0 
                   NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,            NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,              NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,        76,116,54,75,92,41,52,56,108,10,74,66,102,128,118,93,78,119,57,72,87,            55,76,58,88,80,90,75,44,54,8,10,6,77,55,60,       #1
                   NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,            NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,              NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,        70,69,50,40,65,32,66,64,62,86,24,64,99,91,70,67,62,97,55,58,68,                  37,66,56,28,72,26,22,54,30,63,50,45,63,55,70,     #2 
                   NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,            NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,              68,57,98,43,90,133,19,68,59,51,53,45,55,19,30,84,110,73,175,112,162,   8,11,17,64,69,86,34,30,38,49,74,39,17,27,38,32,48,38,37,41,54,                   13,28,26,46,62,87,19,18,39,56,39,5,35,28,43,      #4 
                   NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,            123,34,19,94,109,44,127,213,179,131,138,142,148,118,122,   51,47,54,21,27,29,78,44,25,56,39,37,81,37,36,87,28,16,96,57,48,        35,56,26,21,22,42,86,41,86,NA,NA,NA,82,40,51,41,52,61,11,30,15,                  43,17,25,86,70,66,29,23,34,33,47,32,43,45,46,     #10 
                   17,101,120,10,75,158,117,72,92,135,136,101,103,212,114,  70,36,55,27,20,36,64,72,24,49,28,39,49,41,50,              74,85,56,16,23,9,30,68,37,61,54,51,32,54,29,54,48,33,61,41,57,         19,25,37,62,69,78,13,12,21,NA,NA,NA,107,73,51,57,47,58,54,56,33,                 44,24,30,50,27,41,36,28,40,29,22,24,25,59,32)    #16 

aged_gripStrength$grip_strength <- grip_strength[1:nrow(aged_gripStrength)]

# Reorder columns as specified: animal, sex, virus, timepoint, rota_Trial, rota_data
aged_gripStrength <- aged_gripStrength[, c('animal', 'sex', 'virus','timepoint','trial', 'grip_strength')]

#View(aged_gripStrength)

# drop missed injection animals (94407LL 8.26.25 - still confirming thru IHC, but quite sure)
aged_gripStrength <- aged_gripStrength %>%
  filter(animal != '94407LL')

# export to CSV file
write.csv(aged_gripStrength,file = here::here("Grip Strength", "Dataframes", "grip_raw.csv"))



############################################# SUMMARY Dataframes ############################################################################

# Use my groupSummary FN to get summary stats by animal
grip_summary <- groupSummary(aged_gripStrength, grip_strength, "mean_se", animal, virus, sex, timepoint)
write.csv(grip_summary, file = here::here("Grip Strength", "Dataframes", "grip_animalSummary.csv"))

# summary stats by virus
grip_virusSummary <- groupSummary(grip_summary, mean, "mean_se", virus, timepoint)
grip_virusSummary <- grip_virusSummary %>% mutate(variable = ifelse(variable == 'mean', 'grip_strength', variable)) #rename variable from mean to 'grip strength'
write.csv(grip_virusSummary,file = here::here("Grip Strength", "Dataframes", "grip_virusSummary.csv")) # n is off here ?

# summary stats by sex
grip_sexSummary <- groupSummary(grip_summary, mean, "mean_se", virus,sex,timepoint)
grip_sexSummary <- grip_sexSummary %>% mutate(variable = ifelse(variable == 'mean', 'grip_strength', variable))
write.csv(grip_sexSummary,file = here::here("Grip Strength", "Dataframes", "grip_sexSummary.csv"))


############################################# ANALYSIS ############################################################################

# 1) Grubbs Test : Check for Outliers (rstatix)
identify_outliers(data = grip_summary, variable = "mean")  

# 2) normality check
normResult <- grip_summary %>%
  shapiro_test(mean)  
is_normal <- all(normResult$p > 0.05)
print(normResult)  # not normal -> report p.adj


# 3) 3-Way Repeated-Measures ANOVA (rstatix)
anova_result <- anova_test(
  data= reduced_grip,
  dv = mean,
  wid = animal,
  between = c(sex, virus),
  within = timepoint
)
print(anova_result)

# reduce to only significance
anova_results <- as.data.frame(test$ANOVA) 
anova_sig <- anova_results %>% filter(p<0.05) # filter out ns 

# Add stars 
anova_sig$stars <- case_when(
  anova_sig$p < 0.001 ~ "***", 
  anova_sig$p < 0.01 ~ "**",  
  anova_sig$p < 0.05 ~ "*",   
  TRUE ~ ""  # No stars if p >= 0.05
)
anova_sig

# C) Post-hoc 
# 1) IF NORMAL -> pairwise_t_test w/ Tukey / bonferroni parametric correction
# posthoc <- reduced_grip %>%
#   pairwise_t_test(mean ~ interaction(sex, virus, timepoint), 
#                   p.adjust.method = 'bonferroni') 
# test_type <- "Parametric (Tukey or Bonferroni)"   # confirm that this runs tukey and not bonferroni

# 2) IF NOT NORMAL -> Kruskal-Wallis on grouped data
# main effect : time
kw_result <- grip_virusSummary %>%
  kruskal_test(mean~timepoint)
print(kw_result)

# interaction : time x virus
kw_interaction <- grip_virusSummary %>%
  kruskal_test(mean ~ interaction(virus, timepoint))
print(kw_interaction)

# 3) Dunn posthoc
# main effect : time
dunn_main <- dunn_test(mean ~ timepoint, data = reduced_grip, p.adjust.method = 'holm')
dunn_main

# interaction: create col for Dunn test
reduced_grip$interaction_var <- interaction(reduced_grip$virus, reduced_grip$timepoint)
posthoc <- reduced_grip %>%
  dunn_test(mean ~ interaction_var, p.adjust.method = 'holm')
test_type <- "Non-parametric (Dunn’s test with sidak correction)"
posthoc

# D) extract stats for plotting -> report adjusted p-values always    # this needs to be adj. p values 4.10.25
anova_mainSig <- anova_results[anova_results$Effect == 'timepoint', c('p','stars')] # reduce to sig cols by row (effect)
interact_sigValues <- anova_results[anova_results$Effect == 'virus:timepoint', c('p', 'stars')]




############################################# DATA VIZUALIZATION ############################################################################

# 1) BAR PLOT : Time to Fall

grip_virusSummary %>% filter(!is.na(mean))   # filter NAs before plotting

p <- line_plot(data=grip_virusSummary,
                  y_var = 'mean',
                  y_error = 'se',
                  title = "Grip Strength",
                  y_label = 'Force (g)')
# save plot
ggsave(here::here("Grip Strength", "Graphs", "gripBar.pdf"), 
       plot = p, width = 6, height = 4, dpi = 300)


#### FACETED BOXPLOT: show spread w/in groups over time

facet <- ggplot(grip_summary, aes(x=as.factor(timepoint), y=mean, fill=virus)) +
  geom_boxplot() +
  facet_wrap(~virus) +    #facet plot by virus
  labs(
    title = 'Grip Strength Animal Variability',
    x= 'Time from Injection (weeks)',
    y = 'Force (g)',
    fill = 'Virus Type'
  ) +
  scale_fill_manual(
    values=c('#A3E4D7','#F8C471','#E67E22'),
    labels = c('AAV5-EYFP','AAV5-hTyr-1:100','AAV5-hTy-Full')
  ) +
  scale_x_discrete(breaks = c('0','1','2','4','10','16')      #x-axis labeled as timepoints 
  ) +
  theme_bw() 
facet

#save plot
ggsave(here::here("Grip Strength", "Graphs", "FacetPlot.pdf"),
       plot = facet, width = 6, height = 4, dpi = 300)



#### TRAJECTORY: indiv. change from baseline / each animal 
grip_change <- grip_summary %>% 
  filter(!is.na(mean)) %>% # filter NAs
  group_by(animal, virus) %>%
  mutate(
    baseline =mean[timepoint == 0][1], #explicitly state baseline
    data_change = mean - baseline, #sub baseline at each timepoint
    timepoint = factor(timepoint)) #ensures timepoints print as exact in x-axis

# View(grip_change)

# Create plot
grip_reduced <- grip_change %>% filter(!is.na(data_change)) %>%  #remove NAs (all cohorts without full datasets) 
  mutate(timepoint = as.numeric(as.character(timepoint)))
View(grip_reduced)

table(grip_reduced$timepoint)
quartz()

traj <-  ggplot(grip_reduced, aes(x =timepoint, y = data_change, color = virus, group = animal)) +
  geom_line(aes(linetype = virus), linewidth = 1) +    #Line connecting viral points
  labs(title= 'Grip Change Trajectory',
       x = "Time from Injection (weeks)", 
       y = "Grip Force", 
       color = "Virus Type", 
       linetype = "Virus Type") +
  theme_minimal() +
  scale_color_manual(values=c('#A3E4D7','#F8C471','#E67E22')) +
  scale_x_continuous(breaks = c(0,1,2,4)) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(color = 'black')
  )
print(traj)

ggsave(here::here("Grip Strength", "Graphs", "Trajectory.pdf"),
       plot = traj, width = 6, height = 4, dpi = 300)
