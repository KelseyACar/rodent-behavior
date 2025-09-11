# This Script -------------------------------------------------------------
  # Visualizes and analyzes time series data for any single response variable for 3 factors 
  # In this case it is used to analyse Rotarod data for the NM R01 project. 
  # The three factors are: time, sex, treatment
    ## DEPENDENT VARIABLE: Time-to-Fall (seconds), can convert to RPMS at fall with following equation
  # RPMs at Fall: ((ending speed of full session - starting speed) / full session length) x (seconds at fall) + starting speed
  # RPMS for our study: ((40-4) / 300) X (secs at fall) + 4

# Load Libraries 
  # and custom function script
# Create Data Frames: 
  # Create from raw, manually collected data (rotarod time to fall: 4 trials / animal across 6 timepoints)
  # The NM R01 largely focuses on aged mice, but this script also handles the preliminary data used to write the R01 grant which was done in young animals.
  # At this time (Aug. 2025) the ages are not being compared
# Summary Data Frames: 
  # Drop animals that had missed injections after collecting animal summaries but before group summaries
  # Using custom functions, collect animal summaries (mean time to fall / timepoint / animal), treatment group summaries, and sex by treatment group summaries
# Analysis
  # Check assumptions, reduce files to cols of interest
  # 2 Way Repeated Measures mixed linear regression w/ holms post-hoc correction for multiple comparisons (time, treatment)
  # Adjust to 3 Way for treatment, sex, and time 
  # Analysis was ultimately done by entering the data from these created data frames and summaries into Prism, but outcome confirmed here
# Plots w/ ggplot2 
# Scrap : code not used ultimately
  # intermediate df reductions before data sets were complete
  # analysis when originally doing 2 and 3 way ANOVAs
  # code to extract significance and add to plot (not polished)

# Load Libraries ---------------------------------------------------------------
pacman::p_load('pacman', 'dplyr','ggplot2', 'here', 'FSA', 'rstatix', 'ggpubr', 'tidyr','lme4','lmerTest','emmeans') # auto installs required packages not yet installed
source(here::here('RotarodFunctions_KC.R')) # source in custom functions script 
here::i_am("Rotarod_timeSeries.R")


# Create Data Frames ----------------------------------------------------

#AGED

# define vectors (columns) of metadata
animal <- c('90030R', '90030LL', '102247L', '102247R', '102247LR', #cohort 1
            '90046L', '90046R', '107544L', '107544R', '107544LR', #cohort 2
            '94406R', '94406LR', '94406LL', '94407L', '94407R', '94407LR', '94407LL', #cohort 3
            '90032R', '108807L', '108807R', '107545LLR', '107545LR', '107549R', '107549LL', #cohort 4
            '108829LL', '108829LRR', '108829R', '112935L', '90042LR') #cohort 5

sex <- c( 'M', 'M', 'M', 'M', 'M', #cohort 1
          'F', 'F', 'M', 'M', 'M', #cohort 2
          'M', 'M', 'M', 'F', 'F', 'F', 'F', #cohort 3
          'M','M','M','F','F','F','F', #cohort 4
          'F','F','F','F','F') #cohort 5

virus <- c('AAV5-hTyr-Full','AAV5-hTyr-1:100','AAV5-hTyr-1:100','AAV5-hTyr-Full','AAV5-hTyr-Full', #cohort 1
           'AAV5-hTyr-1:100', 'AAV5-EYFP', 'AAV5-EYFP', 'AAV5-hTyr-1:100', 'AAV5-EYFP', #cohort2
           'AAV5-EYFP', 'AAV5-hTyr-Full', 'AAV5-hTyr-1:100', 'AAV5-hTyr-Full', 'AAV5-hTyr-1:100', 'AAV5-EYFP', 'AAV5-hTyr-Full', #cohort 3
           'AAV5-hTyr-Full', 'AAV5-hTyr-1:100', "AAV5-EYFP", 'AAV5-hTyr-1:100', 'AAV5-hTyr-Full', 'AAV5-EYFP', 'AAV5-hTyr-Full', #cohort 4
           'AAV5-hTyr-1:100','AAV5-EYFP','AAV5-EYFP','AAV5-hTyr-Full','AAV5-hTyr-1:100') #cohort 5 

trial <- c('T1', 'T2', 'T3', 'T4') 
timepoint <- c('0', '1', '2', '4','10','16')

# create a metadata data frame work that integrates the above vectors
# Expand grid so that each timepoint is repeated across all animals first, then trials
aged_rotarod <- expand.grid(animal = animal, timepoint = timepoint, trial = trial)

# map sex and virus to each animal based on their index
aged_rotarod$sex <- sex[match(aged_rotarod$animal, animal)]
aged_rotarod$virus <- virus[match(aged_rotarod$animal, animal)]

# Sort by timepoint, then animal, then trial to ensure correct order in df
aged_rotarod <- aged_rotarod[order(aged_rotarod$timepoint, aged_rotarod$animal, aged_rotarod$trial), ]


# create data vector (all 4 trial values/animal/timepoint)
  # spaces indicate breaks between cohorts, rows indicate time points (Ex. row 1 = time 0: cohort 1 values, then cohort 2, etc. )
time_to_fall <- c(208,167,218,194,105,131,204,94,90,185,236,215,289,292,300,300,135,105,235,249,      136,152,273,165,66,230,123,229,99,131,90,88,99,273,170,133,300,300,295,300,        239,196,108,147,82,124,107,82,109,242,230,265,68,137,192,197,159,213,118,165,109,191,167,96,246,92,286,260,        102,130,99,100,102,125,177,244,77,105,170,134,85,80,73,76,56,97,80,176,170,187,193,169,115,181,251,196,        137,265,266,173,87,63,96,103,200,163,143,117,136,130,198,211,75,261,253,276,   #0, 
                  196,168,248,258,189,264,289,246,199,264,180,300,283,263,265,224,178,130,186,219,    300,300,147,45,51,115,115,169,36,300,300,257,143,300,300,131,300,300,300,235,      244,267,285,300,108,172,205,152,222,300,297,300,215,172,275,283,255,283,300,186,223,202,210,195,241,279,207,85,    224,284,269,300,135,226,172,237,155,168,192,215,80,82,153,137,107,161,92,281,169,151,142,272,199,221,206,208,  300,276,286,300,93,157,113,233,159,203,243,229,130,117,86,179,204,205,234,243, #1,   
                  238,239,295,219,229,214,281,300,225,250,282,269,143,172,235,216,117,137,111,156,    151,215,249,183,205,269,261,177,300,38,214,299,192,294,217,155,300,300,300,300,    195,230,200,278,96,81,97,119,161,279,300,300,107,227,249,228,259,184,188,116,173,237,201,234,119,143,117,234,      196,234,257,222,116,224,265,162,140,220,216,236,115,140,125,300,55,36,71,81,249,281,267,273,126,78,163,101,    275,283,300,293,149,200,188,141,215,236,234,207,126,90,166,251,265,290,300,300,   #2,    
                  19,29,41,62,231,233,184,179,299,286,227,300,236,226,129,235,55,63,100,57,           208,165,240,208,223,275,216,210,300,269,291,262,240,265,196,165,300,300,300,300,   153,272,233,210,37,67,35,32,244,244,298,300,100,101,52,7,205,86,262,273,153,238,211,276,159,229,234,300,           70,79,123,132,159,36,185,33,177,261,292,300,169,98,175,101,76,88,73,69,222,216,248,220,24,86,18,70,            178,198,131,156,156,163,133,179,175,204,206,152,196,196,133,149,211,203,202,179,        #4, 
                  13,28,18,21,208,181,250,207,199,199,249,244,84,149,142,190,12,9,22,29,              80,107,107,142,142,121,140,236,245,160,130,286,131,207,216,235,300,300,300,300,    149,209,201,256,6,43,46,23,258,300,111,81,67,153,128,154,184,254,266,265,195,219,207,176,138,61,126,169,           84,57,78,67,198,166,44,175,198,166,44,175,NA,NA,NA,NA,26,5,36,40,148,156,209,211,62,53,75,69,                  126,184,106,184,157,155,204,173,160,114,76,88,95,95,108,101,212,235,270,263,       #10
                  18,17,32,24,183,125,131,53,214,174,250,232,34,161,155,128,23,37,24,42,              155,193,136,185,225,253,251,183,246,258,295,254,143,80,140,160,262,300,300,300,    104,240,254,259,40,38,17,20,90,48,216,250,49,59,66,42,74,148,186,182,177,138,175,161,145,218,300,300,              66,122,46,116,163,188,218,184,104,197,90,107,NA,NA,NA,NA,11,5,4,5,144,154,91,156,19,60,35,26,                  114,171,127,62,159,162,89,207,173,51,158,95,12,52,79,69,231,216,193,147)                 #16    

# add data vector as a column to metdata data frame
aged_rotarod$time_to_fall <- time_to_fall[1:nrow(aged_rotarod)] # applies time_to_fall vector to aged_rotarod df from index 1 through n rows 

# Reorder columns as specified: animal, sex, virus, timepoint, trial, time_to_fall
aged_rotarod <- aged_rotarod[, c('animal', 'sex', 'virus','timepoint','trial', 'time_to_fall')]

#View(aged_rotarod)

# export to CSV file
write.csv(aged_rotarod, here("Rotarod", "Dataframes", "aged_rotarod.csv"))



#YOUNG : repeat above process
  # (1769 virus first, then 1900 lot)


#animal <- c('94902L', '94902R', '94903L', '94903LR', '99208R', '99194L', '89039R', '89043L', '89043LR','74401L','74401LL','74401LR','74401R','81061L','81061LR','77402R')
#sex <- c('M', 'M', 'F', 'F','F','M', 'F', 'F', 'F','M','M','M','M', 'F','F','F') 
#virus <- c('AAV5-hTyr-Full','AAV5-hTyr-Full','AAV5-hTyr-Full','AAV5-hTyr-Full', 'AAV5-hTyr-1:100','AAV5-hTyr-1:100','AAV5-hTyr-1:100','AAV5-hTyr-1:100','AAV5-hTyr-1:100','AAV5-hTyr-Full','AAV5-hTyr-Full','AAV5-hTyr-Full','AAV5-hTyr-Full','AAV5-hTyr-Full','AAV5-hTyr-Full','AAV5-hTyr-Full')
#timepoint <- c('0', '1', '2', '4', '10')
#trial <- c('T1', 'T2', 'T3', 'T4')


# young_rotarod <- expand.grid(animal = animal, timepoint = timepoint, trial = trial)
# 
# # map sex and virus to each animal based on their index
# young_rotarod$sex <- sex[match(young_rotarod$animal, animal)]
# young_rotarod$virus <- virus[match(young_rotarod$animal, animal)]
# 
# # Sort by timepoint, then animal, then trial to ensure correct order
# young_rotarod <- young_rotarod[order(young_rotarod$timepoint, young_rotarod$animal, young_rotarod$trial), ]
# 
# # add all 4 trial values/animal/timepoint. ex. 90030R time=0 t1-t4, repeats for all animals of all cohorts at time=0. Then add for time=1 starting back at 90030R.  linebreak indicates timepoints
# time_to_fall <- c(233,233,260,271,137,219,114,266,164,206,203,60,300,300,151,229,245,184,300,300,94,139,179,178,159,168,200,196,199,228,249,203,201,139,232,255,NA,NA,NA,NA,220,250,139,264,298,139,300,267,140,300,176,215,300,300,118,300,153,208,252,300,NA,NA,NA,NA,   #0
#                 284,282,267,287,258,290,168,221,225,117,284,225,224,194,300,194,300,300,300,210,300,300,90,148,135,164,171,193,235,176,136,107,237,231,278,292,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,                      #1
#                 97,88,23,69,252,93,110,127,94,112,56,16,94,51,26,17,300,186,288,298,300,163,300,300,232,286,263,241,212,230,249,270,300,300,300,296,143,103,131,138,61,177,39,130,237,250,234,96,109,94,55,92,145,183,160,156,137,185,172,96,111,118,52,83,              #2    
#                 30,34,60,4,132,96,95,124,83,74,35,39,39,54,76,48,131,223,192,206,NA,NA,NA,NA,199,249,228,297,200,260,257,264,232,300,300,217,71,49,63,69,71,73,62,86,45,180,173,87,55,13,23,27,97,120,60,93,115,150,68,9,21,31,53,27,                                    #4
#                 NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,107,258,274,300,NA,NA,NA,NA,256,262,272,300,199,183,181,259,199,244,228,300,33,25,18,21,76,76,120,151,NA,NA,NA,NA,58,80,75,92)                                                                          #10
# 
# young_rotarod$time_to_fall <- time_to_fall[1:nrow(young_rotarod)]
# 
# # Reorder columns as specified: animal, sex, virus, timepoint, rota_Trial, rota_data
# young_rotarod <- young_rotarod[, c('animal', 'sex', 'virus','timepoint','trial', 'time_to_fall')]
# 
# View(young_rotarod)
# 
# # export to CSV file
# write.csv(young_rotarod, here("Rotarod", "Dataframes", "young_rotarod.csv"))







############################################# Summary Data Frames ############################################################################

# Use custom function: groupSummary - get summary stats per animal, then by virus, lastly, by sex and virus

#AGED

# 1) Animal Summary -> get mean of 4 trials/timepoint
aged_animal_summary <- groupSummary(aged_rotarod, time_to_fall, "mean_se", animal, sex, virus, timepoint)
write.csv(aged_animal_summary, here("Rotarod", "Dataframes", "aged_animal_summary.csv"))

# drop animals that are missed injections (verified by IHC / ephys) 
aged_animal_summary <- aged_animal_summary %>%
  filter(animal != '94407LL')

# 2) Virus Summary -> using mean of 4trials/animal calculated in above df, get mean/treatment group/timepoint
aged_virus_summary <- groupSummary(aged_animal_summary, mean, "mean_se", virus, timepoint)
aged_virus_summary <- aged_virus_summary %>% mutate(variable = ifelse(variable == 'mean', 'time_to_fall', variable)) #rename variable from animal mean to 'time_to_fall'
write.csv(aged_virus_summary, file = here::here("Rotarod", "Dataframes", "aged_virus_summary.csv"))

# summary by sex and viruses (interaction)
aged_sex_summary <-groupSummary(aged_animal_summary, mean, "mean_se", sex, virus, timepoint)
aged_sex_summary <- aged_sex_summary %>% mutate(variable = ifelse(variable == 'mean', 'time_to_fall', variable))
write.csv(aged_sex_summary, file = here::here("Rotarod", "Dataframes", "aged_sex_summary.csv"))


# #YOUNG
# #Step 1) animal summary -> get mean of 4 trials / timepoint
# young_animal_summary <- groupSummary(young_rotarod, time_to_fall, "mean_se", animal, sex, virus, timepoint)
# write.csv(young_animal_summary, file = here::here("Rotarod", "Dataframes", "young_animal_summary.csv"))
# 
# # summary by viruses
# young_virus_summary <- groupSummary(young_animal_summary, mean, "mean_se", virus, timepoint)
# young_virus_summary <- young_virus_summary %>% mutate(variable = ifelse(variable == 'mean', 'time_to_fall', variable))
# write.csv(young_virus_summary, file = here::here("Rotarod", "Dataframes", "young_virus_summary.csv"))
# 
# # summary by sex and viruses (interaction)
# young_sex_summary <-groupSummary(young_animal_summary, mean, "mean_se", sex, virus, timepoint)
# young_sex_summary <- young_sex_summary %>% mutate(variable = ifelse(variable == 'mean', 'time_to_fall', variable))
# write.csv(young_sex_summary, file = here::here("Rotarod", "Dataframes", "young_sex_summary.csv"))



############################################# ANALYSIS ############################################################################

# Check 3Way ANOVA Assumptions

# 1) Grubbs Test : Check for Outliers (rstatix)
identify_outliers(data = aged_animal_summary, variable = "mean")  # no outliers

# 2) Reduce df : select columns of interest (1-4 = metdata, 7 = animal trial mean)
reduced_rotarod=aged_animal_summary[,c(1,2,3,4,7)]
#View(reduced_rotarod)

# 3) Normality check
normResult <- reduced_rotarod %>%
  shapiro_test(mean)  # Use y_var_sym instead of y_var
is_normal <- all(normResult$p > 0.05)
print(normResult)                         # not normal -> report p.adj


# 2-way linear fixed effect model (lmer from lme4)

# mixed model: repeated measures
  # virus and timepoint are fixed effect, animal is random
model <- lmer(mean~virus * timepoint + (1 | animal), data = reduced_rotarod) # (1 | animal) accounts for repeated measures of animals

# ANOVA table for fixed effects
anova(model)

# Post-hoc: Dunnett contrast: compare every virus against the control w/in each timepoint
  # emmean = Estimated Marginal Means, (Least-Squares Means)
emm <- emmeans(model, ~ virus | timepoint)
dunnet_results <- contrast(emm, method = 'dunnett', ref = 'AAV5-EYFP', adjust = 'holm')
dunnet_results

# NOTE: these stats are different than reported in Prism where a 2 Way ANOVA was applied with Tukey multiple measures correction


############################################# DATA VIZUALIZATION ############################################################################

####### Make updates to below graphs based on age
## could also join the two dfs together but I'm not sure I want to do that yet.


### 1) Line PLOT : Time to Fall by virus
color_palette <- c('AAV5-EYFP' = '#A3E4D7', 'AAV5-hTyr-1:100' = '#F8C471', 'AAV5-hTyr-Full' = '#E67E22')

#AGED
aged_virus_summary %>% filter(!is.na(mean))   # filter NAs before plotting

aged <- line_plot(data=aged_virus_summary,
                  y_var = 'mean',
                  y_error = 'se',
                  title = "Rotarod",
                  y_label = 'Time to Fall (s)')


aged  
ggsave(here("Rotarod", "Graphs", "agedrotarod.pdf"), 
       plot = aged, width = 6, height = 4, dpi = 300)


# repeat plot by sex
aged_sex_summary %>% filter(!is.na(mean), !is.na(virus))   # filter NAs before plotting

p2 <- line_plot2(data=aged_sex_summary,
                 y_var = 'mean',
                 y_error = 'se',
                 title = 'Rotarod',
                 y_label = 'Time to Fall (s)')
p2
ggsave(here("Rotarod", "Graphs", "agedrotarod_bySex.pdf"), 
       plot = p2, width = 6, height = 4, dpi = 300)



#### FACETED BOXPLOT: show spread w/in groups over time

# AGED
aged_rotarod$timepoint <- factor(aged_rotarod$timepoint, levels=c('0','1','2','4','10'))    # ensure timepoint is a factor w/ ordered levels

facet_box <- ggplot(aged_animal_summary, aes(x=as.factor(timepoint), y=mean, fill=virus)) +
  geom_boxplot() +
  facet_wrap(~virus) +    #facet plot by virus
  labs(
    title = 'Rotarod Animal Variability',
    x= 'Time from Injection (weeks)',
    y = 'Time to Fall (sec)',
    fill = 'Virus Type'
  ) +
  scale_fill_manual(
    values=c('#A3E4D7','#F8C471','#E67E22'),
    labels = c('AAV5-EYFP','AAV5-hTyr-1:100','AAV5-hTy-Full')
  ) +
  scale_x_discrete(breaks = c('0','1','2','4','10','16')      #x-axis labeled as timepoints 
  ) +
  theme_bw() 
facet_box
ggsave(here("Rotarod", "Graphs", "agedFacetPlot.pdf"),
       plot = facet_box, width = 6, height = 4, dpi = 300)



#### TRAJECTORY: indiv. change from baseline / each animal 
#AGED
aged_change <- aged_animal_summary %>% 
  filter(!is.na(mean)) %>% # filter NAs
  group_by(animal, virus) %>%
  mutate(
    baseline =mean[timepoint == 0][1], #explicitly state baseline
    data_change = mean - baseline, #sub baseline at each timepoint
    timepoint = factor(timepoint)) #ensures timepoints print as exact in x-axis

View(aged_change)

# Create plot        
aged_traj <-  ggplot(aged_change, aes(x =timepoint, y = data_change, color = virus, group = animal)) +
  geom_line(linewidth = 1) +    #Line connecting viral points
  labs(title= 'Individual Change Trajectory : Aged',
       x = "Time from Injection (weeks)", 
       y = "Change from Baseline Time to Fall", 
       color = "Virus Type", 
       linetype = "Virus Type") +
  theme_minimal() +
  scale_color_manual(values=c('#A3E4D7','#F8C471','#E67E22')) +
  scale_x_discrete(labels = levels(aged_change$timepoint)) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
aged_traj
ggsave(here("Rotarod", "Graphs", "agedTrajectoryPlot.pdf"),
       plot = aged_traj, width = 6, height = 4, dpi = 300)



# scrap -------------------------------------------------------------------


# filter for any 'NA' in repeated measures (animals which have not yet had injection)
# aged_filtered <- aged_animal_summary %>%
# filter(aged_animal_summary$virus != 'NA')
# View(aged_filtered)

# Pre-group summaries: Reduce df : to only timepoints with data (only need to do this until Cohort 3 wk 16 added)
# check observations for anova (missing observations will break ANOVA so must run an alternative methods until dataset complete)
#table(aged_animal_summary$sex, aged_animal_summary$virus, aged_animal_summary$timepoint) 

#reduced_rotarod <- aged_animal_summary[!aged_animal_summary$timepoint == '16', ]
#reduced_rotarod$timepoint=factor(reduced_rotarod$timepoint) #set time as factor



# IF NOT NORMAL -> Kruskal-Wallis on grouped data
# main effect : time
# kw_result <- aged_virus_summary %>%
#   kruskal_test(mean~timepoint)
# print(kw_result)

# interaction : time x virus
# kw_interaction <- aged_virus_summary %>%
#   kruskal_test(mean ~ interaction(virus, timepoint))
# print(kw_interaction)
# 
# # Run 3-Way Repeated-Measures ANOVA
# anova_result <- anova_test(
#   data= reduced_rotarod,
#   dv = mean,
#   wid = animal,
#   between = c(sex, virus),
#   within = timepoint
# )

# print(anova_result)

# reduce to only significance
# anova_results <- as.data.frame(anova_result$ANOVA) 
# anova_sig <- anova_results %>% filter(p<0.05) # filter out ns 

# Add stars 
# anova_sig$stars <- case_when(
#   anova_sig$p < 0.001 ~ "***", 
#   anova_sig$p < 0.01 ~ "**",  
#   anova_sig$p < 0.05 ~ "*",   
#   TRUE ~ ""  # No stars if p >= 0.05
# )
# anova_sig

# Run Post-hoc: Dunn Test
# posthoc <- reduced_rotarod %>%
#   group_by(timepoint) %>%
#   dunn_test(mean ~ virus, p.adjust.method = "holm")
# posthoc

# Run Post-hoc Test: Dunn
#   # main effect : time
# dunn_main <- dunn_test(mean ~ virus, data = reduced_rotarod, p.adjust.method = 'holm')
# dunn_main

# Interaction: create col for Dunn test
# reduced_rotarod$interaction_var <- interaction(reduced_rotarod$virus, reduced_rotarod$timepoint)
# posthoc <- reduced_rotarod %>%
# dunn_test(mean ~ interaction_var, p.adjust.method = 'holm')
# test_type <- "Non-parametric (Dunn’s test with sidak correction)"
# View(posthoc)

# # D) extract stats for plotting -> report adjusted p-values always    # this needs to be adj. p values 4.10.25
# anova_mainSig <- anova_sig[anova_sig$Effect == 'timepoint', c('p','stars')] # reduce to sig cols by row (effect)
# interact_sigValues <- anova_sig[anova_sig$Effect == 'virus:timepoint', c('p', 'stars')]



# C) Post-hoc 
# 1) IF NORMAL -> pairwise_t_test w/ Tukey / bonferroni parametric correction
# posthoc <- reduced_rotarod %>%
#   pairwise_t_test(mean ~ interaction(sex, virus, timepoint), 
#                   p.adjust.method = 'bonferroni') 
# test_type <- "Parametric (Tukey or Bonferroni)"   # confirm that this runs tukey and not bonferroni

# #add stats to plot
# p <- aged + 
#   annotate("text", x = 0, y = 40, 
#            label = paste("Timepoint: p=", anova_mainSig$p), #plot adj. p-values
#            size = 4, hjust = 0, vjust = 0) +
#   annotate("text", x = 2.4, y = 40, 
#            label = paste(anova_mainSig$stars), #plot stars
#            size = 4, hjust = 0, vjust = 0) +
#   annotate("text", x = 0, y = 25, 
#            label = paste("Virus x Timepoint: p=", interact_sigValues$p), 
#            size = 4, hjust = 0, vjust = 0) +
#   annotate("text", x = 4.0, y = 25, 
#            label = paste(interact_sigValues$stars), 
#            size = 4, hjust = 0, vjust = 0)

