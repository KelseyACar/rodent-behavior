# This Script -------------------------------------------------------------
  # Pulls in residency data from Auto-Track software used with Columbus Instruments open-field locomotion chambers & visualizes animals location as heat maps
# Load Libraries
# Build Dataframe
  # call in pre-made metadata excel sheet of project info: animal ID, time, treatment, etc.
  # call in a pre-made filename sheet
  # loop over raw data files and merge with metadata
# Pre-Processing
  # log transform time spent at each x/y
  # global scale time spent with min/max so that heat map colors represent the same value in all plots
# Heat maps
  # set palette
  # for each treatment type, filter to 1 representative animal and plot a heat map of their residency, faceted by time
  # combine all treatment plots into 1 grand plot if desired

# Load Libraries ----------------------------------------------------------
pacman::p_load('dplyr', 'tidyverse', 'ggplot2', 'here', 'purrr', 'readr', 'viridis', 'dichromat', 'gridExtra') # auto installs required packages not yet installed
here::i_am("AutoTrack_residency.R")
source(here('Loco_functions.R'))

# to collect pack citations, uncomment and run line below:
# cite_packages()

# to collect pack versions, uncomment and run line below:
# report_packages()


# Build Dataframe  --------------------------------------------------------
# call in metadata
metadata <- read.csv(here("metadata", "animal_metadata.csv"), header = TRUE) # animal metadata
file_mapping <- read.csv(here("metadata", "openField_residency_fileName_meta.csv"), header = TRUE) # filenames metadata

# open files by timepoint - make sure these order animal in the same way the animal vector is set up
file_paths <- list(
  timepoint_0 = list.files(path = here("Locomotion", "wk_0", "Residency"), full.names= TRUE),
  timepoint_1 = list.files(path = here("Locomotion", "wk_1", "Residency"), full.names= TRUE),
  timepoint_2 = list.files(path = here("Locomotion", "wk_2", "Residency"), full.names= TRUE))

# Check for trailing whitespace in file_mapping names
any(grepl("\\s+$", file_mapping$file_name))

# Convert file list to df
Residency_files <- tibble(
  file_name = basename(unlist(file_paths)),
  full_path = unlist(file_paths)
)
# View(Residency_files)

# Check for mismatch between file_mapping and Residency_files
#anti_join(file_mapping, Residency_files, by = "file_name")
#Residency_files$file_name <- trimws(Residency_files$file_name)

# merge metadata and file mapping
full_metadata <- file_mapping %>%
  left_join(metadata, by = "AnimalID")
# View(full_metadata)


# merge files with metadata       
residency_df <- full_metadata %>%
  left_join(Residency_files, by = "file_name")
# View(residency_df)


# filer NA for timepoints / animals that don't have files yet (remove at end of project) 4.5.25
residency_df <- residency_df %>%
  filter(!is.na(full_path), file.exists(full_path))


# read and process files into full df
res_process <- function(full_path, AnimalID, Sex, Treatment, Timepoint, ... ) {
  # read files skipping unnecessary rows
  raw <- read_csv(full_path, skip = 24, col_names = FALSE)
  
  # set axes
  y_labels <- as.numeric(unlist(raw[1,-1])) # x-axis 
  x_labels <- as.numeric(raw[-1,1][[1]]) # y-axis
  
  # extract matrix values
  values <- raw[-1,-1] %>% 
    mutate_all(as.numeric) %>%
    as.matrix()
  
  #set df to long format
  df_long <-expand.grid(
    x = x_labels,
    y = y_labels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  df_long$time_spent <- as.vector(values)
  
  # Add metadata 
  df_long <- df_long %>%
    mutate(
      animalID = AnimalID,
      sex = Sex,
      treatment = Treatment,
      timepoint = Timepoint
    )
  return(df_long)
}

# apply function to each file
res_final <- pmap_dfr(residency_df %>% select(full_path, AnimalID, Sex, Treatment, Timepoint),res_process)

# creates an extra y axis of NA values, filter out 
res_final <- res_final %>%
  filter(res_final$y != 'NA')
View(res_final)





# Pre-Processing ----------------------------------------------------------
# log transform, then min/max scale residency plots 
res_scaled <- res_final %>%
  mutate(time_log = log1p(time_spent)) %>% #log 1p (natural log) accounts for data containing zeros
  mutate(
    time_scaled = (time_log - min(time_log)) /    #global scale so a color means the same thing in every graph
      (max(time_log) - min(time_log))
  ) %>% 
  ungroup()

# View(res_scaled)
write.csv(res_scaled, file = here("Locomotion", "Dataframes", "residency_df.csv"))


# Heat Maps ---------------------------------------------------------------
# generate HCL palette (hue, chroma, luminance)
palette <- colorspace::sequential_hcl(5, palette = 'batlow')


# A plot
A <-ggplot(res_scaled %>% filter(animalID == "2AM"), aes(x= x, y= y, fill= time_scaled))+
  geom_tile(color= "white") + 
  facet_grid(~timepoint) +
  scale_fill_gradientn(colours = palette) +
  coord_equal() +
  theme_minimal(base_size = 8) +
  labs(
    title= 'A : 2AM') +  
  theme(
    legend.position = "right", 
    plot.title=element_text(size = 10, hjust =0),
    axis.text.y=element_text(size=6),
    axis.ticks=element_blank(),
    axis.text=element_text(size=6),
    legend.title=element_text(size=6),
    legend.text=element_text(size=4),
    strip.background = element_rect(colour = 'white')
  )
A
ggsave(here("Locomotion", "Graphs", "residency plots","A_residency.pdf"),
       plot = A, width = 6, height = 4, dpi = 300)


# 1:100 plot : don't use 2AM
B <-ggplot(res_scaled %>% filter(animalID == "5BM"), aes(x= x, y= y, fill= time_scaled))+
  geom_tile(color= "white") + 
  facet_grid(~timepoint) +
  scale_fill_gradientn(colours = palette) +
  coord_equal() +
  theme_minimal(base_size = 8) +
  labs(
    title= 'B : 5BM') +  
  theme(
    legend.position = "right", 
    plot.title=element_text(size = 10, hjust =0),
    axis.text.y=element_text(size=6),
    axis.ticks=element_blank(),
    axis.text=element_text(size=6),
    legend.title=element_text(size=6),
    legend.text=element_text(size=4),
    strip.background = element_rect(colour = 'white')
  )
B

ggsave(here("locomotion", "Graphs","residency plots", "B_residency.pdf"),
       plot = B, width = 6, height = 4, dpi = 300)


# C plot
C <-ggplot(res_scaled %>% filter(animalID == "4CF"), aes(x= x, y= y, fill= time_scaled))+
  geom_tile(color= "white") + 
  facet_grid(~timepoint) +
  scale_fill_gradientn(colours = palette, limits = c(0,1)) +
  coord_equal() +
  theme_minimal(base_size = 8) +
  labs(
    title= 'C : 4CF') +  
  theme(
    legend.position = "right", 
    plot.title=element_text(size = 10, hjust =0),
    axis.text.y=element_text(size=6),
    axis.ticks=element_blank(),
    axis.text=element_text(size=6),
    legend.title=element_text(size=6),
    legend.text=element_text(size=4),
    strip.background = element_rect(colour = 'white')
  )
C

ggsave(here("Locomotion", "Graphs","residency plots", "C_residency.pdf"),
       plot = C, width = 6, height = 4, dpi = 300)


# combine to one plot
combined <- grid.arrange(A, B, C)

ggsave(here("Locomotion", "Graphs","residency plots", "residencyHeatMap_repAnimal.pdf"),
       plot = combined, width = 6, height = 4, dpi = 300)





