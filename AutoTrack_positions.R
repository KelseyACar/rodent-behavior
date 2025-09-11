# This Script -------------------------------------------------------------
  # Visualizes Autotrack open-field locomotor position data in 2D and 3D plots over time 
  # requires user to create a metadata file with each file name mapped to the animal, treatment, timepoint, etc. 
# Load Libraries 
# Create Dataframe
  # call in metadata frame (created outside of R)
  # loop over files of each timepoint, skipping unnecessary rows from program (like 21!)
  # read files & merge df
# Pre-process
  # reduce df to time desired for plots
  # convert time from 0.1 seconds to minutes
# 2D plots
  # x & y axes track plots over the first 15 min
  # 1 representative plot / virus type, facet grid by timepoint
# 3D plots
  # 1) Binary resting locations
    # re-read & pre-process full df 
    # set desired column to binary: here I did multiple in case we wanted to plot seperate things, but we stuck with z state of resting (1 = resting, 0 = not)
    # set color gradient for where in the session the resting (z_State == 1) occurred accross full length of the session
    # subset animal & timepoint of interest and plot 
    # plot saves directly to directory subfolder
  # 2) Animated version of the above plots
    # change dot color to stick w/ viral group colors & plot -> screen record to save animation (couldn't figure out export)
    # lose color gradient based on time in these plots
  # 3) Full session trajectory plot ("tornado" plot from Harris's VTA paper in Nature Comms)
    # re-subset the dataframe to select each timepoint, but retain full session length
    # plot that thang


# Load Libraries ---------------------------------------------------------------
pacman::p_load('dplyr', 'tidyverse', 'ggplot2', 'here', 'purrr', 'readr', 'viridis', 'scatterplot3d', 'plotly','png','orca','gridExtra') # auto installs required packages not yet installed
source(here::here('my_functionsKC.R'))
here::i_am("AutoTrack_positions.R")


# Create Dataframe --------------------------------------------------------

# call in metadata
metadata <- read.csv(here::here("metadata", "aged_animal_metadata.csv"), header = TRUE) # animal metadata
file_mapping <- read.csv(here::here("metadata", "openField_positions_fileName_meta.csv"), header = TRUE) # filenames metadata

# open files by timepoint - make sure these order animal in the same way the animal vector is set up
file_paths <- list(
  timepoint_0 = list.files(path = here::here("Locomotion", "Pre-inj. Run 2", "Position"), full.names= TRUE),
  timepoint_1 = list.files(path = here::here("Locomotion", "1wk post", "Position"), full.names= TRUE),
  timepoint_2 = list.files(path = here::here("Locomotion", "2wk post", "Position"), full.names= TRUE),
  timepoint_4 = list.files(path = here::here("Locomotion" ,"4wk post", "Position"), full.names= TRUE),
  timepoint_10 = list.files(path = here::here("Locomotion" ,"10wk post", "Position"), full.names= TRUE),
  timepoint_16 = list.files(path = here::here("Locomotion", "16wk post","Position"), full.names = TRUE)
)
# confirm correct order
# print(file_paths)

# Check for trailing whitespace in file_mapping
any(grepl("\\s+$", file_mapping$file_name))

# Convert file list to df
position_files <- tibble(
  file_name = basename(unlist(file_paths)),
  full_path = unlist(file_paths)
)
#View(position_files)

# Check for mismatch between file_mapping and position_files
anti_join(file_mapping, position_files, by = "file_name")
position_files$file_name <- trimws(position_files$file_name)

# merge metadata and file mapping
full_metadata <- file_mapping %>%
  left_join(metadata, by = "AnimalID")
View(full_metadata)

                                  
# merge files with metadata
position_df <- full_metadata %>%
  left_join(position_files, by = "file_name")
View(position_df)


# filer NA for timepoints / animals that don't have files yet (remove at end of project) 
position_df <- position_df %>%
  filter(!is.na(full_path), file.exists(full_path))

# read and process files into full df
pos_final <- position_df %>%
  mutate(data = pmap(., function(full_path, AnimalID, Sex, Virus, timepoint, ...){ #pmap passes columns in row-by-row manner to maintain metadata 
    
    # read files skipping unnecessary rows
    df <- read.csv(full_path, skip = 21, header = TRUE)
    df <- df %>% filter(Record.Number <= 36001)
    
    return(df)
  })) %>%
  unnest(data) %>%
  select(-file_name, -DOB, -InjectionDate, -full_path, -V.Axis.Rearing, -X)

# clean up
pos_final <- clean_names(pos_final)
pos_final <- pos_final %>% drop_na()

View(pos_final)

write.csv(pos_final, here("Locomotion", "Dataframes", "aged_positions_df.csv"))





################## Pre-process ##################

# reduce df to desired time frame (currently 15min = 900sec)
  # If time adjusted, adjust midpoint and color limits in graphs.
 pos_reduced <- pos_final %>%
  filter(time_s >=0 & time_s <= 900)

# convert time from 0.1 sec to minutes
pos_reduced <- pos_reduced %>%
  mutate(time_as_min =time_s / 60)
View(pos_reduced)



################## 2D plots ##################
# DO NOT use animals : 102247L, 90030LL, 90046R, 107549LL, 94407L, 107549LL (no rearing data), 107545LLR (euthanized at 10wks post-injection), 90030R (no rearing timepoint 16)

# Full
Full <- ggplot(pos_reduced %>% filter(animal_id == '90032R'), aes(x=x_position, y=y_position, color = time_as_min)) +
  geom_path(linewidth = 0.3) +
  facet_grid(~timepoint) +
  scale_colour_gradient2(low = "#fc8d59",mid = "#ffffbf",high = '#91bfdf',
                         midpoint = 7.5,
                         limits = c(0,15)) +
  theme_minimal(base_size = 10) +
  labs(title = 'AAV-hTyr-Full') +
  coord_fixed() + # maintains 1:1 aspect ratio
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
Full
ggsave(here::here("Locomotion", "Graphs", "position plots","2D track plots", "Full_2Dpos_15min.pdf"),
       plot = Full, width = 8, height = 8, dpi = 300)


# Dilute
Dilute <- ggplot(pos_reduced %>% filter(animal_id == '94406LL'), aes(x=x_position, y=y_position, color = time_as_min)) +
  geom_path(linewidth = 0.3) +
  facet_grid(~timepoint) +
  scale_colour_gradient2(low = "#fc8d59",mid = "#ffffbf",high = '#91bfdf',
                         midpoint = 7.5,
                         limits = c(0,15)) +
  theme_minimal(base_size = 10) +
  labs(title = 'AAV-hTyr-1:100') +
  coord_fixed() + # maintains 1:1 aspect ratio
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
Dilute
ggsave(here::here("Locomotion", "Graphs", "position plots", "2D track plots", "Dilute_2Dpos_15min.pdf"),
       plot = Dilute, width = 8, height = 8, dpi = 300)


# EYFP
EYFP <- ggplot(pos_reduced %>% filter(animal_id == '94406R'), aes(x=x_position, y=y_position, color = time_as_min)) +
  geom_path(linewidth = 0.3) +
  facet_grid(~timepoint) +
  scale_colour_gradient2(low = "#fc8d59",mid = "#ffffbf",high = '#91bfdf',
                         midpoint = 7.5,
                         limits = c(0,15)) +
  theme_minimal(base_size = 10) +
  labs(title = 'AAV-eYFP') +
  coord_fixed() + # maintains 1:1 aspect ratio
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
EYFP
ggsave(here::here("Locomotion", "Graphs", "position plots", "2D track plots", "EYFP_2Dpos_15min.pdf"),
       plot = EYFP, width = 8, height = 8, dpi = 300)


# plot together 
combined <- grid.arrange(EYFP, Dilute, Full)
combined
ggsave(here::here("Locomotion", "Graphs","position plots","2D track plots", "allGroups_2Dpos.pdf"),
       plot = combined, width = 6, height = 4, dpi = 300)











################## 3D plots ##################

# 1) 3D plot represent where in the chambers resting is happening w/ scatterplot3d

pos_final <- read.csv(here("Locomotion", "Dataframes", "aged_positions_df.csv"))

# Pre-processing

#convert Z-axes of interest to binary (z = rearing)
pos_3D <- pos_final %>%
  mutate(z_axis_rearing = ifelse(nchar(z_axis_rearing) >0 & !is.na(z_axis_rearing), 1, 0), # if rearing = 1, if not = 0 (in case this were interesting)
         z_state = ifelse(pos_final$state == 'Resting', 1, 0)) # the real unit of interest
#View(pos_3D)


# create color vector to plot points in a gradient over time of session  
# Initialize color vector
pos_3D$rear_color <- "gray80"

# Loop through each unique timepoint
for (animal in unique(pos_3D$animal_id)) {
  for (tp in unique(pos_3D$timepoint)) {
    
    # subset animal and timepoint so gradient is reflective across df
    sub_df <- pos_3D[pos_3D$animal_id == animal & pos_3D$timepoint == tp, ]
   
    # get row indices in full df
    full_idx <- which(pos_3D$animal_id == animal & pos_3D$timepoint == tp)
    
    # subset to just rearing
    rear_idx <- which(sub_df$z_state == 1)
    
    # Skip if no rearing bouts
    if (length(rear_idx) == 0) next
    
    # Rank time_s within the rearing points
    rear_times <- sub_df$time_s[rear_idx]
    rank_time <- rank(rear_times, ties.method = "first")
    
    # Make the palette match the number of rearing points
    rear_colors <- viridis(length(rank_time), option = 'inferno')[rank_time]
    
    # Assign colors to the full df
    pos_3D$rear_color[full_idx[rear_idx]] <- rear_colors
  }
}




# Filter for specific animal / virus type
# DO NOT use animals : 102247L, 90030LL, 90046R, 107549LL, 94407L (no rearing data), 107545LLR (euthanized at 10wks post-injection), 90030R (no rearing timepoint 16)
   # adjust timepoint == X to plot each timepoint, subsetting by time (0-900 = first 15min, 0-3600 = full session)

# full
full <- subset(pos_3D, animal_id == '90032R' & timepoint == 16)

graph_path <- here("Locomotion", "Graphs", "position plots", "3D resting", "3D_resting_Full_wk16.pdf")    
pdf(graph_path, width = 7, height = 7) # save to pdf must come before plot with base R
full_3D <- with(full[full$time_s >= 0 & full$time_s <= 3600, ],  # select session length
                scatterplot3d(x = x_position, y = y_position, z = state == 'Resting', main = "16", color = full$rear_color, pch = 20))

dev.off()


# dilute
dilute <-subset(pos_3D, animal_id =='94406LL' & timepoint == 16)

graph_path <- here("Locomotion", "Graphs", "position plots","3D resting", "3D_resting_Dilute_wk16.pdf") # set subfolder to save plot
pdf(graph_path, width = 7, height = 7) # save to pdf must come before plot with base R
dilute_3D <- with(dilute[dilute$time_s >= 0 & dilute$time_s <= 3600, ], 
                  scatterplot3d(x = x_position, y = y_position, z = state == 'Resting', main = "10", color = dilute$rear_color, pch = 20))
dev.off()


#EYFP
eyfp <- subset(pos_3D, animal_id =='94406R' & timepoint == 16)

graph_path <- here("Locomotion", "Graphs", "position plots","3D resting", "3D_resting_EYFP_wk16.pdf") # set subfolder to save plot
pdf(graph_path, width = 7, height = 7) # save to pdf must come before plot with base R
dilute_3D <- with(eyfp[eyfp$time_s >= 0 & eyfp$time_s <= 3600, ], 
                  scatterplot3d(x = x_position, y = y_position, z = state == 'Resting', main = "10", color = eyfp$rear_color, pch = 20))
dev.off()




# 2) animated version of above with plotly               # screen record to export
eyfp <- subset(pos_3D, animal_id =='107544LR') # re-subset df to contain all timepoints / animal 
dilute <-subset(pos_3D, animal_id =='94406LL')
full <- subset(pos_3D, animal_id == '90032R')

# 3D plot, frames = timepoints
plot_ly(data = eyfp, x = ~x_position, y = ~y_position, z = ~z_state, 
        type = 'scatter3d', 
        mode = 'markers', 
        marker = list(color = 'limegreen'), 
        frame = ~timepoint) %>%
  layout(title = list(text = 'Resting Locations: AAV-eYFP', y = 0.95, x = 0.5))
      



# 3) 3D trajectory plots (R version of Harris's 'tornado plot')

# Full
full2 <- subset(pos_final, animal_id == '90032R' & timepoint == 16)
dilute2 <- subset(pos_final, animal_id =='94406LL' & timepoint == 0)
eyfp2 <- subset(pos_final, animal_id == '94406R' & timepoint ==16)

p <- plot_ly(
  data = eyfp2, # change per group
  x = ~x_position,
  y = ~y_position,
  z = ~time_s, 
  type = 'scatter3d',
  mode = 'lines',
  line = list(
    color = ~time_s,          # optional gradient by time
    colorscale = 'Plasma')   # change to Plasma, Inferno, etc.
) %>%
  layout(
    title = "eYFP: wk16", # change per plot
    scene = list(
      aspectmode = 'manual',
      aspectratio = list(x=1, y=1, z=5), #manually set aspect ratio
      xaxis = list(showgrid = FALSE, zeroline = FALSE, title = 'X'), #remove grids, set titles
      yaxis = list(showgrid = FALSE, title = 'Y'),
      zaxis = list(showgrid = FALSE, title = 'Time'),
      
      camera = list(
        eye = list(x=1.5, y= 1.8, z = 1.4))
      )
)
p

