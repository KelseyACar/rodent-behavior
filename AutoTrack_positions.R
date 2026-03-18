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
source(here('Loco_functions.R'))
here::i_am("AutoTrack_positions.R")

# to collect pack citations, uncomment and run line below:
# cite_packages()

# to collect pack versions, uncomment and run line below:
# report_packages()


# Create Dataframe --------------------------------------------------------

# call in metadata
metadata <- read.csv(here("metadata", "animal_metadata.csv"), header = TRUE) # animal metadata
file_mapping <- read.csv(here("metadata", "openField_positions_fileName_meta.csv"), header = TRUE) # filenames metadata

# open files by timepoint - make sure these order animal in the same way the animal vector is set up
file_paths <- list(
  timepoint_0 = list.files(path = here("Locomotion", "wk_0", "Position"), full.names= TRUE),
  timepoint_1 = list.files(path = here("Locomotion", "wk_1", "Position"), full.names= TRUE),
  timepoint_2 = list.files(path = here("Locomotion", "wk_2", "Position"), full.names= TRUE))

# confirm correct order
# print(file_paths)

# Check for trailing whitespace in file_mapping
any(grepl("\\s+$", file_mapping$file_name))

# Convert file list to df
position_files <- tibble(
  file_name = basename(unlist(file_paths)),
  full_path = unlist(file_paths))


# Check for mismatch between file_mapping and position_files
anti_join(file_mapping, position_files, by = "file_name")
position_files$file_name <- trimws(position_files$file_name)
#View(position_files)

# merge metadata and file mapping
full_metadata <- file_mapping %>%
  left_join(metadata, by = "AnimalID")
#View(full_metadata)

# merge files with metadata
position_df <- full_metadata %>%
  left_join(position_files, by = "file_name")
#View(position_df)



# read and process files into full df
pos_final <- position_df %>%
  mutate(data = pmap(., function(file_name, full_path, AnimalID, Sex, Treatment, timepoint, ...){ #pmap passes columns in row-by-row manner to maintain metadata 
    
    # read files skipping unnecessary rows
    df <- read.csv(full_path, skip = 21, header = TRUE)
    df <- df %>% filter(Record.Number <= 36001)
    
    return(df)
  })) %>%
  unnest(data) %>%
  select(-file_name, -full_path, -V.Axis.Rearing)


# clean up
pos_final <- clean_names(pos_final)
pos_final <- pos_final %>% drop_na()

View(pos_final)

write.csv(pos_final, here("Locomotion", "Dataframes", "positions_df.csv"))





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


# Adjust animal ID and treatment type as needed
A <- ggplot(pos_reduced %>% filter(animal_id == '2AM'), aes(x=x_position, y=y_position, color = time_as_min)) +
  geom_path(linewidth = 0.3) +
  facet_grid(~timepoint) +
  scale_colour_gradient2(low = "#fc8d59",mid = "#ffffbf",high = '#91bfdf',
                         midpoint = 7.5,
                         limits = c(0,15)) +
  theme_minimal(base_size = 10) +
  labs(title = 'Treatment A') +
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
A
ggsave(here("Locomotion", "Graphs", "position plots","2D track plots", "A_2Dpos_15min.pdf"),
       plot = A, width = 8, height = 8, dpi = 300)



# B
B <- ggplot(pos_reduced %>% filter(animal_id == '5BM'), aes(x=x_position, y=y_position, color = time_as_min)) +
  geom_path(linewidth = 0.3) +
  facet_grid(~timepoint) +
  scale_colour_gradient2(low = "#fc8d59",mid = "#ffffbf",high = '#91bfdf',
                         midpoint = 7.5,
                         limits = c(0,15)) +
  theme_minimal(base_size = 10) +
  labs(title = 'Treatment B') +
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
B
ggsave(here("Locomotion", "Graphs", "position plots", "2D track plots", "B_2Dpos_15min.pdf"),
       plot = B, width = 8, height = 8, dpi = 300)


# C
C <- ggplot(pos_reduced %>% filter(animal_id == '4CF'), aes(x=x_position, y=y_position, color = time_as_min)) +
  geom_path(linewidth = 0.3) +
  facet_grid(~timepoint) +
  scale_colour_gradient2(low = "#fc8d59",mid = "#ffffbf",high = '#91bfdf',
                         midpoint = 7.5,
                         limits = c(0,15)) +
  theme_minimal(base_size = 10) +
  labs(title = 'Treatment C') +
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
C
ggsave(here("Locomotion", "Graphs", "position plots", "2D track plots", "C_2Dpos_15min.pdf"),
       plot = C, width = 8, height = 8, dpi = 300)


# plot together 
combined <- grid.arrange(A, B, C)
combined
ggsave(here("Locomotion", "Graphs","position plots","2D track plots", "allGroups_2Dpos.pdf"),
       plot = combined, width = 6, height = 4, dpi = 300)







################## 3D plots ##################

# 1) 3D plot represent where in the chambers resting is happening 

pos_final <- read.csv(here("Locomotion", "Dataframes", "positions_df.csv"))


#convert axes of interest to binary: can replace with any action of interest (Resting, rearing, etc.)
pos_3D <- pos_final %>%
  mutate(z_state = ifelse(pos_final$state == 'Resting', 1, 0))
#View(pos_3D)


# create color vector to plot points in a gradient over time of session  
# Initialize color vector
pos_3D$z_color <- "gray80"

# Loop through each unique timepoint
for (animal in unique(pos_3D$animal_id)) {
  for (tp in unique(pos_3D$timepoint)) {
    
    # subset animal and timepoint so gradient is reflective across df
    sub_df <- pos_3D[pos_3D$animal_id == animal & pos_3D$timepoint == tp, ]
    
    # get row indices in full df
    full_idx <- which(pos_3D$animal_id == animal & pos_3D$timepoint == tp)
    
    # subset to just state of interest ('z_state')
    z_idx <- which(sub_df$z_state == 1)
    
    # Skip if no rearing bouts
    if (length(z_idx) == 0) next
    
    # Rank time_s within the rearing points
    z_times <- sub_df$time_s[z_idx]
    rank_time <- rank(z_times, ties.method = "first")
    
    # Make the palette match the number of rearing points
    z_colors <- viridis(length(rank_time), option = 'inferno')[rank_time]
    
    # Assign colors to the full df
    pos_3D$z_color[full_idx[z_idx]] <- z_colors
  }
}




# Filter for specific animal / treatment type
# adjust timepoint to plot each timepoint and subset by time to capture session length of interest (0-900 = first 15min, 0-3600 = full session)

# Treatment A
Aa <- subset(pos_3D, animal_id == '2AM' & timepoint == 0)

graph_path <- here("Locomotion", "Graphs", "position plots", "3D resting", "3D_resting_A_wk0.pdf")    
pdf(graph_path, width = 7, height = 7) # save to pdf must come before plot with base R
full_3D <- with(Aa[Aa$time_s >= 0 & Aa$time_s <= 3600, ],  # select session length
                scatterplot3d(x = x_position, y = y_position, z = state == 'Resting', main = "0", color = Aa$z_color, pch = 20))
dev.off()


# Treatment B
Ba <-subset(pos_3D, animal_id =='5BM' & timepoint == 0)

graph_path <- here("Locomotion", "Graphs", "position plots","3D resting", "3D_resting_Ba_wk0.pdf") # set subfolder to save plot
pdf(graph_path, width = 7, height = 7) # save to pdf must come before plot with base R
Ba_3D <- with(Ba[Ba$time_s >= 0 & Ba$time_s <= 3600, ], 
                  scatterplot3d(x = x_position, y = y_position, z = state == 'Resting', main = "10", color = Ba$z_color, pch = 20))
dev.off()


# Treatment C
Cc <- subset(pos_3D, animal_id =='4CF' & timepoint == 0)

graph_path <- here("Locomotion", "Graphs", "position plots","3D resting", "3D_resting_Cc_wk0.pdf") # set subfolder to save plot
pdf(graph_path, width = 7, height = 7) # save to pdf must come before plot with base R
a_3D <- with(Cc[Cc$time_s >= 0 & Cc$time_s <= 3600, ], 
                  scatterplot3d(x = x_position, y = y_position, z = state == 'Resting', main = "10", color = Cc$z_color, pch = 20))
dev.off()




# 2) animated version of above with plotly               # screen record to export
A <- subset(pos_3D, animal_id =='2AM') # re-subset df to contain all timepoints / animal 
B <-subset(pos_3D, animal_id =='5BM')
C <- subset(pos_3D, animal_id == '4CF')


# 3D plot, frames = timepoints
plot_ly(data = C, x = ~x_position, y = ~y_position, z = ~z_state, 
        type = 'scatter3d', 
        mode = 'markers', 
        marker = list(color = 'limegreen'), 
        frame = ~timepoint) %>%
  layout(title = list(text = 'Resting Locations: Treatment C', y = 0.95, x = 0.5))




# 3) 3D trajectory plots


A2 <- subset(pos_final, animal_id == '2AM' & timepoint == 0)
B2 <- subset(pos_final, animal_id =='5BM' & timepoint == 0)
C2 <- subset(pos_final, animal_id == '4CF' & timepoint ==0)

p <- plot_ly(
  data = C2, # change per group
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
    title = "C: wk 0", # change per plot
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


