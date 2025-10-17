# Rodent Behavior Analysis Scripts

A collection of R scripts for analyzing laboratory rodent behavior data from common behavioral assays (locomotion, rotarod, gait). These scripts were developed for time-series analysis of motor behavior in mice, with repeated measures across multiple timepoints. Test data is provided to test the application and adjust scripts to fit your lab's behavioral data outputs.

## Quick Start
1. Install required packages: `pacman::p_load(dplyr, tidyverse, ...)`
2. Organize data in folders by timepoint. Adjust read.csv to fit file type.
3. Update metadata vectors in script.
4. Run script section by section
   
## Overview
This repository contains analysis pipelines for three behavioral testing modalities:
- **Open-field locomotion** (Columbus Instruments Opto-varimex Auto-Track System
- **Rotarod performance** (manually collected data, applicable to any manual variable)
- **Gait analysis** (MouseSpecifics^TM^ DigiGait)

Each script handles data import, data frame build, preprocessing, statistical analysis (linear mixed-effects models), and visualization. While developed for a specific experimental design, these scripts can be adapted for similar rodent behavioral studies.

## Requirements
### R Version
R ≥ 4.0.0 recommended

### Required Packages
```r
pacman::p_load('dplyr', 'tidyverse', 'ggplot2', 'here', 'purrr', 'readr', 'viridis',
               'scatterplot3d', 'plotly', 'lme4', 'emmeans', 'writexl', 'rstatix', 'ggpubr',
               'janitor', 'mosaic', 'gridExtra', 'orca', 'png', 'dichromat')
```

The scripts use `pacman::p_load()` which automatically installs missing packages.

## Repository Structure

```
├── metadata/
│   ├── animal_metadata.csv
│   └── openField_positions_fileName_meta.csv
|   └── openField_residencyw_fileName_meta.csv
├── Locomotion/
│   ├── wk_0/
│   ├── wk_1/
│   ├── wk_2/
│   ├── Dataframes/
│   └── Graphs/
├── Digigait/
│   ├── wk_0/
│   ├── wk_1/
│   ├── wk_2/
│   ├── Dataframes/
│   └── Graphs/
├── Rotarod/
│   ├── Dataframes/
│   └── Graphs/
├── AutoTrack_openField.R
├── AutoTrack_positions.R
├── Digigait_timeSeries.R
├── Rotarod_timeSeries.R
├── Loco_functions.R
├── Gait_Functions.R
└── RotarodFunctions.R
```

## Scripts
### 1. AutoTrack_openField.R
Analyzes open-field locomotion data from Columbus Instruments Auto-Track system  (v.5.5.4)

**Input:** CSV files with bin-level data (exported from Auto-Track software)
- Skips first 21 rows of each file (Auto-Track header format)
- Reads from `Bins` subdirectories organized by timepoint
**Notes:**
1. Test data is exported in 5 min bins, this is mutable in software. 
2. This version of Auto-Track occassionaly exports additional rows of data outside of the session length. Must delete these prior to running script or data frame bind will not work.
   
**Outputs:**
- Auto-Track variables: Distance, Resting Time, Ambulatory/Stereotypic behavior, Rearing Events, Time Rearing, Average Speed)
- Time in Center vs. Edge (anxiety-like behavior)
- Summary statistics by treatment group
- Linear mixed-effects model results

**Key visualizations:** Line plots with error bars for all locomotor measures

### 2. AutoTrack_positions.R
Visualizes spatial movement patterns from Auto-Track position data.

**Input:** CSV files with X/Y/Z position data at 0.1s intervals

**Outputs:**
- 2D trajectory plots (colored by time)
- 3D resting location plots
- Interactive 3D animations (plotly)
- Full-session trajectory visualizations

**Note:** Requires metadata file mapping filenames to animal IDs

### 3. Digigait_timeSeries.R
Analyzes gait parameters from MouseSpecifics^TM^ DigiGait system.

**Input:** 
- DigiGait "Indices" export files (combines two-row headers automatically)
- Manually recorded speed data (min/max belt speeds achieved)

**Outputs:**
- Gait indices by limb and belt speed
- Stride parameters, paw angles, stance metrics
- Speed achievement analysis
- Summary statistics and LMER results

**Key features:** 
- Handles multiple belt speeds separately
- Combines fore/hind paw analyses
- Trajectory plots showing individual performance changes

### 4. Rotarod_timeSeries.R
Analyzes rotarod time-to-fall data (adaptable to any single-response variable).

**Input:** Manually entered data (adaptable to any number of trials per animal per timepoint)

**Outputs:**
- Mean time-to-fall per animal/timepoint
- Treatment group comparisons
- Trajectory plots showing individual changes from baseline
- Faceted boxplots showing within-group variability

**Note:** Includes formula to convert time-to-fall to RPM at fall if needed

## Usage

### Basic Workflow

1. **Organize your data** according to the expected folder structure
2. **Update metadata vectors** in each script:
   ```r
   animal <- c('1CM', '2AM', '3AF', '4CF','5BM', '6BF') 
   sex <- c('M', 'M', 'F', 'F', 'M', 'F')  
   treatment <- c('C','A','A','C','B','B') 
   timepoints <- c(0,1,2)
   ```
3. **Source the appropriate functions file** (already at top of each script):
   ```r
   source(here('my_functions.R'))  # for locomotion scripts
   source(here('Gait_Functions.R'))  # for DigiGait
   source(here('RotarodFunctions.R'))  # for rotarod
   ```
4. **Run the script section by section** or in full

### Customization Points

**File paths:** Update these lists to point to your data directories and number of timepoints:
```r
file_paths <- list(
  timepoint_0 = list.files(path = here("YourFolder", "Timepoint0"), full.names= TRUE),
  timepoint_1 = list.files(path = here("YourFolder", "Timepoint1"), full.names= TRUE),
  timepoint_2 = list.files(path = here("YourFolder", "Timepoint2"), full.names= TRUE)
)
```

**Variables of interest:** Modify the `variables` vector for your specific measures:
```r
sum_variables <- c("ambulatory_s", "stereotypic_s", "resting_s", ...)
```

**Statistical model:** Adjust the LMER formula for your experimental design:
```r
model <- lmer(variable ~ treatment * timepoint + (1 | animal), data = df)
```

**Plot aesthetics:** Colors, labels, and themes can be modified in the plotting sections

## Data Format Requirements

### Auto-Track Locomotion
- CSV files with 21-row header (standard Auto-Track export)
- Required columns: Bin #, Zone, Distance (cm), Ambulatory (s), Resting (s), Rearing Events, etc.
- Files organized in separate folders by timepoint

### Auto-Track Positions
- CSV files with Record Number, Time (s), X Position, Y Position, State columns
- Requires separate metadata CSV mapping filenames to animal IDs

### DigiGait
- Standard DigiGait "Indices" export (2-row header format)
- Files organized by timepoint in separate folders
- Speed data entered manually in script

### Rotarod
- Data entered directly as vectors in the script
- Format: All trials for one animal at one timepoint, then next animal, etc.

## Statistical Analysis
All scripts use **linear mixed-effects models (LMER)** from the `lme4` package:
- **Fixed effects:** Treatment, timepoint, and their interaction
- **Random effect:** Animal (accounts for repeated measures)
- **Post-hoc comparisons:** Holm-adjusted pairwise contrasts using `emmeans`

The scripts include:
- Outlier detection (Grubbs test)
- Normality checks (Shapiro-Wilk)
- Model residual plots (commented out by default)
- Export of contrast results to Excel

## Output Files
Each script generates:
- **Dataframes/** - Processed data and summary statistics (CSV)
- **Graphs/** - Publication-ready plots (PDF, 300 dpi)
- **contrasts_results.xlsx** - Statistical comparison results

## Custom Functions
The scripts rely on custom functions defined in separate files:

**loco_functions.R:**
- `dropErrors()` - Drops animal|timepoint for any reason (Ex. equipment errors, study dropouts)
-  `animalSums()` - Summarizes variables per animal
- `timeInZones()` - Calculates time in center/edge
- `plot_openField()` - Creates treatment group line plots

**Gait_Functions.R:**
- `gaitSummaries()` - Summarizes gait indices

**RotarodFunctions.R:**
- `groupSummary()` - Flexible grouped summaries
- `line_plot()`, `line_plot2()` - Treatment/sex line plots

## Notes and Tips
- **File order matters:** Files are read in order from `list.files()`. Verify the order matches your metadata vectors.
- **Recording errors:** Use the `dropErrors()` function to remove specific variables for animals with equipment malfunctions
- **Missing data:** The scripts handle NA values appropriately in summaries and models
- **Bin ranges:** Adjust `bin_range = 1:6` to analyze different session lengths (1:6 = first 30 min for 5-min bins)
- **Visualization time:** For position plots, modify time filters to show desired session segments

## Adapting for Your Study
To use these scripts with your own data:

1. Maintain the same file organization structure or update file paths
2. Update metadata vectors (animal IDs, treatments, sex, timepoints)
3. Verify column names match expected format (or update `expected_cols`)
4. Adjust bin numbers, belt speeds, or trial numbers as needed
5. Modify plotting aesthetics (colors, labels) to match your experimental groups
6. Update statistical model if your design differs (e.g., add covariates, change random effects)

## Citation
If you use these scripts in your research, please cite this repository and acknowledge the software packages used (run `cite_packages()` to generate citations).

## Contact
Kelsey-Carter@omrf.org, [LinkedIn](https://www.linkedin.com/in/kelsey-carter-86427a199)
Feel free to reach out with suggestions, edits, and corrections. This is my first full project and as such, is a learning process. 

