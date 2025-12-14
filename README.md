---
editor_options: 
  markdown: 
    wrap: 72
---

## Data for *Habitat suitability for Sociable Lapwing (Vanellus gregarius) increases across its global range, but populations continue to decline*

## Description of the data and file structure

This dataset includes the following scripts

1)  TimeCalibration_SL_observation_landsat_img.R

This script uses Tasseled Cap indices (predictor variables) from all the
Landsat (4-9) scenes (already extracted) that overlap with each presence
point between 1990 and 2024, and corresponding 10 random points
(pseudoabsences) per presence point. The script then looks for all the
predictor variable values within 45 days before and after the date of
the observation and takes the median value for each predictor value.
Both presence and absence files are merged. Resulting file
'model_df_gee.csv' has a response column (presence and absence) with
temporally matched predictor columns. Raw presence and absence are not
uploaded due to their very large size. 'model_df_gee.csv' is provided
here and can be used directly in the Google Earth Engine (next step).

2)  GEE_randomforest_workflow.java

This script is written in Google Earth Engine and takes the file
generated in step 1) 'model_df_gee.csv' to train a random forest model
to predict monthly habitat suitability. Additionally it takes
'regional_changes.shp' that allows to generate regional changes in
habitat suitability for 15 examples regions through the annual migration
cycle. Due to the large size of 15 sites (and 34 years), statistics is
exported into 5 different files -

-   regional_changes_stats_I.csv
-   regional_changes_stats_II.csv
-   regional_changes_stats_III.csv
-   regional_changes_stats_IV.csv
-   regional_changes_stats_V.csv

For convenience, these files are provided in the GitHub repository.
Nevertheless, the script can be run in GEE to generate these files again
if needed.

3)  Plotting and summarizing results from the Random Forest Model a)
    Plotting_habitat_suitability.R This script takes the annual median
    habitat suitability layers exported from Google Earth Engine and
    summarizes them using BirdLife international home range polygons.
    This is Figure 2 in the manuscript.

<!-- -->

b)  Plotting_regional_trends.R This script takes the regional changes
    statistics files exported from Google Earth Engine and summarizes
    them to generate Figure 4 in the manuscript.

\*Users need a Google Earth Engine account to run this script and need
to upload the 'model_df_gee.csv' file and regional changes file in order
to run the script. \## Code/software

Users need R programming software and an active Google Earth Engine
account to view and process these data. Both scripts are shared. \##
Access information

Other publicly accessible locations of the data:

-   Not applicable

Data was derived from the following sources:

-   Not applicable
