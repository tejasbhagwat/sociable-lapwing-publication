# Data for Agricultural land predicts changes in the habitat suitability for Critically Endangered Sociable Lapwing (Vanellus gregarius) through the annual cycle

Dataset DOI: [10.5061/dryad.6q573n6b6](10.5061/dryad.6q573n6b6)

## Description of the data and file structure

This dataset includes the following scripts

1\) TimeCalibration_SL_observation_landsat_img.R

This script imports 'TSS_soc_lap_allyears_1990-2023_LSAT.csv' which has predictor variables extracted from all the Landsat scenes (4-9) that 
overlap with each presence point between 1990 and 2023, and corresponding 10 random per presence point. 
The script then looks for all the predictor variable values within 45 days before and after the date of the observation and takes the median value for each predictor value. 
Both presence and absence files are merged. Resulting file 'model_df_gee.csv' has a response column (presence and absence) with temporally matched predictor columns. 
Raw presence and absence are not uploaded due to their very large size. 'model_df_gee.csv' is provided here and can be used directly in the Google Earth Engine (next step) 

2\) GEE_randomforest_workflow.java

This script is written in Google Earth Engine and takes the file generated in step 1) 'model_df_gee.csv' to train a random forest model to predict monthly habitat suitability. 
Additionally it takes 'regional_changes.shp' that allows to generate regional changes in habitat suitability for 15 examples regions through the annual migration cycle. 

3\) Plotting and summarizing results from the Random Forest Model 



*Users need a Google Earth Engine account to run this script and need to upload the 'model_df_gee.csv' file and regional changes file in order to run the script. 
## Code/software

Users need R programming software and an active Google Earth Engine account to view and process these data. Both scripts are shared. 
## Access information

Other publicly accessible locations of the data:

* Not applicable

Data was derived from the following sources:

* Not applicable

