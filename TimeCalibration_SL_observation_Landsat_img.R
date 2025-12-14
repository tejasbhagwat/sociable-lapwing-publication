#================================================================================================================
# Time calibration of Sociable Lapwings observation data with Landsat satellite imagery
#================================================================================================================
#----------------------------------------------------------------------------------------------------------------
# Load libraries
library(dplyr); library(ggplot2); library(sf); library(tmap); 
library(lubridate);library(units); library(terra);
library(ggspatial); library(rnaturalearth); library(stringr); 
library(mapview); 
# Install the jsonlite package if not already installed
#if (!require(jsonlite)) install.packages("jsonlite")
#================================================================================================================
#' Loading combined occurrence points (PRESENCE POINTS) from 1990-2024 (movement data are already filtered for LL 
#' and LH), processed by HU to extract all TC indices for all Landsat images available. 
#================================================================================================================
# Load and filter data
SL_lsat_tc = read.csv("TasseledCap_All_Landsat_1990_2024.csv") %>% # TSS_soc_lap_allyears_1990-2023_LSAT.csv (Original file name)
             filter(NDVI != -9999, TCB != -9999, TCG != -9999, TCW != -9999) %>% dplyr::select(-.geo) %>% 
             # Change the obs_date format and change the satellite imagery date column name and format
             mutate(obs_date = as.Date(as.POSIXct(obs_date / 1000, origin = "1970-01-01", tz = "UTC", 
                                                  tryFormats = "%Y-%m-%d"))) %>%
             # Rename the sat_date column and change the format
             dplyr::rename(sat_date = YYYYMMDD) %>% mutate(sat_date = ymd(sat_date))
head(SL_lsat_tc)
#----------------------------------------------------------------------------------------------------------------
#' Match closest Landsat date to observation date +/- 45 days
#' (I am not using the mean value approach right now in order to make simpler to join random points 
#' for the same dates with summarized Landsat values, its going to be difficult to find the exact Landsat
#' scenes random points which are spread within 100-km from presence points. To go back to the old approach
#' refer to the old script 'time_calibration.R' in the Time calibration folder.)
#----------------------------------------------------------------------------------------------------------------
# Check how many total observations are there
length(unique(SL_lsat_tc$unique_id))

# Find the closest available Landsat image values (HU extracted all of them)
SL_lsat_tc_3_months = SL_lsat_tc %>% 
  mutate(date_diff = abs(sat_date - obs_date)) %>%                      # Calculate absolute difference between dates
  group_by(unique_id, obs_date, lat, lon) %>%                           # Group by ID and observation date
  slice_min(date_diff, n = 5, with_ties = F) %>%                        # Select the row with the minimum difference 
  filter(date_diff <= 45) %>%                                           # Filter records with date difference +/- 45 days
  summarise(ndvi_median = median(NDVI),                                 # Calculate median NDVI over 3-months
            tcb_median = median(TCB),                                   # Calculate median TCB over 3-months
            tcg_median = median(TCG),                                   # Calculate median TCG over 3 months
            tcw_median = median(TCW)) %>%                               # Calculate median TCW over 3 months
  ungroup()                                                             # Ungroup the data
  #dplyr::rename(sat_name = system.index) %>%                           # Rename the column that had the Landsat scene name reference
  #mutate(sat_name = str_trim(str_extract(sat_name, "L.*?_\\d{8}")))    # Extract only the scene name
#=================================================================================================================================
#' Processing random points
#' Like occurrence points these are also processed for TC indices. We took 10 points per presence point. This 
#' means there are at least 29,000 in total (2900 presence points * 10 random points). This generated thousands 
#' of rows as all Landsat scenes per random point are extracted. The goal now, is to reduce these files into 
#' a single file with each random point having the closest Landsat image date to the observation date. It shares
#' the observation dates and a unique ID with the presence points. So that should remain consistent. Most likely
#' it will also have the same Landsat scene ID that matches to that of the presence points. But in some cases
#' it maybe that the random point falls outside the extent of the scene from which values for presence point were
#' extracted... 
#=================================================================================================================================
#---------------------------------------------------------------------------------------------------------------------------------
# Load and filter data
#---------------------------------------------------------------------------------------------------------------------------------
SL_lsat_tc_rp = list.files("folder_path/Random_points/", pattern = "\\.csv$", full.names = T) # Not shared due to data size
SL_lsat_tc_rp = SL_lsat_tc_rp %>% lapply(read.csv) %>% bind_rows() %>% # THIS WILL BE A LOT OF DATA... 
  # Remove rows with missing values  
  filter(NDVI != -9999, TCB != -9999, TCG != -9999, TCW != -9999) %>%
  # Change the obs_date format and change the satellite imagery date column name and format
  mutate(obs_date = as.Date(as.POSIXct(obs_date / 1000, origin = "1970-01-01", tz = "UTC", tryFormats = "%Y-%m-%d"))) %>%
  # Rename the sat_date column and change the format
  dplyr::rename(sat_date = YYYYMMDD) %>% mutate(sat_date = ymd(sat_date))

#---------------------------------------------------------------------------------------------------------------------------------
# Match closest 5 Landsat dates to observation date
#---------------------------------------------------------------------------------------------------------------------------------
# Check how many total observations are there
length(unique(SL_lsat_tc_rp$id))
colnames(SL_lsat_tc_rp)  
# Find the closest available Landsat image values (HU extracted all of them)
SL_lsat_tc_rp_3_months = SL_lsat_tc_rp %>% 
  dplyr::select(-lat,-lon,-unique_id) %>% 
  mutate(#sat_name = str_trim(str_extract(sat_name, "L.*?_\\d{8}")),                    # Extract only the scene name
         lat = sapply(.geo, function(geo_json) fromJSON(geo_json)$coordinates[2]),      # Extract latitude from .geo column
         lon = sapply(.geo, function(geo_json) fromJSON(geo_json)$coordinates[1])) %>%  # Extract longitude from .geo column
  dplyr::select(-.geo) %>%
  mutate(date_diff = abs(sat_date - obs_date)) %>%                                      # Calculate absolute difference between dates
  group_by(id, obs_date, lat, lon) %>%                                                  # Columns to keep
  slice_min(date_diff, n = 5, with_ties = F) %>%                                        # Select the row with the minimum difference 
  filter(date_diff <= 45) %>%                                                           # Filter records with date difference less than 45 days
  # Summarise the data
  summarise(ndvi_median = median(NDVI),                                                 # Calculate median NDVI over 3-months
            tcb_median = median(TCB),                                                   # Calculate median TCB over 3-months
            tcg_median = median(TCG),                                                   # Calculate median TCG over 3 months
            tcw_median = median(TCW)) %>%                                               # Calculate median TCW over 3 months
  ungroup()                                                                             # Ungroup the data
                                                                                                   
#---------------------------------------------------------------------------------------------------------------------------------  
# Number of BG points where no imagery was found within 15 days of observation date
length(unique(SL_lsat_tc_rp_15_days$id)) # 

#=================================================================================================================================
#' Combine presence and random points
#' Due to the close vicinity of occurrence points and relatively large 100km buffer around each of them, random points
#' have overlapped with many adjacent buffers. This means random points have significantly less unique ids than 
#' he original data. This is fine, as long as we have roughly x10 total random points relative to the presence points
#=================================================================================================================================
# colnames(SL_lsat_tc_3_months); colnames(SL_lsat_tc_rp_3_months)
# Add presence-absence columns
SL_lsat_tc_3_months = SL_lsat_tc_3_months %>% select(obs_date, ndvi_median, tcb_median, tcg_median, tcw_median, 
                                                     lat, lon) %>% mutate(status = 1)
                                            
SL_lsat_tc_rp_3_months = SL_lsat_tc_rp_3_months %>% select(obs_date, ndvi_median, tcb_median, tcg_median, tcw_median,
                                                           lat, lon) %>% mutate(status = 0)
                                                  
# Combine the two
model_df = bind_rows(SL_lsat_tc_3_months, SL_lsat_tc_rp_3_months)
#=================================================================================================================================
# Write the final data
write.csv(model_df, "model_df_gee.csv", row.names = F)
#=================================================================================================================================