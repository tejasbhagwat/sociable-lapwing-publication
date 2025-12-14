#===============================================================================================================================
# This Script takes the annual habitat suitability layers generated from GEE and calculates the net change in habitat suitability
# over three decades: 1995-2004, 2005-2014, and 2015-2024. It then extracts the mean change for seasonal home ranges defined
# by the Birdlife International range map for the Sociable Lapwing and visualizes these changes using ggplot2. All the files
# required to run this script are available in the Dryad repository
#===============================================================================================================================
#--------------------------------------------------------------------------------------------------------------------------------
# Load necessary libraries
#--------------------------------------------------------------------------------------------------------------------------------
library(terra); 
library(sf); 
library(exactextractr); 
library(dplyr); 
library(tidyr)
library(ggplot2);
#===============================================================================================================================
#--------------------------------------------------------------------------------------------------------------------------------
# Separate Annual cycle by Birdlife range map for Sociable Lapwing
#--------------------------------------------------------------------------------------------------------------------------------
BL_annual_cycle = read_sf("birdlife_home_range/soc_lap_BirdLifeRange.gpkg") %>%
  dplyr::select(OBJECTID, label) %>% 
  st_make_valid() %>% 
  group_by(label) %>% 
  summarise()
#--------------------------------------------------------------------------------------------------------------------------------
#===============================================================================================================================
#' Load HSI median yearly layers exported from GEE. These are not shared due to their big sizes, but can be generated from
#' the GEE scripts provided in the repository
#===============================================================================================================================
#--------------------------------------------------------------------------------------------------------------------------------
# Load in the HSI rasters for years of interest
#--------------------------------------------------------------------------------------------------------------------------------
HSI_95 = rast(paste0(getwd(),"/annual_suitability_maps/HSI_1995.tif")); HSI_04 = rast(paste0(getwd(),"/annual_suitability_maps/HSI_2004.tif"))
HSI_05 = rast(paste0(getwd(),"/annual_suitability_maps/HSI_2005.tif")); HSI_14 = rast(paste0(getwd(),"/annual_suitability_maps/HSI_2014.tif"))
HSI_15 = rast(paste0(getwd(),"/annual_suitability_maps/HSI_2015.tif")); HSI_24 = rast(paste0(getwd(),"/annual_suitability_maps/HSI_2024.tif"))

#--------------------------------------------------------------------------------------------------------------------------------
# Calculate net change in HSI over the decades
#--------------------------------------------------------------------------------------------------------------------------------
HSI_net_change_1995_2004 = ((HSI_04 - HSI_95) / HSI_95)
HSI_net_change_2005_2014 = ((HSI_14 - HSI_05) / HSI_05)
HSI_net_change_2015_2024 = ((HSI_24 - HSI_15) / HSI_15)
#plot(HSI_net_change_2015_2024)
#--------------------------------------------------------------------------------------------------------------------------------
# Combine them 
#--------------------------------------------------------------------------------------------------------------------------------
HSI_net_change_stack = c(HSI_net_change_1995_2004, HSI_net_change_2005_2014, HSI_net_change_2015_2024)
names(HSI_net_change_stack) = c("NC_1995_2004","NC_2005_2014","NC_2015_2024")
#--------------------------------------------------------------------------------------------------------------------------------
# Extract mean change across BL polygons
BL_stat = data.frame(BL_annual_cycle, exact_extract(HSI_net_change_stack, BL_annual_cycle, fun = "mean"))
#--------------------------------------------------------------------------------------------------------------------------------
# Convert to long format for ggplot
BL_stat_long = BL_stat %>%
  st_drop_geometry() %>%
  pivot_longer(cols = starts_with("mean.NC_"), names_to = "year", values_to = "meanNC") %>%
  mutate(year = recode(year,
                       "mean.NC_1995_2004" = "2005",
                       "mean.NC_2005_2014" = "2015",
                       "mean.NC_2015_2024" = "2024"))
#--------------------------------------------------------------------------------------------------------------------------------
# Define custom color palette
custom_palette <- c("#004488", "#BB5566", "#029E73")
#--------------------------------------------------------------------------------------------------------------------------------
# Plot with ggplot
p1 = ggplot(BL_stat_long, aes(x = factor(year), y = meanNC*100, fill = label)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(title = "",
       x = "",
       y = "% Change in habitat suitability",
       fill = "") +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_discrete(labels = c("2005" = "1995-2004",
                              "2015" = "2005-2014",
                              "2024" = "2015-2024")) +
  scale_fill_manual(values = custom_palette) #+

#ggsave("output_file.png", p1, width = 6, height = 4, dpi = 300)
#--------------------------------------------------------------------------------------------------------------------------------