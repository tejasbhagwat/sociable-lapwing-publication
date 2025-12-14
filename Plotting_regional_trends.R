#===============================================================================================================================
#' This script uses regional yearly median habitat suitability statistics exported from GEE to plot trends in habitat suitability 
#' over time for different regions and seasonal home ranges
#===============================================================================================================================
#--------------------------------------------------------------------------------------------------------------------------------
# Load necessary libraries
library(sf); library(dplyr); library(ggplot2); library(tidyverse); library(tidyr); library(stringr)
#--------------------------------------------------------------------------------------------------------------------------------
custom_palette <- c("#004488", "#BB5566", "#029E73")
#--------------------------------------------------------------------------------------------------------------------------------
# Load .csv file (This is from GEE)
regional_stats1 = read.csv("regional_changes/regional_changes_stats_I.csv")[-8]
regional_stats2 = read.csv("regional_changes/regional_changes_stats_II.csv")
regional_stats3 = read.csv("regional_changes/regional_changes_stats_III.csv")
regional_stats4 = read.csv("regional_changes/regional_changes_stats_IV.csv")
regional_stats5 = read.csv("regional_changes/regional_changes_stats_V.csv")
#--------------------------------------------------------------------------------------------------------------------------------
# Combine all dataframes into one
#--------------------------------------------------------------------------------------------------------------------------------
regional_stats = bind_rows(regional_stats1, regional_stats2, regional_stats3, regional_stats4, regional_stats5) %>% 
  select(-c(id, system.index)) %>% rename(median_occurrence_probability = median)
colnames(regional_stats)

#--------------------------------------------------------------------------------------------------------------------------------
# Replace "Tallymarzhan-Turkmenistan-Uzbekistan" region with "Tallymarzhan-Turkmenistan" 
regional_stats$region = gsub("Tallymarzhan-Turkmenistan-Uzbekistan", "Tallymarzhan-Turkmenistan", regional_stats$region)
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# Plot year vs median occurrence probability per ID using ggplot
ggplot(regional_stats, aes(x = year, y = median_occurrence_probability, color = season)) +
  # add a linear regression line
  geom_smooth(method = "loess", se = T, span = 0.75) +
  geom_point(size = 2, alpha = 0.7) +
  # Show all years on the x-axis 
  scale_x_continuous(breaks = seq(1995, 2025, by = 5)) +
  # Add custom color palette
  scale_color_manual(values = c("#BB5566","#004488","#029E73"),
                     labels = c("breed" = "Breeding",
                                "stopover" = "Stopover",
                                "winter" = "Wintering")) +
  labs(
    title = "",
    x = "",
    y = "Median value for habitat suitability"
  ) +
  theme_bw() +
  # Tilt x-axis labels to fit the years 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  theme(legend.position = "bottom") + 
  facet_wrap(~ region, scales = "free_y", ) + 
  # Reduce panel labels 
  theme(strip.text = element_text(size = 9, face = "bold"))
#--------------------------------------------------------------------------------------------------------------------------------
# Save the plot
#ggsave("output_file.png", width = 10, height = 8, dpi = 300)
#--------------------------------------------------------------------------------------------------------------------------------
