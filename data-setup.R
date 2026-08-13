# Libraries
library(dplyr)
library(tidyverse)
library(sf)
library(ggplot2)
library(patchwork)

# File paths 
output_path <- "Z:/01 Mozambique/02 Code/moz-hh-impact-study/Data/Output"
shp_path <- "Z:/01 Mozambique/02 Code/descriptives-flood-poverty/results"

# Import analysis data 
data_for_analysis <- readRDS(paste0(output_path, "/data_for_analysis_v2.rds"))
posto_shp <- st_read(paste0(shp_path, "/moz_admin3_shp.shp"))
district_shp <- st_read(paste0(shp_path, "/moz_admin2_shp.shp"))
aes_2022 <- st_read("Z:/01 Mozambique/02 Code/moz-hh-impact-study/Data/Output/aes_2017_hh_2022/aes_2017_hh_2022.shp")
aes_census_2017 <- st_read("Z:/01 Mozambique/02 Code/moz-hh-impact-study/Data/Output/aes_2017/aes_2017.shp")

# Load the data for analysis:
# - excludes the interviews before cyclones Ana, Dumako, and Gombe
# - has the flood exposure measures for each cyclone (1km, 5km, 10km, 20km)
# - treatment variables is in (binary - 0/1) being affected by any of the three
# - EA information which can be joined with geodata
# - Be careful as this dataset also restricts the sample to before November 2022
# - (but this might not matter if reference period is flood in the last 6 months)
# - Variable named consumption is the World Bank aggregate

data_model <- data_for_analysis$data_model

# Keep only 2022 
data_2022 <- data_model %>%
  filter(
    post == 1
  )

# Add posto and district identifiers 

# First 6 digits = Posto Administrativo
data_2022$PA_ID <- substr(data_2022$BA_ID, 1, 6)

# First 4 digits = District
data_2022$ID_DIST <- substr(data_2022$BA_ID, 1, 4)

# Reorder geographical identifiers to the front
data_2022 <- data_2022 %>%
  relocate(ID_DIST, PA_ID, BA_ID)

### 1. Summary stats ###

# flood_any_5km_no_eloise - treatment variable of interest

# Number of districts per province 
freq_table <- data_2022 %>%
  distinct(province_name, ID_DIST) %>%  
  count(province_name, name = "n_districts")

freq_table

# Nampula: 22 districts, Sofala: 12 districts, Tete: 15 districts,
# Zambezia: 22 districts (71 total)


# Now treatment and control split
eligible_districts <- data_2022 %>%
  group_by(ID_DIST) %>%
  summarise(
    n_groups = n_distinct(flood_any_5km_no_eloise),
    .groups = "drop"
  ) %>%
  filter(n_groups == 2) %>%
  pull(ID_DIST)

# 32 districts left after those with all treatment or all control eliminated
district_treat <- data_2022 %>%
  filter(ID_DIST %in% eligible_districts) %>%
  group_by(ID_DIST, flood_any_5km_no_eloise) %>%
  summarise(
    sample_size = n_distinct(hhid),
    .groups = "drop"
  )

# Dataset with eligible districts
data_eligible <- data_2022 %>%
  filter(ID_DIST %in% eligible_districts)

# Panel 1: all households (treated + control combined)
all_units <- data_eligible %>%
  group_by(ID_DIST) %>%
  summarise(
    sample_size = n_distinct(hhid),
    .groups = "drop"
  )

# Panel 2 & 3: control and treated separately
district_treat <- data_eligible %>%
  group_by(ID_DIST, flood_any_5km_no_eloise) %>%
  summarise(
    sample_size = n_distinct(hhid),
    .groups = "drop"
  )

control_units <- district_treat %>%
  filter(flood_any_5km_no_eloise == 0)

treated_units <- district_treat %>%
  filter(flood_any_5km_no_eloise == 1)


# Plots with sample size 
p1 <- ggplot(all_units,
             aes(x = sample_size)) +
  geom_histogram(
    binwidth = 10,
    fill = "#dceaf6",
    color = "#4c63b6"
  ) +
  labs(
    title = "Treated and control units",
    x = "Sample size",
    y = "Number of districts"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p1
ggsave(
  filename = "./plots/total_sample_dist.png",
  plot = p1
)

p2 <- ggplot(control_units,
             aes(x = sample_size)) +
  geom_histogram(
    binwidth = 10,
    fill = "#ffe85c",
    color = "#ffbf00"
  ) +
  labs(
    title = "Control units",
    x = "Sample size",
    y = "Number of districts"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p2
ggsave(
  filename = "./plots/control_sample_dist.png",
  plot = p2
)


p3 <- ggplot(treated_units,
             aes(x = sample_size)) +
  geom_histogram(
    binwidth = 10,
    fill = "#66c24a",
    color = "#2f8f2f"
  ) +
  labs(
    title = "Treated units",
    x = "Sample size",
    y = "Number of districts"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p3
ggsave(
  filename = "./plots/treat_sample_dist.png",
  plot = p3
)

# Filter the shapefile with only affected districts
adm2_flood <- district_shp %>%
  filter(ID_DIST %in% eligible_districts)

# Save shapefile
st_write(adm2_flood, "./output/adm2_flood.shp",
         append = F)

# Save shapefile with analysis EAs 
eas <- unique(data_eligible$AES)

aes_analysis <- aes_2022 %>%
  filter(AES %in% eas)

st_write(aes_analysis, "./output/aes_analysis.shp",
         append = F)

## Summary stats ##

district_diff <- data_eligible %>%
  group_by(ID_DIST, flood_any_5km_no_eloise) %>%
  summarise(
    mean_consumption = mean(log_cons, na.rm = TRUE),
    sd_consumption = sd(log_cons, na.rm = TRUE),
    n_households = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = flood_any_5km_no_eloise,
    values_from = c(mean_consumption,
                    sd_consumption,
                    n_households)
  ) %>%
  mutate(
    diff = mean_consumption_1 - mean_consumption_0,
    
    # Standard error of difference in means
    se_diff = sqrt(
      (sd_consumption_1^2 / n_households_1) +
        (sd_consumption_0^2 / n_households_0)
    ),
    
    # 95% CI
    ci_lower = diff - 1.96 * se_diff,
    ci_upper = diff + 1.96 * se_diff
  )
mean(district_diff$diff, na.rm = TRUE)

# Export the table with differences in log consumption between T and C 
write.csv(
  district_diff,
  "./output/district_diff.csv",
  row.names = FALSE
)


# Plot the difference and confidence intervals

ggplot(
  district_diff,
  aes(
    y = reorder(ID_DIST, diff),
    x = diff
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "red"
  ) +
  geom_errorbarh(
    aes(xmin = ci_lower, xmax = ci_upper),
    height = 0.2
  ) +
  geom_point(size = 2) +
  labs(
    x = "Difference in mean log consumption (Treated - Control)",
    y = "District",
    title = "District-level treatment effects with 95% confidence intervals"
  ) +
  theme_bw()

# Export the data file for analysis 
saveRDS(
  data_eligible,
  file = "./output/data_eligible.rds"
)

# Find EAs that fall into analysis districts

# Create centroids
ea_centroids <- st_centroid(aes_census_2017)

# Keep EAs whose centroid falls inside an analysis district
eas_in_analysis <- aes_census_2017[
  lengths(st_within(ea_centroids, adm2_flood)) > 0,
]

# Export shapefile
st_write(eas_in_analysis, "./output/population_eas.shp",
         append = F)