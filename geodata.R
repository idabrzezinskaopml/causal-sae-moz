# Prepare population data for causal SAE
# First, assign treatment at EA level (5km buffer)
# You will need rasters with flood extent for Ana, Dumako, Gombe
# Then calculate zonal statistics from WorldPop for the geocovariates 
# And calculate the population for each EA 
# Also the flow accumulation measure for the propensity score
# Remember to also use 

# Libraries
library(sf)
library(terra)
library(exactextractr)
library(dplyr)
library(tidyverse)

# Load population-level EA shapefile (for analysis districts)
pop_ea <- st_read("./output/population_eas.shp")
sample_data <- readRDS("./output/data_eligible.rds")
flood_ana <- rast("./data/flood_ana.tif")
flood_dumako <- rast("./data/flood_dumako.tif")
flood_gombe <- rast("./data/flood_gombe.tif")
flow_acc <- rast("./data/flow_accumulation.tif")


# Load WorldPop covariates
cov_dir <- "C:/Users/idabr/Southampton/causal-sae/data/WorldPop covariates"

tif_files <- list.files(
  cov_dir,
  pattern = "\\.tif$",
  full.names = TRUE
)

cov_stack <- rast(tif_files)

cov_stack
names(cov_stack)

# Check CRS
crs(cov_stack)
st_crs(pop_ea)

# Treat population layer separately - it needs the sum function 
cov_stack_nopop <- cov_stack[[names(cov_stack) != "moz_pop_2022_CN_100m_R2025A_v1"]]
pop_layer <- cov_stack[[names(cov_stack) == "moz_pop_2022_CN_100m_R2025A_v1"]]

# Zonal statistics - mean values of WorldPop covariates
zonal_means <- exact_extract(
  cov_stack_nopop,
  pop_ea,
  "mean",
  progress = TRUE
)

# Add AES identifier
zonal_means <- bind_cols(
  pop_ea_5km %>% st_drop_geometry() %>% select(AES),
  zonal_means
)

# Sum for population (WorldPop)
zonal_sum <- exact_extract(
  pop_layer,
  pop_ea,
  "sum",
  progress = T
)

# Add AES identifier 
zonal_sum <- bind_cols(
  pop_ea %>% st_drop_geometry() %>% select(AES),
  zonal_sum
) 

# Rename second column
names(zonal_sum)[2] <- "population"

# Mean for flow accumulation
zonal_fa <- exact_extract(
  flow_acc,
  pop_ea,
  "mean",
  progress = T
)

# Add AES identifier 
zonal_fa <- bind_cols(
  pop_ea %>% st_drop_geometry() %>% select(AES),
  zonal_fa
) 

# Change name to be consistent with sample data
names(zonal_fa)[2] <- "flow_acc"


# Create a 5km buffer around EAs 
pop_ea_5km <- st_buffer(pop_ea, dist = 5000)


# List flood rasters
flood_rasters <- list(
  ana = flood_ana,
  dumako = flood_dumako,
  gombe = flood_gombe
)

# Calculate zonal max for each raster
zonal_max <- map_dfc(
  flood_rasters,
  exact_extract,
  y = pop_ea_5km, # 5km buffer
  fun = "max"
)

# Recode NA as 0s (means no flood)
zonal_max[is.na(zonal_max)] <- 0

# Check result
summary(zonal_max)

# Add AES identifier 
zonal_max <- bind_cols(
  pop_ea %>% st_drop_geometry() %>% select(AES),
  zonal_max
) 

# Rename so they match sample data
zonal_max <- zonal_max %>%
  rename(
    flood_ana_5km = ana,
    flood_dumako_5km = dumako,
    flood_gombe_5km = gombe
  )

# Left join zonal stats to the households
pop_ea_cov <- pop_ea %>%
  left_join(zonal_means) %>%
  left_join(zonal_sum) %>%
  left_join(zonal_max) %>%
  left_join(zonal_fa)

# Create binary for being affected by at least one flood 
pop_ea_cov <- pop_ea_cov %>%
  mutate(
    flood_any_5km_no_eloise = if_else(
      flood_ana_5km + flood_dumako_5km + flood_gombe_5km > 0,
      1,
      0
    )
  )


### Sample data ###

# Left join zonal stats to the households
sample_data_cov <- sample_data %>%
  left_join(zonal_means) %>%
  left_join(zonal_sum)


# Save
saveRDS(
  sample_data_cov,
  file = "./output/sample_data_cov.rds"
)

saveRDS(
  pop_ea_cov,
  file = "./output/pop_ea_cov.rds"
)

saveRDS(
  pop_ea_5km,
  file= "./output/pop_ea_5km.rds"
)