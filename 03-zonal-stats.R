# Prepare population data for causal SAE
# First, assign treatment at EA level (5km buffer)
# You will need rasters with flood extent for Ana, Dumako, Gombe
# Then calculate zonal statistics from WorldPop for the geocovariates 
# And calculate the population for each EA 
# Also the flow accumulation measure for the propensity score

# Libraries
library(sf)
library(terra)
library(exactextractr)
library(dplyr)
library(tidyverse)


# Function to standardise rasters using z-scores 

std_rast <- function(r) {
  
  # Get mean and sd of the raster
  mean_val <- as.numeric(global(r, fun = "mean", na.rm = TRUE))
  sd_val <- as.numeric(global(r, fun = "sd", na.rm = TRUE))
  
  # Store values of a raster
  vals <- values(r)
  
  # Standardise raster values
  values(r) <- (vals - mean_val) / sd_val
  
  return(r)
}


# Load population-level EA shapefile (for analysis districts)
pop_ea <- st_read("./output/population_eas.shp")
sample_data <- readRDS("./output/data_eligible.rds")
flood_ana <- rast("./data/flood_ana.tif")
flood_dumako <- rast("./data/flood_dumako.tif")
flood_gombe <- rast("./data/flood_gombe.tif")
#flow_acc <- rast("./data/flow_accumulation.tif")
cov_dir <- "./data/Geodata_processed" 

# List all WorldPop covariates
tif_files <- list.files(
  cov_dir,
  pattern = "\\.tif$",
  full.names = TRUE
)

# Stack all geodata
cov_stack <- rast(tif_files)
cov_stack
names(cov_stack)

# Check CRS
crs(cov_stack)
st_crs(pop_ea)

# Treat population layer separately - it needs the sum function 
cov_stack_nopop <- cov_stack[[names(cov_stack) != "moz_pop_2022_CN_100m_R2025A_v1"]]
pop_layer <- cov_stack[[names(cov_stack) == "moz_pop_2022_CN_100m_R2025A_v1"]]

# Standardise all covariates WorldPop
#r_std_list <- lapply(1:nlyr(cov_stack_nopop), function(i) std_rast(cov_stack_nopop[[i]]))

# Combine back to a SpatRaster
#cov_stack_std <- rast(cov_stack_nopop)

# Zonal statistics - mean values of WorldPop covariates
zonal_means <- exact_extract(
  cov_stack_nopop,
  pop_ea,
  "mean",
  progress = TRUE
)

# Add AES identifier
zonal_means <- bind_cols(
  pop_ea %>% st_drop_geometry() %>% select(AES),
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

# # Mean for flow accumulation
# zonal_fa <- exact_extract(
#   flow_acc_std,
#   pop_ea,
#   "mean",
#   progress = T
# )
# 
# # Add AES identifier 
# zonal_fa <- bind_cols(
#   pop_ea %>% st_drop_geometry() %>% select(AES),
#   zonal_fa
# ) 
# 
# # Change name to be consistent with sample data
# names(zonal_fa)[2] <- "flow_acc"

# Create a 5km buffer around EAs 
pop_ea_5km <- pop_ea %>%
  st_transform(32736) %>%  # reproject to meters before calculating buffers
  st_buffer(5000) %>%
  st_transform(st_crs(pop_ea))

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
  left_join(zonal_max) #%>%
  #left_join(zonal_fa)

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
  select(-flow_acc) %>%
  left_join(zonal_means) %>%
  left_join(zonal_sum)# %>%
  #left_join(zonal_fa)


# Save
saveRDS(
  sample_data_cov,
  file = "./output/sample_data_cov_std.rds"
)

saveRDS(
  pop_ea_cov,
  file = "./output/pop_ea_cov_std.rds"
)

saveRDS(
  pop_ea_5km,
  file= "./output/pop_ea_5km.rds"
)