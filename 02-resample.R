library(terra)

# Read a base raster from WorldPop - this will be used to define the national boundary for Mozambique

base_r <- rast("data\\Geodata\\moz_buildings_count_BCB_gl_100m_v1_1.tif")


# See raster parameters 
crs(base_r)
ext(base_r)
res(base_r)

# Set spatial resolution
template_100m <- base_r

# Set input and output directories
input_dir <- "data\\Geodata"
output_dir <- "data\\Geodata_processed"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# Function to check CRS, resample to 100m resolution, and align extent and mask with the national boundary for Mozambique
resample_to_100m <- function(file_path) {
  
  r <- rast(file_path)
  
  # CRS
  if (!same.crs(r, template_100m)) {
    cat("CRS differs. Reprojecting.\n")
    r <- project(r, template_100m, method = "bilinear")
  }
  
  # Extent
  
  r <- crop(r, template_100m)
  
  # Resolution
  if (!all(abs(res(r) - res(template_100m)) < 1e-9)) {
    cat("Resolution differs. Resampling.\n")
    r <- resample(r, template_100m, method = "bilinear")
  }
  
  # Final geometry check (origin/alignment)
  if (!compareGeom(r, template_100m, stopOnError = FALSE)) {
    cat("Grid alignment differs. Resampling.\n")
    r <- resample(r, template_100m, method = "bilinear")
  }
  
  # Mask to get the exact shape of Mozambique
  r <- mask(r, template_100m)
  
  output_path <- file.path(output_dir, basename(file_path))
  writeRaster(r, output_path, overwrite = TRUE)
  
  return(output_path)
}

# List all tif files in the input directory
tif_files <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)

# Process all files
results <- lapply(tif_files, function(f) {
  tryCatch({
    message("\nProcessing: ", basename(f))
    resample_to_100m(f)
  }, error = function(e) {
    message("Error processing ", basename(f), ": ", e$message)
    return(NULL)
  })
})

# Print summary
successful <- sum(!sapply(results, is.null))
message("\nProcessing complete!")
message("Successfully processed ", successful, " out of ", length(tif_files), " files")

