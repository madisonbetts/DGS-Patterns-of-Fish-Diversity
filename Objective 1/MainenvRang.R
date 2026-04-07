# =============================
# INSTALL + LOAD PACKAGES
# =============================
# install.packages(c("sf","dplyr","FedData","elevatr","terra","readr","stringr","readxl"))

library(sf)
library(dplyr)
library(FedData)
library(elevatr)
library(terra)
library(readr)
library(stringr)
library(readxl)

# =============================
# SET WORKING DIRECTORY
# =============================
setwd("/Users/parinaz/Desktop/MR")

# =============================
# LOAD AND CLEAN DATA
# =============================
fish <- read_csv("filtered_data.csv")

fish <- fish %>%
  filter(!is.na(Longitude), !is.na(Latitude))

fish_sf <- st_as_sf(
  fish,
  coords = c("Longitude","Latitude"),
  crs = 4326
)

sf_use_s2(FALSE)

fish_sf <- st_make_valid(fish_sf)
fish_sf <- fish_sf[!st_is_empty(fish_sf), ]

# =============================
# DOWNLOAD NHD FLOWLINES (USA)
# =============================
bbox <- st_bbox(fish_sf)
bbox_poly <- st_as_sfc(bbox)
bbox_poly <- st_sf(geometry = bbox_poly)

nhd_data <- get_nhd(
  template = bbox_poly,
  label = "fish_nhd",
  extraction.dir = tempdir()
)

streams <- nhd_data$Flowline

streams <- st_make_valid(streams)
streams <- streams[!st_is_empty(streams), ]

# Keep real streams only
if ("FType" %in% names(streams)) {
  streams <- streams %>% filter(FType == 460)
}

# =============================
# PROJECT DATA (for distances)
# =============================
fish_proj <- st_transform(fish_sf, 5070)
streams_proj <- st_transform(streams, 5070)

# Ensure valid geometries
fish_proj <- st_make_valid(fish_proj)
streams_proj <- st_make_valid(streams_proj)

# Ensure same CRS
streams_proj <- st_transform(streams_proj, st_crs(fish_proj))

# =============================
# HANDLE STREAM LENGTH SAFELY
# =============================
if (!"lengthkm" %in% names(streams_proj)) {
  streams_proj$lengthkm <- as.numeric(st_length(streams_proj)) / 1000
}

streams_proj <- streams_proj[!is.na(streams_proj$lengthkm), ]

# =============================
# FIND NEAREST STREAM
# =============================
nearest_idx <- st_nearest_feature(fish_proj, streams_proj)

matched_length <- streams_proj$lengthkm[nearest_idx]

print(summary(matched_length))

# Log-transform stream size
fish_proj$stream_size <- log1p(matched_length)

# =============================
# ELEVATION + SLOPE
# =============================
bbox_dem <- st_bbox(fish_sf)
bbox_poly_dem <- st_as_sfc(bbox_dem)
bbox_poly_dem <- st_sf(geometry = bbox_poly_dem)

dem <- tryCatch({
  get_elev_raster(
    locations = bbox_poly_dem,
    z = 11,
    clip = "bbox"
  )
}, error = function(e) {
  message("⚠️ z=11 failed, switching to z=9")
  get_elev_raster(
    locations = bbox_poly_dem,
    z = 9,
    clip = "bbox"
  )
})

dem <- terra::rast(dem)

slope_raster <- terrain(
  dem,
  v = "slope",
  unit = "degrees"
)

# Extract using projected coordinates (IMPORTANT FIX)
fish_proj$elevation <- terra::extract(dem, fish_proj)[,2]
fish_proj$slope <- terra::extract(slope_raster, fish_proj)[,2]

# =============================
# MERGE BACK TO ORIGINAL SF
# =============================
fish_sf$stream_size <- fish_proj$stream_size
fish_sf$elevation <- fish_proj$elevation
fish_sf$slope <- fish_proj$slope

# =============================
# SAVE OUTPUT
# =============================
write.csv(
  st_drop_geometry(fish_sf),
  "fish_env_variables.csv",
  row.names = FALSE
)

st_write(
  fish_sf,
  "fish_env_variables.gpkg",
  delete_dsn = TRUE
)

# =============================
# FINAL CHECKS
# =============================
print(summary(fish_sf$stream_size))
print(summary(fish_sf$elevation))
print(summary(fish_sf$slope))

# =============================
# ADD RANGE SIZE (FishTraits)
# =============================
fish_env <- read_csv("fish_env_variables.csv")

# Clean species names
fish_env <- fish_env %>%
  mutate(Spec_Latin_clean = str_to_lower(str_trim(Spec_Latin_GenDivRange)))

# Load FishTraits
fish_traits <- read_excel("FishTraits_14.3.xls")

fish_traits <- fish_traits %>%
  mutate(ScientificName = str_to_lower(str_trim(paste(GENUS, SPECIES))))

fish_traits_sub <- fish_traits %>%
  select(ScientificName, AREAKM2)

# Merge range size
fish_env <- fish_env %>%
  left_join(fish_traits_sub, by = c("Spec_Latin_clean" = "ScientificName"))

# Check results
print(summary(fish_env$AREAKM2))

# Save final dataset
write_csv(fish_env, "fish_env_with_range.csv")