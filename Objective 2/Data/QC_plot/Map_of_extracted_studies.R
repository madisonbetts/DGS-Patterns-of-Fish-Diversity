# -------------------------------
# Plot all extracted study sites in the lower 48
# colored by Family from Study_metadata.xlsx
# works for repeated species codes like ONCL-1, ONCL-2, ONCL-3
# -------------------------------

library(readxl)
library(dplyr)
library(ggplot2)
library(maps)
library(stringr)

# -------------------------------
# paths
# -------------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"

metadata_file <- file.path(base_dir, "Study_metadata.xlsx")

# -------------------------------
# list valid study folders
# expects names like XXXX-1, XXXX-2, XXXX-3, etc.
# -------------------------------
subdirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
subdirs <- subdirs[grepl("^[A-Za-z0-9]+-[0-9]+$", basename(subdirs))]

# -------------------------------
# helper: identify lon/lat columns in a dataframe
# -------------------------------
find_lon_lat_cols <- function(df) {
  nms <- names(df)
  nms_lower <- tolower(nms)
  
  lon_idx <- grep("(^lon$|^long$|^longitude$|lon|long|longitude)", nms_lower)
  lat_idx <- grep("(^lat$|^latitude$|lat|latitude)", nms_lower)
  
  if (length(lon_idx) == 0 || length(lat_idx) == 0) return(NULL)
  
  out <- df[, c(nms[lon_idx[1]], nms[lat_idx[1]]), drop = FALSE]
  names(out) <- c("lon", "lat")
  
  out <- out %>%
    mutate(
      lon = suppressWarnings(as.numeric(lon)),
      lat = suppressWarnings(as.numeric(lat))
    ) %>%
    filter(!is.na(lon), !is.na(lat)) %>%
    filter(lon >= -130, lon <= -65, lat >= 24, lat <= 50)
  
  if (nrow(out) == 0) return(NULL)
  
  out %>%
    distinct(lon, lat)
}

# -------------------------------
# helper: load .RData from each folder/data/
# and extract best lon/lat dataframe
# preference is given to objects with "coords" in the name
# -------------------------------
extract_site_coords <- function(folder_path) {
  
  folder_name <- basename(folder_path)
  data_dir <- file.path(folder_path, "data")
  
  if (!dir.exists(data_dir)) {
    message("No data/ folder found in: ", folder_name)
    return(NULL)
  }
  
  rdata_files <- list.files(data_dir, pattern = "\\.RData$", full.names = TRUE)
  
  if (length(rdata_files) == 0) {
    message("No .RData file found in: ", folder_name, "/data")
    return(NULL)
  }
  
  tmp_env <- new.env()
  load(rdata_files[1], envir = tmp_env)
  obj_names <- ls(tmp_env)
  
  # prioritize coordinate-like object names
  obj_names <- c(
    obj_names[grepl("coord", tolower(obj_names))],
    obj_names[!grepl("coord", tolower(obj_names))]
  ) %>% unique()
  
  for (nm in obj_names) {
    obj <- get(nm, envir = tmp_env)
    
    if (is.data.frame(obj)) {
      coords <- find_lon_lat_cols(obj)
      if (!is.null(coords)) {
        coords$source_object <- nm
        coords$folder_name <- folder_name
        return(coords)
      }
    }
    
    if (is.list(obj) && !is.data.frame(obj)) {
      sub_nms <- names(obj)
      if (length(sub_nms) > 0) {
        sub_nms <- c(
          sub_nms[grepl("coord", tolower(sub_nms))],
          sub_nms[!grepl("coord", tolower(sub_nms))]
        ) %>% unique()
        
        for (sub_nm in sub_nms) {
          sub_obj <- obj[[sub_nm]]
          if (is.data.frame(sub_obj)) {
            coords <- find_lon_lat_cols(sub_obj)
            if (!is.null(coords)) {
              coords$source_object <- paste0(nm, "$", sub_nm)
              coords$folder_name <- folder_name
              return(coords)
            }
          }
        }
      }
    }
  }
  
  message("No suitable lon/lat dataframe found in: ", folder_name)
  return(NULL)
}

# -------------------------------
# extract coordinates from all study folders
# -------------------------------
site_list <- lapply(subdirs, extract_site_coords)

sites_df <- bind_rows(site_list) %>%
  mutate(
    study_code = folder_name,
    species_code = sub("-[0-9]+$", "", folder_name),
    dataset_num = as.integer(sub("^.*-([0-9]+)$", "\\1", folder_name))
  )

# -------------------------------
# read metadata
# -------------------------------
meta <- read_excel(metadata_file)
names(meta) <- trimws(names(meta))

code_col <- "Study_ID"

if (!code_col %in% names(meta)) {
  stop(paste0("Column '", code_col, "' not found in Study_metadata.xlsx"))
}

if (!"Family" %in% names(meta)) {
  stop("Column 'Family' not found in Study_metadata.xlsx")
}

meta2 <- meta %>%
  mutate(
    study_code = as.character(.data[[code_col]])
  ) %>%
  select(study_code, Family, everything())

# -------------------------------
# join metadata
# -------------------------------
sites_df <- sites_df %>%
  left_join(meta2, by = "study_code")

# -------------------------------
# report unmatched folders
# -------------------------------
unmatched <- sites_df %>%
  filter(is.na(Family)) %>%
  distinct(folder_name)

if (nrow(unmatched) > 0) {
  message("These studies did not match a Family in metadata:")
  print(unmatched)
}

# -------------------------------
# report duplicated species codes
# e.g. ONCL-1, ONCL-2, ONCL-3
# -------------------------------
duplicate_species <- sites_df %>%
  distinct(species_code, study_code) %>%
  count(species_code, name = "n_datasets") %>%
  filter(n_datasets > 1)

if (nrow(duplicate_species) > 0) {
  message("Species codes represented by multiple datasets:")
  print(duplicate_species)
  
  print(
    sites_df %>%
      distinct(species_code, study_code) %>%
      semi_join(duplicate_species, by = "species_code") %>%
      arrange(species_code, study_code)
  )
}

# -------------------------------
# map of lower 48
# -------------------------------
us_map <- map_data("state")

p <- ggplot() +
  geom_polygon(
    data = us_map,
    aes(x = long, y = lat, group = group),
    fill = "grey92",
    color = "grey55",
    linewidth = 0.2
  ) +
  geom_point(
    data = sites_df,
    aes(x = lon, y = lat, color = Family),
    size = 2.8,
    alpha = 0.9
  ) +
  coord_quickmap(
    xlim = c(-125, -66),
    ylim = c(24, 50)
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    color = "Family",
    title = "Extracted study sites across North America."
  )

print(p)

# ggsave(
#   filename = file.path(base_dir, "all_sites_lower48_by_family.png"),
#   plot = p,
#   width = 12,
#   height = 7,
#   dpi = 400
# )