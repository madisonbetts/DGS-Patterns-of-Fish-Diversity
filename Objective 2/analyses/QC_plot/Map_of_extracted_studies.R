# -------------------------------
# Plot all extracted study sites across the U.S. + Canada
# colored by Family from Study_metadata.xlsx
# legend includes dataset counts per family, e.g. Leuciscidae (x23)
# -------------------------------

library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

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
    filter(lon >= -170, lon <= -50, lat >= 20, lat <= 85)
  
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
  left_join(meta2, by = "study_code") %>%
  mutate(Family = ifelse(is.na(Family), "Unknown", Family))

# -------------------------------
# report unmatched folders
# -------------------------------
unmatched <- sites_df %>%
  filter(Family == "Unknown") %>%
  distinct(folder_name)

if (nrow(unmatched) > 0) {
  message("These studies did not match a Family in metadata:")
  print(unmatched)
}

# -------------------------------
# family counts for legend and summary
# dataset-level counts per family
# -------------------------------
family_counts <- sites_df %>%
  distinct(study_code, species_code, Family) %>%
  count(Family, name = "n_datasets") %>%
  arrange(desc(n_datasets), Family)

family_labels <- setNames(
  paste0(family_counts$Family, " (x", family_counts$n_datasets, ")"),
  family_counts$Family
)

# keep legend order consistent with counts
sites_df$Family <- factor(sites_df$Family, levels = family_counts$Family)

# -------------------------------
# convert sites to sf
# -------------------------------
sites_sf <- st_as_sf(
  sites_df,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

# -------------------------------
# get map layers
# country polygons for U.S. + Canada
# state/province internal boundaries
# -------------------------------
countries_sf <- ne_countries(
  country = c("United States of America", "Canada"),
  scale = "medium",
  returnclass = "sf"
)

states_provinces_sf <- ne_states(
  country = c("United States of America", "Canada"),
  returnclass = "sf"
)

# drop AK/HI/territories if you want a cleaner lower-48 + southern Canada map
states_provinces_sf <- states_provinces_sf %>%
  filter(
    !name_en %in% c(
      "Alaska", "Hawaii",
      "Yukon", "Northwest Territories", "Nunavut"
    )
  )

# -------------------------------
# plot
# -------------------------------
p <- ggplot() +
  geom_sf(
    data = countries_sf,
    fill = "grey94",
    color = "grey40",
    linewidth = 0.3
  ) +
  geom_sf(
    data = states_provinces_sf,
    fill = NA,
    color = "grey65",
    linewidth = 0.2
  ) +
  geom_sf(
    data = sites_sf,
    aes(color = Family),
    size = 2.2,
    alpha = 0.9
  ) +
  coord_sf(
    xlim = c(-130, -52),
    ylim = c(23, 60),
    expand = FALSE
  ) +
  scale_color_discrete(
    labels = family_labels
  ) +
  theme_classic() +
  theme(
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9)
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    color = "Family",
    title = "Extracted sites from North American freshwater fish riverscape genetics studies"
  )

print(p)

# -------------------------------
# optional save
# -------------------------------
#ggsave(
#  filename = file.path(base_dir, "all_sites_us_canada_by_family.png"),
#  plot = p,
#  width = 12,
#  height = 8,
#  dpi = 400
#)

# -------------------------------
# summary counts
# -------------------------------
n_sites <- sites_df %>%
  distinct(lon, lat) %>%
  nrow()

n_datasets <- sites_df %>%
  distinct(study_code) %>%
  nrow()

n_species_codes <- sites_df %>%
  distinct(species_code) %>%
  nrow()

cat("\n==================== SUMMARY ====================\n")
cat("Datasets: ", n_datasets, "\n", sep = "")
cat("Species codes: ", n_species_codes, "\n", sep = "")
cat("Unique sites: ", n_sites, "\n", sep = "")
cat("\nDataset counts per family:\n")

for (i in seq_len(nrow(family_counts))) {
  cat("  ", family_counts$Family[i], ": ", family_counts$n_datasets[i], "\n", sep = "")
}

cat("=================================================\n\n")