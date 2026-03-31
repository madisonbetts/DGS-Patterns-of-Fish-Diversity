# -------------------------------
# Plot all extracted study sites across the U.S. + Canada
# colored by Family from Study_metadata.xlsx
# legend entries include dataset counts per family, e.g. Leuciscidae (x23)
# legend title = Dataset (n = total number of datasets)
# adds north arrow + default scale bar + black inset box
# saves as png + pdf
# -------------------------------

library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(grid)

# -------------------------------
# paths
# -------------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"

metadata_file <- file.path(base_dir, "Study_metadata.xlsx")

out_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/analyses/QC_plot"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

png_file <- file.path(out_dir, "all_sites_us_canada_by_family.png")
pdf_file <- file.path(out_dir, "all_sites_us_canada_by_family.pdf")

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
# family counts for legend entries
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

# total datasets for legend title
n_datasets <- sites_df %>%
  distinct(study_code) %>%
  nrow()

legend_title <- paste0("Datasets (n = ", n_datasets, ")")

# keep legend order consistent with counts
sites_df$Family <- factor(sites_df$Family, levels = family_counts$Family)

# common families plotted first; rarer families stay visible on top
sites_df <- sites_df %>%
  left_join(family_counts, by = "Family") %>%
  arrange(desc(n_datasets), Family)

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

states_provinces_sf <- states_provinces_sf %>%
  filter(
    !name_en %in% c(
      "Alaska", "Hawaii",
      "Yukon", "Northwest Territories", "Nunavut"
    )
  )

# -------------------------------
# map extent
# -------------------------------
map_xlim <- c(-130, -52)
map_ylim <- c(23, 60)

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
    size = 0.75,
    alpha = 0.9
  ) +
  coord_sf(
    xlim = map_xlim,
    ylim = map_ylim,
    expand = FALSE
  ) +
  # manual scale bar in lower-left white space
  annotate("rect", xmin = -128.2, xmax = -125.7, ymin = 24.4, ymax = 25.6,
           fill = "black", color = "black") +
  annotate("rect", xmin = -125.7, xmax = -123.2, ymin = 24.4, ymax = 25.6,
           fill = "white", color = "black") +
  annotate("text",
           x = -122.0, y = 25.0,
           label = "1000 km",
           hjust = 0, size = 3.2) +
  # north arrow
  #annotation_north_arrow(
  #  location = "tr",
  #  which_north = "true",
  #  height = unit(0.15, "in"),
  #  width = unit(0.15, "in"),
  #  pad_x = unit(0.10, "in"),
  #  pad_y = unit(0.10, "in"),
  #  style = north_arrow_orienteering(
  #    fill = c("black", "black"))
  #) +
  # north arrow v2
# north arrow
annotation_north_arrow(
  location = "tr",
  which_north = "true",
  height = unit(0.22, "in"),
  width  = unit(0.22, "in"),
  pad_x  = unit(0.10, "in"),
  pad_y  = unit(0.10, "in"),
  style = north_arrow_orienteering(
    fill = c("black", "black"),
    line_col = "black",
    text_col = NA   # suppress built-in N text
  )
) + # north arrow text "N"
  #annotate(
  #  "text",
  #  x = -55.2, y = 56.3,   # tweak if needed
  #  label = "N",
  #  size = 3.2
  #) + 
  # color scale
  scale_color_discrete(
    breaks = family_counts$Family,
    labels = family_labels,
    drop = FALSE,
  ) +
  theme_classic() +
  theme(
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.key.size = unit(0.1, "in")
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    color = legend_title
  )

print(p)

# -------------------------------
# save
# -------------------------------
ggsave(
  filename = png_file,
  plot = p,
  width = 6.5,
  height = 4,
  dpi = 400
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 6.5,
  height = 4,
  device = cairo_pdf
)

# -------------------------------
# summary counts
# -------------------------------
n_sites <- sites_df %>%
  distinct(lon, lat) %>%
  nrow()

n_species_codes <- sites_df %>%
  distinct(species_code) %>%
  nrow()

cat("\n==================== SUMMARY ====================\n")
cat("Datasets: ", n_datasets, "\n", sep = "")
cat("Species codes: ", n_species_codes, "\n", sep = "")
cat("Unique sites: ", n_sites, "\n", sep = "")
cat("PNG saved to: ", png_file, "\n", sep = "")
cat("PDF saved to: ", pdf_file, "\n", sep = "")
cat("=================================================\n\n")
