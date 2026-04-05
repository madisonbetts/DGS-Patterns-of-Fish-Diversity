# -------------------------------
# Plot all extracted study sites across the U.S. + Canada
# colored by Family from Study_metadata.xlsx
# legend moved inside lower-right of the map panel
# unique spp / unique sites moved to upper-left
# keeps panel border, north arrow, and scale bar
# smaller legend + extra left margin so y-axis labels are not clipped
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

png_file <- file.path(out_dir, "all_sites_by_family_inset_legend_small.png")
pdf_file <- file.path(out_dir, "all_sites_by_family_inset_legend_small.pdf")

# -------------------------------
# list valid study folders
# expects names like XXXX-1, XXXX-2, XXXX-3, etc.
# -------------------------------
subdirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
subdirs <- subdirs[grepl("^[A-Za-z0-9]+-[0-9]+$", basename(subdirs))]
subdirs <- sort(subdirs)

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
# helper: choose/load readable RData robustly
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
  
  expected_file <- file.path(data_dir, paste0(folder_name, ".RData"))
  if (file.exists(expected_file)) {
    rdata_files <- c(expected_file, setdiff(rdata_files, expected_file))
  }
  
  for (rf in rdata_files) {
    
    finfo <- file.info(rf)
    if (is.na(finfo$size) || finfo$size == 0) {
      message("Skipping empty .RData in ", folder_name, ": ", basename(rf))
      next
    }
    
    tmp_env <- new.env()
    
    ok <- tryCatch({
      load(rf, envir = tmp_env)
      TRUE
    }, error = function(e) {
      message("Skipping unreadable .RData in ", folder_name, ": ",
              basename(rf), " | ", conditionMessage(e))
      FALSE
    })
    
    if (!ok) next
    
    obj_names <- ls(tmp_env)
    expected_obj <- paste0(gsub("-", "_", folder_name), "_coords")
    
    obj_names <- c(
      expected_obj,
      obj_names[grepl("_coords$", obj_names, ignore.case = TRUE)],
      obj_names[grepl("coord", obj_names, ignore.case = TRUE)],
      obj_names
    ) %>% unique()
    
    obj_names <- obj_names[obj_names %in% ls(tmp_env)]
    
    for (nm in obj_names) {
      obj <- get(nm, envir = tmp_env)
      
      if (is.data.frame(obj)) {
        coords <- find_lon_lat_cols(obj)
        if (!is.null(coords)) {
          coords$source_object <- nm
          coords$folder_name <- folder_name
          coords$rdata_file <- basename(rf)
          return(coords)
        }
      }
      
      if (is.list(obj) && !is.data.frame(obj)) {
        sub_nms <- names(obj)
        if (length(sub_nms) > 0) {
          sub_nms <- c(
            sub_nms[grepl("_coords$", sub_nms, ignore.case = TRUE)],
            sub_nms[grepl("coord", sub_nms, ignore.case = TRUE)],
            sub_nms
          ) %>% unique()
          
          for (sub_nm in sub_nms) {
            sub_obj <- obj[[sub_nm]]
            if (is.data.frame(sub_obj)) {
              coords <- find_lon_lat_cols(sub_obj)
              if (!is.null(coords)) {
                coords$source_object <- paste0(nm, "$", sub_nm)
                coords$folder_name <- folder_name
                coords$rdata_file <- basename(rf)
                return(coords)
              }
            }
          }
        }
      }
    }
  }
  
  message("No suitable readable lon/lat dataframe found in: ", folder_name)
  return(NULL)
}

# -------------------------------
# extract coordinates from all study folders
# -------------------------------
site_list <- lapply(subdirs, extract_site_coords)
site_list <- site_list[!vapply(site_list, is.null, logical(1))]

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
    study_code = trimws(as.character(.data[[code_col]]))
  ) %>%
  select(study_code, Family, everything())

# -------------------------------
# join metadata
# -------------------------------
sites_df <- sites_df %>%
  left_join(meta2, by = "study_code") %>%
  mutate(Family = ifelse(is.na(Family), "Unknown", Family))

family_counts <- sites_df %>%
  distinct(study_code, species_code, Family) %>%
  count(Family, name = "n_datasets") %>%
  arrange(desc(n_datasets), Family)

family_labels <- setNames(
  paste0(family_counts$Family, " (x", family_counts$n_datasets, ")"),
  family_counts$Family
)

n_datasets <- sites_df %>%
  distinct(study_code) %>%
  nrow()

n_sites <- sites_df %>%
  distinct(lon, lat) %>%
  nrow()

n_species_codes <- sites_df %>%
  distinct(species_code) %>%
  nrow()

legend_title <- paste0("Datasets (n = ", n_datasets, ")")

summary_label <- paste(
  paste0("Unique spp: ", n_species_codes),
  paste0("Unique sites: ", n_sites),
  sep = "\n"
)

sites_df$Family <- factor(sites_df$Family, levels = family_counts$Family)

sites_df <- sites_df %>%
  left_join(family_counts, by = "Family") %>%
  arrange(desc(n_datasets), Family)

sites_sf <- st_as_sf(
  sites_df,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

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

map_xlim <- c(-130, -52)
map_ylim <- c(23, 60)

summary_x <- -128.7
summary_y <- 59.2

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
    expand = FALSE,
    clip = "on"
  ) +
  annotate(
    "rect",
    xmin = -128.2, xmax = -125.7, ymin = 24.4, ymax = 25.6,
    fill = "black", color = "black"
  ) +
  annotate(
    "rect",
    xmin = -125.7, xmax = -123.2, ymin = 24.4, ymax = 25.6,
    fill = "white", color = "black"
  ) +
  annotate(
    "text",
    x = -122.0, y = 25.0,
    label = "1000 km",
    hjust = 0, size = 3.0
  ) +
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
      text_col = NA
    )
  ) +
  annotate(
    "text",
    x = summary_x,
    y = summary_y,
    label = summary_label,
    hjust = 0,
    vjust = 1,
    size = 2.9
  ) +
  scale_color_discrete(
    breaks = family_counts$Family,
    labels = family_labels,
    drop = FALSE
  ) +
  guides(
    color = guide_legend(
      ncol = 2,
      byrow = FALSE,
      title.position = "top",
      title.hjust = 0,
      override.aes = list(size = 1.4, alpha = 1)
    )
  ) +
  theme_classic() +
  theme(
    legend.position = c(0.985, 0.035),
    legend.justification = c(1, 0),
    legend.direction = "vertical",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key = element_blank(),
    legend.key.height = unit(0.065, "in"),
    legend.key.width = unit(0.065, "in"),
    legend.spacing.x = unit(0.03, "in"),
    legend.spacing.y = unit(0.008, "in"),
    axis.title = element_text(size = 10.5),
    axis.text = element_text(size = 8.8),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(t = 3, r = 4, b = 3, l = 12, unit = "pt")
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    color = legend_title
  )

print(p)

ggsave(
  filename = png_file,
  plot = p,
  width = 6.5,
  height = 3.75,
  dpi = 400
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 6.5,
  height = 3.75,
  device = cairo_pdf
)

cat("\n==================== SUMMARY ====================\n")
cat("Datasets: ", n_datasets, "\n", sep = "")
cat("Unique spp: ", n_species_codes, "\n", sep = "")
cat("Unique sites: ", n_sites, "\n", sep = "")
cat("PNG saved to: ", png_file, "\n", sep = "")
cat("PDF saved to: ", pdf_file, "\n", sep = "")
cat("=================================================\n\n")
