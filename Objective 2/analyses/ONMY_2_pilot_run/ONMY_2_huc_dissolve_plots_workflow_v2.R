# -------------------------------
# ONMY-2 pilot run:
# 1) load ONMY-2 coordinates
# 2) build polygon around points
# 3) buffer polygon by 10 km
# 4) download WBD HUC12 for buffered extent
# 5) dissolve HUC12s up to HUC10, HUC8, HUC6, HUC4, and HUC2
# 6) plot each HUC level behind points
# note: this version saves plots only
# -------------------------------

library(dplyr)
library(ggplot2)
library(sf)
library(terra)
library(FedData)
library(rnaturalearth)
library(rnaturalearthdata)

# -------------------------------
# paths
# -------------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"

out_dir <- path.expand("~/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/analyses/ONMY_2_pilot_run")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

study_code <- "ONMY-2"

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
    filter(lon >= -170, lon <= -50, lat >= 20, lat <= 85) %>%
    distinct(lon, lat)

  if (nrow(out) == 0) return(NULL)

  out
}

# -------------------------------
# helper: load .RData from folder/data/
# and extract best lon/lat dataframe
# preference to objects with "coord" in the name
# -------------------------------
extract_site_coords <- function(study_code, base_dir) {

  folder_path <- file.path(base_dir, study_code)
  data_dir <- file.path(folder_path, "data")

  if (!dir.exists(data_dir)) {
    stop("No data/ folder found in: ", study_code)
  }

  rdata_files <- list.files(data_dir, pattern = "\\.RData$", full.names = TRUE)

  if (length(rdata_files) == 0) {
    stop("No .RData file found in: ", study_code)
  }

  tmp_env <- new.env()
  load(rdata_files[1], envir = tmp_env)
  obj_names <- ls(tmp_env)

  obj_names <- c(
    obj_names[grepl("coord", tolower(obj_names))],
    obj_names[!grepl("coord", tolower(obj_names))]
  ) |> unique()

  for (nm in obj_names) {
    obj <- get(nm, envir = tmp_env)

    if (is.data.frame(obj)) {
      coords <- find_lon_lat_cols(obj)
      if (!is.null(coords)) {
        coords$study_code <- study_code
        coords$source_object <- nm
        return(coords)
      }
    }

    if (is.list(obj) && !is.data.frame(obj)) {
      sub_nms <- names(obj)

      if (!is.null(sub_nms) && length(sub_nms) > 0) {
        sub_nms <- c(
          sub_nms[grepl("coord", tolower(sub_nms))],
          sub_nms[!grepl("coord", tolower(sub_nms))]
        ) |> unique()

        for (sub_nm in sub_nms) {
          sub_obj <- obj[[sub_nm]]

          if (is.data.frame(sub_obj)) {
            coords <- find_lon_lat_cols(sub_obj)
            if (!is.null(coords)) {
              coords$study_code <- study_code
              coords$source_object <- paste0(nm, "$", sub_nm)
              return(coords)
            }
          }
        }
      }
    }
  }

  stop("No suitable lon/lat dataframe found in: ", study_code)
}

# -------------------------------
# helper: unwrap get_wbd output to sf
# -------------------------------
unwrap_to_sf <- function(x) {

  if (inherits(x, "sf")) {
    return(x)
  }

  if (inherits(x, "sfc")) {
    return(st_as_sf(x))
  }

  if (inherits(x, "SpatVector")) {
    return(st_as_sf(x))
  }

  if (inherits(x, "Spatial")) {
    return(st_as_sf(x))
  }

  if (is.list(x)) {
    for (i in seq_along(x)) {
      out <- try(unwrap_to_sf(x[[i]]), silent = TRUE)
      if (!inherits(out, "try-error")) return(out)
    }
  }

  stop("Could not extract an sf object from get_wbd() output.")
}

# -------------------------------
# helper: find HUC12 code column
# -------------------------------
find_huc12_col <- function(x) {
  nms <- names(x)
  nms_lower <- tolower(nms)

  exact_idx <- which(nms_lower %in% c("huc12", "huc_12", "hu_12"))
  if (length(exact_idx) > 0) return(nms[exact_idx[1]])

  contains_idx <- grep("huc.?12|hu.?12", nms_lower)
  if (length(contains_idx) > 0) return(nms[contains_idx[1]])

  stop("Could not find a HUC12 column in the WBD object.")
}

# -------------------------------
# helper: dissolve HUCs by code length
# -------------------------------
dissolve_huc <- function(wbd_sf, huc12_col, level) {
  code_n <- as.integer(level)

  out <- wbd_sf %>%
    mutate(
      huc12_code = as.character(.data[[huc12_col]]),
      huc_code = substr(huc12_code, 1, code_n)
    ) %>%
    select(huc_code) %>%
    group_by(huc_code) %>%
    summarise(do_union = TRUE, .groups = "drop") %>%
    st_make_valid()

  out$huc_level <- paste0("HUC", level)
  out
}

# -------------------------------
# helper: make map
# -------------------------------
make_huc_plot <- function(huc_sf, sites_sf, states_sf, buffer_sf, xlim, ylim, level) {
  ggplot() +
    geom_sf(
      data = states_sf,
      fill = "grey96",
      color = "grey70",
      linewidth = 0.2
    ) +
    geom_sf(
      data = huc_sf,
      fill = NA,
      color = "grey40",
      linewidth = 0.28
    ) +
    geom_sf(
      data = buffer_sf,
      fill = NA,
      color = "black",
      linewidth = 0.55,
      linetype = "dashed"
    ) +
    geom_sf(
      data = sites_sf,
      color = "firebrick3",
      size = 2.8,
      alpha = 0.95
    ) +
    coord_sf(
      xlim = xlim,
      ylim = ylim,
      expand = FALSE
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14)
    ) +
    labs(
      x = "Longitude",
      y = "Latitude",
      title = paste0("ONMY-2: ", level, " boundaries within 10 km buffered study extent")
    )
}

# -------------------------------
# extract coordinates
# -------------------------------
sites_df <- extract_site_coords(study_code, base_dir)

if (nrow(sites_df) == 0) {
  stop("No site coordinates were found for ", study_code)
}

print(sites_df)

# -------------------------------
# points to sf
# -------------------------------
sites_sf <- st_as_sf(
  sites_df,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

# -------------------------------
# build buffered polygon around sites
# convex hull in projected CRS so 10 km buffer is in meters
# -------------------------------
sites_aea <- st_transform(sites_sf, 5070)

onmy_poly_aea <- sites_aea %>%
  st_union() %>%
  st_convex_hull() %>%
  st_as_sf() %>%
  st_set_crs(5070)

onmy_poly_buffer_aea <- st_buffer(onmy_poly_aea, dist = 10000)

onmy_poly <- st_transform(onmy_poly_aea, 4326)
onmy_poly_buffer <- st_transform(onmy_poly_buffer_aea, 4326)

# -------------------------------
# dynamic plot extent from buffered polygon bbox
# -------------------------------
bb <- st_bbox(onmy_poly_buffer)
xlim <- c(bb["xmin"], bb["xmax"])
ylim <- c(bb["ymin"], bb["ymax"])

# -------------------------------
# download WBD HUC12 using buffered polygon as template
# -------------------------------
wbd_obj <- get_wbd(
  template = onmy_poly_buffer,
  label = "ONMY_2_pilot_run",
  extraction.dir = file.path(out_dir, "FedData_wbd"),
  force.redo = FALSE
)

wbd_sf <- unwrap_to_sf(wbd_obj)
wbd_sf <- st_make_valid(wbd_sf)

# assign/match CRS before intersection
if (is.na(st_crs(wbd_sf))) {
  st_crs(wbd_sf) <- st_crs(onmy_poly_buffer)
} else if (st_crs(wbd_sf) != st_crs(onmy_poly_buffer)) {
  wbd_sf <- st_transform(wbd_sf, st_crs(onmy_poly_buffer))
}

# keep only polygons intersecting the buffered polygon
wbd_sf <- wbd_sf[st_intersects(wbd_sf, onmy_poly_buffer, sparse = FALSE)[, 1], ]

# -------------------------------
# identify HUC12 code column and dissolve upward
# -------------------------------
huc12_col <- find_huc12_col(wbd_sf)

huc12_sf <- wbd_sf %>%
  mutate(
    huc_code = as.character(.data[[huc12_col]]),
    huc_level = "HUC12"
  ) %>%
  select(huc_code, huc_level) %>%
  st_make_valid()

huc10_sf <- dissolve_huc(wbd_sf, huc12_col, 10)
huc8_sf  <- dissolve_huc(wbd_sf, huc12_col, 8)
huc6_sf  <- dissolve_huc(wbd_sf, huc12_col, 6)
huc4_sf  <- dissolve_huc(wbd_sf, huc12_col, 4)
huc2_sf  <- dissolve_huc(wbd_sf, huc12_col, 2)

cat("\nHUC polygon counts within buffered extent:\n")
cat("HUC12:", nrow(huc12_sf), "\n")
cat("HUC10:", nrow(huc10_sf), "\n")
cat("HUC8: ", nrow(huc8_sf),  "\n")
cat("HUC6: ", nrow(huc6_sf),  "\n")
cat("HUC4: ", nrow(huc4_sf),  "\n")
cat("HUC2: ", nrow(huc2_sf),  "\n\n")

# -------------------------------
# get U.S. states for background
# -------------------------------
states_sf <- ne_states(
  country = "United States of America",
  returnclass = "sf"
)

states_sf <- states_sf %>%
  filter(!name_en %in% c("Alaska", "Hawaii", "Puerto Rico"))

if (st_crs(states_sf) != st_crs(onmy_poly_buffer)) {
  states_sf <- st_transform(states_sf, st_crs(onmy_poly_buffer))
}

# -------------------------------
# make plots
# -------------------------------
p_huc12 <- make_huc_plot(huc12_sf, sites_sf, states_sf, onmy_poly_buffer, xlim, ylim, "HUC12")
p_huc10 <- make_huc_plot(huc10_sf, sites_sf, states_sf, onmy_poly_buffer, xlim, ylim, "HUC10")
p_huc8  <- make_huc_plot(huc8_sf,  sites_sf, states_sf, onmy_poly_buffer, xlim, ylim, "HUC8")
p_huc6  <- make_huc_plot(huc6_sf,  sites_sf, states_sf, onmy_poly_buffer, xlim, ylim, "HUC6")
p_huc4  <- make_huc_plot(huc4_sf,  sites_sf, states_sf, onmy_poly_buffer, xlim, ylim, "HUC4")
p_huc2  <- make_huc_plot(huc2_sf,  sites_sf, states_sf, onmy_poly_buffer, xlim, ylim, "HUC2")

print(p_huc12)
print(p_huc10)
print(p_huc8)
print(p_huc6)
print(p_huc4)
print(p_huc2)

# -------------------------------
# save plots only
# -------------------------------
ggsave(
  filename = file.path(out_dir, "ONMY_2_wbd_huc12_map.png"),
  plot = p_huc12,
  width = 11,
  height = 7,
  dpi = 400
)

ggsave(
  filename = file.path(out_dir, "ONMY_2_wbd_huc10_map.png"),
  plot = p_huc10,
  width = 11,
  height = 7,
  dpi = 400
)

ggsave(
  filename = file.path(out_dir, "ONMY_2_wbd_huc8_map.png"),
  plot = p_huc8,
  width = 11,
  height = 7,
  dpi = 400
)

ggsave(
  filename = file.path(out_dir, "ONMY_2_wbd_huc6_map.png"),
  plot = p_huc6,
  width = 11,
  height = 7,
  dpi = 400
)

ggsave(
  filename = file.path(out_dir, "ONMY_2_wbd_huc4_map.png"),
  plot = p_huc4,
  width = 11,
  height = 7,
  dpi = 400
)

ggsave(
  filename = file.path(out_dir, "ONMY_2_wbd_huc2_map.png"),
  plot = p_huc2,
  width = 11,
  height = 7,
  dpi = 400
)
