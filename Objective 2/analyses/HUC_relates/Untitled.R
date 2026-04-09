# ============================================================
# Link HUC2-HUC12 to all Objective 2 study sites
# Exclude sites outside the conterminous US from HUC outputs
# Canonical output = long table
# Convenience output = wide table (rows = HUC level, cols = site)
# Also writes excluded-site QA tables
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(nhdplusTools)
  library(tigris)
})

options(tigris_use_cache = TRUE)

# ------------------------------------------------------------
# 1) paths
# ------------------------------------------------------------
base_data_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/data"

huc_cache_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/geodata/HUC"

dir.create(huc_cache_dir, recursive = TRUE, showWarnings = FALSE)

wbd_rds        <- file.path(huc_cache_dir, "WBD_HUC12_CONUS_cached.rds")
wbd_gpkg       <- file.path(huc_cache_dir, "WBD_HUC12_CONUS_cached.gpkg")
conus_rds      <- file.path(huc_cache_dir, "CONUS_states_union_cached.rds")
summary_csv    <- file.path(huc_cache_dir, "HUC_linking_summary.csv")

# ------------------------------------------------------------
# 2) helper functions
# ------------------------------------------------------------
get_study_code <- function(rdata_file) {
  tools::file_path_sans_ext(basename(rdata_file))
}

standardize_coords <- function(df, study_code) {
  df <- as.data.frame(df)
  
  if ("site_id" %in% names(df) && !"site" %in% names(df)) {
    df <- dplyr::rename(df, site = site_id)
  }
  
  req <- c("site", "lat", "lon")
  if (!all(req %in% names(df))) {
    stop(
      paste0(
        study_code,
        ": coords object must contain site (or site_id), lat, lon. Found columns: ",
        paste(names(df), collapse = ", ")
      )
    )
  }
  
  df %>%
    dplyr::mutate(
      site = as.character(site),
      lat  = as.numeric(lat),
      lon  = as.numeric(lon)
    ) %>%
    dplyr::filter(!is.na(site), !is.na(lat), !is.na(lon)) %>%
    dplyr::distinct(site, .keep_all = TRUE)
}

find_coords_object <- function(env, study_code) {
  obj_names <- ls(env)
  
  is_coords_df <- function(x) {
    obj <- get(x, envir = env)
    if (!is.data.frame(obj)) return(FALSE)
    nms <- names(obj)
    has_site <- any(c("site", "site_id") %in% nms)
    has_latlon <- all(c("lat", "lon") %in% nms)
    has_site && has_latlon
  }
  
  coord_named <- obj_names[stringr::str_detect(tolower(obj_names), "coord")]
  coord_named_hits <- coord_named[vapply(coord_named, is_coords_df, logical(1))]
  if (length(coord_named_hits) >= 1) return(coord_named_hits[1])
  
  generic_hits <- obj_names[vapply(obj_names, is_coords_df, logical(1))]
  if (length(generic_hits) >= 1) return(generic_hits[1])
  
  stop(paste0(study_code, ": could not identify coords object in .RData"))
}

make_huc_wide <- function(site_huc_long) {
  site_huc_long %>%
    dplyr::select(site, HUC2, HUC4, HUC6, HUC8, HUC10, HUC12) %>%
    tidyr::pivot_longer(
      cols = c(HUC2, HUC4, HUC6, HUC8, HUC10, HUC12),
      names_to = "huc_level",
      values_to = "huc_code"
    ) %>%
    dplyr::mutate(
      huc_level = factor(
        huc_level,
        levels = c("HUC2", "HUC4", "HUC6", "HUC8", "HUC10", "HUC12")
      )
    ) %>%
    dplyr::arrange(huc_level) %>%
    tidyr::pivot_wider(
      names_from = site,
      values_from = huc_code
    ) %>%
    as.data.frame()
}

write_if_possible <- function(x, path) {
  tryCatch({
    utils::write.csv(x, path, row.names = FALSE)
  }, error = function(e) {
    message("Could not write: ", path)
    message("  ", e$message)
  })
}

# ------------------------------------------------------------
# 3) build/read a CONUS polygon mask
# ------------------------------------------------------------
# This mask is used only to decide whether a site is in the
# conterminous U.S. Sites outside it (e.g., Canada) get excluded
# from HUC outputs before any nearest/edge logic is applied.

if (file.exists(conus_rds)) {
  conus_poly <- readRDS(conus_rds)
} else {
  conus_states <- tigris::states(cb = TRUE, year = 2024) %>%
    dplyr::filter(!STUSPS %in% c("AK", "HI", "PR", "VI", "GU", "MP", "AS"))
  
  conus_poly <- conus_states %>%
    sf::st_make_valid() %>%
    sf::st_union() %>%
    sf::st_as_sf()
  
  saveRDS(conus_poly, conus_rds)
}

# ------------------------------------------------------------
# 4) get/cache WBD HUC12 polygons
# ------------------------------------------------------------
# WBD is a national U.S. hydrologic unit dataset, and the
# HUC12 polygons carry HUC2-HUC12 attributes. The helper
# download_wbd() can fail, so use a guarded download and stop
# cleanly if it returns NULL. :contentReference[oaicite:2]{index=2}

if (file.exists(wbd_rds)) {
  message("Reading cached CONUS WBD from RDS...")
  wbd <- readRDS(wbd_rds)
} else {
  message("No cached WBD found. Trying nhdplusTools::download_wbd() ...")
  
  wbd_path <- tryCatch(
    nhdplusTools::download_wbd(),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  
  if (is.null(wbd_path)) {
    stop(
      paste(
        "download_wbd() failed and returned NULL.",
        "Do not continue to st_read().",
        "Either rerun later, update nhdplusTools/GDAL, or manually download",
        "national WBD from The National Map and point st_read() to that file."
      )
    )
  }
  
  message("Reading downloaded WBD...")
  wbd_raw <- sf::st_read(wbd_path, quiet = TRUE)
  
  keep_cols <- intersect(
    c("HUC2", "HUC4", "HUC6", "HUC8", "HUC10", "HUC12", "geometry"),
    names(wbd_raw)
  )
  
  if (length(setdiff(c("HUC2", "HUC4", "HUC6", "HUC8", "HUC10", "HUC12"), keep_cols)) > 0) {
    stop("Downloaded WBD is missing one or more required HUC fields.")
  }
  
  wbd_raw <- wbd_raw[, keep_cols]
  wbd_raw <- sf::st_make_valid(wbd_raw)
  
  # crop to CONUS so Canadian points cannot be snapped into U.S. HUCs
  if (sf::st_crs(conus_poly) != sf::st_crs(wbd_raw)) {
    conus_poly_use <- sf::st_transform(conus_poly, sf::st_crs(wbd_raw))
  } else {
    conus_poly_use <- conus_poly
  }
  
  wbd <- suppressWarnings(sf::st_intersection(wbd_raw, conus_poly_use))
  wbd <- sf::st_make_valid(wbd)
  
  saveRDS(wbd, wbd_rds)
  tryCatch({
    sf::st_write(wbd, wbd_gpkg, delete_dsn = TRUE, quiet = TRUE)
  }, error = function(e) {
    message("Could not write gpkg cache; continuing with RDS only.")
  })
}

# ------------------------------------------------------------
# 5) find all study .RData files
# ------------------------------------------------------------
rdata_files <- list.files(
  path = base_data_dir,
  pattern = "\\.RData$",
  recursive = TRUE,
  full.names = TRUE
)

rdata_files <- rdata_files[stringr::str_detect(rdata_files, "/data/[^/]+\\.RData$")]

if (length(rdata_files) == 0) {
  stop("No .RData files found under base_data_dir.")
}

message("Found ", length(rdata_files), " .RData files.")

# ------------------------------------------------------------
# 6) process each study
# ------------------------------------------------------------
summary_list <- vector("list", length(rdata_files))

for (i in seq_along(rdata_files)) {
  rdata_file <- rdata_files[i]
  study_code <- get_study_code(rdata_file)
  
  message("\n------------------------------")
  message("Processing ", i, " / ", length(rdata_files), ": ", study_code)
  message("------------------------------")
  
  study_env <- new.env(parent = emptyenv())
  
  status <- tryCatch({
    
    load(rdata_file, envir = study_env)
    
    coords_obj_name <- find_coords_object(study_env, study_code)
    coords_df <- get(coords_obj_name, envir = study_env)
    coords_df <- standardize_coords(coords_df, study_code)
    
    pts <- sf::st_as_sf(
      coords_df,
      coords = c("lon", "lat"),
      crs = 4326,
      remove = FALSE
    )
    
    # transform mask and WBD to point CRS as needed
    conus_use <- if (sf::st_crs(conus_poly) != sf::st_crs(pts)) {
      sf::st_transform(conus_poly, sf::st_crs(pts))
    } else {
      conus_poly
    }
    
    wbd_use <- if (sf::st_crs(wbd) != sf::st_crs(pts)) {
      sf::st_transform(wbd, sf::st_crs(pts))
    } else {
      wbd
    }
    
    # classify points inside/outside CONUS first
    in_conus_mat <- sf::st_within(pts, conus_use, sparse = FALSE)
    in_conus <- as.logical(in_conus_mat[, 1])
    
    pts_in  <- pts[in_conus, ]
    pts_out <- pts[!in_conus, ]
    
    # join only CONUS points to HUC polygons
    # no global nearest fallback for excluded points
    joined_in <- sf::st_join(
      pts_in,
      wbd_use[, c("HUC2", "HUC4", "HUC6", "HUC8", "HUC10", "HUC12")],
      join = sf::st_within,
      left = TRUE
    )
    
    # optional tiny edge-case fix ONLY for points already in CONUS
    missing_idx <- which(is.na(joined_in$HUC12))
    if (length(missing_idx) > 0) {
      nearest_ids <- sf::st_nearest_feature(joined_in[missing_idx, ], wbd_use)
      joined_in[missing_idx, c("HUC2", "HUC4", "HUC6", "HUC8", "HUC10", "HUC12")] <-
        sf::st_drop_geometry(wbd_use[nearest_ids, c("HUC2", "HUC4", "HUC6", "HUC8", "HUC10", "HUC12")])
    }
    
    # canonical long output = only CONUS sites retained
    site_huc_long <- joined_in %>%
      sf::st_drop_geometry() %>%
      dplyr::select(site, lat, lon, HUC2, HUC4, HUC6, HUC8, HUC10, HUC12)
    
    # excluded sites QA table
    excluded_sites <- pts_out %>%
      sf::st_drop_geometry() %>%
      dplyr::select(site, lat, lon) %>%
      dplyr::mutate(exclusion_reason = "Outside CONUS (e.g., Canada)")
    
    # wide output built only from retained CONUS sites
    site_huc_wide <- make_huc_wide(site_huc_long)
    
    obj_prefix <- gsub("-", "_", study_code)
    
    long_name     <- paste0(obj_prefix, "_site_huc_long")
    wide_name     <- paste0(obj_prefix, "_site_huc_wide")
    excluded_name <- paste0(obj_prefix, "_site_huc_excluded")
    
    assign(long_name, site_huc_long, envir = study_env)
    assign(wide_name, site_huc_wide, envir = study_env)
    assign(excluded_name, excluded_sites, envir = study_env)
    
    out_dir <- dirname(rdata_file)
    
    write_if_possible(
      site_huc_long,
      file.path(out_dir, paste0(study_code, "_site_huc_long.csv"))
    )
    
    write_if_possible(
      site_huc_wide,
      file.path(out_dir, paste0(study_code, "_site_huc_wide.csv"))
    )
    
    write_if_possible(
      excluded_sites,
      file.path(out_dir, paste0(study_code, "_site_huc_excluded.csv"))
    )
    
    save(list = ls(study_env), file = rdata_file, envir = study_env)
    
    list(
      study_code = study_code,
      rdata_file = rdata_file,
      coords_object = coords_obj_name,
      n_total_sites = nrow(coords_df),
      n_sites_retained_conus = nrow(site_huc_long),
      n_sites_excluded_outside_conus = nrow(excluded_sites),
      n_missing_huc12_after_join = sum(is.na(site_huc_long$HUC12)),
      status = "OK"
    )
    
  }, error = function(e) {
    message("ERROR in ", study_code, ": ", e$message)
    
    list(
      study_code = study_code,
      rdata_file = rdata_file,
      coords_object = NA_character_,
      n_total_sites = NA_integer_,
      n_sites_retained_conus = NA_integer_,
      n_sites_excluded_outside_conus = NA_integer_,
      n_missing_huc12_after_join = NA_integer_,
      status = paste("ERROR:", e$message)
    )
  })
  
  summary_list[[i]] <- status
}

# ------------------------------------------------------------
# 7) write overall summary
# ------------------------------------------------------------
summary_df <- dplyr::bind_rows(summary_list)
write_if_possible(summary_df, summary_csv)

message("\nDone.")
message("Summary written to: ", summary_csv)