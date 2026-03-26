# -------------------------------
# ONMY-2 pairwise-eligible site composite by HUC level
# Sites in the same eligible HUC group at a focal scale share a color
# Grey transparent points are not pairwise-eligible at that scale
# Saves one composite plot only; does not print to RStudio
# -------------------------------

library(dplyr)
library(ggplot2)
library(sf)
library(terra)
library(FedData)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)

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
  
  if (inherits(x, "sf")) return(x)
  if (inherits(x, "sfc")) return(st_as_sf(x))
  if (inherits(x, "SpatVector")) return(st_as_sf(x))
  if (inherits(x, "Spatial")) return(st_as_sf(x))
  
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
# helper: classify points into pairwise-eligible HUC groups
# eligible = at least 2 sites in same focal HUC unit
# -------------------------------
prep_sites_for_plot <- function(sites_sf, huc_sf, level_label) {
  
  join_df <- st_join(
    sites_sf,
    huc_sf %>% select(huc_code),
    join = st_intersects,
    left = TRUE
  ) %>%
    st_drop_geometry() %>%
    mutate(site_id = seq_len(nrow(sites_sf)))
  
  huc_counts <- join_df %>%
    count(huc_code, name = "n_sites_in_huc")
  
  site_status <- join_df %>%
    left_join(huc_counts, by = "huc_code") %>%
    mutate(
      n_sites_in_huc = ifelse(is.na(n_sites_in_huc), 0L, n_sites_in_huc),
      eligible = n_sites_in_huc >= 2,
      plot_group = ifelse(eligible, as.character(huc_code), "unused")
    ) %>%
    select(site_id, huc_code, n_sites_in_huc, eligible, plot_group)
  
  out <- sites_sf %>%
    mutate(site_id = seq_len(nrow(sites_sf))) %>%
    left_join(site_status, by = "site_id") %>%
    mutate(
      huc_level = level_label,
      eligible = ifelse(is.na(eligible), FALSE, eligible),
      plot_group = ifelse(is.na(plot_group), "unused", plot_group),
      n_sites_in_huc = ifelse(is.na(n_sites_in_huc), 0L, n_sites_in_huc)
    )
  
  out
}

# -------------------------------
# helper: qualitative palette
# -------------------------------
get_qual_palette <- function(n) {
  base_pal <- c(
    "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
    "#e6ab02", "#a6761d", "#1f78b4", "#b2df8a", "#fb9a99",
    "#cab2d6", "#fdbf6f", "#6a3d9a", "#33a02c", "#ff7f00",
    "#b15928", "#8dd3c7", "#80b1d3", "#bebada", "#fb8072"
  )
  
  if (n <= length(base_pal)) {
    base_pal[seq_len(n)]
  } else {
    grDevices::colorRampPalette(base_pal)(n)
  }
}

# -------------------------------
# helper: make one panel
# -------------------------------
make_unique_pair_plot <- function(huc_sf, sites_plot_sf, states_sf, buffer_sf, xlim, ylim, level_label) {
  
  used_sf <- sites_plot_sf %>% filter(eligible)
  unused_sf <- sites_plot_sf %>% filter(!eligible)
  
  used_groups <- sort(unique(used_sf$plot_group))
  
  if (length(used_groups) > 0) {
    pal <- get_qual_palette(length(used_groups))
    names(pal) <- used_groups
  } else {
    pal <- NULL
  }
  
  p <- ggplot() +
    geom_sf(
      data = states_sf,
      fill = "grey96",
      color = "grey75",
      linewidth = 0.18
    ) +
    geom_sf(
      data = huc_sf,
      fill = NA,
      color = "grey50",
      linewidth = 0.20
    ) +
    geom_sf(
      data = buffer_sf,
      fill = NA,
      color = "black",
      linewidth = 0.45,
      linetype = "dashed"
    ) +
    geom_sf(
      data = unused_sf,
      color = "grey60",
      size = 2.2,
      alpha = 0.22
    )
  
  if (nrow(used_sf) > 0) {
    p <- p +
      geom_sf(
        data = used_sf,
        aes(color = plot_group),
        size = 2.35,
        alpha = 0.95
      ) +
      scale_color_manual(values = pal, guide = "none")
  }
  
  p +
    coord_sf(
      xlim = xlim,
      ylim = ylim,
      expand = FALSE
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      axis.title = element_blank(),
      axis.text = element_text(size = 8),
      axis.ticks = element_line(linewidth = 0.25),
      plot.margin = margin(3, 3, 3, 3)
    ) +
    labs(title = level_label)
}

# -------------------------------
# extract coordinates
# -------------------------------
sites_df <- extract_site_coords(study_code, base_dir)

if (nrow(sites_df) == 0) {
  stop("No site coordinates were found for ", study_code)
}

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

if (is.na(st_crs(wbd_sf))) {
  st_crs(wbd_sf) <- st_crs(onmy_poly_buffer)
} else if (st_crs(wbd_sf) != st_crs(onmy_poly_buffer)) {
  wbd_sf <- st_transform(wbd_sf, st_crs(onmy_poly_buffer))
}

wbd_sf <- wbd_sf[st_intersects(wbd_sf, onmy_poly_buffer, sparse = FALSE)[, 1], ]

# -------------------------------
# identify HUC12 code column and dissolve upward
# -------------------------------
huc12_col <- find_huc12_col(wbd_sf)

huc12_sf <- wbd_sf %>%
  mutate(huc_code = as.character(.data[[huc12_col]])) %>%
  select(huc_code) %>%
  st_make_valid()

huc10_sf <- dissolve_huc(wbd_sf, huc12_col, 10)
huc8_sf  <- dissolve_huc(wbd_sf, huc12_col, 8)
huc6_sf  <- dissolve_huc(wbd_sf, huc12_col, 6)
huc4_sf  <- dissolve_huc(wbd_sf, huc12_col, 4)
huc2_sf  <- dissolve_huc(wbd_sf, huc12_col, 2)

# -------------------------------
# prep sites for plotting by HUC scale
# -------------------------------
sites_huc12_plot <- prep_sites_for_plot(sites_sf, huc12_sf, "HUC12")
sites_huc10_plot <- prep_sites_for_plot(sites_sf, huc10_sf, "HUC10")
sites_huc8_plot  <- prep_sites_for_plot(sites_sf, huc8_sf,  "HUC8")
sites_huc6_plot  <- prep_sites_for_plot(sites_sf, huc6_sf,  "HUC6")
sites_huc4_plot  <- prep_sites_for_plot(sites_sf, huc4_sf,  "HUC4")
sites_huc2_plot  <- prep_sites_for_plot(sites_sf, huc2_sf,  "HUC2")

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
# make panels
# -------------------------------
p12 <- make_unique_pair_plot(huc12_sf, sites_huc12_plot, states_sf, onmy_poly_buffer, xlim, ylim, "HUC12")
p10 <- make_unique_pair_plot(huc10_sf, sites_huc10_plot, states_sf, onmy_poly_buffer, xlim, ylim, "HUC10")
p8  <- make_unique_pair_plot(huc8_sf,  sites_huc8_plot,  states_sf, onmy_poly_buffer, xlim, ylim, "HUC8")
p6  <- make_unique_pair_plot(huc6_sf,  sites_huc6_plot,  states_sf, onmy_poly_buffer, xlim, ylim, "HUC6")
p4  <- make_unique_pair_plot(huc4_sf,  sites_huc4_plot,  states_sf, onmy_poly_buffer, xlim, ylim, "HUC4")
p2  <- make_unique_pair_plot(huc2_sf,  sites_huc2_plot,  states_sf, onmy_poly_buffer, xlim, ylim, "HUC2")

# -------------------------------
# composite plot
# -------------------------------
p_combo <- (p12 | p10 | p8) / (p6 | p4 | p2) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2),
    panel.spacing = unit(0.25, "lines")
  )

# -------------------------------
# save only
# -------------------------------
ggsave(
  filename = file.path(out_dir, "ONMY_2_pairwise_sites_composite_unique_colors.png"),
  plot = p_combo,
  width = 14,
  height = 7.5,
  dpi = 400,
  bg = "white"
)



# -------------------------------
# ONMY-2 pilot run:
# multiscale HUC grouping maps + simple pairwise regressions
# -------------------------------

library(dplyr)
library(ggplot2)
library(sf)
library(terra)
library(FedData)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(scales)

# -------------------------------
# paths
# -------------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"

onmy_rdata <- file.path(base_dir, "ONMY-2", "data", "ONMY-2.RData")

out_dir <- path.expand("~/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/analyses/ONMY_2_pilot_run")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# load ONMY-2 data reproducibly
# -------------------------------
load(onmy_rdata)

gen_dist_mat <- as.matrix(ONMY_2_fst)
gen_dist_mat[gen_dist_mat < 0] <- 0

sites_df <- ONMY_2_coords %>%
  dplyr::select(site_id, lat, lon) %>%
  mutate(
    site_id = as.integer(site_id),
    lat = as.numeric(lat),
    lon = as.numeric(lon)
  ) %>%
  arrange(site_id)

if (nrow(sites_df) != nrow(gen_dist_mat)) {
  stop("Number of rows in ONMY_2_coords does not match dimensions of ONMY_2_fst.")
}

if (nrow(gen_dist_mat) != ncol(gen_dist_mat)) {
  stop("ONMY_2_fst is not square.")
}

# enforce site IDs as row/col names if missing
if (is.null(rownames(gen_dist_mat))) rownames(gen_dist_mat) <- sites_df$site_id
if (is.null(colnames(gen_dist_mat))) colnames(gen_dist_mat) <- sites_df$site_id

# reorder if possible to match site_id
if (all(as.character(sites_df$site_id) %in% rownames(gen_dist_mat))) {
  gen_dist_mat <- gen_dist_mat[
    as.character(sites_df$site_id),
    as.character(sites_df$site_id)
  ]
}

# Euclidean geographic distance matrix for pilot run
geo_dist_mat <- as.matrix(
  dist(as.matrix(sites_df[, c("lon", "lat")]))
)

rownames(geo_dist_mat) <- sites_df$site_id
colnames(geo_dist_mat) <- sites_df$site_id

# -------------------------------
# helpers
# -------------------------------
unwrap_to_sf <- function(x) {
  if (inherits(x, "sf")) return(x)
  if (inherits(x, "sfc")) return(st_as_sf(x))
  if (inherits(x, "SpatVector")) return(st_as_sf(x))
  if (inherits(x, "Spatial")) return(st_as_sf(x))
  if (is.list(x)) {
    for (i in seq_along(x)) {
      out <- try(unwrap_to_sf(x[[i]]), silent = TRUE)
      if (!inherits(out, "try-error")) return(out)
    }
  }
  stop("Could not extract an sf object from get_wbd() output.")
}

find_huc12_col <- function(x) {
  nms <- names(x)
  nms_lower <- tolower(nms)
  
  exact_idx <- which(nms_lower %in% c("huc12", "huc_12", "hu_12"))
  if (length(exact_idx) > 0) return(nms[exact_idx[1]])
  
  contains_idx <- grep("huc.?12|hu.?12", nms_lower)
  if (length(contains_idx) > 0) return(nms[contains_idx[1]])
  
  stop("Could not find a HUC12 column in the WBD object.")
}

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

prep_sites_for_plot <- function(sites_sf, huc_sf, level_label) {
  join_df <- st_join(
    sites_sf,
    huc_sf %>% select(huc_code),
    join = st_intersects,
    left = TRUE
  ) %>%
    st_drop_geometry() %>%
    mutate(site_row = seq_len(nrow(sites_sf)))
  
  huc_counts <- join_df %>%
    count(huc_code, name = "n_sites_in_huc")
  
  site_status <- join_df %>%
    left_join(huc_counts, by = "huc_code") %>%
    mutate(
      n_sites_in_huc = ifelse(is.na(n_sites_in_huc), 0L, n_sites_in_huc),
      eligible = n_sites_in_huc >= 2,
      plot_group = ifelse(eligible, as.character(huc_code), "unused")
    ) %>%
    select(site_row, huc_code, n_sites_in_huc, eligible, plot_group)
  
  out <- sites_sf %>%
    mutate(site_row = seq_len(nrow(.))) %>%
    left_join(site_status, by = "site_row") %>%
    mutate(
      huc_level = level_label,
      eligible = ifelse(is.na(eligible), FALSE, eligible),
      plot_group = ifelse(is.na(plot_group), "unused", plot_group),
      n_sites_in_huc = ifelse(is.na(n_sites_in_huc), 0L, n_sites_in_huc)
    )
  
  out
}

get_qual_palette <- function(n) {
  base_pal <- c(
    "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
    "#e6ab02", "#a6761d", "#1f78b4", "#b2df8a", "#fb9a99",
    "#cab2d6", "#fdbf6f", "#6a3d9a", "#33a02c", "#ff7f00",
    "#b15928", "#8dd3c7", "#80b1d3", "#bebada", "#fb8072"
  )
  
  if (n <= length(base_pal)) {
    base_pal[seq_len(n)]
  } else {
    grDevices::colorRampPalette(base_pal)(n)
  }
}

get_pairwise_df <- function(sites_plot_sf, gen_mat, geo_mat) {
  sites_plot_sf <- sites_plot_sf %>%
    st_drop_geometry() %>%
    mutate(site_row = seq_len(nrow(.)))
  
  df_list <- list()
  
  groups <- unique(sites_plot_sf$plot_group)
  groups <- groups[groups != "unused"]
  
  for (g in groups) {
    idx <- which(sites_plot_sf$plot_group == g)
    
    if (length(idx) < 2) next
    
    combs <- t(combn(idx, 2))
    
    tmp <- data.frame(
      i = combs[, 1],
      j = combs[, 2],
      group = g
    )
    
    tmp$gen <- mapply(function(a, b) gen_mat[a, b], tmp$i, tmp$j)
    tmp$geo <- mapply(function(a, b) geo_mat[a, b], tmp$i, tmp$j)
    
    df_list[[g]] <- tmp
  }
  
  bind_rows(df_list)
}

make_unique_pair_map <- function(huc_sf, sites_plot_sf, states_sf, buffer_sf, xlim, ylim, level_label) {
  used_sf <- sites_plot_sf %>% filter(eligible)
  unused_sf <- sites_plot_sf %>% filter(!eligible)
  
  used_groups <- sort(unique(used_sf$plot_group))
  
  if (length(used_groups) > 0) {
    pal <- get_qual_palette(length(used_groups))
    names(pal) <- used_groups
  } else {
    pal <- NULL
  }
  
  p <- ggplot() +
    geom_sf(
      data = states_sf,
      fill = "grey96",
      color = "grey75",
      linewidth = 0.18
    ) +
    geom_sf(
      data = huc_sf,
      fill = NA,
      color = "grey50",
      linewidth = 0.20
    ) +
    geom_sf(
      data = buffer_sf,
      fill = NA,
      color = "black",
      linewidth = 0.45,
      linetype = "dashed"
    ) +
    geom_sf(
      data = unused_sf,
      shape = 21,
      fill = alpha("grey70", 0.35),
      color = "black",
      size = 2.3,
      stroke = 0.3
    )
  
  if (nrow(used_sf) > 0) {
    p <- p +
      geom_sf(
        data = used_sf,
        aes(color = plot_group),
        size = 2.35,
        alpha = 0.95
      ) +
      scale_color_manual(values = pal, guide = "none")
  }
  
  p +
    coord_sf(
      xlim = xlim,
      ylim = ylim,
      expand = FALSE
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 11, hjust = 0.5),
      axis.title = element_blank(),
      axis.text = element_text(size = 8),
      axis.ticks = element_line(linewidth = 0.25),
      plot.margin = margin(3, 3, 3, 3)
    ) +
    labs(title = level_label)
}

make_scatter_plot <- function(pair_df, level_label, palette_named) {
  if (nrow(pair_df) == 0) {
    return(
      ggplot() +
        theme_void() +
        labs(title = level_label)
    )
  }
  
  fit <- lm(gen ~ geo, data = pair_df)
  
  r_val <- suppressWarnings(cor(pair_df$gen, pair_df$geo, use = "complete.obs"))
  p_val <- summary(fit)$coefficients[2, 4]
  
  label_txt <- paste0(
    "r = ", round(r_val, 3),
    "\n",
    "p = ", signif(p_val, 3)
  )
  
  ggplot(pair_df, aes(x = geo, y = gen, color = group)) +
    geom_point(size = 2.0, alpha = 0.9) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.55) +
    scale_color_manual(values = palette_named, guide = "none") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 11, hjust = 0.5),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      plot.margin = margin(3, 3, 3, 3)
    ) +
    labs(
      title = level_label,
      x = "Geographic distance",
      y = "Genetic distance"
    ) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = label_txt,
      hjust = 1.05, vjust = 1.3,
      size = 3.2
    )
}

make_panel <- function(huc_sf, sites_plot_sf, level_label,
                       states_sf, buffer_sf, xlim, ylim,
                       gen_mat, geo_mat) {
  used_groups <- sort(unique((sites_plot_sf %>% filter(eligible))$plot_group))
  
  if (length(used_groups) > 0) {
    pal <- get_qual_palette(length(used_groups))
    names(pal) <- used_groups
  } else {
    pal <- c()
  }
  
  map_plot <- make_unique_pair_map(
    huc_sf = huc_sf,
    sites_plot_sf = sites_plot_sf,
    states_sf = states_sf,
    buffer_sf = buffer_sf,
    xlim = xlim,
    ylim = ylim,
    level_label = level_label
  )
  
  pair_df <- get_pairwise_df(sites_plot_sf, gen_mat, geo_mat)
  
  scatter_plot <- make_scatter_plot(
    pair_df = pair_df,
    level_label = level_label,
    palette_named = pal
  )
  
  map_plot | scatter_plot
}

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

if (is.na(st_crs(wbd_sf))) {
  st_crs(wbd_sf) <- st_crs(onmy_poly_buffer)
} else if (st_crs(wbd_sf) != st_crs(onmy_poly_buffer)) {
  wbd_sf <- st_transform(wbd_sf, st_crs(onmy_poly_buffer))
}

wbd_sf <- wbd_sf[st_intersects(wbd_sf, onmy_poly_buffer, sparse = FALSE)[, 1], ]

# -------------------------------
# identify HUC12 code column and dissolve upward
# -------------------------------
huc12_col <- find_huc12_col(wbd_sf)

huc12_sf <- wbd_sf %>%
  mutate(huc_code = as.character(.data[[huc12_col]])) %>%
  select(huc_code) %>%
  st_make_valid()

huc10_sf <- dissolve_huc(wbd_sf, huc12_col, 10)
huc8_sf  <- dissolve_huc(wbd_sf, huc12_col, 8)
huc6_sf  <- dissolve_huc(wbd_sf, huc12_col, 6)
huc4_sf  <- dissolve_huc(wbd_sf, huc12_col, 4)
huc2_sf  <- dissolve_huc(wbd_sf, huc12_col, 2)

# -------------------------------
# classify sites for plotting by scale
# -------------------------------
sites_huc12_plot <- prep_sites_for_plot(sites_sf, huc12_sf, "HUC12")
sites_huc10_plot <- prep_sites_for_plot(sites_sf, huc10_sf, "HUC10")
sites_huc8_plot  <- prep_sites_for_plot(sites_sf, huc8_sf,  "HUC8")
sites_huc6_plot  <- prep_sites_for_plot(sites_sf, huc6_sf,  "HUC6")
sites_huc4_plot  <- prep_sites_for_plot(sites_sf, huc4_sf,  "HUC4")
sites_huc2_plot  <- prep_sites_for_plot(sites_sf, huc2_sf,  "HUC2")

# -------------------------------
# states background
# -------------------------------
states_sf <- ne_states(
  country = "United States of America",
  returnclass = "sf"
) %>%
  filter(!name_en %in% c("Alaska", "Hawaii", "Puerto Rico"))

if (st_crs(states_sf) != st_crs(onmy_poly_buffer)) {
  states_sf <- st_transform(states_sf, st_crs(onmy_poly_buffer))
}

# -------------------------------
# build paired panels
# -------------------------------
panel12 <- make_panel(huc12_sf, sites_huc12_plot, "HUC12",
                      states_sf, onmy_poly_buffer, xlim, ylim,
                      gen_dist_mat, geo_dist_mat)

panel10 <- make_panel(huc10_sf, sites_huc10_plot, "HUC10",
                      states_sf, onmy_poly_buffer, xlim, ylim,
                      gen_dist_mat, geo_dist_mat)

panel8 <- make_panel(huc8_sf, sites_huc8_plot, "HUC8",
                     states_sf, onmy_poly_buffer, xlim, ylim,
                     gen_dist_mat, geo_dist_mat)

panel6 <- make_panel(huc6_sf, sites_huc6_plot, "HUC6",
                     states_sf, onmy_poly_buffer, xlim, ylim,
                     gen_dist_mat, geo_dist_mat)

panel4 <- make_panel(huc4_sf, sites_huc4_plot, "HUC4",
                     states_sf, onmy_poly_buffer, xlim, ylim,
                     gen_dist_mat, geo_dist_mat)

panel2 <- make_panel(huc2_sf, sites_huc2_plot, "HUC2",
                     states_sf, onmy_poly_buffer, xlim, ylim,
                     gen_dist_mat, geo_dist_mat)

# -------------------------------
# final composite
# -------------------------------
p_final <- (panel12 / panel10 / panel8) | (panel6 / panel4 / panel2)

# -------------------------------
# save only
# -------------------------------
ggsave(
  filename = file.path(out_dir, "ONMY_2_multiscale_map_scatter_composite.png"),
  plot = p_final,
  width = 20,
  height = 14,
  dpi = 400,
  bg = "white"
)