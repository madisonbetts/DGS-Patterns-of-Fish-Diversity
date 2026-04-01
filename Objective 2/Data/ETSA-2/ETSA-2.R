# ========================================
# ETSA-2
# Etheostoma sagitta
# Kentucky Arrow Darter
#
# Builds:
#   ETSA_2_coords
#   ETSA_2_fst
#   ETSA_2_rivdists
#
# Saves to:
#   ETSA-2/data/ETSA-2.RData
#
# Plots in RStudio:
#   - site map
#   - IBD using in-river distance only
# ========================================

library(readxl)
library(dplyr)
library(ggplot2)
library(geosphere)
library(maps)

study_code <- "ETSA-2"

base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETSA-2"
culley_dir <- file.path(base_dir, "Culley_etal_2025_KAD_LE")
data_dir <- file.path(base_dir, "data")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

excel_file <- file.path(culley_dir, "Culley_etal_2025_KAD_LE_Sites_Distance_Fst.xlsx")
stopifnot(file.exists(excel_file))

# -----------------------------
# 1) read sheets
# -----------------------------
sheet_names <- excel_sheets(excel_file)
stopifnot(all(c("Site Information", "Distance Matrix (meters)", "Pairwise FST") %in% sheet_names))

site_info <- read_excel(excel_file, sheet = "Site Information")
dist_raw  <- read_excel(excel_file, sheet = "Distance Matrix (meters)")
fst_raw   <- read_excel(excel_file, sheet = "Pairwise FST")

# -----------------------------
# 2) clean site information
# -----------------------------
names(site_info) <- trimws(names(site_info))

site_name_col <- names(site_info)[tolower(names(site_info)) %in% c("site name", "site_name")]
id_col        <- names(site_info)[tolower(names(site_info)) %in% c("id", "site", "site id", "site_id")]
lat_col       <- names(site_info)[grepl("^lat$|latitude", tolower(names(site_info)))]
lon_col       <- names(site_info)[grepl("^long$|^lon$|longitude", tolower(names(site_info)))]

stopifnot(length(site_name_col) >= 1, length(id_col) >= 1, length(lat_col) >= 1, length(lon_col) >= 1)

site_name_col <- site_name_col[1]
id_col        <- id_col[1]
lat_col       <- lat_col[1]
lon_col       <- lon_col[1]

site_info <- site_info %>%
  transmute(
    site_name = as.character(.data[[site_name_col]]),
    site_orig = as.integer(.data[[id_col]]),
    lat = as.numeric(.data[[lat_col]]),
    lon = as.numeric(.data[[lon_col]])
  ) %>%
  filter(!is.na(site_orig)) %>%
  arrange(site_orig)

# -----------------------------
# 3) helper for square matrix sheets
# -----------------------------
clean_square_sheet <- function(x) {
  xdf <- as.data.frame(x, stringsAsFactors = FALSE)
  names(xdf)[1] <- "row_id"

  row_ids <- suppressWarnings(as.integer(xdf$row_id))
  col_ids <- suppressWarnings(as.integer(names(xdf)[-1]))

  keep_cols <- !is.na(col_ids)
  col_ids <- col_ids[keep_cols]

  vals <- as.matrix(xdf[, c(TRUE, keep_cols), drop = FALSE][, -1, drop = FALSE])
  mode(vals) <- "numeric"

  vals <- vals[!is.na(row_ids), , drop = FALSE]
  row_ids <- row_ids[!is.na(row_ids)]

  rownames(vals) <- row_ids
  colnames(vals) <- col_ids

  stopifnot(nrow(vals) == ncol(vals))
  vals
}

# -----------------------------
# 4) clean FST matrix
# -----------------------------
fst_mat <- clean_square_sheet(fst_raw)
stopifnot(setequal(as.integer(rownames(fst_mat)), site_info$site_orig))
stopifnot(setequal(as.integer(colnames(fst_mat)), site_info$site_orig))

fst_mat <- fst_mat[as.character(site_info$site_orig), as.character(site_info$site_orig), drop = FALSE]
fst_mat[lower.tri(fst_mat)] <- t(fst_mat)[lower.tri(fst_mat)]
fst_mat[fst_mat < 0] <- 0
diag(fst_mat) <- 0

# -----------------------------
# 5) clean in-river distance matrix
# -----------------------------
rivdist_mat <- clean_square_sheet(dist_raw)
stopifnot(setequal(as.integer(rownames(rivdist_mat)), site_info$site_orig))
stopifnot(setequal(as.integer(colnames(rivdist_mat)), site_info$site_orig))

rivdist_mat <- rivdist_mat[as.character(site_info$site_orig), as.character(site_info$site_orig), drop = FALSE]
rivdist_mat[lower.tri(rivdist_mat)] <- t(rivdist_mat)[lower.tri(rivdist_mat)]
diag(rivdist_mat) <- 0

# -----------------------------
# 6) build final objects
# -----------------------------
ETSA_2_coords <- data.frame(
  site = as.character(seq_len(nrow(site_info))),
  lat = site_info$lat,
  lon = site_info$lon,
  stringsAsFactors = FALSE
)

site_lookup <- data.frame(
  site = ETSA_2_coords$site,
  site_orig = site_info$site_orig,
  site_name = site_info$site_name,
  stringsAsFactors = FALSE
)

ETSA_2_fst <- fst_mat
ETSA_2_rivdists <- rivdist_mat

rownames(ETSA_2_fst) <- ETSA_2_coords$site
colnames(ETSA_2_fst) <- ETSA_2_coords$site
rownames(ETSA_2_rivdists) <- ETSA_2_coords$site
colnames(ETSA_2_rivdists) <- ETSA_2_coords$site

# -----------------------------
# 7) pairwise dataframe for in-river IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1 = rownames(ETSA_2_fst)[row(ETSA_2_fst)[upper.tri(ETSA_2_fst)]],
  site2 = colnames(ETSA_2_fst)[col(ETSA_2_fst)[upper.tri(ETSA_2_fst)]],
  fst = ETSA_2_fst[upper.tri(ETSA_2_fst)],
  rivdist_km = ETSA_2_rivdists[upper.tri(ETSA_2_rivdists)] / 1000,
  stringsAsFactors = FALSE
)

ibd_df$site1_name <- site_lookup$site_name[match(ibd_df$site1, site_lookup$site)]
ibd_df$site2_name <- site_lookup$site_name[match(ibd_df$site2, site_lookup$site)]

# -----------------------------
# 8) map
# -----------------------------
state_map <- map_data("state")

lon_rng <- range(ETSA_2_coords$lon, na.rm = TRUE)
lat_rng <- range(ETSA_2_coords$lat, na.rm = TRUE)

x_pad <- max(0.20, diff(lon_rng) * 0.12)
y_pad <- max(0.15, diff(lat_rng) * 0.12)

xlim_use <- c(lon_rng[1] - x_pad, lon_rng[2] + x_pad)
ylim_use <- c(lat_rng[1] - y_pad, lat_rng[2] + y_pad)

plot_sites <- merge(ETSA_2_coords, site_lookup, by = "site", sort = FALSE)

map_plot <- ggplot() +
  geom_polygon(
    data = state_map,
    aes(x = long, y = lat, group = group),
    fill = "grey96",
    color = "grey55",
    linewidth = 0.25
  ) +
  geom_point(
    data = plot_sites,
    aes(x = lon, y = lat),
    size = 3
  ) +
  geom_text(
    data = plot_sites,
    aes(x = lon, y = lat, label = paste0(site, ". ", site_name)),
    nudge_y = 0.012,
    size = 3.2
  ) +
  coord_fixed(xlim = xlim_use, ylim = ylim_use) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude", title = "ETSA-2 sampling locations")

print(map_plot)

ggsave(
  filename = file.path(base_dir, "ETSA-2_map.png"),
  plot = map_plot,
  width = 7.5,
  height = 5.75,
  dpi = 300
)

# -----------------------------
# 9) IBD plot using in-river distance only
# plotted in RStudio; not saved
# -----------------------------
ibd_plot_river <- ggplot(ibd_df, aes(x = rivdist_km, y = fst)) +
  geom_point(size = 2.7, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "In-river distance (km)",
    y = expression(F[ST]),
    title = "ETSA-2 isolation by distance (in-river)"
  )

print(ibd_plot_river)

# -----------------------------
# 10) checks + save
# -----------------------------
stopifnot(identical(rownames(ETSA_2_fst), ETSA_2_coords$site))
stopifnot(identical(colnames(ETSA_2_fst), ETSA_2_coords$site))
stopifnot(identical(rownames(ETSA_2_rivdists), ETSA_2_coords$site))
stopifnot(identical(colnames(ETSA_2_rivdists), ETSA_2_coords$site))
stopifnot(isTRUE(all.equal(ETSA_2_fst, t(ETSA_2_fst))))
stopifnot(isTRUE(all.equal(ETSA_2_rivdists, t(ETSA_2_rivdists))))

ETSA_2_rivdists = (ETSA_2_rivdists / 1000)


save(
  ETSA_2_fst,
  ETSA_2_coords,
  ETSA_2_rivdists,
  file = file.path(data_dir, "ETSA-2.RData")
)
