# -----------------------------
# SAFO-8
# White, Hanks, and Wagner (2020)
# Brook trout in the Loyalsock Creek watershed, Pennsylvania
#
# IMPORTANT
# - This workflow uses the provided D.rds file directly for the published
#   33 x 33 pairwise FST matrix from the BGR example files.
# - The site order is assumed to follow the order of D.rds and the site
#   sequence shown in Figure 4 / associated BGR materials.
# - Coordinates below are best-available approximate sampling centroids
#   reconstructed from the labeled map in Figure 4A and watershed geography.
#   They are not exact GPS points.
# - If exact site coordinates are recovered later, update ONLY the coords
#   object; the FST matrix order should stay the same.
# -----------------------------

library(dplyr)
library(ggplot2)
library(geosphere)
library(maps)

# -----------------------------
# 0) paths
# -----------------------------
study_code <- "SAFO-8"

base_dir <- file.path(
  "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026",
  "Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data",
  study_code
)

data_dir <- file.path(base_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# directory containing the supplied BGR files
# put D.rds here, or edit this path to where your copy lives
source_dir <- base_dir
d_file <- file.path(source_dir, "D.rds")

if (!file.exists(d_file)) {
  stop("D.rds not found. Put the supplied D.rds in ", source_dir,
       " or edit `source_dir` in this script.")
}

# -----------------------------
# 1) read published FST matrix
# -----------------------------
fst_raw <- readRDS(d_file)

stopifnot(is.matrix(fst_raw) || is.data.frame(fst_raw))
fst_raw <- as.matrix(fst_raw)

if (nrow(fst_raw) != 33 || ncol(fst_raw) != 33) {
  stop("Expected a 33 x 33 FST matrix in D.rds; got ",
       nrow(fst_raw), " x ", ncol(fst_raw))
}

# -----------------------------
# 2) site names in D.rds order
# this follows the supplied BGR example ordering
# -----------------------------
site_lookup <- tibble::tribble(
  ~site, ~site_name,
  1, "FLAG",
  2, "CONK",
  3, "MILA",
  4, "BEAR",
  5, "USCO",
  6, "DSCO",
  7, "DPOL",
  8, "UPOL",
  9, "DSHA",
  10, "USHA",
  11, "DOUB",
  12, "DSEA",
  13, "USEA",
  14, "YELL",
  15, "ROCK",
  16, "LEV",
  17, "STRB",
  18, "SCAR",
  19, "UNT",
  20, "LICK",
  21, "USWE",
  22, "DSWE",
  23, "DRHO",
  24, "SWAM",
  25, "MIHI",
  26, "HUCK",
  27, "BRUN",
  28, "GRAN",
  29, "SNAK",
  30, "SSR",
  31, "RED",
  32, "DSLB",
  33, "JACO"
)

# -----------------------------
# 3) approximate coordinates
# best-available map-digitized centroids from Figure 4A
# -----------------------------
SAFO_8_coords <- tibble::tribble(
  ~site, ~site_name, ~lat,    ~lon,
  1, "FLAG", 41.559, -76.739,
  2, "CONK", 41.485, -76.915,
  3, "MILA", 41.510, -76.941,
  4, "BEAR", 41.523, -76.978,
  5, "USCO", 41.504, -77.013,
  6, "DSCO", 41.497, -77.020,
  7, "DPOL", 41.492, -77.028,
  8, "UPOL", 41.484, -77.034,
  9, "DSHA", 41.486, -77.061,
  10, "USHA", 41.475, -77.074,
  11, "DOUB", 41.495, -77.126,
  12, "DSEA", 41.487, -77.111,
  13, "USEA", 41.478, -77.105,
  14, "YELL", 41.540, -77.081,
  15, "ROCK", 41.533, -77.016,
  16, "LEV",  41.578, -77.084,
  17, "STRB", 41.567, -77.077,
  18, "SCAR", 41.514, -77.154,
  19, "UNT",  41.593, -77.192,
  20, "LICK", 41.544, -77.225,
  21, "USWE", 41.563, -77.437,
  22, "DSWE", 41.533, -77.470,
  23, "DRHO", 41.516, -77.326,
  24, "SWAM", 41.492, -77.377,
  25, "MIHI", 41.471, -77.414,
  26, "HUCK", 41.417, -77.240,
  27, "BRUN", 41.400, -77.136,
  28, "GRAN", 41.357, -77.288,
  29, "SNAK", 41.362, -77.345,
  30, "SSR",  41.343, -77.385,
  31, "RED",  41.330, -77.414,
  32, "DSLB", 41.315, -77.451,
  33, "JACO", 41.354, -77.628
)

stopifnot(nrow(SAFO_8_coords) == nrow(fst_raw))

# -----------------------------
# 4) final FST object
# use numeric row/col names to match workflow
# -----------------------------
SAFO_8_fst <- fst_raw
SAFO_8_fst[SAFO_8_fst < 0] <- 0
diag(SAFO_8_fst) <- 0
rownames(SAFO_8_fst) <- SAFO_8_coords$site
colnames(SAFO_8_fst) <- SAFO_8_coords$site

SAFO_8_coords <- SAFO_8_coords %>%
  select(site, site_name, lat, lon)

# -----------------------------
# 5) map of sampling locations
# -----------------------------
map_df <- map_data("state")

xpad <- 0.18
ypad <- 0.10
xlim_use <- range(SAFO_8_coords$lon) + c(-xpad, xpad)
ylim_use <- range(SAFO_8_coords$lat) + c(-ypad, ypad)

ggplot() +
  geom_polygon(
    data = map_df,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey60",
    linewidth = 0.2
  ) +
  geom_point(
    data = SAFO_8_coords,
    aes(x = lon, y = lat),
    size = 2.2
  ) +
  geom_text(
    data = SAFO_8_coords,
    aes(x = lon, y = lat, label = site_name),
    nudge_y = 0.01,
    size = 2.6
  ) +
  coord_quickmap(
    xlim = xlim_use,
    ylim = ylim_use,
    expand = FALSE
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude"
  )

# -----------------------------
# 6) IBD plot
# straight-line geographic distance as a quick check
# -----------------------------
coord_mat <- as.matrix(SAFO_8_coords[, c("lon", "lat")])
geo_dist_km <- geosphere::distm(coord_mat, fun = geosphere::distHaversine) / 1000

upper_idx <- upper.tri(SAFO_8_fst)

ibd_df <- data.frame(
  fst = SAFO_8_fst[upper_idx],
  dist_km = geo_dist_km[upper_idx]
)

ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Great-circle distance (km)",
    y = expression(F[ST])
  )

# -----------------------------
# 7) save RData
# -----------------------------
save(
  SAFO_8_fst,
  SAFO_8_coords,
  file = file.path(data_dir, "SAFO-8.RData")
)

# -----------------------------
# 8) optional assignment to global env
# -----------------------------
assign("SAFO_8_fst", SAFO_8_fst, envir = .GlobalEnv)
assign("SAFO_8_coords", SAFO_8_coords, envir = .GlobalEnv)

cat("
==================== SUMMARY ====================
")
cat("Study code: ", study_code, "
", sep = "")
cat("Sites: ", nrow(SAFO_8_coords), "
", sep = "")
cat("FST matrix dimensions: ", nrow(SAFO_8_fst), " x ", ncol(SAFO_8_fst), "
", sep = "")
cat("RData saved to: ", file.path(data_dir, "SAFO-8.RData"), "
", sep = "")
cat("=================================================

")
