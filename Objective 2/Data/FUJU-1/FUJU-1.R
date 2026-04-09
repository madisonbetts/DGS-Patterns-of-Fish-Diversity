# ============================================================
# FUJU-1
# Barrens Topminnow (Fundulus julisia)
# Hurt et al. 2017
#
# Clean workflow:
# - hardcodes the published 10 x 10 microsatellite pairwise FST matrix
# - uses best-available approximate site coordinates
# - plots sites + quick IBD in RStudio
# - saves FUJU_1_fst and FUJU_1_coords to data/FUJU-1.RData
#
# IMPORTANT
# - Coordinates are inferred / approximate, not field GPS points.
# - FST values are transcribed from the published pairwise matrix.
# - Site order follows the paper:
#     1 McMahan Creek
#     2 Pedigo Highway
#     3 Pedigo Farm
#     4 Benedict Spring
#     5 Pond Spring
#     6 Clayborne Spring
#     7 Collier Spring
#     8 Short Spring
#     9 Merkle Big Spring
#    10 Farris Spring
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(geosphere)
})

# -----------------------------
# 0) paths
# -----------------------------
study_code <- "FUJU-1"

base_dir <- file.path(
  "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026",
  "Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data",
  study_code
)

data_dir <- file.path(base_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) site metadata in paper order
# -----------------------------
site_lookup <- tibble::tribble(
  ~site, ~site_name,
   1, "McMahan Creek",
   2, "Pedigo Highway",
   3, "Pedigo Farm",
   4, "Benedict Spring",
   5, "Pond Spring",
   6, "Clayborne Spring",
   7, "Collier Spring",
   8, "Short Spring",
   9, "Merkle Big Spring",
  10, "Farris Spring"
)

# -----------------------------
# 2) hardcoded FST matrix
# values from Hurt et al. 2017 Table 4
#
# Published order in Table 4:
# Benedict, Pond, McMahan, Pedigo H., Pedigo F., Clayborne,
# Collier, Farris, Merkle, Short
#
# Workflow order here:
# McMahan, Pedigo H., Pedigo F., Benedict, Pond, Clayborne,
# Collier, Short, Merkle, Farris
# -----------------------------
FUJU_1_fst <- matrix(c(
  # row 1 = McMahan Creek
  0.0000, 0.2792, 0.2939, 0.6634, 0.4673, 0.6676, 0.3393, 0.2552, 0.4104, 0.3919,
  # row 2 = Pedigo Highway
  0.2792, 0.0000, 0.0744, 0.5469, 0.3976, 0.5294, 0.0978, 0.1662, 0.3556, 0.3339,
  # row 3 = Pedigo Farm
  0.2939, 0.0744, 0.0000, 0.5513, 0.3401, 0.5420, 0.1222, 0.2233, 0.2977, 0.2775,
  # row 4 = Benedict Spring
  0.6634, 0.5469, 0.5513, 0.0000, 0.5365, 0.0007, 0.5899, 0.2573, 0.5720, 0.5873,
  # row 5 = Pond Spring
  0.4673, 0.3976, 0.3401, 0.5365, 0.0000, 0.5271, 0.4142, 0.3727, 0.0639, 0.0618,
  # row 6 = Clayborne Spring
  0.6676, 0.5294, 0.5420, 0.0007, 0.5271, 0.0000, 0.5771, 0.2560, 0.5559, 0.5713,
  # row 7 = Collier Spring
  0.3393, 0.0978, 0.1222, 0.5899, 0.4142, 0.5771, 0.0000, 0.2309, 0.3685, 0.3548,
  # row 8 = Short Spring
  0.2552, 0.1662, 0.2233, 0.2573, 0.3727, 0.2560, 0.2309, 0.0000, 0.3372, 0.3262,
  # row 9 = Merkle Big Spring
  0.4104, 0.3556, 0.2977, 0.5720, 0.0639, 0.5559, 0.3685, 0.3372, 0.0000, 0.0193,
  # row 10 = Farris Spring
  0.3919, 0.3339, 0.2775, 0.5873, 0.0618, 0.5713, 0.3548, 0.3262, 0.0193, 0.0000
), nrow = 10, byrow = TRUE)

FUJU_1_fst[is.na(FUJU_1_fst)] <- 0
FUJU_1_fst[FUJU_1_fst < 0] <- 0
diag(FUJU_1_fst) <- 0
rownames(FUJU_1_fst) <- colnames(FUJU_1_fst) <- as.character(site_lookup$site)

# -----------------------------
# 3) approximate coordinates
# best-available inferred locations
# -----------------------------
FUJU_1_coords_full <- tibble::tribble(
  ~site, ~site_name,            ~lat,     ~lon,
   1,    "McMahan Creek",       35.6940, -86.0610,
   2,    "Pedigo Highway",      35.7030, -86.0594,
   3,    "Pedigo Farm",         35.7015, -86.0625,
   4,    "Benedict Spring",     35.5635, -86.0215,
   5,    "Pond Spring",         35.2650, -86.2450,
   6,    "Clayborne Spring",    35.5585, -86.0285,
   7,    "Collier Spring",      35.3950, -86.4970,
   8,    "Short Spring",        35.3570, -86.1900,
   9,    "Merkle Big Spring",   35.2070, -86.1650,
  10,    "Farris Spring",       35.1292, -86.1345
)

stopifnot(nrow(FUJU_1_coords_full) == 10)

FUJU_1_coords <- FUJU_1_coords_full %>%
  select(site, lat, lon)

# -----------------------------
# 4) map of sampling locations
# -----------------------------
map_df <- map_data("state")

xpad <- 0.18
ypad <- 0.10

ggplot() +
  geom_polygon(
    data = map_df,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = FUJU_1_coords_full,
    aes(x = lon, y = lat),
    size = 2.5
  ) +
  geom_text(
    data = FUJU_1_coords_full,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.01,
    size = 3
  ) +
  coord_quickmap(
    xlim = range(FUJU_1_coords_full$lon) + c(-xpad, xpad),
    ylim = range(FUJU_1_coords_full$lat) + c(-ypad, ypad),
    expand = FALSE
  ) +
  theme_classic() +
  labs(
    title = "FUJU-1 sampling locations",
    x = "Longitude",
    y = "Latitude"
  )

# -----------------------------
# 5) IBD plot
# straight-line distance quick check
# -----------------------------
coord_mat <- as.matrix(FUJU_1_coords[, c("lon", "lat")])
geo_dist_km <- geosphere::distm(coord_mat, fun = geosphere::distHaversine) / 1000

upper_idx <- upper.tri(FUJU_1_fst)

ibd_df <- data.frame(
  fst = FUJU_1_fst[upper_idx],
  dist_km = geo_dist_km[upper_idx]
)

ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    title = "FUJU-1 IBD quick-check",
    x = "Great-circle distance (km)",
    y = expression(F[ST])
  )

# -----------------------------
# 6) save RData
# -----------------------------
save(
  FUJU_1_fst,
  FUJU_1_coords,
  file = file.path(data_dir, "FUJU-1.RData")
)

assign("FUJU_1_fst", FUJU_1_fst, envir = .GlobalEnv)
assign("FUJU_1_coords", FUJU_1_coords, envir = .GlobalEnv)

cat("\n==================== SUMMARY ====================\n")
cat("Study code: ", study_code, "\n", sep = "")
cat("Sites: ", nrow(FUJU_1_coords), "\n", sep = "")
cat("FST matrix dimensions: ", nrow(FUJU_1_fst), " x ", ncol(FUJU_1_fst), "\n", sep = "")
cat("RData saved to: ", file.path(data_dir, "FUJU-1.RData"), "\n", sep = "")
cat("=================================================\n\n")
