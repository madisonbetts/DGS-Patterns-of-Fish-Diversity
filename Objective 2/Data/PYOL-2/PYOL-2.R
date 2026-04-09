# ============================================================
# PYOL-2 | Flathead catfish
# Pylodictis olivaris
# Waraniak et al. 2024 J. Fish Biology
#
# Objective 2 extraction workflow
# - HUC-8 / basin-level coordinate proxies hard-coded from Figure 1
# - pairwise FST matrix hard-coded as best estimates from Figure 3a
#   heatmap + text constraints in the paper
# - plots shown in RStudio only; not saved
# - saves PYOL_2_fst and PYOL_2_coords to data/PYOL-2.RData
#
# IMPORTANT:
# The paper does not tabulate the full HUC-8 FST matrix numerically.
# These FST values are best-fit estimates reconstructed from Figure 3a
# and the textual ranges reported in the Results.
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(geosphere)
})

study_code <- "PYOL-2"
study_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/PYOL-2"
out_dir <- file.path(study_dir, "data")

if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

# -----------------------------
# 1) site metadata
# HUC-8 / basin-level coordinate proxies estimated from Figure 1
# and site descriptions in the paper
# -----------------------------
PYOL_2_coords_full <- data.frame(
  site_id = 1:9,
  site_name = c(
    "Delaware River",
    "Backyard Ponds",
    "Schuylkill",
    "Lower Susquehanna",
    "Lower Susquehanna-Penns",
    "Lower Susquehanna-Swatara",
    "Lower Juniata",
    "Upper Susquehanna-Lackawanna",
    "Upper Susquehanna-Tunkhannock"
  ),
  locality_label = c(
    "Single Delaware River sample near Trenton / lower reach proxy",
    "Private pond near Collegeville, Pennsylvania",
    "Schuylkill River centroid of sampled reach (73–113 km upstream of mouth)",
    "Lower Susquehanna centroid of sampled reach",
    "Penns Creek confluence / Lower Susquehanna-Penns proxy",
    "Swatara Creek confluence / Lower Susquehanna-Swatara proxy",
    "Lower Juniata River centroid of sampled reach",
    "Upper Susquehanna-Lackawanna basin centroid proxy",
    "Upper Susquehanna-Tunkhannock basin centroid proxy"
  ),
  lat = c(
    40.2200,
    40.1700,
    40.1200,
    39.6200,
    40.5400,
    40.3700,
    40.5100,
    41.2100,
    41.4700
  ),
  lon = c(
    -74.9300,
    -75.4500,
    -75.4300,
    -76.1200,
    -76.9300,
    -76.6000,
    -77.1900,
    -75.8600,
    -75.7800
  ),
  stringsAsFactors = FALSE
)

PYOL_2_coords <- PYOL_2_coords_full %>%
  select(site_id, lat, lon)

write.csv(
  PYOL_2_coords_full,
  file = file.path(study_dir, "PYOL-2_coords.csv"),
  row.names = FALSE
)

# -----------------------------
# 2) pairwise FST matrix
# best-fit reconstruction from Figure 3a heatmap plus text constraints
#
# text constraints used:
# - Susquehanna vs Schuylkill: FST = 0.155–0.180
# - Susquehanna vs Delaware: FST = 0.209–0.256
# - Backyard Ponds vs Schuylkill: FST = 0.161
# - Backyard Ponds vs Susquehanna: FST = -0.038 to 0.010,
#   truncated here to 0 where needed
# - within Susquehanna basins: FST = 0.002–0.033
# -----------------------------
PYOL_2_fst <- matrix(0, nrow = 9, ncol = 9)
rownames(PYOL_2_fst) <- as.character(1:9)
colnames(PYOL_2_fst) <- as.character(1:9)

# Delaware River
PYOL_2_fst[1,2] <- 0.180
PYOL_2_fst[1,3] <- 0.170
PYOL_2_fst[1,4] <- 0.220
PYOL_2_fst[1,5] <- 0.230
PYOL_2_fst[1,6] <- 0.220
PYOL_2_fst[1,7] <- 0.240
PYOL_2_fst[1,8] <- 0.250
PYOL_2_fst[1,9] <- 0.240

# Backyard Ponds
PYOL_2_fst[2,3] <- 0.161
PYOL_2_fst[2,4] <- 0.005
PYOL_2_fst[2,5] <- 0.000
PYOL_2_fst[2,6] <- 0.008
PYOL_2_fst[2,7] <- 0.010
PYOL_2_fst[2,8] <- 0.006
PYOL_2_fst[2,9] <- 0.004

# Schuylkill
PYOL_2_fst[3,4] <- 0.155
PYOL_2_fst[3,5] <- 0.160
PYOL_2_fst[3,6] <- 0.158
PYOL_2_fst[3,7] <- 0.170
PYOL_2_fst[3,8] <- 0.180
PYOL_2_fst[3,9] <- 0.175

# Susquehanna / Juniata basins
PYOL_2_fst[4,5] <- 0.003
PYOL_2_fst[4,6] <- 0.008
PYOL_2_fst[4,7] <- 0.015
PYOL_2_fst[4,8] <- 0.025
PYOL_2_fst[4,9] <- 0.030

PYOL_2_fst[5,6] <- 0.005
PYOL_2_fst[5,7] <- 0.010
PYOL_2_fst[5,8] <- 0.020
PYOL_2_fst[5,9] <- 0.025

PYOL_2_fst[6,7] <- 0.012
PYOL_2_fst[6,8] <- 0.018
PYOL_2_fst[6,9] <- 0.022

PYOL_2_fst[7,8] <- 0.028
PYOL_2_fst[7,9] <- 0.032

PYOL_2_fst[8,9] <- 0.010

PYOL_2_fst <- PYOL_2_fst + t(PYOL_2_fst)
diag(PYOL_2_fst) <- 0
PYOL_2_fst[is.na(PYOL_2_fst)] <- 0
PYOL_2_fst[PYOL_2_fst < 0] <- 0

stopifnot(isTRUE(all.equal(PYOL_2_fst, t(PYOL_2_fst))))
stopifnot(identical(rownames(PYOL_2_fst), as.character(PYOL_2_coords$site_id)))

# -----------------------------
# 3) map plot
# include USA and Canada and zoom to point extent
# -----------------------------
world_df <- ggplot2::map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

plot_df <- PYOL_2_coords_full %>%
  mutate(label = paste0(site_id, " ", site_name))

x_pad <- max(0.8, diff(range(plot_df$lon)) * 0.18)
y_pad <- max(0.8, diff(range(plot_df$lat)) * 0.18)

p_map <- ggplot() +
  geom_polygon(
    data = world_df,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = plot_df,
    aes(x = lon, y = lat),
    size = 2.2
  ) +
  geom_text(
    data = plot_df,
    aes(x = lon, y = lat, label = label),
    nudge_y = 0.08,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(min(plot_df$lon) - x_pad, max(plot_df$lon) + x_pad),
    ylim = c(min(plot_df$lat) - y_pad, max(plot_df$lat) + y_pad)
  ) +
  theme_bw() +
  labs(
    title = "PYOL-2 sampling units",
    x = "Longitude",
    y = "Latitude"
  )

print(p_map)

# -----------------------------
# 4) IBD plot
# straight-line distance for QC plotting only
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(PYOL_2_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  fst = PYOL_2_fst[upper.tri(PYOL_2_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "PYOL-2 IBD plot",
    x = "Geographic distance (km)",
    y = "Pairwise FST"
  ) +
  theme_bw()

print(p_ibd)

# -----------------------------
# 5) save outputs
# -----------------------------
save(
  PYOL_2_fst,
  PYOL_2_coords,
  file = file.path(out_dir, "PYOL-2.RData")
)
