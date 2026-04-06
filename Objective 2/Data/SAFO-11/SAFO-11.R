# ============================================================
# SAFO-11 | Brook trout
# Salvelinus fontinalis
# Annett et al. 2012, Transactions of the American Fisheries Society
#
# Objective 2 workflow
# - pairwise Weir and Cockerham FST values transcribed from Table 4
# - coordinates hardcoded from Table 1
# - plots shown in RStudio only; not saved
# - saves SAFO_11_fst and SAFO_11_coords to data/SAFO-11.RData
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(geosphere)
})

study_code <- "SAFO-11"
study_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SAFO-11"
out_dir <- file.path(study_dir, "data")

if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

# -----------------------------
# 1) site metadata
# order follows Tables 1 and 4:
# SA, MA, QU, RB, HA, CO
# coordinates taken directly from Table 1
# lat/lon are listed in decimal minutes in the paper and converted here
# to decimal degrees
# -----------------------------
SAFO_11_coords_check <- data.frame(
  site = 1:6,
  site_code = c("SA", "MA", "QU", "RB", "HA", "CO"),
  site_name = c(
    "Santuit River",
    "Mashpee River",
    "Quashnet River",
    "Red Brook",
    "Sandwich Hatchery",
    "Connetquot River"
  ),
  lat = c(
    41 + 37.672/60,
    41 + 37.300/60,
    41 + 35.533/60,
    41 + 45.915/60,
    41 + 45.159/60,
    40 + 45.783/60
  ),
  lon = c(
    -(70 + 27.062/60),
    -(70 + 28.823/60),
    -(70 + 30.463/60),
    -(70 + 38.035/60),
    -(70 + 29.381/60),
    -(73 + 9.166/60)
  ),
  stringsAsFactors = FALSE
)

SAFO_11_coords <- SAFO_11_coords_check %>%
  transmute(site = site, lat = lat, lon = lon)

# -----------------------------
# 2) pairwise FST matrix
# values transcribed from Table 4 (below diagonal in paper)
# order:
# SA, MA, QU, RB, HA, CO
# -----------------------------
site_ids <- as.character(SAFO_11_coords$site)

SAFO_11_fst <- matrix(0, nrow = 6, ncol = 6)
rownames(SAFO_11_fst) <- site_ids
colnames(SAFO_11_fst) <- site_ids

lower_vals <- list(
  c(2, 1, 0.105),
  c(3, 1, 0.215), c(3, 2, 0.137),
  c(4, 1, 0.159), c(4, 2, 0.129), c(4, 3, 0.105),
  c(5, 1, 0.183), c(5, 2, 0.183), c(5, 3, 0.228), c(5, 4, 0.160),
  c(6, 1, 0.152), c(6, 2, 0.130), c(6, 3, 0.168), c(6, 4, 0.050), c(6, 5, 0.149)
)

for (x in lower_vals) {
  i <- x[1]; j <- x[2]; v <- x[3]
  SAFO_11_fst[i, j] <- v
  SAFO_11_fst[j, i] <- v
}

diag(SAFO_11_fst) <- 0
SAFO_11_fst[is.na(SAFO_11_fst)] <- 0
SAFO_11_fst[SAFO_11_fst < 0] <- 0

stopifnot(identical(as.character(SAFO_11_coords$site), rownames(SAFO_11_fst)))

# -----------------------------
# 3) map plot
# include USA and Canada, then zoom to point extent
# -----------------------------
world_df <- ggplot2::map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

plot_df <- SAFO_11_coords_check %>%
  mutate(label = paste0(site, " ", site_code))

x_pad <- max(0.8, diff(range(plot_df$lon)) * 0.12)
y_pad <- max(0.5, diff(range(plot_df$lat)) * 0.12)

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
    size = 2
  ) +
  geom_text(
    data = plot_df,
    aes(x = lon, y = lat, label = label),
    nudge_y = 0.05,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(min(plot_df$lon) - x_pad, max(plot_df$lon) + x_pad),
    ylim = c(min(plot_df$lat) - y_pad, max(plot_df$lat) + y_pad)
  ) +
  labs(
    title = "SAFO-11 sampling sites",
    subtitle = "Coordinates from Table 1",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw()

print(p_map)

# -----------------------------
# 4) IBD plot
# straight-line distance for QC only
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(SAFO_11_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  fst = SAFO_11_fst[upper.tri(SAFO_11_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "SAFO-11 IBD plot",
    x = "Geographic distance (km)",
    y = "Pairwise FST"
  ) +
  theme_bw()

print(p_ibd)

# -----------------------------
# 5) save outputs
# -----------------------------
save(
  SAFO_11_fst,
  SAFO_11_coords,
  file = file.path(out_dir, "SAFO-11.RData")
)

write.csv(
  SAFO_11_coords_check,
  file = file.path(out_dir, "SAFO-11_coords_check.csv"),
  row.names = FALSE
)
