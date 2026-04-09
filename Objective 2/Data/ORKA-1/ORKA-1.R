# ============================================================
# ORKA-1 | Umpqua Chub
# Oregonichthys kalawatseti
# O'Malley et al. 2013. Transactions of the American Fisheries Society
#
# Objective 2 extraction workflow
# - pairwise FST matrix hard-coded from Table 3
# - representative coordinates hard-coded from user-provided XYs
# - plots shown in RStudio only; not saved
# - saves ORKA_1_fst and ORKA_1_coords to data/ORKA-1.RData
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(geosphere)
})

study_code <- "ORKA-1"
study_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ORKA-1"
out_dir <- file.path(study_dir, "data")

if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

# -----------------------------
# 1) site metadata
# representative coordinates hard-coded from user-provided XYs
# site order matches the FST matrix
# -----------------------------
coords_lookup <- data.frame(
  site = 1:6,
  site_name = c(
    "Smith River",
    "Elk Creek",
    "Calapooya Creek",
    "Olalla Creek",
    "Cow Creek",
    "South Umpqua River"
  ),
  locality_label = c(
    "Smith River centroid",
    "Elk Creek near Hardscrabble Creek",
    "Calapooya Creek centroid",
    "Olalla Creek centroid",
    "Cow Creek centroid",
    "South Umpqua River / 3-C Rock reach centroid"
  ),
  lat = c(
    43.7786469,
    43.6426277,
    43.4395618,
    43.1079612,
    42.8022295,
    42.9341111
  ),
  lon = c(
    -123.9930866,
    -123.2290462,
    -123.2847024,
    -123.5132418,
    -123.5904268,
    -122.991295
  ),
  stringsAsFactors = FALSE
)

ORKA_1_coords <- coords_lookup %>%
  transmute(
    site = site,
    lat = lat,
    lon = lon
  )

stopifnot(identical(ORKA_1_coords$site, 1:6))

# -----------------------------
# 2) pairwise FST matrix
# hard-coded from O'Malley et al. 2013 Table 3
# site order:
# 1 Smith River
# 2 Elk Creek
# 3 Calapooya Creek
# 4 Olalla Creek
# 5 Cow Creek
# 6 South Umpqua River
# -----------------------------
ORKA_1_fst <- matrix(0, nrow = 6, ncol = 6)
rownames(ORKA_1_fst) <- as.character(1:6)
colnames(ORKA_1_fst) <- as.character(1:6)

lower_vals <- list(
  c(2, 1, 0.120),
  c(3, 1, 0.070), c(3, 2, 0.042),
  c(4, 1, 0.083), c(4, 2, 0.040), c(4, 3, 0.010),
  c(5, 1, 0.121), c(5, 2, 0.069), c(5, 3, 0.053), c(5, 4, 0.055),
  c(6, 1, 0.105), c(6, 2, 0.051), c(6, 3, 0.030), c(6, 4, 0.028), c(6, 5, 0.040)
)

for (x in lower_vals) {
  i <- x[1]
  j <- x[2]
  v <- x[3]
  ORKA_1_fst[i, j] <- v
  ORKA_1_fst[j, i] <- v
}

diag(ORKA_1_fst) <- 0
ORKA_1_fst[is.na(ORKA_1_fst)] <- 0
ORKA_1_fst[ORKA_1_fst < 0] <- 0

stopifnot(identical(rownames(ORKA_1_fst), as.character(ORKA_1_coords$site)))

# -----------------------------
# 3) map plot
# show US and Canada, zoom to point extent
# -----------------------------
world_df <- ggplot2::map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

plot_df <- coords_lookup %>%
  mutate(label = paste0(site, " ", site_name))

x_pad <- max(0.5, diff(range(plot_df$lon)) * 0.12)
y_pad <- max(0.4, diff(range(plot_df$lat)) * 0.12)

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
    nudge_y = 0.03,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(min(plot_df$lon) - x_pad, max(plot_df$lon) + x_pad),
    ylim = c(min(plot_df$lat) - y_pad, max(plot_df$lat) + y_pad)
  ) +
  labs(
    title = "ORKA-1 sampling sites",
    subtitle = "Representative coordinates hard-coded from user-provided XYs",
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
  x = as.matrix(ORKA_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  fst = ORKA_1_fst[upper.tri(ORKA_1_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "ORKA-1 IBD plot",
    x = "Geographic distance (km)",
    y = "Pairwise FST"
  ) +
  theme_bw()

print(p_ibd)

# -----------------------------
# 5) save outputs
# save the matrix + coords df in Objective 2 format
# -----------------------------
save(
  ORKA_1_fst,
  ORKA_1_coords,
  file = file.path(out_dir, "ORKA-1.RData")
)

# optional debug export of site names / labels
# write.csv(coords_lookup, file.path(out_dir, "ORKA-1_coords_lookup_debug.csv"), row.names = FALSE)
