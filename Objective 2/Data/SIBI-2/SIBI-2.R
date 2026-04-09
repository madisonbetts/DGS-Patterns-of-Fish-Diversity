# -----------------------------
# SIBI-1
# Siphateles thalassinus / tui chub paper
# Full 8-basin microsatellite FST matrix from Table 3
# Outputs: SIBI_1_fst, SIBI_1_coords
# Save: SIBI-1.RData in data subfolder
# -----------------------------

library(geosphere)
library(ggplot2)
library(ggrepel)
library(tigris)
library(dplyr)

options(tigris_use_cache = TRUE)

# -----------------------------
# paths
# -----------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SIBI-1"
data_dir <- file.path(base_dir, "data")

# -----------------------------
# 1) coordinate dataframe
# Representative coordinates for the 8 basin groups used in Table 3.
# These are basin anchors / centroids based on the named sampling
# localities in Table 1.
# -----------------------------
SIBI_1_coords <- data.frame(
  site_num = 1:8,
  basin_code = c("CH", "WV", "GL", "PR", "LA", "SL", "CV", "GV"),
  basin = c(
    "Cow Head",
    "Warner",
    "Goose Lake",
    "Pit River",
    "Lake Abert",
    "Summer Lake",
    "Catlow Valley",
    "Guano Valley"
  ),
  site_name = c(
    "Cow Head Slough",
    "Twenty Mile Slough",
    "Goose Lake basin centroid",
    "Pit River / Big Sage centroid",
    "Crooked Creek",
    "Summer Lake group centroid",
    "Catlow group centroid",
    "Fish Creek"
  ),
  lon = c(
    -120.9810,
    -119.9950,
    -120.4710,
    -120.3100,
    -120.7610,
    -120.8450,
    -118.9520,
    -119.2570
  ),
  lat = c(
    41.5280,
    42.1840,
    41.8390,
    41.2700,
    42.6400,
    43.0250,
    42.2650,
    41.9950
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) pairwise FST matrix
# Table 3 values below the diagonal in basin order:
# CH, WV, GL, PR, LA, SL, CV, GV
# -----------------------------
SIBI_1_fst <- matrix(0, nrow = 8, ncol = 8)
rownames(SIBI_1_fst) <- colnames(SIBI_1_fst) <- as.character(1:8)

SIBI_1_fst["2", "1"] <- 0.13

SIBI_1_fst["3", "1"] <- 0.21
SIBI_1_fst["3", "2"] <- 0.10

SIBI_1_fst["4", "1"] <- 0.23
SIBI_1_fst["4", "2"] <- 0.12
SIBI_1_fst["4", "3"] <- 0.03

SIBI_1_fst["5", "1"] <- 0.22
SIBI_1_fst["5", "2"] <- 0.17
SIBI_1_fst["5", "3"] <- 0.13
SIBI_1_fst["5", "4"] <- 0.18

SIBI_1_fst["6", "1"] <- 0.26
SIBI_1_fst["6", "2"] <- 0.18
SIBI_1_fst["6", "3"] <- 0.15
SIBI_1_fst["6", "4"] <- 0.19
SIBI_1_fst["6", "5"] <- 0.18

SIBI_1_fst["7", "1"] <- 0.20
SIBI_1_fst["7", "2"] <- 0.19
SIBI_1_fst["7", "3"] <- 0.24
SIBI_1_fst["7", "4"] <- 0.27
SIBI_1_fst["7", "5"] <- 0.25
SIBI_1_fst["7", "6"] <- 0.27

SIBI_1_fst["8", "1"] <- 0.19
SIBI_1_fst["8", "2"] <- 0.20
SIBI_1_fst["8", "3"] <- 0.24
SIBI_1_fst["8", "4"] <- 0.28
SIBI_1_fst["8", "5"] <- 0.23
SIBI_1_fst["8", "6"] <- 0.24
SIBI_1_fst["8", "7"] <- 0.18

SIBI_1_fst[upper.tri(SIBI_1_fst)] <- t(SIBI_1_fst)[upper.tri(SIBI_1_fst)]
SIBI_1_fst[SIBI_1_fst < 0] <- 0

# -----------------------------
# 3) map of sampling locations
# Plot study extent with full US states for context
# -----------------------------
# -----------------------------
# 3) map of sampling locations
# -----------------------------
library(tigris)
library(sf)

states_sf <- tigris::states(cb = TRUE, year = 2022) |>
  sf::st_transform(4326)

x_pad <- 2.5
y_pad <- 2.0

xlim_use <- range(SIBI_1_coords$lon) + c(-x_pad, x_pad)
ylim_use <- range(SIBI_1_coords$lat) + c(-y_pad, y_pad)

p_map <- ggplot() +
  geom_sf(
    data = states_sf,
    fill = "gray95",
    color = "gray70",
    linewidth = 0.25
  ) +
  geom_point(
    data = SIBI_1_coords,
    aes(x = lon, y = lat),
    size = 3
  ) +
  geom_text_repel(
    data = SIBI_1_coords,
    aes(x = lon, y = lat, label = paste0(site_num, ": ", basin_code)),
    size = 3.2,
    max.overlaps = 100
  ) +
  coord_sf(
    xlim = xlim_use,
    ylim = ylim_use,
    expand = FALSE
  ) +
  theme_classic() +
  labs(
    title = "SIBI-1 sampling basins",
    x = "Longitude",
    y = "Latitude"
  )

print(p_map)
# -----------------------------
# 4) geographic distance matrix
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = SIBI_1_coords[, c("lon", "lat")],
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- colnames(geo_dist_km) <- as.character(SIBI_1_coords$site_num)

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(SIBI_1_fst)[row(SIBI_1_fst)[upper.tri(SIBI_1_fst)]],
  site2   = colnames(SIBI_1_fst)[col(SIBI_1_fst)[upper.tri(SIBI_1_fst)]],
  fst     = SIBI_1_fst[upper.tri(SIBI_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 6) IBD plot
# -----------------------------
p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    title = "SIBI-1 IBD",
    x = "Geographic distance (km)",
    y = expression(F[ST])
  )

print(p_ibd)

# -----------------------------
# 7) save objects
# -----------------------------
save(
  SIBI_1_fst,
  SIBI_1_coords,
  file = file.path(data_dir, "SIBI-1.RData")
)
