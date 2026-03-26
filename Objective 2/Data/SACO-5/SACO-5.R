# ============================================================
# SACO-5
# Salvelinus confluentus (bull trout)
# DeHaan et al. 2011. Northwest Science 85:463-475.
#
# Objective 2 workflow:
#   1) build pairwise FST matrix from published Table 4
#   2) build best-available site coordinates for Elwha River sections
#   3) plot sampling sites with U.S. + Canada context
#   4) plot IBD
#   5) save SACO-5.RData in /data subfolder
#
# NOTE:
# This paper does not provide a raw genotype table in the article.
# Pairwise FST values below were transcribed directly from Table 4.
# Coordinates are best-available approximate centroids for the four
# Elwha River sections shown in Figure 1 and described in the text:
#   lower = below Elwha Dam
#   middle = between Elwha Dam and Glines Canyon Dam
#   upper = above Glines Canyon Dam to ~rkm 44
#   headwaters = ~rkm 58-65, upstream of Carlson Canyon
#
# Coordinate rationale:
#   - river mouth from USGS near-mouth station
#   - Elwha Dam is ~4.9 miles upstream from the mouth
#   - Glines Canyon Dam is at river kilometer 21.6
#   - Carlson Canyon is ~33 miles upstream from the mouth
#   - section centroids were tightened against Figure 1 in the paper
#
# These XYs are appropriate for broad-scale Objective 2 IBD work,
# but they should still be treated as approximations rather than
# exact capture coordinates.
# ============================================================

rm(list = ls())

# -----------------------------
# packages
# -----------------------------
libs <- c("ggplot2", "dplyr", "maps", "geosphere")

for (p in libs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Package not installed: ", p)
  }
}

library(ggplot2)
library(dplyr)
library(maps)
library(geosphere)

# -----------------------------
# paths
# -----------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"
study_dir <- file.path(base_dir, "SACO-5")
data_dir  <- file.path(study_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) site coordinates
# best-available approximate centroids for the four Elwha sections
# ordered downstream -> upstream to match Table 4 interpretation
# -----------------------------
SACO_5_coords <- data.frame(
  site = as.character(1:4),
  site_name = c(
    "Lower Elwha",
    "Middle Elwha",
    "Upper Elwha",
    "Elwha Headwaters"
  ),
  lat = c(
    48.1180,  # below Elwha Dam
    48.0505,  # between Elwha Dam and Glines Canyon Dam
    47.9365,  # above Glines Canyon Dam to ~rkm 44
    47.8320   # ~rkm 58-65, upstream of Carlson Canyon
  ),
  lon = c(
    -123.5580,
    -123.5765,
    -123.6480,
    -123.7210
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) pairwise FST matrix
# transcribed from Table 4 in DeHaan et al. 2011
# lower triangle values:
#   headwaters-lower  = 0.072
#   headwaters-middle = 0.088
#   headwaters-upper  = 0.072
#   lower-middle      = 0.004
#   lower-upper       = 0.004
#   middle-upper      = 0.007
# negative values are set to 0 by workflow rule
# -----------------------------
SACO_5_fst <- matrix(0, nrow = 4, ncol = 4)
rownames(SACO_5_fst) <- colnames(SACO_5_fst) <- SACO_5_coords$site

SACO_5_fst["1", "2"] <- 0.004
SACO_5_fst["1", "3"] <- 0.004
SACO_5_fst["1", "4"] <- 0.072
SACO_5_fst["2", "3"] <- 0.007
SACO_5_fst["2", "4"] <- 0.088
SACO_5_fst["3", "4"] <- 0.072

SACO_5_fst[lower.tri(SACO_5_fst)] <- t(SACO_5_fst)[lower.tri(SACO_5_fst)]
diag(SACO_5_fst) <- 0
SACO_5_fst[SACO_5_fst < 0] <- 0

# -----------------------------
# 3) geographic distance matrix
# Euclidean / great-circle distance in km
# -----------------------------
geo_dist_m <- geosphere::distm(
  x = SACO_5_coords[, c("lon", "lat")],
  fun = geosphere::distHaversine
)
geo_dist_km <- geo_dist_m / 1000
rownames(geo_dist_km) <- colnames(geo_dist_km) <- SACO_5_coords$site

# -----------------------------
# 4) pairwise dataframe for IBD
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(SACO_5_fst)[row(SACO_5_fst)[upper.tri(SACO_5_fst)]],
  site2   = colnames(SACO_5_fst)[col(SACO_5_fst)[upper.tri(SACO_5_fst)]],
  fst     = SACO_5_fst[upper.tri(SACO_5_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
) %>%
  left_join(
    SACO_5_coords %>% select(site, site1_name = site_name),
    by = c("site1" = "site")
  ) %>%
  left_join(
    SACO_5_coords %>% select(site, site2_name = site_name),
    by = c("site2" = "site")
  )

# -----------------------------
# 5) map of sites
# -----------------------------
us_map <- map_data("state")
can_map <- map_data("world") %>%
  filter(region == "Canada")

p_map <- ggplot() +
  geom_polygon(
    data = can_map,
    aes(x = long, y = lat, group = group),
    fill = "grey92",
    color = "grey55",
    linewidth = 0.2
  ) +
  geom_polygon(
    data = us_map,
    aes(x = long, y = lat, group = group),
    fill = "grey96",
    color = "grey55",
    linewidth = 0.2
  ) +
  geom_point(
    data = SACO_5_coords,
    aes(x = lon, y = lat),
    size = 3
  ) +
  geom_text(
    data = SACO_5_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.08,
    size = 4
  ) +
  coord_fixed(1.3, xlim = c(-130, -116), ylim = c(45, 50.5)) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "SACO-5 bull trout sampling sections",
    subtitle = "Elwha River, Washington"
  )

# -----------------------------
# 6) IBD plot
# -----------------------------
p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "Isolation by distance for SACO-5"
  )

print(p_map)
print(p_ibd)

# -----------------------------
# 7) save outputs
# -----------------------------
SACO_5_coords <- SACO_5_coords[, -2] # remove unecssary cols 

# export fst and xys
save(
  SACO_5_coords,
  SACO_5_fst,
  #geo_dist_km,
  #ibd_df,
  file = file.path(data_dir, "SACO-5.RData")
)
