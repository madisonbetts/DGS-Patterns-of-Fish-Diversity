# ------------------------------
# LEAU-1 redbreast sunfish
# site coordinates + FST matrix + IBD plot
# ------------------------------

library(ggplot2)
library(geosphere)
library(dplyr)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/LEAU-1/"

sites_path <- file.path(save_dir, "redbreast_sunfish_sites.csv")
fst_path   <- file.path(save_dir, "redbreast_sunfish_fst_matrix.csv")

# ------------------------------
# 0) read data
# ------------------------------
LEAU_1_coords <- read.csv(sites_path, stringsAsFactors = FALSE)
LEAU_1_fst <- as.matrix(read.csv(fst_path, row.names = 1, check.names = FALSE))

# ------------------------------
# 1) optional: drop S7
# thesis says only 1 individual remained there
# ------------------------------
drop_sites <- c("S7")

LEAU_1_coords <- LEAU_1_coords %>%
  filter(!site %in% drop_sites)

LEAU_1_fst <- LEAU_1_fst[
  rownames(LEAU_1_fst) %in% LEAU_1_coords$site,
  colnames(LEAU_1_fst) %in% LEAU_1_coords$site
]

# reorder to match site dataframe
LEAU_1_fst <- LEAU_1_fst[LEAU_1_coords$site, LEAU_1_coords$site]
diag(LEAU_1_fst) <- 0

# ------------------------------
# clean FST matrix
# set negatives to 0 and fix diagonal
# ------------------------------
LEAU_1_fst[LEAU_1_fst < 0] <- 0
diag(LEAU_1_fst) <- 0

# keep only site + coordinates for final object
LEAU_1_coords <- LEAU_1_coords[, c("site", "lat", "lon")]

# ------------------------------
# 2) map of sampling locations
# ------------------------------
ggplot(LEAU_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.05, size = 3.5) +
  coord_equal() +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "LEAU-1 sampling locations"
  )

# ------------------------------
# 3) geographic distance matrix
# straight-line distance in km
# ------------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(LEAU_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- LEAU_1_coords$site
colnames(geo_dist_km) <- LEAU_1_coords$site

# ------------------------------
# 4) pairwise dataframe for IBD plot
# ------------------------------
ibd_df <- data.frame(
  site1 = rownames(LEAU_1_fst)[row(LEAU_1_fst)[upper.tri(LEAU_1_fst)]],
  site2 = colnames(LEAU_1_fst)[col(LEAU_1_fst)[upper.tri(LEAU_1_fst)]],
  distance_km = geo_dist_km[upper.tri(geo_dist_km)],
  fst = LEAU_1_fst[upper.tri(LEAU_1_fst)]
)

# ------------------------------
# 5) IBD plot
# ------------------------------
ggplot(ibd_df, aes(x = distance_km, y = fst)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(
    x = "Geographic distance (km)",
    y = expression(F[ST]),
    title = "LEAU-1 isolation by distance"
  )

# ------------------------------
# 6) quick checks
# ------------------------------
stopifnot(identical(rownames(LEAU_1_fst), LEAU_1_coords$site))
stopifnot(identical(colnames(LEAU_1_fst), LEAU_1_coords$site))
stopifnot(isTRUE(all.equal(LEAU_1_fst, t(LEAU_1_fst))))

# ------------------------------
# 7) save RData
# ------------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  LEAU_1_fst,
  LEAU_1_coords,
  file = file.path(out_dir, "LEAU-1.RData")
)