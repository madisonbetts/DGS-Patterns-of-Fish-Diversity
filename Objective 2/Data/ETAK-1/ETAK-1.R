
# -----------------------------
# ETAK-1 Bluemask Darter
# tributary coordinates + FST matrix + IBD plot
# Taylor et al. 2021, Molecular Ecology
# Pairwise FST from Supplemental Table S1 (100% individuals-per-locus branch)
# -----------------------------

library(ggplot2)
library(geosphere)
library(maps)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETAK-1"

# -----------------------------
# 0) site coordinates
# tributary-level representative coordinates
# these are basin/tributary centroids matched to the coarse scale of the FST data
# order follows Supplemental Table S1
# 1 = Rocky River
# 2 = Cane Creek
# 3 = Caney Fork
# 4 = Collins River
# -----------------------------
ETAK_1_coords <- data.frame(
  site = as.character(1:4),
  #site_name = c("Rocky River", "Cane Creek", "Caney Fork", "Collins River"),
  lat = c(
    35.731,
    35.745,
    35.805,
    35.687
  ),
  lon = c(
    -85.620,
    -85.600,
    -85.530,
    -85.790
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map of sampling tributaries
# -----------------------------
states_map <- map_data("state")
tn_map <- subset(states_map, region %in% c("tennessee"))

map_xlim <- c(-86.1, -85.2)
map_ylim <- c(35.5, 36.0)

ggplot() +
  geom_polygon(
    data = tn_map,
    aes(x = long, y = lat, group = group),
    fill = "grey92", color = "grey40", linewidth = 0.3
  ) +
  geom_point(
    data = ETAK_1_coords,
    aes(x = lon, y = lat),
    size = 2.7
  ) +
  geom_text(
    data = ETAK_1_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.02,
    size = 3.3
  ) +
  coord_fixed(xlim = map_xlim, ylim = map_ylim) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "ETAK-1 sampling tributaries"
  )

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(ETAK_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- ETAK_1_coords$site
colnames(geo_dist_km) <- ETAK_1_coords$site

# -----------------------------
# 3) pairwise FST matrix
# Supplemental Table S1, 100% individuals-per-locus branch
# negative values -> 0 (none here)
# -----------------------------
ETAK_1_fst <- matrix(
  c(
    0.00, 0.11, 0.09, 0.20,
    0.11, 0.00, 0.05, 0.29,
    0.09, 0.05, 0.00, 0.27,
    0.20, 0.29, 0.27, 0.00
  ),
  nrow = 4,
  byrow = TRUE
)

rownames(ETAK_1_fst) <- colnames(ETAK_1_fst) <- ETAK_1_coords$site
ETAK_1_fst[ETAK_1_fst < 0] <- 0

# -----------------------------
# 4) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(ETAK_1_fst)[row(ETAK_1_fst)[upper.tri(ETAK_1_fst)]],
  site2   = colnames(ETAK_1_fst)[col(ETAK_1_fst)[upper.tri(ETAK_1_fst)]],
  fst     = ETAK_1_fst[upper.tri(ETAK_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 5) IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "ETAK-1 isolation by distance"
  )

# -----------------------------
# 6) quick checks
# -----------------------------
stopifnot(identical(rownames(ETAK_1_fst), ETAK_1_coords$site))
stopifnot(identical(colnames(ETAK_1_fst), ETAK_1_coords$site))
stopifnot(isTRUE(all.equal(ETAK_1_fst, t(ETAK_1_fst))))

# -----------------------------
# 7) save RData
# -----------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  ETAK_1_fst,
  ETAK_1_coords,
  file = file.path(out_dir, "ETAK-1.RData")
)
