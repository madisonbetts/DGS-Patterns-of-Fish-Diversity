# -----------------------------
# ICPU-2 channel catfish
# site coordinates + FST matrix + IBD plot
# -----------------------------

library(ggplot2)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ICPU-2"

# -----------------------------
# 0) site coordinates
# GPS coordinates converted to decimal degrees
# -----------------------------
ICPU_2_coords <- data.frame(
  site = c("Ot1", "Ot2", "Ot3", "Ot4", "Ot5", "Mis", "Mad"),
  lat = c(
    45.462889,  # N45°27'46.4"
    45.513222,  # N45°30'47.6"
    45.496667,  # N45°29'48.0"
    45.518889,  # N45°31'08.0"
    45.447528,  # N45°26'51.0"
    45.421306,  # N45°25'14.7"
    45.442478   # N45°26'32.92"
  ),
  lon = c(
    -76.386944, # W076°23'13.0"
    -76.504778, # W076°30'17.2"
    -76.444722, # W076°26'47.0"
    -76.539444, # W076°32'22.0"
    -76.317833, # W076°19'04.2"
    -76.261111, # W076°15'40.0"
    -76.348472  # W076°20'54.9"
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map of sampling locations
# -----------------------------
ggplot(ICPU_2_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.01, size = 4) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "ICPU-2 sampling locations"
  )

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(ICPU_2_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- ICPU_2_coords$site
colnames(geo_dist_km) <- ICPU_2_coords$site

# -----------------------------
# 3) pairwise FST matrix
# Table 2 caption:
# FST values are BELOW the diagonal
# p-values are ABOVE the diagonal
# -----------------------------
ICPU_2_fst <- matrix(
  0,
  nrow = nrow(ICPU_2_coords),
  ncol = nrow(ICPU_2_coords),
  dimnames = list(ICPU_2_coords$site, ICPU_2_coords$site)
)

ICPU_2_fst["Ot2", "Ot1"] <- 0.6985

ICPU_2_fst["Ot3", c("Ot1", "Ot2")] <- c(0.3248, 0.6091)

ICPU_2_fst["Ot4", c("Ot1", "Ot2", "Ot3")] <- c(0.4167, 0.5661, 0.6578)

ICPU_2_fst["Ot5", c("Ot1", "Ot2", "Ot3", "Ot4")] <- c(0.0316, 0.2007, 0.1882, 0.1180)

ICPU_2_fst["Mis", c("Ot1", "Ot2", "Ot3", "Ot4", "Ot5")] <- c(0.0154, 0.0984, 0.1371, 0.0547, 0.8446)

ICPU_2_fst["Mad", c("Ot1", "Ot2", "Ot3", "Ot4", "Ot5", "Mis")] <- c(0.4425, 0.5765, 0.4293, 0.3446, 0.2254, 0.2165)

# reflect lower triangle to upper triangle
ICPU_2_fst[upper.tri(ICPU_2_fst)] <- t(ICPU_2_fst)[upper.tri(ICPU_2_fst)]
diag(ICPU_2_fst) <- 0
ICPU_2_fst[ICPU_2_fst < 0] <- 0

# -----------------------------
# 4) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(ICPU_2_fst)[row(ICPU_2_fst)[upper.tri(ICPU_2_fst)]],
  site2   = colnames(ICPU_2_fst)[col(ICPU_2_fst)[upper.tri(ICPU_2_fst)]],
  fst     = ICPU_2_fst[upper.tri(ICPU_2_fst)],
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
    title = "ICPU-2 isolation by distance"
  )

# -----------------------------
# 6) quick checks
# -----------------------------
stopifnot(identical(rownames(ICPU_2_fst), ICPU_2_coords$site))
stopifnot(identical(colnames(ICPU_2_fst), ICPU_2_coords$site))
stopifnot(isTRUE(all.equal(ICPU_2_fst, t(ICPU_2_fst))))

# -----------------------------
# 7) save RData
# -----------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  ICPU_2_fst,
  ICPU_2_coords,
  file = file.path(out_dir, "ICPU-2.RData")
)