# -----------------------------
# NOGI-1
# -----------------------------

library(ggplot2)
library(geosphere)
library(maps)

save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/NOGI-1"

# -----------------------------
# 1) site coordinates
# population centroids digitized from Fig. 1 + supplement
# -----------------------------
NOGI_1_coords <- data.frame(
  site = as.character(1:5),
  lat = c(35.82, 35.74, 35.68, 35.61, 35.55),
  lon = c(-85.52, -85.60, -85.66, -85.72, -85.80),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) FST matrix (MOESM Table S2)
# -----------------------------
NOGI_1_fst <- matrix(c(
  0.00, 0.08, 0.12, 0.18, 0.25,
  0.08, 0.00, 0.05, 0.14, 0.22,
  0.12, 0.05, 0.00, 0.10, 0.20,
  0.18, 0.14, 0.10, 0.00, 0.15,
  0.25, 0.22, 0.20, 0.15, 0.00
), nrow = 5, byrow = TRUE)

rownames(NOGI_1_fst) <- colnames(NOGI_1_fst) <- NOGI_1_coords$site
NOGI_1_fst[NOGI_1_fst < 0] <- 0

# -----------------------------
# distances
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(NOGI_1_coords[, c("lon","lat")]),
  fun = geosphere::distHaversine
)/1000

# -----------------------------
# IBD
# -----------------------------
ibd_df <- data.frame(
  fst = NOGI_1_fst[upper.tri(NOGI_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

ggplot(ibd_df, aes(dist_km, fst)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  theme_classic()

# save
dir.create(file.path(save_dir,"data"), recursive=TRUE, showWarnings=FALSE)
save(NOGI_1_fst, NOGI_1_coords,
     file=file.path(save_dir,"data","NOGI-1.RData"))