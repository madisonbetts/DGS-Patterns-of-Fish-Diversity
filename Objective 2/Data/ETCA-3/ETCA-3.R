############################
## ETCA-3 Rainbow darter
## Microsatellite θST (Table IVc)
## Outputs: ETCA_3_fst, ETCA_3_coords
############################

library(geosphere)
library(ggplot2)

## -------------------------
## 1) site coordinates (Table I)
## -------------------------
## -------------------------
## 1) site coordinates (Table I — full precision)
## -------------------------
ETCA_3_coords <- data.frame(
  site_id = 1:6,
  #site = c(
  #  "Blanchard River",
  #  "Cuyahoga River",
  #  "Chagrin River",
  #  "Grand River",
  #  "Little Miami River",
  #  "Big Darby Creek"
  #),
  lat = c(
    41.036229,
    41.159442,
    41.488884,
    41.742064,
    39.782865,
    39.947505
  ),
  lon = c(
    -83.576722,
    -81.606751,
    -81.396439,
    -81.047707,
    -83.876253,
    -83.236055
  )
)

## -------------------------
## 2) FST matrix (microsats Table IVc)
## order:
## 1 = Blanchard R.
## 2 = Cuyahoga R.
## 3 = Chagrin R.
## 4 = Grand R.
## 5 = Little Miami R.
## 6 = Big Darby Ck.
## -------------------------
ETCA_3_fst <- matrix(0, nrow = 6, ncol = 6)

## lower.tri() fills by column in R:
## (2,1), (3,1), (4,1), (5,1), (6,1),
## (3,2), (4,2), (5,2), (6,2),
## (4,3), (5,3), (6,3),
## (5,4), (6,4),
## (6,5)

ETCA_3_fst[lower.tri(ETCA_3_fst)] <- c(
  0.046, 0.063, 0.029, 0.055, 0.074,
  0.049, 0.031, 0.046, 0.080,
  0.033, 0.088, 0.118,
  0.052, 0.088,
  0.027
)

ETCA_3_fst <- ETCA_3_fst + t(ETCA_3_fst)
diag(ETCA_3_fst) <- 0

rownames(ETCA_3_fst) <- 1:6
colnames(ETCA_3_fst) <- 1:6

ETCA_3_fst[ETCA_3_fst < 0] <- 0
ETCA_3_fst

## -------------------------
## 3) geographic distance (km)
## -------------------------
coords <- as.matrix(ETCA_3_coords[, c("lon", "lat")])

geo_dist_km <- distm(coords, fun = distHaversine) / 1000
rownames(geo_dist_km) <- 1:6
colnames(geo_dist_km) <- 1:6

## -------------------------
## 4) pairwise dataframe for IBD
## -------------------------
ibd_df <- data.frame(
  site1   = rownames(ETCA_3_fst)[row(ETCA_3_fst)[upper.tri(ETCA_3_fst)]],
  site2   = colnames(ETCA_3_fst)[col(ETCA_3_fst)[upper.tri(ETCA_3_fst)]],
  fst     = ETCA_3_fst[upper.tri(ETCA_3_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

## -------------------------
## 5) site map
## -------------------------
ggplot(ETCA_3_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site_id), vjust = -1) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "ETCA-3 Sampling Sites"
  )

## -------------------------
## 6) IBD plot
## -------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "ETCA-3 Isolation by Distance"
  )


## -------------------------
## 7) save objects
## -------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETCA-3/data"

dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

save(
  ETCA_3_fst,
  ETCA_3_coords,
  file = file.path(save_dir, "ETCA_3.RData")
)


