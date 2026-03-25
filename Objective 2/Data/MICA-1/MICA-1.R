############################
## MICA-1
## Shoal Bass (Micropterus cataractae)
## Taylor et al. 2018 North American Journal of Fisheries Management 38:549-564
## DOI: 10.1002/nafm.10048
##
## Inputs taken from paper:
## - Table 1: 13 named sampling sites
## - Figure 1 / Figure 4: site arrangement in the ACF basin
## - Table 7: pairwise FST among 13 sites
##
## Notes on coordinates:
## The paper does not provide site lat/lon directly, so coordinates below are
## inferred from the named reaches using the best matching public gage / locality
## anchors for each site.
##
## Workflow:
## 1) build coordinate dataframe in matrix order
## 2) build symmetrical FST matrix
## 3) set negative FST to 0
## 4) make site map
## 5) compute Euclidean geographic distances
## 6) make IBD plot
## 7) save only MICA_1_fst and MICA_1_coords to .RData
############################

# -----------------------------
# libraries
# -----------------------------
library(ggplot2)
library(maps)
library(geosphere)

# -----------------------------
# paths
# -----------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/MICA-1"
data_dir <- file.path(base_dir, "data")
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) coordinate dataframe
# matrix order follows Table 7:
# 1 Upper Chattahoochee River
# 2 Chestatee River
# 3 Big Creek
# 4 Chattahoochee River below Morgan Falls Dam
# 5 Little Uchee Creek
# 6 Upper Flint River above Fall Line
# 7 Big Lazer Creek
# 8 Upper Flint River below Fall Line
# 9 Lower Flint River below Lake Blackshear
# 10 Lower Flint River below Lake Worth
# 11 Lower Flint River
# 12 Ichawaynochaway Creek
# 13 Chipola River
# -----------------------------
MICA_1_coords <- data.frame(
  site_id = 1:13,
  lat = c(
    34.5407222,  # 1 Chattahoochee River near Cornelia, GA
    34.5282778,  # 2 Chestatee River near Dahlonega, GA
    34.0064889,  # 3 Big Creek at Riverside Rd near Roswell, GA
    33.9683000,  # 4 Chattahoochee River below Morgan Falls Dam, GA
    32.3518130,  # 5 Little Uchee Creek, AL
    32.8385000,  # 6 Flint River near Thomaston, GA
    32.7424167,  # 7 Big Lazer Creek at GA 41 near Talbotton, GA
    32.7221944,  # 8 Flint River at US 19 near Carsonville, GA
    31.9643389,  # 9 Flint River at Cobb / below Lake Blackshear, GA
    31.6005310,  # 10 Flint River below Lake Worth Dam, GA
    31.3069444,  # 11 Flint River at Newton, GA
    31.3827778,  # 12 Ichawaynochaway Creek at Milford, GA
    30.5340832   # 13 Chipola River near Altha, FL
  ),
  lon = c(
    -83.6227750,
    -83.9402222,
    -84.3499258,
    -84.3828000,
    -85.0724335,
    -84.4242500,
    -84.5555000,
    -84.2324167,
    -83.9326762,
    -84.1396380,
    -84.3388889,
    -84.5463889,
    -85.1651996
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) pairwise FST matrix from Table 7
# values below diagonal in paper
# negatives retained initially, then set to 0 below
# -----------------------------
n_sites <- nrow(MICA_1_coords)

MICA_1_fst <- matrix(NA_real_, n_sites, n_sites)
diag(MICA_1_fst) <- 0

# fill lower triangle from Table 7
MICA_1_fst[2,1]   <-  0.03384
MICA_1_fst[3,1]   <-  0.02320
MICA_1_fst[3,2]   <-  0.06987
MICA_1_fst[4,1]   <-  0.02192
MICA_1_fst[4,2]   <-  0.01638
MICA_1_fst[4,3]   <-  0.03870
MICA_1_fst[5,1]   <-  0.30018
MICA_1_fst[5,2]   <-  0.31559
MICA_1_fst[5,3]   <-  0.35637
MICA_1_fst[5,4]   <-  0.29983
MICA_1_fst[6,1]   <-  0.03198
MICA_1_fst[6,2]   <-  0.05069
MICA_1_fst[6,3]   <-  0.03631
MICA_1_fst[6,4]   <-  0.01246
MICA_1_fst[6,5]   <-  0.30621
MICA_1_fst[7,1]   <-  0.04290
MICA_1_fst[7,2]   <-  0.03967
MICA_1_fst[7,3]   <-  0.06998
MICA_1_fst[7,4]   <-  0.02994
MICA_1_fst[7,5]   <-  0.28072
MICA_1_fst[7,6]   <-  0.04112
MICA_1_fst[8,1]   <-  0.03024
MICA_1_fst[8,2]   <-  0.02511
MICA_1_fst[8,3]   <-  0.04092
MICA_1_fst[8,4]   <- -0.00627
MICA_1_fst[8,5]   <-  0.28492
MICA_1_fst[8,6]   <-  0.00278
MICA_1_fst[8,7]   <-  0.02918
MICA_1_fst[9,1]   <-  0.04094
MICA_1_fst[9,2]   <-  0.05044
MICA_1_fst[9,3]   <-  0.03205
MICA_1_fst[9,4]   <-  0.00798
MICA_1_fst[9,5]   <-  0.26841
MICA_1_fst[9,6]   <-  0.02312
MICA_1_fst[9,7]   <-  0.05114
MICA_1_fst[9,8]   <-  0.00447
MICA_1_fst[10,1]  <-  0.03925
MICA_1_fst[10,2]  <-  0.06033
MICA_1_fst[10,3]  <-  0.02136
MICA_1_fst[10,4]  <-  0.01147
MICA_1_fst[10,5]  <-  0.31081
MICA_1_fst[10,6]  <-  0.02619
MICA_1_fst[10,7]  <-  0.06781
MICA_1_fst[10,8]  <-  0.01653
MICA_1_fst[10,9]  <-  0.00359
MICA_1_fst[11,1]  <-  0.05256
MICA_1_fst[11,2]  <-  0.06244
MICA_1_fst[11,3]  <-  0.06678
MICA_1_fst[11,4]  <-  0.01213
MICA_1_fst[11,5]  <-  0.27974
MICA_1_fst[11,6]  <-  0.01536
MICA_1_fst[11,7]  <-  0.04756
MICA_1_fst[11,8]  <-  0.01207
MICA_1_fst[11,9]  <-  0.01778
MICA_1_fst[11,10] <-  0.01960
MICA_1_fst[12,1]  <-  0.03618
MICA_1_fst[12,2]  <-  0.04757
MICA_1_fst[12,3]  <-  0.04677
MICA_1_fst[12,4]  <-  0.00143
MICA_1_fst[12,5]  <-  0.28500
MICA_1_fst[12,6]  <-  0.01513
MICA_1_fst[12,7]  <-  0.04821
MICA_1_fst[12,8]  <-  0.00488
MICA_1_fst[12,9]  <-  0.01072
MICA_1_fst[12,10] <-  0.00814
MICA_1_fst[12,11] <-  0.00672
MICA_1_fst[13,1]  <-  0.12232
MICA_1_fst[13,2]  <-  0.12416
MICA_1_fst[13,3]  <-  0.13356
MICA_1_fst[13,4]  <-  0.10779
MICA_1_fst[13,5]  <-  0.19874
MICA_1_fst[13,6]  <-  0.12395
MICA_1_fst[13,7]  <-  0.07330
MICA_1_fst[13,8]  <-  0.08866
MICA_1_fst[13,9]  <-  0.13851
MICA_1_fst[13,10] <-  0.14254
MICA_1_fst[13,11] <-  0.12035
MICA_1_fst[13,12] <-  0.13370

# make symmetrical
MICA_1_fst[upper.tri(MICA_1_fst)] <- t(MICA_1_fst)[upper.tri(MICA_1_fst)]

# set negatives to 0 per workflow
MICA_1_fst[MICA_1_fst < 0] <- 0

# use numeric site IDs as row/col names
rownames(MICA_1_fst) <- MICA_1_coords$site_id
colnames(MICA_1_fst) <- MICA_1_coords$site_id

# -----------------------------
# 3) site map
# -----------------------------
states_map <- map_data("state")

ggplot() +
  geom_polygon(
    data = states_map,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = MICA_1_coords,
    aes(x = lon, y = lat),
    size = 3
  ) +
  geom_text(
    data = MICA_1_coords,
    aes(x = lon, y = lat, label = site_id),
    nudge_y = 0.12,
    size = 3
  ) +
  coord_fixed(
    1.2,
    xlim = c(-86.5, -83.0),
    ylim = c(30.0, 35.0)
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "MICA-1 sampling sites"
  )


#---------------------------
# 4) geographic distance matrix
# -----------------------------
mica_1_geo_dist_m <- geosphere::distm(
  x = MICA_1_coords[, c("lon", "lat")],
  fun = geosphere::distHaversine
)

mica_1_geo_dist_km <- mica_1_geo_dist_m / 1000
rownames(mica_1_geo_dist_km) <- MICA_1_coords$site_id
colnames(mica_1_geo_dist_km) <- MICA_1_coords$site_id

# -----------------------------
# 5) pairwise dataframe for IBD
# -----------------------------
mica_1_ibd_df <- data.frame(
  site1   = rownames(MICA_1_fst)[row(MICA_1_fst)[upper.tri(MICA_1_fst)]],
  site2   = colnames(MICA_1_fst)[col(MICA_1_fst)[upper.tri(MICA_1_fst)]],
  fst     = MICA_1_fst[upper.tri(MICA_1_fst)],
  dist_km = mica_1_geo_dist_km[upper.tri(mica_1_geo_dist_km)]
)

# -----------------------------
# 6) IBD plot
# -----------------------------
ggplot(mica_1_ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 12) +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "MICA-1 isolation by distance"
  )


# -----------------------------
# 7) save outputs
# save only fst and coords as requested
# -----------------------------
save(
  MICA_1_fst,
  MICA_1_coords,
  file = file.path(data_dir, "MICA-1.RData")
)
