# -----------------------------
# COGR-1 Shoshone sculpin
# Campbell et al. 2023
# Table 3 pairwise FST + Figure 1 site coordinates
# -----------------------------

library(dplyr)
library(geosphere)
library(ggplot2)
library(maps)

# -----------------------------
# paths
# -----------------------------
setwd("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/COGR-1")

# -----------------------------
# sites
# site coordinates were georeferenced from Fig. 1 and aligned
# to named springs/creeks in the Hagerman Valley, Idaho.
# final object keeps only site, lat, lon per workflow.
# -----------------------------
COGR_1_sites <- data.frame(
  site = 1:20,
  lat = c(
    42.8240, 42.7930, 42.7810, 42.7470, 42.7190,
    42.7010, 42.6630, 42.5930, 42.5870, 42.5810,
    42.5750, 42.5470, 42.5180, 42.5000, 42.4570,
    42.4460, 42.4470, 42.4290, 42.3920, 42.3410
  ),
  lon = c(
    -114.8750, -114.8500, -114.8440, -114.8200, -114.8200,
    -114.8020, -114.7740, -114.8070, -114.7950, -114.7580,
    -114.7520, -114.7270, -114.7260, -114.7260, -114.7020,
    -114.6970, -114.6690, -114.6940, -114.6900, -114.6600
  )
)

# -----------------------------
# FST matrix (Table 3)
# -----------------------------
COGR_1_fst <- matrix(0, 20, 20)
rownames(COGR_1_fst) <- colnames(COGR_1_fst) <- 1:20

COGR_1_fst[2,1] <- 0.01

COGR_1_fst[3,1] <- 0.06; COGR_1_fst[3,2] <- 0.04

COGR_1_fst[4,1] <- 0.21; COGR_1_fst[4,2] <- 0.24; COGR_1_fst[4,3] <- 0.23

COGR_1_fst[5,1] <- 0.20; COGR_1_fst[5,2] <- 0.23; COGR_1_fst[5,3] <- 0.22; COGR_1_fst[5,4] <- 0.02

COGR_1_fst[6,1] <- 0.26; COGR_1_fst[6,2] <- 0.27; COGR_1_fst[6,3] <- 0.27; COGR_1_fst[6,4] <- 0.17; COGR_1_fst[6,5] <- 0.15

COGR_1_fst[7,1] <- 0.25; COGR_1_fst[7,2] <- 0.27; COGR_1_fst[7,3] <- 0.27; COGR_1_fst[7,4] <- 0.20; COGR_1_fst[7,5] <- 0.17; COGR_1_fst[7,6] <- 0.01

COGR_1_fst[8,1] <- 0.12; COGR_1_fst[8,2] <- 0.14; COGR_1_fst[8,3] <- 0.19; COGR_1_fst[8,4] <- 0.26; COGR_1_fst[8,5] <- 0.22; COGR_1_fst[8,6] <- 0.21; COGR_1_fst[8,7] <- 0.21

COGR_1_fst[9,1] <- 0.26; COGR_1_fst[9,2] <- 0.30; COGR_1_fst[9,3] <- 0.29; COGR_1_fst[9,4] <- 0.11; COGR_1_fst[9,5] <- 0.12; COGR_1_fst[9,6] <- 0.24; COGR_1_fst[9,7] <- 0.28; COGR_1_fst[9,8] <- 0.28

COGR_1_fst[10,1] <- 0.23; COGR_1_fst[10,2] <- 0.25; COGR_1_fst[10,3] <- 0.25; COGR_1_fst[10,4] <- 0.23; COGR_1_fst[10,5] <- 0.19; COGR_1_fst[10,6] <- 0.12; COGR_1_fst[10,7] <- 0.14; COGR_1_fst[10,8] <- 0.12; COGR_1_fst[10,9] <- 0.26

COGR_1_fst[11,1] <- 0.29; COGR_1_fst[11,2] <- 0.32; COGR_1_fst[11,3] <- 0.31; COGR_1_fst[11,4] <- 0.11; COGR_1_fst[11,5] <- 0.13; COGR_1_fst[11,6] <- 0.26; COGR_1_fst[11,7] <- 0.29; COGR_1_fst[11,8] <- 0.30; COGR_1_fst[11,9] <- 0.00; COGR_1_fst[11,10] <- 0.28

COGR_1_fst[12,1] <- 0.30; COGR_1_fst[12,2] <- 0.34; COGR_1_fst[12,3] <- 0.33; COGR_1_fst[12,4] <- 0.12; COGR_1_fst[12,5] <- 0.14; COGR_1_fst[12,6] <- 0.27; COGR_1_fst[12,7] <- 0.31; COGR_1_fst[12,8] <- 0.32; COGR_1_fst[12,9] <- 0.00; COGR_1_fst[12,10] <- 0.30; COGR_1_fst[12,11] <- 0.00

COGR_1_fst[13,1] <- 0.30; COGR_1_fst[13,2] <- 0.33; COGR_1_fst[13,3] <- 0.33; COGR_1_fst[13,4] <- 0.12; COGR_1_fst[13,5] <- 0.15; COGR_1_fst[13,6] <- 0.28; COGR_1_fst[13,7] <- 0.32; COGR_1_fst[13,8] <- 0.31; COGR_1_fst[13,9] <- 0.00; COGR_1_fst[13,10] <- 0.30; COGR_1_fst[13,11] <- 0.00; COGR_1_fst[13,12] <- 0.00

COGR_1_fst[14,1] <- 0.30; COGR_1_fst[14,2] <- 0.33; COGR_1_fst[14,3] <- 0.33; COGR_1_fst[14,4] <- 0.13; COGR_1_fst[14,5] <- 0.15; COGR_1_fst[14,6] <- 0.29; COGR_1_fst[14,7] <- 0.32; COGR_1_fst[14,8] <- 0.32; COGR_1_fst[14,9] <- 0.01; COGR_1_fst[14,10] <- 0.30; COGR_1_fst[14,11] <- 0.00; COGR_1_fst[14,12] <- 0.00; COGR_1_fst[14,13] <- 0.00

COGR_1_fst[15,1] <- 0.35; COGR_1_fst[15,2] <- 0.38; COGR_1_fst[15,3] <- 0.40; COGR_1_fst[15,4] <- 0.18; COGR_1_fst[15,5] <- 0.21; COGR_1_fst[15,6] <- 0.39; COGR_1_fst[15,7] <- 0.41; COGR_1_fst[15,8] <- 0.37; COGR_1_fst[15,9] <- 0.09; COGR_1_fst[15,10] <- 0.37; COGR_1_fst[15,11] <- 0.10; COGR_1_fst[15,12] <- 0.09; COGR_1_fst[15,13] <- 0.06; COGR_1_fst[15,14] <- 0.07

COGR_1_fst[16,1] <- 0.32; COGR_1_fst[16,2] <- 0.36; COGR_1_fst[16,3] <- 0.37; COGR_1_fst[16,4] <- 0.16; COGR_1_fst[16,5] <- 0.18; COGR_1_fst[16,6] <- 0.36; COGR_1_fst[16,7] <- 0.40; COGR_1_fst[16,8] <- 0.35; COGR_1_fst[16,9] <- 0.07; COGR_1_fst[16,10] <- 0.35; COGR_1_fst[16,11] <- 0.08; COGR_1_fst[16,12] <- 0.08; COGR_1_fst[16,13] <- 0.05; COGR_1_fst[16,14] <- 0.06; COGR_1_fst[16,15] <- 0.00

COGR_1_fst[17,1] <- 0.37; COGR_1_fst[17,2] <- 0.42; COGR_1_fst[17,3] <- 0.43; COGR_1_fst[17,4] <- 0.26; COGR_1_fst[17,5] <- 0.27; COGR_1_fst[17,6] <- 0.44; COGR_1_fst[17,7] <- 0.46; COGR_1_fst[17,8] <- 0.41; COGR_1_fst[17,9] <- 0.17; COGR_1_fst[17,10] <- 0.42; COGR_1_fst[17,11] <- 0.18; COGR_1_fst[17,12] <- 0.17; COGR_1_fst[17,13] <- 0.16; COGR_1_fst[17,14] <- 0.17; COGR_1_fst[17,15] <- 0.13; COGR_1_fst[17,16] <- 0.09

COGR_1_fst[18,1] <- 0.32; COGR_1_fst[18,2] <- 0.35; COGR_1_fst[18,3] <- 0.36; COGR_1_fst[18,4] <- 0.16; COGR_1_fst[18,5] <- 0.19; COGR_1_fst[18,6] <- 0.35; COGR_1_fst[18,7] <- 0.39; COGR_1_fst[18,8] <- 0.35; COGR_1_fst[18,9] <- 0.07; COGR_1_fst[18,10] <- 0.34; COGR_1_fst[18,11] <- 0.09; COGR_1_fst[18,12] <- 0.08; COGR_1_fst[18,13] <- 0.06; COGR_1_fst[18,14] <- 0.07; COGR_1_fst[18,15] <- 0.00; COGR_1_fst[18,16] <- 0.00; COGR_1_fst[18,17] <- 0.10

COGR_1_fst[19,1] <- 0.51; COGR_1_fst[19,2] <- 0.52; COGR_1_fst[19,3] <- 0.53; COGR_1_fst[19,4] <- 0.30; COGR_1_fst[19,5] <- 0.33; COGR_1_fst[19,6] <- 0.51; COGR_1_fst[19,7] <- 0.51; COGR_1_fst[19,8] <- 0.51; COGR_1_fst[19,9] <- 0.20; COGR_1_fst[19,10] <- 0.51; COGR_1_fst[19,11] <- 0.20; COGR_1_fst[19,12] <- 0.18; COGR_1_fst[19,13] <- 0.17; COGR_1_fst[19,14] <- 0.18; COGR_1_fst[19,15] <- 0.15; COGR_1_fst[19,16] <- 0.14; COGR_1_fst[19,17] <- 0.25; COGR_1_fst[19,18] <- 0.13

COGR_1_fst[20,1] <- 0.60; COGR_1_fst[20,2] <- 0.59; COGR_1_fst[20,3] <- 0.62; COGR_1_fst[20,4] <- 0.38; COGR_1_fst[20,5] <- 0.39; COGR_1_fst[20,6] <- 0.60; COGR_1_fst[20,7] <- 0.55; COGR_1_fst[20,8] <- 0.59; COGR_1_fst[20,9] <- 0.29; COGR_1_fst[20,10] <- 0.59; COGR_1_fst[20,11] <- 0.31; COGR_1_fst[20,12] <- 0.29; COGR_1_fst[20,13] <- 0.26; COGR_1_fst[20,14] <- 0.27; COGR_1_fst[20,15] <- 0.24; COGR_1_fst[20,16] <- 0.20; COGR_1_fst[20,17] <- 0.36; COGR_1_fst[20,18] <- 0.19; COGR_1_fst[20,19] <- 0.15

COGR_1_fst[upper.tri(COGR_1_fst)] <- t(COGR_1_fst)[upper.tri(COGR_1_fst)]
COGR_1_fst[COGR_1_fst < 0] <- 0
diag(COGR_1_fst) <- 0

# -----------------------------
# pairwise geographic distances for IBD
# -----------------------------
coords <- COGR_1_sites %>% select(lon, lat)

geo_dist_km <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
rownames(geo_dist_km) <- colnames(geo_dist_km) <- COGR_1_sites$site

ibd_df <- data.frame(
  site1   = rownames(COGR_1_fst)[row(COGR_1_fst)[upper.tri(COGR_1_fst)]],
  site2   = colnames(COGR_1_fst)[col(COGR_1_fst)[upper.tri(COGR_1_fst)]],
  fst     = COGR_1_fst[upper.tri(COGR_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# map of sampling sites
# -----------------------------
world_map <- map_data("world")
usa_map   <- subset(world_map, region %in% c("USA", "Canada", "Mexico"))

map_df <- COGR_1_sites %>%
  mutate(label = as.character(site))

p_map <- ggplot() +
  geom_polygon(
    data = usa_map,
    aes(x = long, y = lat, group = group),
    fill = "grey92",
    color = "grey55",
    linewidth = 0.2
  ) +
  geom_point(
    data = map_df,
    aes(x = lon, y = lat),
    color = "red3",
    size = 3
  ) +
  geom_text(
    data = map_df,
    aes(x = lon, y = lat, label = label),
    nudge_y = 0.01,
    size = 3.5
  ) +
  coord_fixed(
    ratio = 1.3,
    xlim = range(map_df$lon) + c(-0.03, 0.03),
    ylim = range(map_df$lat) + c(-0.03, 0.03)
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "COGR-1 sampling sites"
  )

# -----------------------------
# IBD plot
# -----------------------------
p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Geographic distance (km)",
    y = expression(F[ST]),
    title = "COGR-1 isolation by distance"
  )

# -----------------------------
# save RData
# -----------------------------
save(
  COGR_1_fst,
  COGR_1_sites,
  file = "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/COGR-1/data/COGR-1.RData"
)

# -----------------------------
# print outputs
# -----------------------------
print(COGR_1_sites)
print(round(COGR_1_fst, 4))
print(head(ibd_df))
print(p_map)
print(p_ibd)
