# -----------------------------
# ONCL-3 Lahontan cutthroat trout
# Nielsen and Sage 2002
# Table 5 pairwise FST + deep-searched site anchors
# -----------------------------

library(dplyr)
library(geosphere)
library(ggplot2)
library(maps)

# -----------------------------
# paths
# -----------------------------
setwd("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ONCL-3")

# -----------------------------
# sites (site, lat, lon ONLY)
# best-available named waterbody anchors matched to paper sites
# note: these are feature anchors for named lakes/streams or source drainages,
# not exact collection GPS locations from the study.
# -----------------------------
ONCL_3_sites <- data.frame(
  site = 1:13,
  lat = c(
    40.9170,  # 1 Pilot Peak wild - Utah (Morrison Creek / Pilot Peak drainage)
    38.6382,  # 2 Slinkard Creek
    38.8070,  # 3 East Carson River
    39.5196,  # 4 Macklin Creek
    38.6513,  # 5 Heenan Creek
    41.5153,  # 6 Summit Lake
    39.6171,  # 7 Edwards Creek
    39.4067,  # 8 Independence Lake
    40.0670,  # 9 Pyramid Lake Hatchery (Pyramid Lake anchor)
    40.0675,  # 10 Pilot Peak Hatchery (Pyramid Lake / LNFH stock anchor)
    38.5890,  # 11 Four Mile Canyon Creek
    41.2700,  # 12 Frazer Creek
    41.7320   # 13 West Marys River
  ),
  lon = c(
    -114.0470,
    -119.5077,
    -119.7270,
    -120.6255,
    -119.6541,
    -119.0632,
    -117.6718,
    -120.2537,
    -119.5650,
    -119.5550,
    -119.7020,
    -116.8000,
    -115.5100
  )
)

# -----------------------------
# FST matrix (Table 5; above diagonal)
# -----------------------------
ONCL_3_fst <- matrix(0, 13, 13)
rownames(ONCL_3_fst) <- colnames(ONCL_3_fst) <- 1:13

ONCL_3_fst[1,2]  <- 0.453
ONCL_3_fst[1,3]  <- 0.235
ONCL_3_fst[1,4]  <- 0.207
ONCL_3_fst[1,5]  <- 0.211
ONCL_3_fst[1,6]  <- 0.349
ONCL_3_fst[1,7]  <- 0.450
ONCL_3_fst[1,8]  <- 0.278
ONCL_3_fst[1,9]  <- 0.218
ONCL_3_fst[1,10] <- 0.121
ONCL_3_fst[1,11] <- 0.431
ONCL_3_fst[1,12] <- 0.157
ONCL_3_fst[1,13] <- 0.151

ONCL_3_fst[2,3]  <- 0.339
ONCL_3_fst[2,4]  <- 0.385
ONCL_3_fst[2,5]  <- 0.401
ONCL_3_fst[2,6]  <- 0.470
ONCL_3_fst[2,7]  <- 0.661
ONCL_3_fst[2,8]  <- 0.454
ONCL_3_fst[2,9]  <- 0.458
ONCL_3_fst[2,10] <- 0.319
ONCL_3_fst[2,11] <- 0.624
ONCL_3_fst[2,12] <- 0.374
ONCL_3_fst[2,13] <- 0.326

ONCL_3_fst[3,4]  <- 0.047
ONCL_3_fst[3,5]  <- 0.079
ONCL_3_fst[3,6]  <- 0.183
ONCL_3_fst[3,7]  <- 0.541
ONCL_3_fst[3,8]  <- 0.041
ONCL_3_fst[3,9]  <- 0.196
ONCL_3_fst[3,10] <- 0.089
ONCL_3_fst[3,11] <- 0.335
ONCL_3_fst[3,12] <- 0.106
ONCL_3_fst[3,13] <- 0.105

ONCL_3_fst[4,5]  <- 0.207
ONCL_3_fst[4,6]  <- 0.352
ONCL_3_fst[4,7]  <- 0.423
ONCL_3_fst[4,8]  <- 0.318
ONCL_3_fst[4,9]  <- 0.233
ONCL_3_fst[4,10] <- 0.246
ONCL_3_fst[4,11] <- 0.451
ONCL_3_fst[4,12] <- 0.089
ONCL_3_fst[4,13] <- 0.111

ONCL_3_fst[5,6]  <- 0.243
ONCL_3_fst[5,7]  <- 0.488
ONCL_3_fst[5,8]  <- 0.216
ONCL_3_fst[5,9]  <- 0.133
ONCL_3_fst[5,10] <- 0.144
ONCL_3_fst[5,11] <- 0.467
ONCL_3_fst[5,12] <- 0.160
ONCL_3_fst[5,13] <- 0.135

ONCL_3_fst[6,7]  <- 0.550
ONCL_3_fst[6,8]  <- 0.343
ONCL_3_fst[6,9]  <- 0.259
ONCL_3_fst[6,10] <- 0.310
ONCL_3_fst[6,11] <- 0.586
ONCL_3_fst[6,12] <- 0.336
ONCL_3_fst[6,13] <- 0.260

ONCL_3_fst[7,8]  <- 0.494
ONCL_3_fst[7,9]  <- 0.575
ONCL_3_fst[7,10] <- 0.402
ONCL_3_fst[7,11] <- 0.695
ONCL_3_fst[7,12] <- 0.522
ONCL_3_fst[7,13] <- 0.561

ONCL_3_fst[8,9]  <- 0.022
ONCL_3_fst[8,10] <- 0.276
ONCL_3_fst[8,11] <- 0.498
ONCL_3_fst[8,12] <- 0.238
ONCL_3_fst[8,13] <- 0.184

ONCL_3_fst[9,10] <- 0.186
ONCL_3_fst[9,11] <- 0.484
ONCL_3_fst[9,12] <- 0.263
ONCL_3_fst[9,13] <- 0.170

ONCL_3_fst[10,11] <- 0.461
ONCL_3_fst[10,12] <- 0.051
ONCL_3_fst[10,13] <- 0.070

ONCL_3_fst[11,12] <- 0.426
ONCL_3_fst[11,13] <- 0.415

ONCL_3_fst[12,13] <- 0.115

# mirror + clean
ONCL_3_fst[lower.tri(ONCL_3_fst)] <- t(ONCL_3_fst)[lower.tri(ONCL_3_fst)]
ONCL_3_fst[ONCL_3_fst < 0] <- 0
diag(ONCL_3_fst) <- 0

# -----------------------------
# pairwise geographic distances for IBD
# -----------------------------
coords <- ONCL_3_sites %>% select(lon, lat)
geo_dist_km <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
rownames(geo_dist_km) <- colnames(geo_dist_km) <- ONCL_3_sites$site

ibd_df <- data.frame(
  site1   = rownames(ONCL_3_fst)[row(ONCL_3_fst)[upper.tri(ONCL_3_fst)]],
  site2   = colnames(ONCL_3_fst)[col(ONCL_3_fst)[upper.tri(ONCL_3_fst)]],
  fst     = ONCL_3_fst[upper.tri(ONCL_3_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# map of sampling sites
# -----------------------------
world_map <- map_data("world")
usa_map   <- subset(world_map, region %in% c("USA", "Canada", "Mexico"))

map_df <- ONCL_3_sites %>%
  mutate(label = as.character(site))

x_pad <- diff(range(map_df$lon)) * 0.15
y_pad <- diff(range(map_df$lat)) * 0.15

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
    nudge_y = 0.12,
    size = 3.5
  ) +
  coord_fixed(
    ratio = 1.2,
    xlim = range(map_df$lon) + c(-x_pad, x_pad),
    ylim = range(map_df$lat) + c(-y_pad, y_pad)
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "ONCL-3 sampling sites"
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
    title = "ONCL-3 isolation by distance"
  )

# -----------------------------
# save RData
# -----------------------------
save(
  ONCL_3_fst,
  ONCL_3_sites,
  file = "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ONCL-3/data/ONCL-3.RData"
)

# -----------------------------
# print outputs
# -----------------------------
print(ONCL_3_sites)
print(round(ONCL_3_fst, 4))
print(head(ibd_df))
print(p_map)
print(p_ibd)
