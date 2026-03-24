# -----------------------------
# ONMY-3 California golden trout
# Cordes et al. 2006
# Table 5 pairwise FST + Figure 1 sites
# -----------------------------

library(dplyr)
library(geosphere)
library(ggplot2)
library(maps)

# -----------------------------
# paths
# -----------------------------
setwd("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ONMY-3")

# -----------------------------
# sites (site, lat, lon ONLY)
# site order follows Table 5:
# VC, GC, LWM, JL, CSL, CL2, CL4, MJC, MBS, USS, SFK
# named-feature anchors were used where available;
# stream localities without exact sample GPS were placed
# at the named feature / best-available confluence anchor.
# -----------------------------
ONMY_3_sites <- data.frame(
  site = 1:11,
  lat = c(
    36.3655,    # 1 VC  Volcanic Creek
    36.36605,   # 2 GC  Groundhog Meadows / Groundhog Creek anchor
    36.3702162, # 3 LWM Little Whitney Meadow
    36.4330905, # 4 JL  Johnson Lake
    36.4573567, # 5 CSL Chicken Spring Lake
    36.4893625, # 6 CL2 Cottonwood Lake Number Two
    36.49898,   # 7 CL4 Cottonwood Lake Number Four
    36.3738273, # 8 MJC Middle Johnson Creek / Johnson Creek anchor
    36.4116,    # 9 MBS Mouth of Barigan Stringer / Barigan Stringer anchor
    36.434381,  # 10 USS Upper Stokes Stringer / Stokes Stringer anchor
    36.3480556  # 11 SFK South Fork Kern River above Ramshaw barrier
  ),
  lon = c(
    -118.3562,    # 1 VC
    -118.30787,   # 2 GC
    -118.346199,  # 3 LWM
    -118.3371076, # 4 JL
    -118.2272578, # 5 CSL
    -118.2155585, # 6 CL2
    -118.22969,   # 7 CL4
    -118.3450879, # 8 MJC
    -118.2776,    # 9 MBS
    -118.2611984, # 10 USS
    -118.24795    # 11 SFK
  )
)

# -----------------------------
# FST matrix (Table 5)
# order: VC, GC, LWM, JL, CSL, CL2, CL4, MJC, MBS, USS, SFK
# negative values set to 0
# -----------------------------
ONMY_3_fst <- matrix(0, 11, 11)
rownames(ONMY_3_fst) <- colnames(ONMY_3_fst) <- 1:11

ONMY_3_fst[2,1]  <- 0.05
ONMY_3_fst[3,1]  <- 0.01; ONMY_3_fst[3,2]  <- 0.00
ONMY_3_fst[4,1]  <- 0.20; ONMY_3_fst[4,2]  <- 0.15; ONMY_3_fst[4,3]  <- 0.15
ONMY_3_fst[5,1]  <- 0.36; ONMY_3_fst[5,2]  <- 0.37; ONMY_3_fst[5,3]  <- 0.38; ONMY_3_fst[5,4]  <- 0.28
ONMY_3_fst[6,1]  <- 0.18; ONMY_3_fst[6,2]  <- 0.12; ONMY_3_fst[6,3]  <- 0.09; ONMY_3_fst[6,4]  <- 0.01; ONMY_3_fst[6,5]  <- 0.24
ONMY_3_fst[7,1]  <- 0.19; ONMY_3_fst[7,2]  <- 0.12; ONMY_3_fst[7,3]  <- 0.09; ONMY_3_fst[7,4]  <- 0.00; ONMY_3_fst[7,5]  <- 0.25; ONMY_3_fst[7,6]  <- 0.01
ONMY_3_fst[8,1]  <- 0.08; ONMY_3_fst[8,2]  <- 0.02; ONMY_3_fst[8,3]  <- -0.02; ONMY_3_fst[8,4]  <- 0.20; ONMY_3_fst[8,5]  <- 0.41; ONMY_3_fst[8,6]  <- 0.19; ONMY_3_fst[8,7]  <- 0.20
ONMY_3_fst[9,1]  <- 0.08; ONMY_3_fst[9,2]  <- -0.01; ONMY_3_fst[9,3]  <- -0.04; ONMY_3_fst[9,4]  <- 0.12; ONMY_3_fst[9,5]  <- 0.32; ONMY_3_fst[9,6]  <- 0.11; ONMY_3_fst[9,7]  <- 0.13; ONMY_3_fst[9,8]  <- 0.05
ONMY_3_fst[10,1] <- 0.06; ONMY_3_fst[10,2] <- -0.05; ONMY_3_fst[10,3] <- -0.07; ONMY_3_fst[10,4] <- 0.14; ONMY_3_fst[10,5] <- 0.33; ONMY_3_fst[10,6] <- 0.12; ONMY_3_fst[10,7] <- 0.15; ONMY_3_fst[10,8] <- 0.03; ONMY_3_fst[10,9] <- 0.02
ONMY_3_fst[11,1] <- 0.04; ONMY_3_fst[11,2] <- 0.01; ONMY_3_fst[11,3] <- -0.02; ONMY_3_fst[11,4] <- 0.05; ONMY_3_fst[11,5] <- 0.27; ONMY_3_fst[11,6] <- 0.02; ONMY_3_fst[11,7] <- 0.04; ONMY_3_fst[11,8] <- 0.10; ONMY_3_fst[11,9] <- 0.02; ONMY_3_fst[11,10] <- 0.03

ONMY_3_fst[upper.tri(ONMY_3_fst)] <- t(ONMY_3_fst)[upper.tri(ONMY_3_fst)]
ONMY_3_fst[ONMY_3_fst < 0] <- 0
diag(ONMY_3_fst) <- 0

# -----------------------------
# pairwise geographic distances for IBD
# -----------------------------
coords <- ONMY_3_sites %>% select(lon, lat)
geo_dist_km <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
rownames(geo_dist_km) <- colnames(geo_dist_km) <- ONMY_3_sites$site

ibd_df <- data.frame(
  site1   = rownames(ONMY_3_fst)[row(ONMY_3_fst)[upper.tri(ONMY_3_fst)]],
  site2   = colnames(ONMY_3_fst)[col(ONMY_3_fst)[upper.tri(ONMY_3_fst)]],
  fst     = ONMY_3_fst[upper.tri(ONMY_3_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# map of sampling sites
# -----------------------------
world_map <- map_data("world")
usa_map   <- subset(world_map, region %in% c("USA", "Mexico"))

map_df <- ONMY_3_sites %>%
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
    xlim = range(map_df$lon) + c(-0.05, 0.05),
    ylim = range(map_df$lat) + c(-0.05, 0.05)
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "ONMY-3 sampling sites"
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
    title = "ONMY-3 isolation by distance"
  )

# -----------------------------
# save RData
# -----------------------------
save(
  ONMY_3_fst,
  ONMY_3_sites,
  file = "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ONMY-3/data/ONMY-3.RData"
)

# -----------------------------
# print outputs
# -----------------------------
print(ONMY_3_sites)
print(round(ONMY_3_fst, 4))
print(head(ibd_df))
print(p_map)
print(p_ibd)
