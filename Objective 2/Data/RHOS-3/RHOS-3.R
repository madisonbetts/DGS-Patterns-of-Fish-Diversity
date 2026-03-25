# -----------------------------
# RHOS-3 Rhinichthys osculus
# Pilger et al. 2017
# Upper Gila River, New Mexico
#
# Workflow:
# - site coordinates are taken directly from Table S1
# - pairwise values in Table S4 are Slatkin's linearized FST
# - raw FST is back-calculated as linearized_FST / (1 + linearized_FST)
# - the matrix is forced to be symmetric
# - any negative FST values are set to 0
# -----------------------------
library(geosphere)
library(ggplot2)
library(maps)

# -----------------------------
# paths
# -----------------------------
out_rdata <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/RHOS-3/data/RHOS-3.RData"

# -----------------------------
# 1) sites in Table S4 order
# original_site gives the site number from Pilger et al. Table S1 / Fig. 1
# site is the renumbered Objective 2 site ID used in the matrix
# -----------------------------
RHOS_3_coords <- data.frame(
  site = 1:8,
  #original_site = c(16, 15, 14, 13, 12, 11, 9, 2),
  #stream = c(
  #  "Upper East Fork Gila River",
  #  "Black Canyon",
  #  "West Fork Gila River",
  #  "Middle Fork Gila River",
  #  "West Fork Gila River at Heart Bar Game Reserve",
  #  "Little Creek",
  #  "Gila River below confluence with East Fork",
  #  "Blue Creek"
  #),
  #classification = c("Tributary", "Tributary", "Tributary", "Tributary", "Mainstem", "Tributary", "Mainstem", "Tributary"),
  lat = c(33.30011181, 33.184035, 33.22968549, 33.22604743, 33.2141233, 33.19301909, 33.17707829, 32.833278),
  lon = c(-108.1229968, -108.033769, -108.2611784, -108.2412087, -108.2289433, -108.2253339, -108.2074147, -108.827975),
  #elev_m = c(1878, 2059, 1738, 1725, 1711, 1729, 1691, 1422),
  #river_km_from_site1 = c(197.0, 183.0, 170.0, 169.0, 166.0, 160.0, 159.0, 43.0),
  stringsAsFactors = FALSE
)

stopifnot(nrow(RHOS_3_coords) == 8)
stopifnot(identical(RHOS_3_coords$site, 1:8))

# -----------------------------
# 2) pairwise FST matrix
# back-calculated from Slatkin's linearized FST in Table S4
# FST = linearized_FST / (1 + linearized_FST)
# -----------------------------
RHOS_3_fst <- matrix(0, nrow = 8, ncol = 8)
rownames(RHOS_3_fst) <- colnames(RHOS_3_fst) <- as.character(1:8)

RHOS_3_fst[1,2] <- 0.0000000000
RHOS_3_fst[1,3] <- 0.0243902439
RHOS_3_fst[1,4] <- 0.0108803165
RHOS_3_fst[1,5] <- 0.0128331688
RHOS_3_fst[1,6] <- 0.0099009901
RHOS_3_fst[1,7] <- 0.0338164251
RHOS_3_fst[1,8] <- 0.1371872304
RHOS_3_fst[2,3] <- 0.0147783251
RHOS_3_fst[2,4] <- 0.0118577075
RHOS_3_fst[2,5] <- 0.0019960080
RHOS_3_fst[2,6] <- 0.0089197225
RHOS_3_fst[2,7] <- 0.0243902439
RHOS_3_fst[2,8] <- 0.1212653779
RHOS_3_fst[3,4] <- 0.0000000000
RHOS_3_fst[3,5] <- 0.0118577075
RHOS_3_fst[3,6] <- 0.0234375000
RHOS_3_fst[3,7] <- 0.0099009901
RHOS_3_fst[3,8] <- 0.1341991342
RHOS_3_fst[4,5] <- 0.0059642147
RHOS_3_fst[4,6] <- 0.0128331688
RHOS_3_fst[4,7] <- 0.0029910269
RHOS_3_fst[4,8] <- 0.1289198606
RHOS_3_fst[5,6] <- 0.0000000000
RHOS_3_fst[5,7] <- 0.0089197225
RHOS_3_fst[5,8] <- 0.1063449508
RHOS_3_fst[6,7] <- 0.0157480315
RHOS_3_fst[6,8] <- 0.1235758107
RHOS_3_fst[7,8] <- 0.1150442478

# make symmetric and clean
RHOS_3_fst[lower.tri(RHOS_3_fst)] <- t(RHOS_3_fst)[lower.tri(RHOS_3_fst)]
RHOS_3_fst[RHOS_3_fst < 0] <- 0
diag(RHOS_3_fst) <- 0

stopifnot(identical(dim(RHOS_3_fst), c(8L, 8L)))
stopifnot(isTRUE(all.equal(RHOS_3_fst, t(RHOS_3_fst), check.attributes = FALSE)))
stopifnot(identical(as.character(RHOS_3_coords$site), rownames(RHOS_3_fst)))

# -----------------------------
# 3) geographic distances (km)
# -----------------------------
coords <- RHOS_3_coords[, c("lon", "lat")]
geo_dist_km <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
rownames(geo_dist_km) <- colnames(geo_dist_km) <- RHOS_3_coords$site

# -----------------------------
# 4) IBD dataframe
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(RHOS_3_fst)[row(RHOS_3_fst)[upper.tri(RHOS_3_fst)]],
  site2   = colnames(RHOS_3_fst)[col(RHOS_3_fst)[upper.tri(RHOS_3_fst)]],
  fst     = RHOS_3_fst[upper.tri(RHOS_3_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 5) map of sampling sites
# -----------------------------
world_map <- map_data("world")
usa_map   <- subset(world_map, region == "USA")

map_df <- RHOS_3_coords
map_df$label <- as.character(map_df$site)

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
    nudge_y = 0.02,
    size = 3.3
  ) +
  coord_fixed(
    ratio = 1.3,
    xlim = range(map_df$lon) + c(-0.2, 0.2),
    ylim = range(map_df$lat) + c(-0.15, 0.15)
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "RHOS-3 sampling sites"
  )

# -----------------------------
# 6) IBD plot
# -----------------------------
p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Geographic distance (km)",
    y = expression(F[ST]),
    title = "RHOS-3 isolation by distance"
  )

# -----------------------------
# 7) save RData
# -----------------------------
save(
  RHOS_3_fst,
  RHOS_3_coords,
  file = out_rdata
)

# -----------------------------
# 8) print outputs
# -----------------------------
print(RHOS_3_coords)
print(round(RHOS_3_fst, 4))
print(head(ibd_df))
print(p_map)
print(p_ibd)
