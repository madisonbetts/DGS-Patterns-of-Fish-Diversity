# -----------------------------
# PACL-1 Pantosteus clarkii
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
out_rdata <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/PACL-1/data/PACL-1.RData"

# -----------------------------
# 1) sites in Table S4 order
# original_site gives the site number from Pilger et al. Table S1 / Fig. 1
# site is the renumbered Objective 2 site ID used in the matrix
# -----------------------------
PACL_1_coords <- data.frame(
  site = 1:12,
  #original_site = c(16, 15, 14, 13, 12, 11, 10, 9, 6, 5, 4, 3),
  #stream = c(
  #  "Upper East Fork Gila River",
  #  "Black Canyon",
  #  "West Fork Gila River",
  #  "Middle Fork Gila River",
  #  "West Fork Gila River at Heart Bar Game Reserve",
  #  "Little Creek",
  #  "Lower East Fork Gila River",
  #  "Gila River below confluence with East Fork",
  #  "Gila River at Gila Farm",
  #  "Gila River at Riverside",
  #  "Gila River at Bird Area",
  #  "Gila River at confluence with Ash Canyon"
#  ),
#  classification = c("Tributary", "Tributary", "Tributary", "Tributary", "Mainstem", "Tributary", "Tributary", "Mainstem", "Mainstem", "Mainstem", "Mainstem", "Mainstem"),
  lat = c(33.30011181, 33.184035, 33.22968549, 33.22604743, 33.2141233, 33.19301909, 33.177019, 33.17707829, 33.03671285, 32.9388458, 32.84751, 32.713959),
  lon = c(-108.1229968, -108.033769, -108.2611784, -108.2412087, -108.2289433, -108.2253339, -108.199432, -108.2074147, -108.5328707, -108.6056954, -108.593688, -108.71031),
 # elev_m = c(1878, 2059, 1738, 1725, 1711, 1729, 1694, 1691, 1412, 1360, 1331, 1240),
#  river_km_from_site1 = c(197.0, 183.0, 170.0, 169.0, 166.0, 160.0, 155.0, 159.0, 90.0, 75.0, 60.0, 30.0),
  stringsAsFactors = FALSE
)

stopifnot(nrow(PACL_1_coords) == 12)
stopifnot(identical(PACL_1_coords$site, 1:12))

# -----------------------------
# 2) pairwise FST matrix
# back-calculated from Slatkin's linearized FST in Table S4
# FST = linearized_FST / (1 + linearized_FST)
# -----------------------------
PACL_1_fst <- matrix(0, nrow = 12, ncol = 12)
rownames(PACL_1_fst) <- colnames(PACL_1_fst) <- as.character(1:12)

PACL_1_fst[1,2] <- 0.0049751244
PACL_1_fst[1,3] <- 0.0029910269
PACL_1_fst[1,4] <- 0.0019960080
PACL_1_fst[1,5] <- 0.0009990010
PACL_1_fst[1,6] <- 0.0176817289
PACL_1_fst[1,7] <- 0.0029910269
PACL_1_fst[1,8] <- 0.0029910269
PACL_1_fst[1,9] <- 0.0215264188
PACL_1_fst[1,10] <- 0.0186457311
PACL_1_fst[1,11] <- 0.0118577075
PACL_1_fst[1,12] <- 0.0196078431
PACL_1_fst[2,3] <- 0.0079365079
PACL_1_fst[2,4] <- 0.0039840637
PACL_1_fst[2,5] <- 0.0049751244
PACL_1_fst[2,6] <- 0.0167158309
PACL_1_fst[2,7] <- 0.0069513406
PACL_1_fst[2,8] <- 0.0069513406
PACL_1_fst[2,9] <- 0.0366088632
PACL_1_fst[2,10] <- 0.0291262136
PACL_1_fst[2,11] <- 0.0224828935
PACL_1_fst[2,12] <- 0.0186457311
PACL_1_fst[3,4] <- 0.0009990010
PACL_1_fst[3,5] <- 0.0000000000
PACL_1_fst[3,6] <- 0.0167158309
PACL_1_fst[3,7] <- 0.0069513406
PACL_1_fst[3,8] <- 0.0039840637
PACL_1_fst[3,9] <- 0.0224828935
PACL_1_fst[3,10] <- 0.0205680705
PACL_1_fst[3,11] <- 0.0138067061
PACL_1_fst[3,12] <- 0.0118577075
PACL_1_fst[4,5] <- 0.0000000000
PACL_1_fst[4,6] <- 0.0089197225
PACL_1_fst[4,7] <- 0.0000000000
PACL_1_fst[4,8] <- 0.0000000000
PACL_1_fst[4,9] <- 0.0196078431
PACL_1_fst[4,10] <- 0.0176817289
PACL_1_fst[4,11] <- 0.0147783251
PACL_1_fst[4,12] <- 0.0157480315
PACL_1_fst[5,6] <- 0.0108803165
PACL_1_fst[5,7] <- 0.0029910269
PACL_1_fst[5,8] <- 0.0000000000
PACL_1_fst[5,9] <- 0.0243902439
PACL_1_fst[5,10] <- 0.0147783251
PACL_1_fst[5,11] <- 0.0128331688
PACL_1_fst[5,12] <- 0.0147783251
PACL_1_fst[6,7] <- 0.0167158309
PACL_1_fst[6,8] <- 0.0079365079
PACL_1_fst[6,9] <- 0.0356798457
PACL_1_fst[6,10] <- 0.0272373541
PACL_1_fst[6,11] <- 0.0300678952
PACL_1_fst[6,12] <- 0.0281827017
PACL_1_fst[7,8] <- 0.0029910269
PACL_1_fst[7,9] <- 0.0272373541
PACL_1_fst[7,10] <- 0.0215264188
PACL_1_fst[7,11] <- 0.0157480315
PACL_1_fst[7,12] <- 0.0186457311
PACL_1_fst[8,9] <- 0.0234375000
PACL_1_fst[8,10] <- 0.0215264188
PACL_1_fst[8,11] <- 0.0128331688
PACL_1_fst[8,12] <- 0.0157480315
PACL_1_fst[9,10] <- 0.0176817289
PACL_1_fst[9,11] <- 0.0224828935
PACL_1_fst[9,12] <- 0.0224828935
PACL_1_fst[10,11] <- 0.0118577075
PACL_1_fst[10,12] <- 0.0128331688
PACL_1_fst[11,12] <- 0.0000000000

# make symmetric and clean
PACL_1_fst[lower.tri(PACL_1_fst)] <- t(PACL_1_fst)[lower.tri(PACL_1_fst)]
PACL_1_fst[PACL_1_fst < 0] <- 0
diag(PACL_1_fst) <- 0

stopifnot(identical(dim(PACL_1_fst), c(12L, 12L)))
stopifnot(isTRUE(all.equal(PACL_1_fst, t(PACL_1_fst), check.attributes = FALSE)))
stopifnot(identical(as.character(PACL_1_coords$site), rownames(PACL_1_fst)))

# -----------------------------
# 3) geographic distances (km)
# -----------------------------
coords <- PACL_1_coords[, c("lon", "lat")]
geo_dist_km <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
rownames(geo_dist_km) <- colnames(geo_dist_km) <- PACL_1_coords$site

# -----------------------------
# 4) IBD dataframe
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(PACL_1_fst)[row(PACL_1_fst)[upper.tri(PACL_1_fst)]],
  site2   = colnames(PACL_1_fst)[col(PACL_1_fst)[upper.tri(PACL_1_fst)]],
  fst     = PACL_1_fst[upper.tri(PACL_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 5) map of sampling sites
# -----------------------------
world_map <- map_data("world")
usa_map   <- subset(world_map, region == "USA")

map_df <- PACL_1_coords
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
    title = "PACL-1 sampling sites"
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
    title = "PACL-1 isolation by distance"
  )

# -----------------------------
# 7) save RData
# -----------------------------
save(
  PACL_1_fst,
  PACL_1_coords,
  file = out_rdata
)

# -----------------------------
# 8) print outputs
# -----------------------------
print(PACL_1_coords)
print(round(PACL_1_fst, 4))
print(head(ibd_df))
print(p_map)
print(p_ibd)
