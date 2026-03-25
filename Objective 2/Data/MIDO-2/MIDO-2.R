# -----------------------------
# MIDO-2 Micropterus dolomieu
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
out_rdata <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/MIDO-2/data/MIDO-2.RData"

# -----------------------------
# 1) sites in Table S4 order
# original_site gives the site number from Pilger et al. Table S1 / Fig. 1
# site is the renumbered Objective 2 site ID used in the matrix
# -----------------------------
MIDO_2_coords <- data.frame(
  site = 1:6,
  original_site = c(13, 12, 10, 9, 6, 5),
  stream = c(
    "Middle Fork Gila River",
    "West Fork Gila River at Heart Bar Game Reserve",
    "Lower East Fork Gila River",
    "Gila River below confluence with East Fork",
    "Gila River at Gila Farm",
    "Gila River at Riverside"
  ),
  classification = c("Tributary", "Mainstem", "Tributary", "Mainstem", "Mainstem", "Mainstem"),
  lat = c(33.22604743, 33.2141233, 33.177019, 33.17707829, 33.03671285, 32.9388458),
  lon = c(-108.2412087, -108.2289433, -108.199432, -108.2074147, -108.5328707, -108.6056954),
  elev_m = c(1725, 1711, 1694, 1691, 1412, 1360),
  river_km_from_site1 = c(169.0, 166.0, 155.0, 159.0, 90.0, 75.0),
  stringsAsFactors = FALSE
)

stopifnot(nrow(MIDO_2_coords) == 6)
stopifnot(identical(MIDO_2_coords$site, 1:6))

# -----------------------------
# 2) pairwise FST matrix
# back-calculated from Slatkin's linearized FST in Table S4
# FST = linearized_FST / (1 + linearized_FST)
# -----------------------------
MIDO_2_fst <- matrix(0, nrow = 6, ncol = 6)
rownames(MIDO_2_fst) <- colnames(MIDO_2_fst) <- as.character(1:6)

MIDO_2_fst[1,2] <- 0.0157480315
MIDO_2_fst[1,3] <- 0.0328820116
MIDO_2_fst[1,4] <- 0.0196078431
MIDO_2_fst[1,5] <- 0.0224828935
MIDO_2_fst[1,6] <- 0.0393852065
MIDO_2_fst[2,3] <- 0.0430622010
MIDO_2_fst[2,4] <- 0.0224828935
MIDO_2_fst[2,5] <- 0.0272373541
MIDO_2_fst[2,6] <- 0.0328820116
MIDO_2_fst[3,4] <- 0.0029910269
MIDO_2_fst[3,5] <- 0.0186457311
MIDO_2_fst[3,6] <- 0.0557129367
MIDO_2_fst[4,5] <- 0.0205680705
MIDO_2_fst[4,6] <- 0.0512333966
MIDO_2_fst[5,6] <- 0.0467111535

# make symmetric and clean
MIDO_2_fst[lower.tri(MIDO_2_fst)] <- t(MIDO_2_fst)[lower.tri(MIDO_2_fst)]
MIDO_2_fst[MIDO_2_fst < 0] <- 0
diag(MIDO_2_fst) <- 0

stopifnot(identical(dim(MIDO_2_fst), c(6L, 6L)))
stopifnot(isTRUE(all.equal(MIDO_2_fst, t(MIDO_2_fst), check.attributes = FALSE)))
stopifnot(identical(as.character(MIDO_2_coords$site), rownames(MIDO_2_fst)))

# -----------------------------
# 3) geographic distances (km)
# -----------------------------
coords <- MIDO_2_coords[, c("lon", "lat")]
geo_dist_km <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
rownames(geo_dist_km) <- colnames(geo_dist_km) <- MIDO_2_coords$site

# -----------------------------
# 4) IBD dataframe
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(MIDO_2_fst)[row(MIDO_2_fst)[upper.tri(MIDO_2_fst)]],
  site2   = colnames(MIDO_2_fst)[col(MIDO_2_fst)[upper.tri(MIDO_2_fst)]],
  fst     = MIDO_2_fst[upper.tri(MIDO_2_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 5) map of sampling sites
# -----------------------------
world_map <- map_data("world")
usa_map   <- subset(world_map, region == "USA")

map_df <- MIDO_2_coords
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
    title = "MIDO-2 sampling sites"
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
    title = "MIDO-2 isolation by distance"
  )

# -----------------------------
# 7) save RData
# -----------------------------
save(
  MIDO_2_fst,
  MIDO_2_coords,
  file = out_rdata
)

# -----------------------------
# 8) print outputs
# -----------------------------
print(MIDO_2_coords)
print(round(MIDO_2_fst, 4))
print(head(ibd_df))
print(p_map)
print(p_ibd)
