# -----------------------------
# AMNA-1 Ameiurus natalis
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
out_rdata <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/AMNA-1/data/AMNA-1.RData"

# -----------------------------
# 1) sites in Table S4 order
# original_site gives the site number from Pilger et al. Table S1 / Fig. 1
# site is the renumbered Objective 2 site ID used in the matrix
# -----------------------------
AMNA_1_coords <- data.frame(
  site = 1:5,
  #original_site = c(16, 13, 12, 10, 9),
  #stream = c(
  #  "Upper East Fork Gila River",
  #  "Middle Fork Gila River",
  #  "West Fork Gila River at Heart Bar Game Reserve",
  #  "Lower East Fork Gila River",
  #  "Gila River below confluence with East Fork"
  #),
  #classification = c("Tributary", "Tributary", "Mainstem", "Tributary", "Mainstem"),
  lat = c(33.30011181, 33.22604743, 33.2141233, 33.177019, 33.17707829),
  lon = c(-108.1229968, -108.2412087, -108.2289433, -108.199432, -108.2074147),
  #elev_m = c(1878, 1725, 1711, 1694, 1691),
  #river_km_from_site1 = c(197.0, 169.0, 166.0, 155.0, 159.0),
  stringsAsFactors = FALSE
)

stopifnot(nrow(AMNA_1_coords) == 5)
stopifnot(identical(AMNA_1_coords$site, 1:5))

# -----------------------------
# 2) pairwise FST matrix
# back-calculated from Slatkin's linearized FST in Table S4
# FST = linearized_FST / (1 + linearized_FST)
# -----------------------------
AMNA_1_fst <- matrix(0, nrow = 5, ncol = 5)
rownames(AMNA_1_fst) <- colnames(AMNA_1_fst) <- as.character(1:5)

AMNA_1_fst[1,2] <- 0.0000000000
AMNA_1_fst[1,3] <- 0.0128331688
AMNA_1_fst[1,4] <- 0.0059642147
AMNA_1_fst[1,5] <- 0.0049751244
AMNA_1_fst[2,3] <- 0.0147783251
AMNA_1_fst[2,4] <- 0.0019960080
AMNA_1_fst[2,5] <- 0.0019960080
AMNA_1_fst[3,4] <- 0.0029910269
AMNA_1_fst[3,5] <- 0.0000000000
AMNA_1_fst[4,5] <- 0.0000000000

# make symmetric and clean
AMNA_1_fst[lower.tri(AMNA_1_fst)] <- t(AMNA_1_fst)[lower.tri(AMNA_1_fst)]
AMNA_1_fst[AMNA_1_fst < 0] <- 0
diag(AMNA_1_fst) <- 0

stopifnot(identical(dim(AMNA_1_fst), c(5L, 5L)))
stopifnot(isTRUE(all.equal(AMNA_1_fst, t(AMNA_1_fst), check.attributes = FALSE)))
stopifnot(identical(as.character(AMNA_1_coords$site), rownames(AMNA_1_fst)))

# -----------------------------
# 3) geographic distances (km)
# -----------------------------
coords <- AMNA_1_coords[, c("lon", "lat")]
geo_dist_km <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
rownames(geo_dist_km) <- colnames(geo_dist_km) <- AMNA_1_coords$site

# -----------------------------
# 4) IBD dataframe
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(AMNA_1_fst)[row(AMNA_1_fst)[upper.tri(AMNA_1_fst)]],
  site2   = colnames(AMNA_1_fst)[col(AMNA_1_fst)[upper.tri(AMNA_1_fst)]],
  fst     = AMNA_1_fst[upper.tri(AMNA_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 5) map of sampling sites
# -----------------------------
world_map <- map_data("world")
usa_map   <- subset(world_map, region == "USA")

map_df <- AMNA_1_coords
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
    title = "AMNA-1 sampling sites"
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
    title = "AMNA-1 isolation by distance"
  )

# -----------------------------
# 7) save RData
# -----------------------------
save(
  AMNA_1_fst,
  AMNA_1_coords,
  file = out_rdata
)

# -----------------------------
# 8) print outputs
# -----------------------------
print(AMNA_1_coords)
print(round(AMNA_1_fst, 4))
print(head(ibd_df))
print(p_map)
print(p_ibd)
