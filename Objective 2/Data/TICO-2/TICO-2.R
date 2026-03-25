# -----------------------------
# TICO-2 Tiaroga cobitis
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
out_rdata <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/TICO-2/data/TICO-2.RData"

# -----------------------------
# 1) sites in Table S4 order
# original_site gives the site number from Pilger et al. Table S1 / Fig. 1
# site is the renumbered Objective 2 site ID used in the matrix
# -----------------------------
TICO_2_coords <- data.frame(
  site = 1:4,
  #original_site = c(9, 5, 4, 3),
  #stream = c(
  #  "Gila River below confluence with East Fork",
  #  "Gila River at Riverside",
  #  "Gila River at Bird Area",
  #  "Gila River at confluence with Ash Canyon"
  #),
  #classification = c("Mainstem", "Mainstem", "Mainstem", "Mainstem"),
  lat = c(33.17707829, 32.9388458, 32.84751, 32.713959),
  lon = c(-108.2074147, -108.6056954, -108.593688, -108.71031),
#  elev_m = c(1691, 1360, 1331, 1240),
#river_km_from_site1 = c(159.0, 75.0, 60.0, 30.0),
  stringsAsFactors = FALSE
)

stopifnot(nrow(TICO_2_coords) == 4)
stopifnot(identical(TICO_2_coords$site, 1:4))

# -----------------------------
# 2) pairwise FST matrix
# back-calculated from Slatkin's linearized FST in Table S4
# FST = linearized_FST / (1 + linearized_FST)
# -----------------------------
TICO_2_fst <- matrix(0, nrow = 4, ncol = 4)
rownames(TICO_2_fst) <- colnames(TICO_2_fst) <- as.character(1:4)

TICO_2_fst[1,2] <- 0.0108803165
TICO_2_fst[1,3] <- 0.0147783251
TICO_2_fst[1,4] <- 0.0147783251
TICO_2_fst[2,3] <- 0.0039840637
TICO_2_fst[2,4] <- 0.0049751244
TICO_2_fst[3,4] <- 0.0049751244

# make symmetric and clean
TICO_2_fst[lower.tri(TICO_2_fst)] <- t(TICO_2_fst)[lower.tri(TICO_2_fst)]
TICO_2_fst[TICO_2_fst < 0] <- 0
diag(TICO_2_fst) <- 0

stopifnot(identical(dim(TICO_2_fst), c(4L, 4L)))
stopifnot(isTRUE(all.equal(TICO_2_fst, t(TICO_2_fst), check.attributes = FALSE)))
stopifnot(identical(as.character(TICO_2_coords$site), rownames(TICO_2_fst)))

# -----------------------------
# 3) geographic distances (km)
# -----------------------------
coords <- TICO_2_coords[, c("lon", "lat")]
geo_dist_km <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
rownames(geo_dist_km) <- colnames(geo_dist_km) <- TICO_2_coords$site

# -----------------------------
# 4) IBD dataframe
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(TICO_2_fst)[row(TICO_2_fst)[upper.tri(TICO_2_fst)]],
  site2   = colnames(TICO_2_fst)[col(TICO_2_fst)[upper.tri(TICO_2_fst)]],
  fst     = TICO_2_fst[upper.tri(TICO_2_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 5) map of sampling sites
# -----------------------------
world_map <- map_data("world")
usa_map   <- subset(world_map, region == "USA")

map_df <- TICO_2_coords
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
    title = "TICO-2 sampling sites"
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
    title = "TICO-2 isolation by distance"
  )

# -----------------------------
# 7) save RData
# -----------------------------
save(
  TICO_2_fst,
  TICO_2_coords,
  file = out_rdata
)

# -----------------------------
# 8) print outputs
# -----------------------------
print(TICO_2_coords)
print(round(TICO_2_fst, 4))
print(head(ibd_df))
print(p_map)
print(p_ibd)
