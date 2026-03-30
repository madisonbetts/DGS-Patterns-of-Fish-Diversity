#-----------------
# MOSA-1
# Morone saxatilis
# striped bass
#-----------------

# updated using uploaded Figure 1 coordinates from LeBlanc et al. (2020)
# BD-MICHI / Bras d'Or-Miramichi removed from sites and FST calculations
# Hudson and Delaware pooled across years to match Figure 1 site-level locations

# -----------------------------
# retained site IDs
# original Figure 1 site numbers in parentheses
# -----------------------------
# 1  Mira River (2)
# 2  Shubenacadie River (3)
# 3  Saint John River (4)
# 4  Kennebec River (5)
# 5  Hudson River (6)
# 6  Delaware River (7)
# 7  Upper Chesapeake Bay (8)
# 8  Potomac River (9)
# 9  Rappahannock River (10)
# 10 James River (11)
# 11 Choptank River (12)
# 12 Nanticoke River (13)
# 13 Roanoke River (14)
# 14 Cape Fear River (15)

# -----------------------------
# coordinates
# digitized from uploaded Figure 1 CSV
# -----------------------------
MOSA_1_coords <- data.frame(
  site = 1:14,
  lat = c(
    45.7648906, 45.0437827, 45.2930963, 44.4798915, 40.7762392, 39.2978876, 39.5288732, 38.0664161, 37.6298330, 36.9628787, 38.6406970, 38.2430349, 35.3449110, 34.3690049
  ),
  lon = c(
    -60.0778873, -63.9253196, -66.4432289, -68.2552946, -74.0792068, -75.3669779, -75.9829762, -76.5114190, -76.5295119, -76.3944973, -76.1579983, -75.8065188, -76.2082210, -77.8017184
  )
)

# -----------------------------
# pairwise FST matrix
# BD-MICHI removed
# HUD2012 + HUD2014 pooled to Hudson River
# DEL2012 + DEL2014 pooled to Delaware River
# -----------------------------
n_sites <- 14
MOSA_1_fst <- matrix(c(
  0.0000, 0.1640, 0.1340, 0.1320, 0.1305, 0.1315, 0.1330, 0.1270, 0.1300, 0.1310, 0.1340, 0.1370, 0.1360, 0.1400, 0.1640, 0.0000, 0.1270, 0.1490, 0.1545, 0.1570, 0.1610, 0.1540, 0.1550, 0.1560, 0.1610, 0.1610, 0.1550, 0.1610, 0.1340, 0.1270, 0.0000, 0.0910, 0.0930, 0.0940, 0.0970, 0.0930, 0.0940, 0.0980, 0.0970, 0.1000, 0.0990, 0.1010, 0.1320, 0.1490, 0.0910, 0.0000, 0.0060, 0.0130, 0.0140, 0.0130, 0.0150, 0.0130, 0.0200, 0.0220, 0.0280, 0.0310, 0.1305, 0.1545, 0.0930, 0.0060, 0.0000, 0.0145, 0.0150, 0.0125, 0.0160, 0.0130, 0.0210, 0.0225, 0.0270, 0.0280, 0.1315, 0.1570, 0.0940, 0.0130, 0.0145, 0.0000, 0.0000, 0.0010, 0.0005, 0.0030, 0.0025, 0.0035, 0.0220, 0.0260, 0.1330, 0.1610, 0.0970, 0.0140, 0.0150, 0.0000, 0.0000, 0.0020, 0.0000, 0.0040, 0.0020, 0.0020, 0.0250, 0.0280, 0.1270, 0.1540, 0.0930, 0.0130, 0.0125, 0.0010, 0.0020, 0.0000, 0.0020, 0.0010, 0.0080, 0.0080, 0.0210, 0.0250, 0.1300, 0.1550, 0.0940, 0.0150, 0.0160, 0.0005, 0.0000, 0.0020, 0.0000, 0.0040, 0.0010, 0.0020, 0.0240, 0.0280, 0.1310, 0.1560, 0.0980, 0.0130, 0.0130, 0.0030, 0.0040, 0.0010, 0.0040, 0.0000, 0.0100, 0.0110, 0.0200, 0.0250, 0.1340, 0.1610, 0.0970, 0.0200, 0.0210, 0.0025, 0.0020, 0.0080, 0.0010, 0.0100, 0.0000, 0.0000, 0.0290, 0.0340, 0.1370, 0.1610, 0.1000, 0.0220, 0.0225, 0.0035, 0.0020, 0.0080, 0.0020, 0.0110, 0.0000, 0.0000, 0.0290, 0.0350, 0.1360, 0.1550, 0.0990, 0.0280, 0.0270, 0.0220, 0.0250, 0.0210, 0.0240, 0.0200, 0.0290, 0.0290, 0.0000, 0.0040, 0.1400, 0.1610, 0.1010, 0.0310, 0.0280, 0.0260, 0.0280, 0.0250, 0.0280, 0.0250, 0.0340, 0.0350, 0.0040, 0.0000
), nrow = 14, byrow = TRUE)

rownames(MOSA_1_fst) <- colnames(MOSA_1_fst) <- as.character(1:14)
diag(MOSA_1_fst) <- 0
MOSA_1_fst[MOSA_1_fst < 0] <- 0

# -----------------------------
# plotting packages
# -----------------------------
library(ggplot2)
library(dplyr)
library(geosphere)

# -----------------------------
# site labels for plotting only
# not written into MOSA_1_coords
# -----------------------------
plot_labs <- c(
  "1. Mira River",
  "2. Shubenacadie River",
  "3. Saint John River",
  "4. Kennebec River",
  "5. Hudson River",
  "6. Delaware River",
  "7. Upper Chesapeake Bay",
  "8. Potomac River",
  "9. Rappahannock River",
  "10. James River",
  "11. Choptank River",
  "12. Nanticoke River",
  "13. Roanoke River",
  "14. Cape Fear River"
)

plot_df <- MOSA_1_coords
plot_df$label <- plot_labs

# -----------------------------
# base map: US + Canada
# -----------------------------
usa <- ggplot2::map_data("world", region = "USA")
can <- ggplot2::map_data("world", region = "Canada")
world_map <- bind_rows(usa, can)

# -----------------------------
# site map
# -----------------------------
ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "grey90", color = "grey50", linewidth = 0.2
  ) +
  geom_point(
    data = plot_df,
    aes(x = lon, y = lat),
    size = 3
  ) +
  geom_text(
    data = plot_df,
    aes(x = lon, y = lat, label = label),
    hjust = 0, nudge_x = 0.25, size = 3
  ) +
  coord_fixed(1.25, xlim = c(-80, -58), ylim = c(33, 47)) +
  theme_bw() +
  labs(
    title = "MOSA-1 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )

# -----------------------------
# pairwise geographic distance
# Euclidean / great-circle distance in km
# -----------------------------
coords_mat <- as.matrix(MOSA_1_coords[, c("lon", "lat")])
geo_dist_km <- geosphere::distm(coords_mat, fun = geosphere::distGeo) / 1000
rownames(geo_dist_km) <- colnames(geo_dist_km) <- as.character(MOSA_1_coords$site)

# -----------------------------
# pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(MOSA_1_fst)[row(MOSA_1_fst)[upper.tri(MOSA_1_fst)]],
  site2   = colnames(MOSA_1_fst)[col(MOSA_1_fst)[upper.tri(MOSA_1_fst)]],
  fst     = MOSA_1_fst[upper.tri(MOSA_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(
    title = "MOSA-1 isolation by distance",
    x = "Geographic distance (km)",
    y = expression(F[ST])
  )

# -----------------------------
# save RData
# -----------------------------
save(
  MOSA_1_fst,
  MOSA_1_coords,
  file = file.path("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/MOSA-1/data/MOSA-1.RData")
)
