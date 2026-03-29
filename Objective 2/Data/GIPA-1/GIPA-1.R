#-----------------
# GIPA-1
# Gila pandora
# Rio Grande chub
#-----------------

# -----------------------------
# site IDs from the study
# UTM coordinates from Table 1 converted from Zone 13N to lat/lon
# -----------------------------
# 1  Alamosa Creek
# 2  Jemez River
# 3  East Fork Jemez River
# 4  Rio Guadalupe
# 5  Rio Cebolla
# 6  Rio de las Vacas
# 7  Rito Penas Negras
# 8  El Rito
# 9  Rio Ojo Caliente
# 10 Rio Vallecitos
# 11 Rio Pueblo de Taos
# 12 Cieneguilla Creek
# 13 Pecos River
# 14 Upper Rio Bonito
# 15 Rio Penasco

# -----------------------------
# coordinates
# converted from UTM Zone 13N listed in Table 1
# -----------------------------
library(sf)

utm_df <- data.frame(
  site_id = 1:15,
  easting = c(
    261152, 347507, 356527, 338560, 345507,
    338464, 339664, 389424, 406630, 403895,
    434728, 476304, 439327, 446819, 477332
  ),
  northing = c(
    3717396, 3961478, 3965771, 3965077, 3971767,
    3977630, 3982014, 4027618, 4022484, 4033935,
    4021645, 4038428, 3937179, 3703745, 3640086
  )
)

utm_sf <- st_as_sf(utm_df, coords = c("easting", "northing"), crs = 32613)
lonlat <- st_transform(utm_sf, 4326)
lonlat_coords <- st_coordinates(lonlat)

GIPA_1_coords <- data.frame(
  site_id = utm_df$site_id,
  lat = lonlat_coords[, "Y"],
  lon = lonlat_coords[, "X"]
)

# -----------------------------
# pairwise FST matrix
# transcribed from Table 2
# lower triangle shown in the paper and rebuilt here
# as a full symmetric matrix
# -----------------------------
n_sites <- 15
GIPA_1_fst <- matrix(0, nrow = n_sites, ncol = n_sites)
rownames(GIPA_1_fst) <- colnames(GIPA_1_fst) <- as.character(1:15)

# site order in matrix matches Table 2:
# 1  Alamosa Creek
# 2  Jemez River
# 3  East Fork Jemez River
# 4  Rio Guadalupe
# 5  Rio Cebolla
# 6  Rio de las Vacas
# 7  Rito Penas Negras
# 8  El Rito
# 9  Rio Ojo Caliente
# 10 Rio Vallecitos
# 11 Rio Pueblo de Taos
# 12 Cieneguilla Creek
# 13 Pecos River
# 14 Upper Rio Bonito
# 15 Rio Penasco

lt <- list(
  `2`  = c(0.0800),
  `3`  = c(0.0820, 0.0270),
  `4`  = c(0.0760, 0.0260, 0.0400),
  `5`  = c(0.1040, 0.0630, 0.0730, 0.0680),
  `6`  = c(0.0750, 0.0360, 0.0380, 0.0040, 0.0640),
  `7`  = c(0.0800, 0.0250, 0.0370, -0.0020, 0.0650, 0.0070),
  `8`  = c(0.0950, 0.0700, 0.0630, 0.0680, 0.0880, 0.0590, 0.0680),
  `9`  = c(0.0690, 0.0310, 0.0260, 0.0310, 0.0530, 0.0190, 0.0340, 0.0400),
  `10` = c(0.0710, 0.0500, 0.0370, 0.0450, 0.0780, 0.0310, 0.0510, 0.0570, 0.0130),
  `11` = c(0.0770, 0.0350, 0.0430, 0.0390, 0.0500, 0.0290, 0.0370, 0.0600, 0.0290, 0.0410),
  `12` = c(0.0810, 0.0550, 0.0520, 0.0520, 0.0810, 0.0460, 0.0550, 0.0550, 0.0380, 0.0500, 0.0470),
  `13` = c(0.0740, 0.0530, 0.0570, 0.0460, 0.0730, 0.0360, 0.0400, 0.0650, 0.0320, 0.0510, 0.0500, 0.0580),
  `14` = c(0.1460, 0.1380, 0.1520, 0.1470, 0.1710, 0.1410, 0.1430, 0.1640, 0.1270, 0.1560, 0.1450, 0.1640, 0.0990),
  `15` = c(0.1310, 0.0830, 0.0970, 0.0980, 0.1080, 0.0900, 0.0940, 0.0950, 0.0580, 0.1000, 0.0770, 0.1010, 0.0700, 0.1650)
)

for (i in 2:n_sites) {
  GIPA_1_fst[i, 1:(i - 1)] <- lt[[as.character(i)]]
}

GIPA_1_fst <- GIPA_1_fst + t(GIPA_1_fst)
GIPA_1_fst[GIPA_1_fst < 0] <- 0
diag(GIPA_1_fst) <- 0

# -----------------------------
# plotting packages
# -----------------------------
library(ggplot2)
library(dplyr)
library(geosphere)

# -----------------------------
# site labels for plotting only
# not written into GIPA_1_coords
# -----------------------------
plot_labs <- c(
  "1. Alamosa Creek",
  "2. Jemez River",
  "3. East Fork Jemez River",
  "4. Rio Guadalupe",
  "5. Rio Cebolla",
  "6. Rio de las Vacas",
  "7. Rito Penas Negras",
  "8. El Rito",
  "9. Rio Ojo Caliente",
  "10. Rio Vallecitos",
  "11. Rio Pueblo de Taos",
  "12. Cieneguilla Creek",
  "13. Pecos River",
  "14. Upper Rio Bonito",
  "15. Rio Penasco"
)

plot_df <- GIPA_1_coords
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
    hjust = 0, nudge_x = 0.08, size = 3
  ) +
  coord_fixed(1.25, xlim = c(-108.6, -104.0), ylim = c(32.5, 37.0)) +
  theme_bw() +
  labs(
    title = "GIPA-1 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )

# -----------------------------
# pairwise geographic distance
# Euclidean / great-circle distance in km
# -----------------------------
coords_mat <- as.matrix(GIPA_1_coords[, c("lon", "lat")])
geo_dist_km <- geosphere::distm(coords_mat, fun = geosphere::distGeo) / 1000
rownames(geo_dist_km) <- colnames(geo_dist_km) <- as.character(GIPA_1_coords$site_id)

# -----------------------------
# pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(GIPA_1_fst)[row(GIPA_1_fst)[upper.tri(GIPA_1_fst)]],
  site2   = colnames(GIPA_1_fst)[col(GIPA_1_fst)[upper.tri(GIPA_1_fst)]],
  fst     = GIPA_1_fst[upper.tri(GIPA_1_fst)],
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
    title = "GIPA-1 isolation by distance",
    x = "Geographic distance (km)",
    y = expression(F[ST])
  )

# -----------------------------
# save RData
# -----------------------------
save(
  GIPA_1_fst,
  GIPA_1_coords,
  file = file.path("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/GIPA-1/data/GIPA-1.RData")
)
