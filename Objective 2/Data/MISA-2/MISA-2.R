# -----------------------------
# MISA-2
# Micropterus salmoides
# statewide Alabama location-level SNP FST
# hard-coded from Appendix 8
# -----------------------------

library(ggplot2)
library(geosphere)
library(maps)

# -----------------------------
# 0) paths
# -----------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/MISA-2"
data_dir <- file.path(base_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) site coordinates
# best-available approximate site centroids
# site order follows Figure 1 / Appendix 8
# -----------------------------
MISA_2_coords <- data.frame(
  site_id = 1:29,
  site_name = c(
    "Lake Harding",
    "Lake Eufaula",
    "Harris Reservoir",
    "Lake Martin",
    "Yates Reservoir",
    "Weiss Reservoir",
    "Neely Henry Reservoir",
    "Logan Martin Reservoir",
    "Lay Lake",
    "Lake Guntersville",
    "Wheeler Reservoir",
    "Wilson Reservoir",
    "Pickwick Reservoir",
    "Bear Creek Reservoir",
    "Jones Bluff",
    "Miller's Ferry",
    "Claiborne",
    "Lewis Smith Reservoir",
    "Lake Tuscaloosa",
    "Sipsey River",
    "Demopolis",
    "Big Bayou Canot",
    "Crab Creek",
    "Tensaw Lake",
    "D'Olive Bay",
    "Dog River",
    "Fowl River",
    "Fish River",
    "Styx River"
  ),
  lat = c(
    32.6860, 31.8982, 33.2834, 32.7262, 32.5741,
    34.1990, 33.8889, 33.5305, 33.0188, 34.4231,
    34.6695, 34.8006, 34.9023, 34.3800, 32.3242,
    32.1001, 31.5468, 33.9417, 33.3427, 33.2571,
    32.5212, 30.7857, 30.7310, 31.0480, 30.6383,
    30.5649, 30.4305, 30.4602, 30.5183
  ),
  lon = c(
    -85.1576, -85.1602, -85.6095, -85.9458, -85.8891,
    -85.6161, -86.0476, -86.1882, -86.5520, -86.3922,
    -87.0472, -87.6259, -88.0387, -87.9500, -86.7843,
    -87.3978, -87.5125, -87.1083, -87.5561, -87.7764,
    -87.8792, -88.0128, -87.9756, -87.8875, -87.9162,
    -88.0883, -88.1358, -87.8039, -87.4628
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) helper: fill symmetric matrix
# from lower-triangle values entered row-wise
# -----------------------------
fill_sym_from_lower <- function(pops, vals, diag_val = 0) {
  n <- length(pops)
  stopifnot(length(vals) == n * (n - 1) / 2)
  
  mat <- matrix(
    0,
    nrow = n,
    ncol = n,
    dimnames = list(pops, pops)
  )
  
  k <- 1
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      mat[i, j] <- vals[k]
      mat[j, i] <- vals[k]
      k <- k + 1
    }
  }
  
  diag(mat) <- diag_val
  mat
}

# -----------------------------
# 3) pairwise FST matrix
# lower triangle from Appendix 8
# row-wise by site order 1:29
# negatives set to 0 after reconstruction
# -----------------------------
fst_vals <- c(
  
  # row 2: Lake Eufaula
  0.0327,
  
  # row 3: Harris Reservoir
  0.0440, 0.1190,
  
  # row 4: Lake Martin
  0.0622, 0.1329, 0.0327,
  
  # row 5: Yates Reservoir
  0.0340, 0.1007, 0.0195, 0.0217,
  
  # row 6: Weiss Reservoir
  0.1719, 0.2874, 0.0629, 0.0964, 0.1044,
  
  # row 7: Neely Henry Reservoir
  0.1790, 0.2906, 0.0679, 0.0908, 0.1008, 0.0019,
  
  # row 8: Logan Martin Reservoir
  0.1360, 0.2425, 0.0337, 0.0609, 0.0599, 0.0128, 0.0097,
  
  # row 9: Lay Lake
  0.0612, 0.1302, 0.0155, 0.0455, 0.0250, 0.0795, 0.0713, 0.0452,
  
  # row 10: Lake Guntersville
  0.1555, 0.2476, 0.0673, 0.0862, 0.1116, 0.0436, 0.0549, 0.0608, 0.1085,
  
  # row 11: Wheeler Reservoir
  0.4040, 0.5177, 0.2797, 0.3103, 0.3598, 0.1635, 0.1851, 0.2262, 0.3208, 0.0756,
  
  # row 12: Wilson Reservoir
  0.3496, 0.4679, 0.2273, 0.2631, 0.3090, 0.1240, 0.1449, 0.1806, 0.2672, 0.0619, 0.0074,
  
  # row 13: Pickwick Reservoir
  0.3806, 0.4993, 0.2595, 0.2928, 0.3437, 0.1567, 0.1771, 0.2146, 0.2987, 0.0754, -0.0007, 0.0073,
  
  # row 14: Bear Creek Reservoir
  0.2142, 0.3297, 0.0898, 0.1202, 0.1395, 0.0161, 0.0185, 0.0358, 0.1066, 0.0381, 0.1207, 0.0860, 0.1125,
  
  # row 15: Jones Bluff
  0.2791, 0.4024, 0.1487, 0.1960, 0.2000, 0.0575, 0.0435, 0.0656, 0.1279, 0.1183, 0.2474, 0.2024, 0.2481, 0.0504,
  
  # row 16: Miller's Ferry
  0.3456, 0.4727, 0.2126, 0.2659, 0.2747, 0.0898, 0.0745, 0.1108, 0.1788, 0.1455, 0.2742, 0.2381, 0.2842, 0.0839, 0.0052,
  
  # row 17: Claiborne
  0.3121, 0.4331, 0.1867, 0.2329, 0.2429, 0.0823, 0.0686, 0.1016, 0.1598, 0.1460, 0.2712, 0.2315, 0.2782, 0.0775, 0.0139, -0.0010,
  
  # row 18: Lewis Smith Reservoir
  0.3082, 0.4326, 0.1705, 0.2105, 0.2274, 0.0587, 0.0440, 0.0708, 0.1589, 0.1107, 0.2211, 0.1871, 0.2204, 0.0485, 0.0280, 0.0360, 0.0371,
  
  # row 19: Lake Tuscaloosa
  0.2527, 0.3687, 0.1316, 0.1686, 0.1809, 0.0465, 0.0329, 0.0635, 0.1192, 0.0954, 0.2117, 0.1651, 0.2045, 0.0312, 0.0215, 0.0313, 0.0229, 0.0217,
  
  # row 20: Sipsey River
  0.3767, 0.5004, 0.2329, 0.2843, 0.2946, 0.1006, 0.0799, 0.1261, 0.1983, 0.1464, 0.2682, 0.2342, 0.2723, 0.0801, 0.0160, 0.0165, 0.0127, 0.0361, 0.0238,
  
  # row 21: Demopolis
  0.3467, 0.4630, 0.1983, 0.2519, 0.2570, 0.0753, 0.0589, 0.0969, 0.1681, 0.1292, 0.2457, 0.2055, 0.2394, 0.0631, 0.0096, 0.0101, 0.0088, 0.0246, 0.0074, 0.0085,
  
  # row 22: Big Bayou Canot
  0.3334, 0.4540, 0.2058, 0.2583, 0.2761, 0.1048, 0.0916, 0.1243, 0.1758, 0.1528, 0.2801, 0.2392, 0.2847, 0.0983, 0.0349, 0.0224, 0.0180, 0.0567, 0.0383, 0.0327, 0.0169,
  
  # row 23: Crab Creek
  0.3306, 0.4529, 0.2074, 0.2587, 0.2717, 0.1087, 0.0943, 0.1243, 0.1754, 0.1537, 0.2877, 0.2488, 0.2966, 0.0987, 0.0243, 0.0196, 0.0362, 0.0507, 0.0554, 0.0381, 0.0247, 0.0070,
  
  # row 24: Tensaw Lake
  0.3185, 0.4390, 0.1915, 0.2380, 0.2482, 0.0922, 0.0747, 0.1073, 0.1628, 0.1399, 0.2673, 0.2241, 0.2700, 0.0809, 0.0164, 0.0070, 0.0025, 0.0350, 0.0222, 0.0169, 0.0067, 0.0109, 0.0182,
  
  # row 25: D'Olive Bay
  0.3016, 0.4194, 0.1768, 0.2324, 0.2380, 0.0925, 0.0797, 0.1032, 0.1486, 0.1440, 0.2838, 0.2385, 0.2866, 0.0860, 0.0216, 0.0167, 0.0116, 0.0430, 0.0279, 0.0276, 0.0119, 0.0036, 0.0063, -0.0010,
  
  # row 26: Dog River
  0.2747, 0.3818, 0.1511, 0.2131, 0.2109, 0.0891, 0.0783, 0.0889, 0.1304, 0.1419, 0.2926, 0.2377, 0.2904, 0.0827, 0.0330, 0.0446, 0.0399, 0.0684, 0.0397, 0.0519, 0.0239, 0.0298, 0.0457, 0.0256, 0.0100,
  
  # row 27: Fowl River
  0.2010, 0.2980, 0.0851, 0.1478, 0.1362, 0.0723, 0.0617, 0.0521, 0.0668, 0.1266, 0.3135, 0.2563, 0.3060, 0.0880, 0.0607, 0.1009, 0.0944, 0.0926, 0.0692, 0.1129, 0.0709, 0.0938, 0.0989, 0.0876, 0.0635, 0.0295,
  
  # row 28: Fish River
  0.1141, 0.2095, 0.0349, 0.0591, 0.0455, 0.0384, 0.0167, 0.0063, 0.0117, 0.0883, 0.3144, 0.2586, 0.3170, 0.0577, 0.0628, 0.1347, 0.1338, 0.0987, 0.0602, 0.1372, 0.0869, 0.1523, 0.1550, 0.1089, 0.0919, 0.0456, 0.0084,
  
  # row 29: Styx River
  0.1768, 0.2683, 0.0737, 0.1224, 0.1128, 0.0610, 0.0535, 0.0443, 0.0647, 0.1058, 0.2944, 0.2389, 0.2835, 0.0757, 0.0681, 0.1096, 0.1016, 0.0894, 0.0672, 0.1269, 0.0754, 0.1103, 0.1158, 0.0809, 0.0750, 0.0432, 0.0266, -0.0120
)

MISA_2_fst <- fill_sym_from_lower(
  pops = as.character(MISA_2_coords$site_id),
  vals = fst_vals,
  diag_val = 0
)

# user workflow: set negative FST to 0
MISA_2_fst[MISA_2_fst < 0] <- 0

# -----------------------------
# 4) checks
# -----------------------------
stopifnot(nrow(MISA_2_fst) == 29)
stopifnot(ncol(MISA_2_fst) == 29)
stopifnot(identical(rownames(MISA_2_fst), as.character(MISA_2_coords$site_id)))
stopifnot(identical(colnames(MISA_2_fst), as.character(MISA_2_coords$site_id)))
stopifnot(isTRUE(all.equal(MISA_2_fst, t(MISA_2_fst))))

# -----------------------------
# 5) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(MISA_2_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- MISA_2_coords$site_id
colnames(geo_dist_km) <- MISA_2_coords$site_id

# -----------------------------
# 6) IBD dataframe
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(MISA_2_fst)[row(MISA_2_fst)[upper.tri(MISA_2_fst)]],
  site2   = colnames(MISA_2_fst)[col(MISA_2_fst)[upper.tri(MISA_2_fst)]],
  fst     = MISA_2_fst[upper.tri(MISA_2_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 7) map
# -----------------------------
usa <- map_data("state")
can <- map_data("world", region = "Canada")

map_plot <- ggplot() +
  geom_polygon(
    data = usa,
    aes(x = long, y = lat, group = group),
    fill = "grey92", color = "black", linewidth = 0.2
  ) +
  geom_polygon(
    data = can,
    aes(x = long, y = lat, group = group),
    fill = "grey86", color = "black", linewidth = 0.2
  ) +
  geom_point(
    data = MISA_2_coords,
    aes(x = lon, y = lat),
    size = 2.2
  ) +
  geom_text(
    data = MISA_2_coords,
    aes(x = lon, y = lat, label = site_id),
    nudge_y = 0.08,
    size = 3
  ) +
  coord_fixed(
    xlim = c(min(MISA_2_coords$lon) - 1, max(MISA_2_coords$lon) + 1),
    ylim = c(min(MISA_2_coords$lat) - 1, max(MISA_2_coords$lat) + 1)
  ) +
  theme_classic() +
  labs(
    title = "MISA-2 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )

print(map_plot)

# -----------------------------
# 8) IBD plot
# -----------------------------
ibd_plot <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "MISA-2 isolation by distance"
  )

print(ibd_plot)

# -----------------------------
# 9) assign objects and save
# -----------------------------
assign("MISA_2_fst", MISA_2_fst, envir = .GlobalEnv)
assign("MISA_2_coords", MISA_2_coords, envir = .GlobalEnv)

save(
  MISA_2_fst,
  MISA_2_coords,
  file = file.path(data_dir, "MISA-2.RData")
)

write.csv(
  MISA_2_coords,
  file = file.path(data_dir, "MISA-2_coords.csv"),
  row.names = FALSE
)