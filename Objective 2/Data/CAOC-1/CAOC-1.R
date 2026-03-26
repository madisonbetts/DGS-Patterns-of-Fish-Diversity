# -----------------------------
# CAOC-1 Sacramento sucker
# Catostomus occidentalis
# site coordinates + FST matrix + IBD plot
# -----------------------------

library(ggplot2)
library(geosphere)
library(maps)

save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/CAOC-1"

# -----------------------------
# 0) site coordinates
# best-available coordinates inferred from Figure 2 + named streams
# original sample order from paper:
# 10 Ash Valley (field microps, genetically occidentalis)
# 11 Ash Valley
# 12 Willow (Ash)
# 13 Dent
# 14 Drews
# 15 Hay
# 16 Auger
# 17 Cottonwood
# 18 Cox
# 19 Dry Creek
# 20 Willow (G)
# 21 Mainstem
# 22 NorthFork
# 23 Parker
# 24 South Fork
# 25 Stones Canyon
# 26 Turner-low
# 27 Turner-mid-b
# -----------------------------
CAOC_1_coords <- data.frame(
  site = as.character(1:18),
  #site_name = c("Ash Valley_m", "Ash Valley_o", "Willow_Ash", "Dent", "Drews", "Hay", "Auger", "Cottonwood", "Cox", "Dry Creek", "Willow_G", "Mainstem", "NorthFork", "Parker", "South Fork", "Stones Canyon", "Turner-low", "Turner-mid-b"),
  lat = c(41.131, 41.131, 41.081, 42.037, 42.012, 41.973, 42.028, 42.020, 41.949, 41.899, 41.201, 41.281, 41.224, 41.173, 41.141, 41.215, 41.224, 41.223),
  lon = c(-121.127, -121.127, -121.190, -120.901, -120.928, -121.010, -120.867, -120.913, -120.979, -120.904, -121.214, -120.995, -120.996, -120.862, -120.930, -121.122, -121.125, -121.131),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map of sampling locations
# -----------------------------
states_map <- map_data("state")
ca_or_map <- subset(states_map, region %in% c("california", "oregon"))

map_xlim <- c(-121.35, -120.75)
map_ylim <- c(41.00, 42.10)

ggplot() +
  geom_polygon(
    data = ca_or_map,
    aes(x = long, y = lat, group = group),
    fill = "grey92", color = "grey40", linewidth = 0.3
  ) +
  geom_point(
    data = CAOC_1_coords,
    aes(x = lon, y = lat),
    size = 2.5
  ) +
  geom_text(
    data = CAOC_1_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.015,
    size = 3.0
  ) +
  coord_fixed(xlim = map_xlim, ylim = map_ylim) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "CAOC-1 sampling sites"
  )

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(CAOC_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- CAOC_1_coords$site
colnames(geo_dist_km) <- CAOC_1_coords$site

# -----------------------------
# 3) pairwise FST matrix
# Appendix S1; above-diagonal FST values
# CAOC-1 includes samples 10-27, because sample 10 was
# considered C. occidentalis in population-level analyses
# -----------------------------
CAOC_1_fst <- matrix(c(
  0.000, 0.000, 0.128, 0.168, 0.143, 0.156, 0.138, 0.205, 0.211, 0.165, 0.223, 0.046, 0.022, 0.046, 0.055, 0.125, 0.041, 0.055,
  0.000, 0.000, 0.125, 0.163, 0.147, 0.164, 0.152, 0.193, 0.218, 0.175, 0.205, 0.057, 0.038, 0.049, 0.063, 0.128, 0.053, 0.068,
  0.128, 0.125, 0.000, 0.142, 0.122, 0.168, 0.150, 0.169, 0.193, 0.160, 0.196, 0.083, 0.072, 0.095, 0.099, 0.185, 0.078, 0.110,
  0.168, 0.163, 0.142, 0.000, 0.002, 0.017, 0.065, 0.035, 0.089, 0.027, 0.170, 0.112, 0.161, 0.120, 0.141, 0.236, 0.109, 0.109,
  0.143, 0.147, 0.122, 0.002, 0.000, 0.021, 0.064, 0.025, 0.074, 0.016, 0.170, 0.100, 0.098, 0.113, 0.120, 0.177, 0.098, 0.117,
  0.156, 0.164, 0.168, 0.017, 0.021, 0.000, 0.074, 0.026, 0.074, 0.034, 0.153, 0.129, 0.161, 0.136, 0.139, 0.199, 0.125, 0.136,
  0.138, 0.152, 0.150, 0.065, 0.064, 0.074, 0.000, 0.074, 0.072, 0.032, 0.165, 0.104, 0.088, 0.099, 0.111, 0.163, 0.095, 0.098,
  0.205, 0.193, 0.169, 0.035, 0.025, 0.026, 0.074, 0.000, 0.085, 0.038, 0.198, 0.147, 0.185, 0.175, 0.161, 0.227, 0.142, 0.175,
  0.211, 0.218, 0.193, 0.089, 0.074, 0.074, 0.072, 0.085, 0.000, 0.054, 0.180, 0.174, 0.169, 0.183, 0.190, 0.252, 0.164, 0.185,
  0.165, 0.175, 0.160, 0.027, 0.016, 0.034, 0.032, 0.038, 0.054, 0.000, 0.178, 0.127, 0.102, 0.124, 0.143, 0.198, 0.120, 0.130,
  0.223, 0.205, 0.196, 0.170, 0.170, 0.153, 0.165, 0.198, 0.180, 0.178, 0.000, 0.186, 0.286, 0.196, 0.206, 0.303, 0.175, 0.200,
  0.046, 0.057, 0.083, 0.112, 0.100, 0.129, 0.104, 0.147, 0.174, 0.127, 0.186, 0.000, 0.000, 0.014, 0.012, 0.061, 0.001, 0.019,
  0.022, 0.038, 0.072, 0.161, 0.098, 0.161, 0.088, 0.185, 0.169, 0.102, 0.286, 0.000, 0.000, 0.000, 0.000, 0.079, 0.000, 0.017,
  0.046, 0.049, 0.095, 0.120, 0.113, 0.136, 0.099, 0.175, 0.183, 0.124, 0.196, 0.014, 0.000, 0.000, 0.030, 0.115, 0.008, 0.015,
  0.055, 0.063, 0.099, 0.141, 0.120, 0.139, 0.111, 0.161, 0.190, 0.143, 0.206, 0.012, 0.000, 0.030, 0.000, 0.066, 0.026, 0.051,
  0.125, 0.128, 0.185, 0.236, 0.177, 0.199, 0.163, 0.227, 0.252, 0.198, 0.303, 0.061, 0.079, 0.115, 0.066, 0.000, 0.071, 0.116,
  0.041, 0.053, 0.078, 0.109, 0.098, 0.125, 0.095, 0.142, 0.164, 0.120, 0.175, 0.001, 0.000, 0.008, 0.026, 0.071, 0.000, 0.009,
  0.055, 0.068, 0.110, 0.109, 0.117, 0.136, 0.098, 0.175, 0.185, 0.130, 0.200, 0.019, 0.017, 0.015, 0.051, 0.116, 0.009, 0.000
), nrow = 18, byrow = TRUE)

rownames(CAOC_1_fst) <- colnames(CAOC_1_fst) <- CAOC_1_coords$site
CAOC_1_fst[CAOC_1_fst < 0] <- 0
diag(CAOC_1_fst) <- 0

# -----------------------------
# 4) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(CAOC_1_fst)[row(CAOC_1_fst)[upper.tri(CAOC_1_fst)]],
  site2   = colnames(CAOC_1_fst)[col(CAOC_1_fst)[upper.tri(CAOC_1_fst)]],
  fst     = CAOC_1_fst[upper.tri(CAOC_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 5) IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "CAOC-1 isolation by distance"
  )

# -----------------------------
# 6) quick checks
# -----------------------------
stopifnot(identical(rownames(CAOC_1_fst), CAOC_1_coords$site))
stopifnot(identical(colnames(CAOC_1_fst), CAOC_1_coords$site))
stopifnot(isTRUE(all.equal(CAOC_1_fst, t(CAOC_1_fst))))

# -----------------------------
# 7) save RData
# -----------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  CAOC_1_fst,
  CAOC_1_coords,
  file = file.path(out_dir, "CAOC-1.RData")
)
