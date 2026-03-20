# -----------------------------
# LARI-1 western brook lamprey
# site coordinates + FST matrix + IBD plot
# -----------------------------

library(ggplot2)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/LARI-1"

# -----------------------------
# 0) site coordinates
# (ONLY sites 1–17 used in FST table)
# -----------------------------
LARI_1_coords <- data.frame(
  site = as.character(1:17),
  lat = c(
    46.306, 46.33, 46.05, 46.22, 45.925,
    45.795, 45.924, 45.935, 45.94, 45.907,
    45.853, 45.742, 45.747, 45.764, 45.57,
    46.027, 45.94
  ),
  lon = c(
    -123.799, -123.912, -123.723, -123.345, -122.924,
    -122.922, -122.497, -122.377, -122.402, -122.382,
    -122.636, -122.545, -122.477, -122.432, -122.317,
    -121.575, -118.384
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map
# -----------------------------
ggplot(LARI_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 2.5) +
  geom_text(aes(label = site), size = 3, nudge_y = 0.05) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "LARI-1 sampling locations"
  )

# -----------------------------
# 2) geographic distance matrix
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(LARI_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- LARI_1_coords$site
colnames(geo_dist_km) <- LARI_1_coords$site

# -----------------------------
# 3) FST matrix
# NOTE:
# values are ABOVE diagonal in paper
# we will input upper triangle and mirror
# -----------------------------
LARI_1_fst <- matrix(
  0,
  nrow = 17,
  ncol = 17,
  dimnames = list(LARI_1_coords$site, LARI_1_coords$site)
)

# ---- row 1 ----
LARI_1_fst["1", 2:17] <- c(
  0.0199, 0.1823, 0.1021, 0.1727, 0.2381,
  0.4379, 0.3138, 0.4032, 0.3997, 0.2606,
  0.3098, 0.2486, 0.2266, 0.1941, 0.0912, 0.3088
)

# ---- row 2 ----
LARI_1_fst["2", 3:17] <- c(
  0.1886, 0.1086, 0.1430, 0.2295, 0.4610,
  0.3216, 0.4251, 0.4273, 0.2534, 0.3121,
  0.2604, 0.2307, 0.1823, 0.1049, 0.3239
)

# ---- row 3 ----
LARI_1_fst["3", 4:17] <- c(
  0.0927, 0.1470, 0.1969, 0.4893, 0.3594,
  0.4523, 0.4585, 0.2033, 0.2620, 0.2448,
  0.2075, 0.1973, 0.2548, 0.3884
)

# ---- row 4 ----
LARI_1_fst["4", 5:17] <- c(
  0.1248, 0.1427, 0.4147, 0.2941, 0.3767,
  0.3849, 0.1781, 0.2055, 0.1932, 0.1547,
  0.1523, 0.1693, 0.4037
)

# ---- row 5 ----
LARI_1_fst["5", 6:17] <- c(
  0.1932, 0.4244, 0.3033, 0.3857, 0.3998,
  0.2495, 0.3026, 0.2708, 0.2275, 0.2004,
  0.1956, 0.3573
)

# ---- row 6 ----
LARI_1_fst["6", 7:17] <- c(
  0.5104, 0.4008, 0.4768, 0.4808, 0.3264,
  0.3585, 0.3504, 0.2998, 0.2771, 0.2246, 0.4070
)

# ---- row 7 ----
LARI_1_fst["7", 8:17] <- c(
  0.2687, -0.0026, 0.0081, 0.4627, 0.5131,
  0.4811, 0.4423, 0.4127, 0.4663, 0.7117
)

# ---- row 8 ----
LARI_1_fst["8", 9:17] <- c(
  0.2407, 0.2117, 0.2638, 0.2972, 0.2539,
  0.2187, 0.2056, 0.3003, 0.5666
)

# ---- row 9 ----
LARI_1_fst["9", 10:17] <- c(
  0.0065, 0.4294, 0.4774, 0.4432,
  0.4051, 0.3776, 0.4329, 0.6754
)

# ---- row 10 ----
LARI_1_fst["10", 11:17] <- c(
  0.4204, 0.4669, 0.4268, 0.3876,
  0.3626, 0.4213, 0.6749
)

# ---- row 11 ----
LARI_1_fst["11", 12:17] <- c(
  0.0334, 0.0484, 0.0359, 0.0324,
  0.2783, 0.4857
)

# ---- row 12 ----
LARI_1_fst["12", 13:17] <- c(
  0.0167, 0.0244, 0.0832, 0.3275, 0.5615
)

# ---- row 13 ----
LARI_1_fst["13", 14:17] <- c(
  0.0021, 0.0520, 0.2585, 0.5200
)

# ---- row 14 ----
LARI_1_fst["14", 15:17] <- c(
  0.0261, 0.2337, 0.4982
)

# ---- row 15 ----
LARI_1_fst["15", 16:17] <- c(
  0.1879, 0.4584
)

# ---- row 16 ----
LARI_1_fst["16", 17] <- 0.2916

# -----------------------------
# mirror + clean
# -----------------------------
LARI_1_fst[lower.tri(LARI_1_fst)] <- t(LARI_1_fst)[lower.tri(LARI_1_fst)]
diag(LARI_1_fst) <- 0
LARI_1_fst[LARI_1_fst < 0] <- 0  # remove negative estimate

# -----------------------------
# 4) IBD dataframe
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(LARI_1_fst)[row(LARI_1_fst)[upper.tri(LARI_1_fst)]],
  site2   = colnames(LARI_1_fst)[col(LARI_1_fst)[upper.tri(LARI_1_fst)]],
  fst     = LARI_1_fst[upper.tri(LARI_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 5) IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "LARI-1 isolation by distance"
  )

# -----------------------------
# 6) save (only fst + coords)
# -----------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  LARI_1_fst,
  LARI_1_coords,
  file = file.path(out_dir, "LARI-1.RData")
)