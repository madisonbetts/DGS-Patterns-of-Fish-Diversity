# -----------------------------
# COBA-1 mottled sculpin
# Cottus bairdi
# Lamphere & Blum 2012
# Table 1 site coordinates + Table 2 site-level pairwise FST
# -----------------------------

library(ggplot2)
library(ggrepel)
library(dplyr)
library(maps)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/COBA-1"
data_dir <- file.path(save_dir, "data")

# -----------------------------
# 1) site coordinates
# exact values from Table 1
# note: the published table headers appear swapped,
# but these values clearly correspond to lat ~35 N and lon ~83 W
# -----------------------------
COBA_1_coords <- data.frame(
  site = as.character(1:8),
  lat = c(
    35.042833,  # A = 35 02 34.20
    35.044597,  # B = 35 02 40.55
    35.052006,  # C = 35 03 07.22
    35.056900,  # D = 35 03 24.84
    35.059606,  # E = 35 03 34.58
    35.067678,  # F = 35 04 03.64
    35.072747,  # G = 35 04 21.89
    35.075242   # H = 35 04 30.87
  ),
  lon = c(
    -83.508408, # A = 83 30 30.27 W
    -83.510108, # B = 83 30 36.39 W
    -83.513450, # C = 83 30 48.42 W
    -83.512247, # D = 83 30 44.09 W
    -83.516503, # E = 83 30 59.41 W
    -83.520906, # F = 83 31 15.26 W
    -83.528503, # G = 83 31 42.61 W
    -83.530278  # H = 83 31 49.00 W
  ),
  stringsAsFactors = FALSE
)

# plotting helper with original paper site letters
COBA_1_coords_plot <- data.frame(
  site = as.character(1:8),
  site_code = LETTERS[1:8],
  lat = COBA_1_coords$lat,
  lon = COBA_1_coords$lon,
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) helper: fill symmetric matrix from lower triangle
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
# Table 2 values entered from the LOWER triangle row-wise
# site order in the paper: A, B, C, D, E, F, G, H
# -----------------------------
fst_vals <- c(
  # B
  0.003,
  
  # C
  0.008, 0.006,
  
  # D
  0.018, 0.013, 0.000,
  
  # E
  0.021, 0.004, 0.000, 0.003,
  
  # F
  0.039, 0.025, 0.027, 0.018, 0.009,
  
  # G
  0.071, 0.055, 0.056, 0.045, 0.039, 0.007,
  
  # H
  0.059, 0.051, 0.053, 0.043, 0.041, 0.014, 0.005
)

COBA_1_fst <- fill_sym_from_lower(
  pops = COBA_1_coords$site,
  vals = fst_vals,
  diag_val = 0
)

COBA_1_fst[COBA_1_fst < 0] <- 0
diag(COBA_1_fst) <- 0

# -----------------------------
# 4) inspect outputs
# -----------------------------
print(COBA_1_coords)
print(round(COBA_1_fst, 3))

cat("\nSite key:\n")
cat("1 = A\n")
cat("2 = B\n")
cat("3 = C\n")
cat("4 = D\n")
cat("5 = E\n")
cat("6 = F\n")
cat("7 = G\n")
cat("8 = H\n")

# -----------------------------
# 5) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(COBA_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- COBA_1_coords$site
colnames(geo_dist_km) <- COBA_1_coords$site

# -----------------------------
# 6) map of sampling locations
# US + Canada, but zoomed to the study area
# -----------------------------
world_map <- map_data("world") |>
  dplyr::filter(region %in% c("USA", "Canada"))

states_map <- map_data("state")

x_pad <- 1.2
y_pad <- 0.8

xlim_use <- range(COBA_1_coords_plot$lon) + c(-x_pad, x_pad)
ylim_use <- range(COBA_1_coords_plot$lat) + c(-y_pad, y_pad)

p_map <- ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "gray98",
    color = "gray80",
    linewidth = 0.2
  ) +
  geom_polygon(
    data = states_map,
    aes(x = long, y = lat, group = group),
    fill = "gray95",
    color = "gray70",
    linewidth = 0.2
  ) +
  geom_point(
    data = COBA_1_coords_plot,
    aes(x = lon, y = lat),
    size = 2.8
  ) +
  geom_text_repel(
    data = COBA_1_coords_plot,
    aes(x = lon, y = lat, label = paste0(site, ": ", site_code)),
    size = 3.2,
    max.overlaps = 100
  ) +
  coord_quickmap(
    xlim = xlim_use,
    ylim = ylim_use,
    expand = FALSE
  ) +
  theme_classic() +
  labs(
    title = "COBA-1 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )

print(p_map)

# -----------------------------
# 7) IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(COBA_1_fst)[row(COBA_1_fst)[upper.tri(COBA_1_fst)]],
  site2   = colnames(COBA_1_fst)[col(COBA_1_fst)[upper.tri(COBA_1_fst)]],
  fst     = COBA_1_fst[upper.tri(COBA_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    title = "COBA-1 isolation by distance",
    x = "Euclidean distance among site centroids (km)",
    y = expression(F[ST])
  )

print(p_ibd)

# -----------------------------
# 8) save
# -----------------------------
save(
  COBA_1_fst,
  COBA_1_coords,
  file = file.path(data_dir, "COBA-1.RData")
)