# -----------------------------
# ARIN-1 Sacramento perch
# Archoplites interruptus
# Coen et al. 2021
# Full metadata set: Cluster A + Benton + Cluster B
# Pairwise microsatellite FST extracted from Fig. 2 heatmap
# Coordinates inferred from Fig. 1 map + named waterbodies / localities
# -----------------------------

library(ggplot2)
library(ggrepel)
library(dplyr)
library(maps)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ARIN-1"
data_dir <- file.path(save_dir, "data")

# -----------------------------
# 1) plotting metadata
# numeric site IDs match Fig. 2 order exactly
# -----------------------------
ARIN_1_site_key <- data.frame(
  site = as.character(1:15),
  code = c(
    "ABBL", "GRYL", "JEWL", "ALMR", "CLLR", "WSVR", "LWSL", "PYRL", "SNWR",
    "BNT23", "BNT4", "CROL", "SINL", "BSCR", "BRPT"
  ),
  site_name = c(
    "Abbotts Lagoon",
    "Gray Lodge Wildlife Area",
    "Jewel Lake",
    "Lake Almanor",
    "Clear Lake Reservoir",
    "West Valley Reservoir",
    "Little Washoe Lake",
    "Pyramid Lake",
    "Stillwater National Wildlife Refuge",
    "Benton Ponds 2 and 3",
    "Benton Pond 4",
    "Crowley Lake",
    "Sindicich Lagoon",
    "Biscar Reservoir",
    "Bridgeport Reservoir"
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) coordinates object for saving
# final saved df has only site, lat, lon
# XYs are best-inference centroids because the paper
# does not publish a coordinate table
# -----------------------------
ARIN_1_coords <- data.frame(
  site = as.character(1:15),
  lat = c(
    38.1269,   # ABBL
    39.3269,   # GRYL
    37.9077,   # JEWL
    40.2500,   # ALMR
    41.9261,   # CLLR
    41.2004,   # WSVR
    39.2500,   # LWSL
    40.0800,   # PYRL
    39.5000,   # SNWR
    37.8191,   # BNT23
    37.8241,   # BNT4
    37.6174,   # CROL
    37.9828,   # SINL
    40.5497,   # BSCR
    38.2888    # BRPT
  ),
  lon = c(
    -122.9500, # ABBL
    -121.8104, # GRYL
    -122.2390, # JEWL
    -121.1500, # ALMR
    -121.0760, # CLLR
    -120.3997, # WSVR
    -119.8000, # LWSL
    -119.6000, # PYRL
    -118.9500, # SNWR
    -118.4765, # BNT23
    -118.4705, # BNT4
    -118.7397, # CROL
    -122.1326, # SINL
    -120.3366, # BSCR
    -119.2312  # BRPT
  ),
  stringsAsFactors = FALSE
)

# merge for plotting labels
ARIN_1_coords_plot <- left_join(ARIN_1_coords, ARIN_1_site_key, by = "site")

# -----------------------------
# 3) helper function
# fill symmetric matrix from upper-triangle values entered row-wise
# -----------------------------
fill_sym_from_upper <- function(pops, vals, diag_val = 0) {
  n <- length(pops)
  stopifnot(length(vals) == n * (n - 1) / 2)
  
  mat <- matrix(
    0,
    nrow = n,
    ncol = n,
    dimnames = list(pops, pops)
  )
  
  k <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      mat[i, j] <- vals[k]
      mat[j, i] <- vals[k]
      k <- k + 1
    }
  }
  
  diag(mat) <- diag_val
  mat
}

# -----------------------------
# 4) pairwise FST matrix
# values transcribed from Fig. 2 heatmap row-by-row
# order:
# ABBL, GRYL, JEWL, ALMR, CLLR, WSVR, LWSL, PYRL, SNWR,
# BNT23, BNT4, CROL, SINL, BSCR, BRPT
# -----------------------------
fst_vals <- c(
  # ABBL vs others
  0.099, 0.054, 0.057, 0.048, 0.073, 0.067, 0.068, 0.086, 0.102, 0.102, 0.134, 0.090, 0.112, 0.170,
  
  # GRYL vs others
  0.035, 0.076, 0.091, 0.102, 0.133, 0.126, 0.176, 0.145, 0.159, 0.221, 0.200, 0.198, 0.293,
  
  # JEWL vs others
  0.048, 0.064, 0.065, 0.135, 0.095, 0.142, 0.113, 0.116, 0.181, 0.196, 0.178, 0.258,
  
  # ALMR vs others
  0.023, 0.033, 0.092, 0.079, 0.100, 0.097, 0.108, 0.166, 0.138, 0.162, 0.215,
  
  # CLLR vs others
  0.038, 0.137, 0.090, 0.121, 0.096, 0.106, 0.167, 0.172, 0.159, 0.237,
  
  # WSVR vs others
  0.108, 0.088, 0.116, 0.097, 0.096, 0.178, 0.158, 0.159, 0.226,
  
  # LWSL vs others
  0.093, 0.127, 0.143, 0.147, 0.223, 0.270, 0.245, 0.320,
  
  # PYRL vs others
  0.067, 0.106, 0.104, 0.153, 0.136, 0.128, 0.226,
  
  # SNWR vs others
  0.151, 0.146, 0.231, 0.235, 0.228, 0.314,
  
  # BNT23 vs others
  0.042, 0.153, 0.140, 0.136, 0.213,
  
  # BNT4 vs others
  0.171, 0.150, 0.136, 0.232,
  
  # CROL vs others
  0.044, 0.044, 0.088,
  
  # SINL vs others
  0.063, 0.126,
  
  # BSCR vs BRPT
  0.093
)

ARIN_1_fst <- fill_sym_from_upper(
  pops = ARIN_1_coords$site,
  vals = fst_vals,
  diag_val = 0
)

ARIN_1_fst[ARIN_1_fst < 0] <- 0
diag(ARIN_1_fst) <- 0

# -----------------------------
# 5) inspect outputs
# -----------------------------
print(ARIN_1_coords)
print(round(ARIN_1_fst, 3))

cat("\nSite key:\n")
for (i in seq_len(nrow(ARIN_1_site_key))) {
  cat(
    ARIN_1_site_key$site[i], " = ",
    ARIN_1_site_key$code[i], " = ",
    ARIN_1_site_key$site_name[i], "\n",
    sep = ""
  )
}

# -----------------------------
# 6) map with US + Canada
# -----------------------------
world_map <- map_data("world") |>
  dplyr::filter(region %in% c("USA", "Canada"))

states_map <- map_data("state")

x_pad <- 4
y_pad <- 3

xlim_use <- range(ARIN_1_coords_plot$lon) + c(-x_pad, x_pad)
ylim_use <- range(ARIN_1_coords_plot$lat) + c(-y_pad, y_pad)

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
    data = ARIN_1_coords_plot,
    aes(x = lon, y = lat),
    size = 2.8
  ) +
  geom_text_repel(
    data = ARIN_1_coords_plot,
    aes(x = lon, y = lat, label = paste0(site, ": ", code)),
    size = 3.1,
    max.overlaps = 100
  ) +
  coord_quickmap(
    xlim = xlim_use,
    ylim = ylim_use,
    expand = FALSE
  ) +
  theme_classic() +
  labs(
    title = "ARIN-1 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )

print(p_map)

# -----------------------------
# 7) IBD plot
# Euclidean distance among site centroids
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(ARIN_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- ARIN_1_coords$site
colnames(geo_dist_km) <- ARIN_1_coords$site

ibd_df <- data.frame(
  site1   = rownames(ARIN_1_fst)[row(ARIN_1_fst)[upper.tri(ARIN_1_fst)]],
  site2   = colnames(ARIN_1_fst)[col(ARIN_1_fst)[upper.tri(ARIN_1_fst)]],
  fst     = ARIN_1_fst[upper.tri(ARIN_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    title = "ARIN-1 isolation by distance",
    x = "Euclidean distance among site centroids (km)",
    y = expression(F[ST])
  )

print(p_ibd)

# -----------------------------
# 8) save
# -----------------------------
save(
  ARIN_1_fst,
  ARIN_1_coords,
  file = file.path(data_dir, "ARIN-1.RData")
)