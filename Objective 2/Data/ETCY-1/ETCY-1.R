# ========================================
# ETCY-1
# Etheostoma cyanoprosopum
# Blueface darter
# Extraction from Fluker et al. 2019
# Pairwise FST from Table 4
# Coordinates from Table 1
# ========================================

# -----------------------------
# 0) setup
# -----------------------------
library(ggplot2)
library(geosphere)
library(maps)

study_code <- "ETCY-1"

base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETCY-1"
data_dir <- file.path(base_dir, "data")

dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) site coordinates
# exact coordinates from Table 1
# site order matches Table 4
# 1  Little Bear Creek
# 2  Turkey Creek
# 3  Bear Creek
# 4  Hubbard Creek above Kinlock Falls (9a)
# 5  Hubbard Creek below Kinlock Falls (9b)
# -----------------------------
ETCY_1_coords <- data.frame(
  site_id = as.character(1:5),
  lat = c(
    34.41825,
    34.36372,
    34.32992,
    34.30768,
    34.30875
  ),
  lon = c(
    -87.60293,
    -87.60761,
    -87.57977,
    -87.50340,
    -87.50225
  ),
  stringsAsFactors = FALSE
)

site_lookup <- data.frame(
  site_id = as.character(1:5),
  site_code = c("1", "2", "3", "9a", "9b"),
  site_name = c(
    "Little Bear Creek",
    "Turkey Creek",
    "Bear Creek",
    "Hubbard Creek above falls",
    "Hubbard Creek below falls"
  ),
  n_ind = c(25, 13, 25, 25, 25),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) helper
# fill symmetric matrix from lower triangle values
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
# transcribed from Table 4
# 9a vs 9b reported as 0.004 and non-significant
# -----------------------------
fst_vals <- c(
  0.115,
  0.060, 0.059,
  0.161, 0.131, 0.092,
  0.133, 0.132, 0.062, 0.004
)

ETCY_1_fst <- fill_sym_from_lower(
  pops = ETCY_1_coords$site_id,
  vals = fst_vals,
  diag_val = 0
)

ETCY_1_fst[ETCY_1_fst < 0] <- 0
diag(ETCY_1_fst) <- 0

# -----------------------------
# 4) geographic distance matrix
# Euclidean great-circle distances in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(ETCY_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- ETCY_1_coords$site_id
colnames(geo_dist_km) <- ETCY_1_coords$site_id

# -----------------------------
# 5) IBD dataframe
# -----------------------------
ibd_df <- data.frame(
  site1 = rownames(ETCY_1_fst)[row(ETCY_1_fst)[upper.tri(ETCY_1_fst)]],
  site2 = colnames(ETCY_1_fst)[col(ETCY_1_fst)[upper.tri(ETCY_1_fst)]],
  fst = ETCY_1_fst[upper.tri(ETCY_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  stringsAsFactors = FALSE
)

ibd_df$site1_name <- site_lookup$site_name[match(ibd_df$site1, site_lookup$site_id)]
ibd_df$site2_name <- site_lookup$site_name[match(ibd_df$site2, site_lookup$site_id)]

# -----------------------------
# 6) sampling map
# tight regional extent around points
# -----------------------------
state_map <- map_data("state")
plot_sites <- merge(ETCY_1_coords, site_lookup, by = "site_id", sort = FALSE)

lon_rng <- range(plot_sites$lon, na.rm = TRUE)
lat_rng <- range(plot_sites$lat, na.rm = TRUE)

x_pad <- max(0.12, diff(lon_rng) * 0.30)
y_pad <- max(0.10, diff(lat_rng) * 0.35)

xlim_use <- c(lon_rng[1] - x_pad, lon_rng[2] + x_pad)
ylim_use <- c(lat_rng[1] - y_pad, lat_rng[2] + y_pad)

map_plot <- ggplot() +
  geom_polygon(
    data = state_map,
    aes(x = long, y = lat, group = group),
    fill = "grey96",
    color = "grey55",
    linewidth = 0.25
  ) +
  geom_point(
    data = plot_sites,
    aes(x = lon, y = lat),
    size = 3
  ) +
  geom_text(
    data = plot_sites,
    aes(x = lon, y = lat, label = paste0(site_id, ". ", site_code)),
    nudge_y = 0.015,
    size = 3.5
  ) +
  coord_fixed(
    xlim = xlim_use,
    ylim = ylim_use
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "ETCY-1 sampling locations"
  )

print(map_plot)

ggsave(
  filename = file.path(base_dir, "ETCY-1_map.png"),
  plot = map_plot,
  width = 7,
  height = 5.25,
  dpi = 300
)

# -----------------------------
# 7) IBD plot
# -----------------------------
ibd_plot <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "ETCY-1 isolation by distance"
  )

print(ibd_plot)

ggsave(
  filename = file.path(base_dir, "ETCY-1_IBD.png"),
  plot = ibd_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# -----------------------------
# 8) checks
# -----------------------------
stopifnot(identical(rownames(ETCY_1_fst), ETCY_1_coords$site_id))
stopifnot(identical(colnames(ETCY_1_fst), ETCY_1_coords$site_id))
stopifnot(isTRUE(all.equal(ETCY_1_fst, t(ETCY_1_fst))))
stopifnot(isTRUE(all.equal(geo_dist_km, t(geo_dist_km))))

# -----------------------------
# 9) save
# -----------------------------
save(
  ETCY_1_fst,
  ETCY_1_coords,
  file = file.path(data_dir, "ETCY-1.RData")
)
