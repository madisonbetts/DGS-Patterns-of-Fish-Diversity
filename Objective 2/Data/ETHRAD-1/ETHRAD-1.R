# ========================================
# ETHRAD-1
# Etheostoma radiosum
# Thesis-based extraction workflow
# Pairwise FST transcribed from thesis Table 3
# Site coordinates hardcoded from user-updated ETHRAD-1_sites.csv
# ========================================

# -----------------------------
# 0) setup
# -----------------------------
library(ggplot2)
library(geosphere)
library(maps)

study_code <- "ETHRAD-1"

base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETHRAD-1"
data_dir <- file.path(base_dir, "data")

dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) site coordinates
# hardcoded from user-updated ETHRAD-1_sites.csv
# site order matches thesis Table 3 exactly
# -----------------------------
ETHRAD_1_coords <- data.frame(
  site_id = c("1", "2", "3", "4", "5", "6", "7", "8"),
  lat = c(33.8184, 33.7833, 33.8525, 33.8025, 33.774, 34.0855, 34.2237, 34.3826),
  lon = c(-96.5567, -96.5, -95.5444, -96.6647, -96.553, -96.7319, -95.4148, -96.1067),
  stringsAsFactors = FALSE
)

site_lookup <- data.frame(
  site_id = c("1", "2", "3", "4", "5", "6", "7", "8"),
  site_code = c("TX01", "TX02", "TX03", "TX04", "TX05", "OK01", "OK02", "OK03"),
  site_name = c("Lower Shawnee Creek", "Duck Creek", "Sanders Creek", "Little Mineral Creek", "Upper Shawnee Creek", "Glasses Creek", "Rock Creek", "Sandy Creek"),
  n_ind = c(30, 30, 13, 28, 30, 29, 30, 30),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) helper function
# lower triangle entered row-wise
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
# transcribed from thesis Table 3
# -----------------------------
fst_vals <- c(
  0.100,
  0.123, 0.025,
  0.153, 0.040, 0.050,
  0.108, 0.023, 0.054, 0.037,
  0.152, 0.068, 0.097, 0.121, 0.035,
  0.071, 0.071, 0.071, 0.095, 0.070, 0.118,
  0.155, 0.037, 0.023, 0.036, 0.058, 0.129, 0.088
)

ETHRAD_1_fst <- fill_sym_from_lower(
  pops = ETHRAD_1_coords$site_id,
  vals = fst_vals,
  diag_val = 0
)

ETHRAD_1_fst[ETHRAD_1_fst < 0] <- 0
diag(ETHRAD_1_fst) <- 0

# -----------------------------
# 4) geographic distance matrix
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(ETHRAD_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- ETHRAD_1_coords$site_id
colnames(geo_dist_km) <- ETHRAD_1_coords$site_id

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1 = rownames(ETHRAD_1_fst)[row(ETHRAD_1_fst)[upper.tri(ETHRAD_1_fst)]],
  site2 = colnames(ETHRAD_1_fst)[col(ETHRAD_1_fst)[upper.tri(ETHRAD_1_fst)]],
  fst = ETHRAD_1_fst[upper.tri(ETHRAD_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  stringsAsFactors = FALSE
)

ibd_df$site1_name <- site_lookup$site_name[match(ibd_df$site1, site_lookup$site_id)]
ibd_df$site2_name <- site_lookup$site_name[match(ibd_df$site2, site_lookup$site_id)]

# -----------------------------
# 6) sampling map
# regional extent around the points with state outlines
# -----------------------------
state_map <- map_data("state")
plot_sites <- merge(ETHRAD_1_coords, site_lookup, by = "site_id", sort = FALSE)

lon_rng <- range(plot_sites$lon, na.rm = TRUE)
lat_rng <- range(plot_sites$lat, na.rm = TRUE)

x_pad <- max(0.35, diff(lon_rng) * 0.22)
y_pad <- max(0.25, diff(lat_rng) * 0.22)

xlim_use <- c(lon_rng[1] - x_pad, lon_rng[2] + x_pad)
ylim_use <- c(lat_rng[1] - y_pad, lat_rng[2] + y_pad)

map_plot <- ggplot() +
  geom_polygon(
    data = state_map,
    aes(x = long, y = lat, group = group),
    fill = "grey96",
    color = "grey55",
    linewidth = 0.30
  ) +
  geom_point(
    data = plot_sites,
    aes(x = lon, y = lat),
    size = 3
  ) +
  geom_text(
    data = plot_sites,
    aes(x = lon, y = lat, label = paste0(site_id, ". ", site_code)),
    nudge_y = 0.05,
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
    title = "ETHRAD-1 sampling locations"
  )

print(map_plot)

ggsave(
  filename = file.path(base_dir, "ETHRAD-1_map.png"),
  plot = map_plot,
  width = 7.5,
  height = 5.5,
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
    title = "ETHRAD-1 isolation by distance"
  )

print(ibd_plot)

ggsave(
  filename = file.path(base_dir, "ETHRAD-1_IBD.png"),
  plot = ibd_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# -----------------------------
# 8) quick checks
# -----------------------------
stopifnot(identical(rownames(ETHRAD_1_fst), ETHRAD_1_coords$site_id))
stopifnot(identical(colnames(ETHRAD_1_fst), ETHRAD_1_coords$site_id))
stopifnot(isTRUE(all.equal(ETHRAD_1_fst, t(ETHRAD_1_fst))))
stopifnot(isTRUE(all.equal(geo_dist_km, t(geo_dist_km))))

# -----------------------------
# 9) save RData
# -----------------------------
save(
  ETHRAD_1_fst,
  ETHRAD_1_coords,
  file = file.path(data_dir, "ETHRAD-1.RData")
)
