# ========================================
# ETHZON-1
# Etheostoma zonistium
# Bandfin darter
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

study_code <- "ETHZON-1"

base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETHZON-1"
data_dir <- file.path(base_dir, "data")

dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) site coordinates
# exact coordinates from Table 1
# site order matches Table 4
# 11 Bear Creek
# 12 Big Sandy River
# 13 Birdsong Creek
# 14 Hatchie River / West Fork Spring Creek
# -----------------------------
ETHZON_1_coords <- data.frame(
  site_id = as.character(1:4),
  lat = c(
    34.83091,
    35.80539,
    35.91634,
    35.08152
  ),
  lon = c(
    -88.05857,
    -88.33439,
    -88.11398,
    -89.11192
  ),
  stringsAsFactors = FALSE
)

site_lookup <- data.frame(
  site_id = as.character(1:4),
  site_code = c("11", "12", "13", "14"),
  site_name = c(
    "Bear Creek",
    "Big Sandy River",
    "Birdsong Creek",
    "Hatchie River"
  ),
  locality = c(
    "Bear Creek",
    "Big Sandy River",
    "Birdsong Creek",
    "West Fork Spring Creek"
  ),
  n_ind = c(14, 23, 23, 23),
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
# -----------------------------
fst_vals <- c(
  0.095,
  0.102, 0.113,
  0.140, 0.180, 0.160
)

ETHZON_1_fst <- fill_sym_from_lower(
  pops = ETHZON_1_coords$site_id,
  vals = fst_vals,
  diag_val = 0
)

ETHZON_1_fst[ETHZON_1_fst < 0] <- 0
diag(ETHZON_1_fst) <- 0

# -----------------------------
# 4) geographic distance matrix
# Euclidean great-circle distances in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(ETHZON_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- ETHZON_1_coords$site_id
colnames(geo_dist_km) <- ETHZON_1_coords$site_id

# -----------------------------
# 5) IBD dataframe
# -----------------------------
ibd_df <- data.frame(
  site1 = rownames(ETHZON_1_fst)[row(ETHZON_1_fst)[upper.tri(ETHZON_1_fst)]],
  site2 = colnames(ETHZON_1_fst)[col(ETHZON_1_fst)[upper.tri(ETHZON_1_fst)]],
  fst = ETHZON_1_fst[upper.tri(ETHZON_1_fst)],
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
plot_sites <- merge(ETHZON_1_coords, site_lookup, by = "site_id", sort = FALSE)

lon_rng <- range(plot_sites$lon, na.rm = TRUE)
lat_rng <- range(plot_sites$lat, na.rm = TRUE)

x_pad <- max(0.18, diff(lon_rng) * 0.20)
y_pad <- max(0.18, diff(lat_rng) * 0.20)

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
    title = "ETHZON-1 sampling locations"
  )

print(map_plot)

ggsave(
  filename = file.path(base_dir, "ETHZON-1_map.png"),
  plot = map_plot,
  width = 7.25,
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
    title = "ETHZON-1 isolation by distance"
  )

print(ibd_plot)

ggsave(
  filename = file.path(base_dir, "ETHZON-1_IBD.png"),
  plot = ibd_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# -----------------------------
# 8) checks
# -----------------------------
stopifnot(identical(rownames(ETHZON_1_fst), ETHZON_1_coords$site_id))
stopifnot(identical(colnames(ETHZON_1_fst), ETHZON_1_coords$site_id))
stopifnot(isTRUE(all.equal(ETHZON_1_fst, t(ETHZON_1_fst))))
stopifnot(isTRUE(all.equal(geo_dist_km, t(geo_dist_km))))

# -----------------------------
# 9) save
# -----------------------------
save(
  ETHZON_1_fst,
  ETHZON_1_coords,
  file = file.path(data_dir, "ETHZON-1.RData")
)
