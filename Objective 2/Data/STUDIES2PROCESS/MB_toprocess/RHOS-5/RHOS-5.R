# ========================================
# RHOS-5
# Rhinichthys osculus
# Speckled dace
#
# Su et al. 2022 / dissertation chapter
# Population genomic analysis of the Speckled Dace species complex
#
# IMPORTANT:
# The published pairwise FST table available from the supplement
# is clearly legible for the 10 focal California / California-adjacent
# population units:
#   Santa Ana
#   Amargosa
#   Owens
#   Long Valley
#   Lahontan
#   Sacramento
#   Central California
#   Klamath
#   Butte Lake
#   Warner
#
# The supplement also includes a second panel for some non-California
# comparisons, but those values are not all legible in the uploaded
# article view. To avoid inventing values, this workflow builds the
# exact 10 x 10 published matrix that is fully recoverable from
# Supplemental Table S.2 panel A.
#
# Coordinates were NOT grabbed from figures for the final dataset.
# Exact GPS coordinates are available in Supplemental Table S.1.
# Because FST is reported at the population-unit level rather than
# at all 38 exact collection localities, the coordinates below are
# centroids calculated from the exact locality GPS coordinates in
# Supplemental Table S.1.
#
# Figure 1 was used only as a QC reference for spatial sanity.
# ========================================

library(ggplot2)
library(geosphere)
library(maps)

study_code <- "RHOS-5"

base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/RHOS-5"
data_dir <- file.path(base_dir, "data")

dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) exact localities from Supplemental Table S.1
# used to derive population-unit centroids
# -----------------------------
site_tbl <- data.frame(
  pop_unit = c(
    "Lahontan", "Lahontan", "Lahontan",
    "Warner", "Warner", "Warner", "Warner",
    "Butte Lake",
    "Klamath", "Klamath",
    "Sacramento", "Sacramento", "Sacramento", "Sacramento", "Sacramento", "Sacramento",
    "Santa Ana", "Santa Ana",
    "Central California", "Central California", "Central California",
    "Amargosa", "Amargosa", "Amargosa", "Amargosa", "Amargosa", "Amargosa",
    "Owens", "Owens", "Owens",
    "Long Valley"
  ),
  locality = c(
    "Poore Creek", "Susie Creek", "Martis Creek",
    "Hart Lake Honey Creek", "Clover Creek", "Twenty Mile Creek", "Foskett Spring",
    "Butte Creek",
    "Lone Pine Bar", "Klamath River near Happy Camp",
    "Ash Creek", "Bear Creek", "Thomas Creek", "Howard Creek", "Quartz Creek", "Dye Creek",
    "Cajon Creek", "Lytle Creek",
    "Brizziolari Creek", "Sisquoc River", "San Lorenzo River",
    "Beatty Oasis", "Valley-CoffR-Spring", "Amargosa Canyon", "Bradford Spring 1", "Bradford Spring 2", "Roger's Spring",
    "Rock Creek Ditch", "Matlick Ditch", "Fish Slough",
    "Whitmore Hot Springs"
  ),
  lat = c(
    38.34278, 40.81223, 39.34555,
    42.42111, 42.47500, 42.11384, 42.06944,
    40.56472,
    41.54883, 41.767682,
    41.16111, 41.11433, 42.30416, 42.30138, 42.29222, 40.07075,
    34.27475, 34.23899,
    35.28000, 34.83086, 36.99305,
    36.91398, 37.04113, 35.84956, 36.40214, 36.40124, 36.47915,
    37.45932, 37.37602, 37.47551,
    37.62977
  ),
  lon = c(
    -119.52881, -116.04667, -120.61777,
    -120.11000, -120.18000, -119.93161, -119.83764,
    -121.29138,
    -123.52476, -123.403174,
    -120.83000, -121.55361, -120.53250, -120.73160, -120.74527, -122.11222,
    -117.45261, -117.49892,
    -120.66889, -120.16056, -122.03277,
    -116.75551, -116.71248, -116.23078, -116.30236, -116.30304, -116.32646,
    -118.59280, -118.42450, -118.40095,
    -118.81017
  ),
  stringsAsFactors = FALSE
)

fst_order <- c(
  "Santa Ana",
  "Amargosa",
  "Owens",
  "Long Valley",
  "Lahontan",
  "Sacramento",
  "Central California",
  "Klamath",
  "Butte Lake",
  "Warner"
)

coords_centroid <- aggregate(
  cbind(lat, lon) ~ pop_unit,
  data = site_tbl,
  FUN = mean
)

coords_centroid <- coords_centroid[match(fst_order, coords_centroid$pop_unit), ]

RHOS_5_coords <- data.frame(
  site = as.character(seq_len(nrow(coords_centroid))),
  lat = coords_centroid$lat,
  lon = coords_centroid$lon,
  stringsAsFactors = FALSE
)

site_lookup <- data.frame(
  site = RHOS_5_coords$site,
  pop_unit = fst_order,
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) helper
# -----------------------------
fill_sym_from_lower <- function(pops, vals, diag_val = 0) {
  n <- length(pops)
  stopifnot(length(vals) == n * (n - 1) / 2)

  mat <- matrix(0, nrow = n, ncol = n, dimnames = list(pops, pops))
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
# exact values from Supplemental Table S.2 panel A
# order:
# Santa Ana, Amargosa, Owens, Long Valley, Lahontan,
# Sacramento, Central California, Klamath, Butte Lake, Warner
# -----------------------------
fst_vals <- c(
  0.683038,
  0.655648, 0.163662,
  0.495095, 0.382165, 0.300142,
  0.549713, 0.328260, 0.246113, 0.259542,
  0.645057, 0.555592, 0.511999, 0.443290, 0.412116,
  0.753472, 0.620699, 0.586420, 0.488401, 0.463870, 0.194572,
  0.549713, 0.571066, 0.535016, 0.456907, 0.428349, 0.165698, 0.326771,
  0.679071, 0.576948, 0.533901, 0.454898, 0.424014, 0.083905, 0.228633, 0.217210,
  0.679669, 0.559563, 0.507790, 0.418762, 0.388788, 0.288543, 0.429950, 0.319681, 0.329080
)

RHOS_5_fst <- fill_sym_from_lower(
  pops = RHOS_5_coords$site,
  vals = fst_vals,
  diag_val = 0
)

RHOS_5_fst[RHOS_5_fst < 0] <- 0
diag(RHOS_5_fst) <- 0

# -----------------------------
# 4) geographic distance matrix
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(RHOS_5_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- RHOS_5_coords$site
colnames(geo_dist_km) <- RHOS_5_coords$site

# -----------------------------
# 5) pairwise dataframe for IBD
# -----------------------------
ibd_df <- data.frame(
  site1 = rownames(RHOS_5_fst)[row(RHOS_5_fst)[upper.tri(RHOS_5_fst)]],
  site2 = colnames(RHOS_5_fst)[col(RHOS_5_fst)[upper.tri(RHOS_5_fst)]],
  fst = RHOS_5_fst[upper.tri(RHOS_5_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  stringsAsFactors = FALSE
)

ibd_df$site1_name <- site_lookup$pop_unit[match(ibd_df$site1, site_lookup$site)]
ibd_df$site2_name <- site_lookup$pop_unit[match(ibd_df$site2, site_lookup$site)]

# -----------------------------
# 6) sampling map
# -----------------------------
state_map <- map_data("state")
plot_sites <- merge(RHOS_5_coords, site_lookup, by = "site", sort = FALSE)

lon_rng <- range(plot_sites$lon, na.rm = TRUE)
lat_rng <- range(plot_sites$lat, na.rm = TRUE)

x_pad <- max(1.2, diff(lon_rng) * 0.08)
y_pad <- max(0.8, diff(lat_rng) * 0.08)

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
    size = 2.8
  ) +
  geom_text(
    data = plot_sites,
    aes(x = lon, y = lat, label = paste0(site, ". ", pop_unit)),
    nudge_y = 0.25,
    size = 3.0
  ) +
  coord_fixed(
    xlim = xlim_use,
    ylim = ylim_use
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "RHOS-5 sampling locations"
  )

print(map_plot)

ggsave(
  filename = file.path(base_dir, "RHOS-5_map.png"),
  plot = map_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# -----------------------------
# 7) IBD plot
# -----------------------------
ibd_plot <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.7, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "RHOS-5 isolation by distance"
  )

print(ibd_plot)

ggsave(
  filename = file.path(base_dir, "RHOS-5_IBD.png"),
  plot = ibd_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# -----------------------------
# 8) checks + save
# -----------------------------
stopifnot(identical(rownames(RHOS_5_fst), RHOS_5_coords$site))
stopifnot(identical(colnames(RHOS_5_fst), RHOS_5_coords$site))
stopifnot(isTRUE(all.equal(RHOS_5_fst, t(RHOS_5_fst))))
stopifnot(isTRUE(all.equal(geo_dist_km, t(geo_dist_km))))

save(
  RHOS_5_fst,
  RHOS_5_coords,
  file = file.path(data_dir, "RHOS-5.RData")
)
