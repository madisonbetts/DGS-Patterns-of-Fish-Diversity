# -----------------------------
# RHOS-2 speckled dace (Rhinichthys osculus complex)
# Oregon Great Basin populations from Hoekzema & Sidlauskas 2014
# Population coordinates + pairwise FST matrix + IBD plot
# -----------------------------

library(ggplot2)
library(ggrepel)
library(geosphere)
library(sf)
library(tigris)

options(tigris_use_cache = TRUE)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/RHOS-2"
data_dir <- file.path(save_dir, "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 0) site / population coordinates
# -----------------------------
# Table 2 reports pairwise FST among 7 populations inferred by STRUCTURE:
# 1 = Malheur
# 2 = Stinking Lake Spring
# 3 = Silver Lake
# 4 = Goose Lake
# 5 = Warner
# 6 = Foskett Spring
# 7 = Lake Abert
#
# Exact x/y were not provided for the basin-level populations in the paper.
# Because the genetic matrix is population-level rather than stream-level,
# coordinates below are representative population centroids / anchors for IBD.
# Foskett Spring uses the reported spring coordinates.
# Stinking Lake Spring uses the published specimen/locality coordinate.
# The remaining populations use basin-centered approximations informed by the
# study map and named stream sets in Table A.1.

RHOS_2_coords <- data.frame(
  site = as.character(1:7),
  population = c(
    "Malheur",
    "Stinking Lake Spring",
    "Silver Lake",
    "Goose Lake",
    "Warner",
    "Foskett Spring",
    "Lake Abert"
  ),
  lat = c(
    43.1600,   # Malheur basin stream-population centroid (approx.)
    43.3266,   # Stinking Lake Spring
    43.0200,   # Silver Lake basin centroid (approx.)
    42.1200,   # Goose Lake basin centroid / Drews Creek area anchor (approx.)
    42.4200,   # Warner basin centroid between Twentymile + Snyder (approx.)
    42.06975,  # Foskett Spring
    42.5600    # Lake Abert basin centroid (approx.)
  ),
  lon = c(
    -118.9800,  # Malheur basin stream-population centroid (approx.)
    -119.3663,  # Stinking Lake Spring
    -121.1800,  # Silver Lake basin centroid (approx.)
    -120.6200,  # Goose Lake basin centroid / Drews Creek area anchor (approx.)
    -120.0500,  # Warner basin centroid (approx.)
    -119.839666,# Foskett Spring
    -120.3300   # Lake Abert basin centroid (approx.)
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) relevant state(s) map + sampling locations
# -----------------------------
or_state <- tigris::states(year = 2023, class = "sf")
or_state <- or_state[or_state$STUSPS == "OR", ]

RHOS_2_coords_sf <- st_as_sf(
  RHOS_2_coords,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

p_map <- ggplot() +
  geom_sf(data = or_state, fill = "grey97", color = "grey40", linewidth = 0.4) +
  geom_sf(data = RHOS_2_coords_sf, size = 2.8) +
  ggrepel::geom_text_repel(
    data = RHOS_2_coords,
    aes(x = lon, y = lat, label = paste0(site, ": ", population)),
    size = 3.4,
    min.segment.length = 0,
    seed = 1
  ) +
  coord_sf(
    xlim = c(-122.2, -117.8),
    ylim = c(41.7, 44.0),
    expand = FALSE
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "RHOS-2 sampling populations"
  )

print(p_map)

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(RHOS_2_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- RHOS_2_coords$site
colnames(geo_dist_km) <- RHOS_2_coords$site

# -----------------------------
# 3) helper function:
# takes a vector of lower-triangle values (row-wise order)
# and reconstructs a full symmetric matrix
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
# 4) pairwise FST matrix
# Values transcribed from Table 2 (below diagonal)
# Order:
# 1 Malheur
# 2 Stinking Lake Spring
# 3 Silver Lake
# 4 Goose Lake
# 5 Warner
# 6 Foskett Spring
# 7 Lake Abert
# -----------------------------
fst_vals <- c(
  # 2 vs 1
  0.109,

  # 3 vs 1:2
  0.088, 0.081,

  # 4 vs 1:3
  0.069, 0.086, 0.004,

  # 5 vs 1:4
  0.097, 0.090, 0.076, 0.067,

  # 6 vs 1:5
  0.114, 0.109, 0.097, 0.083, 0.060,

  # 7 vs 1:6
  0.117, 0.101, 0.066, 0.061, 0.071, 0.180
)

RHOS_2_fst <- fill_sym_from_lower(
  pops = RHOS_2_coords$site,
  vals = fst_vals,
  diag_val = 0
)

# always make fst < 0 equal to 0
RHOS_2_fst[RHOS_2_fst < 0] <- 0

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
RHOS_2_ibd_df <- data.frame(
  site1 = rownames(RHOS_2_fst)[row(RHOS_2_fst)[upper.tri(RHOS_2_fst)]],
  site2 = colnames(RHOS_2_fst)[col(RHOS_2_fst)[upper.tri(RHOS_2_fst)]],
  fst = RHOS_2_fst[upper.tri(RHOS_2_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  stringsAsFactors = FALSE
)

RHOS_2_ibd_df$pop1 <- RHOS_2_coords$population[match(RHOS_2_ibd_df$site1, RHOS_2_coords$site)]
RHOS_2_ibd_df$pop2 <- RHOS_2_coords$population[match(RHOS_2_ibd_df$site2, RHOS_2_coords$site)]

# -----------------------------
# 6) IBD plot
# -----------------------------
p_ibd <- ggplot(RHOS_2_ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "RHOS-2 isolation by distance"
  )

print(p_ibd)

# -----------------------------
# 7) quick checks
# -----------------------------
stopifnot(identical(rownames(RHOS_2_fst), RHOS_2_coords$site))
stopifnot(identical(colnames(RHOS_2_fst), RHOS_2_coords$site))
stopifnot(isTRUE(all.equal(RHOS_2_fst, t(RHOS_2_fst))))

# -----------------------------
# 8) save objects
# -----------------------------

RHOS_2_coords <- RHOS_2_coords[, c("site", "lat", "lon")] # remove pop col
save(
  RHOS_2_fst,
  RHOS_2_coords,
  file = file.path(data_dir, "RHOS-2.RData")
)
