# -----------------------------
# ETSW-2 gulf darter
# Fluker et al. 2010 Conserv Genet 11:2267-2279
# site coordinates + FST matrix + IBD plot + save RData
# -----------------------------

library(ggplot2)
library(geosphere)
library(sf)
library(tigris)

options(tigris_use_cache = TRUE)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETSW-1"
data_dir <- file.path(save_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 0) site coordinates
# numeric site IDs are the matrix row/col names
# the paper did not provide exact XYs for the three Black Warrior
# E. swaini localities used in the microsatellite analysis, so these
# are representative stream coordinates inferred from the locality
# names and mapped positions in the paper.
# -----------------------------
ETSW_1_coords <- data.frame(
  site = as.character(1:3),
  site_code = c("WC", "WF", "MC"),
  site_name = c(
    "Walker County Shoal Creek",
    "Wolf Creek",
    "Mill Creek"
  ),
  lat = c(
    33.53261,   # Walker County Shoal Creek
    33.70880,   # representative Wolf Creek point near Oakman
    33.17472    # representative Mill Creek point near Tuscaloosa/Northport
  ),
  lon = c(
    -87.26333,
    -87.47770,
    -87.56472
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map of sampling locations
# -----------------------------
al_sf <- tigris::states(cb = TRUE, year = 2024, class = "sf")
al_sf <- al_sf[al_sf$STUSPS == "AL", ]

xpad <- 0.20
ypad <- 0.15

ETSW_1_map <- ggplot() +
  geom_sf(data = al_sf, fill = "grey95", color = "grey40", linewidth = 0.4) +
  geom_point(
    data = ETSW_1_coords,
    aes(x = lon, y = lat),
    size = 3
  ) +
  geom_text(
    data = ETSW_1_coords,
    aes(x = lon, y = lat, label = site_name),
    nudge_y = 0.02,
    size = 3.5
  ) +
  coord_sf(
    xlim = range(ETSW_1_coords$lon) + c(-xpad, xpad),
    ylim = range(ETSW_1_coords$lat) + c(-ypad, ypad),
    expand = FALSE
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "ETSW-2 sampling locations"
  )

print(ETSW_1_map)

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
ETSW_1_geo_dist_km <- geosphere::distm(
  as.matrix(ETSW_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(ETSW_1_geo_dist_km) <- ETSW_1_coords$site
colnames(ETSW_1_geo_dist_km) <- ETSW_1_coords$site

# -----------------------------
# 3) helper function
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
# Table 4 from Fluker et al. 2010
# subset order: WC, WF, MC
# any negative values are forced to 0
# -----------------------------
ETSW_1_fst_vals <- c(
  0.058,
  0.083, 0.075
)

ETSW_1_fst <- fill_sym_from_lower(
  pops = ETSW_1_coords$site,
  vals = ETSW_1_fst_vals,
  diag_val = 0
)

ETSW_1_fst[ETSW_1_fst < 0] <- 0

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ETSW_1_ibd_df <- data.frame(
  site1 = rownames(ETSW_1_fst)[row(ETSW_1_fst)[upper.tri(ETSW_1_fst)]],
  site2 = colnames(ETSW_1_fst)[col(ETSW_1_fst)[upper.tri(ETSW_1_fst)]],
  fst = ETSW_1_fst[upper.tri(ETSW_1_fst)],
  dist_km = ETSW_1_geo_dist_km[upper.tri(ETSW_1_geo_dist_km)],
  stringsAsFactors = FALSE
)

ETSW_1_ibd_df$site1_name <- ETSW_1_coords$site_name[match(ETSW_1_ibd_df$site1, ETSW_1_coords$site)]
ETSW_1_ibd_df$site2_name <- ETSW_1_coords$site_name[match(ETSW_1_ibd_df$site2, ETSW_1_coords$site)]

# -----------------------------
# 6) IBD plot
# -----------------------------
ETSW_1_ibd_plot <- ggplot(ETSW_1_ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "ETSW-2 isolation by distance"
  )

print(ETSW_1_ibd_plot)

# -----------------------------
# 7) quick checks
# -----------------------------
stopifnot(identical(rownames(ETSW_1_fst), ETSW_1_coords$site))
stopifnot(identical(colnames(ETSW_1_fst), ETSW_1_coords$site))
stopifnot(isTRUE(all.equal(ETSW_1_fst, t(ETSW_1_fst))))
stopifnot(isTRUE(all.equal(ETSW_1_geo_dist_km, t(ETSW_1_geo_dist_km))))

# -----------------------------
# 8) save RData
# -----------------------------
save(
  ETSW_1_coords,
  ETSW_1_fst,
  file = file.path(data_dir, "ETSW-1.RData")
)
