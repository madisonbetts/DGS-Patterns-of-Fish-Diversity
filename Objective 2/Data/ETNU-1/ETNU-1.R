# -----------------------------
# ETNU-1 watercress darter
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
save_dir  <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETNU-1"
data_dir  <- file.path(save_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 0) site coordinates
# numeric site IDs are the matrix row/col names
# lat/lon were inferred from the study map plus web sources
# when exact spring-head coordinates were not given in the paper
# -----------------------------
ETNU_1_coords <- data.frame(
  site = as.character(1:5),
  site_code = c("GS", "TH", "SS", "RS", "TS"),
  site_name = c(
    "Glenn Spring",
    "Thomas Spring",
    "Seven Springs",
    "Roebuck Spring",
    "Tapawingo Spring"
  ),
  lat = c(
    33.38590,   # Glenn Spring; approximate from Thomas/Glenn Spring map
    33.37545,   # Thomas Spring / Watercress Darter NWR
    33.47206,   # Seven Springs EcoScape / Faith Apostolic Church area
    33.57700,   # Roebuck Spring / Vacca Campus-Roebuck Hawkins Park area
    33.69972    # Tapawingo Spring run at Turkey Creek confluence
  ),
  lon = c(
    -86.95250,
    -86.96222,
    -86.87015,
    -86.70350,
    -86.67555
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

ETNU_1_map <- ggplot() +
  geom_sf(data = al_sf, fill = "grey95", color = "grey40", linewidth = 0.4) +
  geom_point(
    data = ETNU_1_coords,
    aes(x = lon, y = lat),
    size = 3
  ) +
  geom_text(
    data = ETNU_1_coords,
    aes(x = lon, y = lat, label = site_name),
    nudge_y = 0.02,
    size = 3.5
  ) +
  coord_sf(
    xlim = range(ETNU_1_coords$lon) + c(-xpad, xpad),
    ylim = range(ETNU_1_coords$lat) + c(-ypad, ypad),
    expand = FALSE
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "ETNU-1 sampling locations"
  )

print(ETNU_1_map)

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
ETNU_1_geo_dist_km <- geosphere::distm(
  as.matrix(ETNU_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(ETNU_1_geo_dist_km) <- ETNU_1_coords$site
colnames(ETNU_1_geo_dist_km) <- ETNU_1_coords$site

# -----------------------------
# 3) helper function
# takes lower triangle values row-wise and rebuilds a
# full symmetric matrix
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
# order: GS, TH, SS, RS, TS
# any negative values are forced to 0
# -----------------------------
ETNU_1_fst_vals <- c(
  0.051,
  0.208, 0.221,
  0.261, 0.269, 0.274,
  0.268, 0.277, 0.286, 0.000
)

ETNU_1_fst <- fill_sym_from_lower(
  pops = ETNU_1_coords$site,
  vals = ETNU_1_fst_vals,
  diag_val = 0
)

ETNU_1_fst[ETNU_1_fst < 0] <- 0

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ETNU_1_ibd_df <- data.frame(
  site1 = rownames(ETNU_1_fst)[row(ETNU_1_fst)[upper.tri(ETNU_1_fst)]],
  site2 = colnames(ETNU_1_fst)[col(ETNU_1_fst)[upper.tri(ETNU_1_fst)]],
  fst = ETNU_1_fst[upper.tri(ETNU_1_fst)],
  dist_km = ETNU_1_geo_dist_km[upper.tri(ETNU_1_geo_dist_km)],
  stringsAsFactors = FALSE
)

ETNU_1_ibd_df$site1_name <- ETNU_1_coords$site_name[match(ETNU_1_ibd_df$site1, ETNU_1_coords$site)]
ETNU_1_ibd_df$site2_name <- ETNU_1_coords$site_name[match(ETNU_1_ibd_df$site2, ETNU_1_coords$site)]

# -----------------------------
# 6) IBD plot
# -----------------------------
ETNU_1_ibd_plot <- ggplot(ETNU_1_ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "ETNU-1 isolation by distance"
  )

print(ETNU_1_ibd_plot)

# -----------------------------
# 7) quick checks
# -----------------------------
stopifnot(identical(rownames(ETNU_1_fst), ETNU_1_coords$site))
stopifnot(identical(colnames(ETNU_1_fst), ETNU_1_coords$site))
stopifnot(isTRUE(all.equal(ETNU_1_fst, t(ETNU_1_fst))))
stopifnot(isTRUE(all.equal(ETNU_1_geo_dist_km, t(ETNU_1_geo_dist_km))))

# -----------------------------
# 8) save RData
# -----------------------------
save(
  ETNU_1_coords,
  ETNU_1_fst,
  #ETNU_1_geo_dist_km,
  #ETNU_1_ibd_df,
  file = file.path(data_dir, "ETNU-1.RData")
)
