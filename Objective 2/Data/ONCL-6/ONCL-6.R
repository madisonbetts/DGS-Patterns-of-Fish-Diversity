# -----------------------------
# ONCL-6: westslope cutthroat trout
# Oncorhynchus clarkii lewisi
# provisional workflow from pairwise FST heatmap
# -----------------------------

library(ggplot2)
library(geosphere)
library(maps)

# -----------------------------
# 0) paths
# -----------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ONCL-6"
data_dir <- file.path(base_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) population order
# from the supplied FST heatmap
# -----------------------------
pop_codes <- c(
  "USKC_2019",
  "COR_2019",
  "COY_2019",
  "FVC_2019",
  "CAB_2021",
  "NFC_2017",
  "NFC_2020",
  "SHK_2021",
  "SWC_2021"
)

# -----------------------------
# 2) coordinates
# provisional basin/stream anchors linked to code labels
# final coords df only has: site, lat, lon
#
# site order matches pop_codes exactly
# 1 USKC_2019
# 2 COR_2019
# 3 COY_2019
# 4 FVC_2019
# 5 CAB_2021
# 6 NFC_2017
# 7 NFC_2020
# 8 SHK_2021
# 9 SWC_2021
# -----------------------------
ONCL_6_coords <- data.frame(
  site = as.character(1:9),
  lat = c(
    47.50,  # USKC
    47.67,  # COR
    47.55,  # COY
    48.10,  # FVC
    48.30,  # CAB
    46.55,  # NFC_2017
    46.60,  # NFC_2020
    47.10,  # SHK_2021
    47.75   # SWC_2021
  ),
  lon = c(
    -116.30,
    -116.75,
    -116.60,
    -114.20,
    -115.70,
    -115.30,
    -115.20,
    -115.80,
    -113.90
  ),
  stringsAsFactors = FALSE
)

# optional lookup table for your own editing
ONCL_6_lookup <- data.frame(
  site = as.character(1:9),
  pop_code = pop_codes,
  stringsAsFactors = FALSE
)

# -----------------------------
# 3) helper: fill symmetric matrix from lower triangle
# row-wise order:
# row2 has 1 value, row3 has 2 values, etc.
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
# transcribed from supplied heatmap, row by row
# lower triangle in population order above
# -----------------------------
fst_vals <- c(
  # row 2: COR_2019
  0.013,
  
  # row 3: COY_2019
  0.024, 0.023,
  
  # row 4: FVC_2019
  0.013, 0.036, 0.024,
  
  # row 5: CAB_2021
  0.028, 0.063, 0.062, 0.053,
  
  # row 6: NFC_2017
  0.062, 0.037, 0.040, 0.039, 0.030,
  
  # row 7: NFC_2020
  0.039, 0.055, 0.056, 0.054, 0.046, 0.023,
  
  # row 8: SHK_2021
  0.057, 0.065, 0.056, 0.062, 0.059, 0.048, 0.064,
  
  # row 9: SWC_2021
  0.046, 0.044, 0.045, 0.052, 0.062, 0.049, 0.052, 0.068
)

ONCL_6_fst <- fill_sym_from_lower(
  pops = ONCL_6_coords$site,
  vals = fst_vals,
  diag_val = 0
)

ONCL_6_fst[ONCL_6_fst < 0] <- 0

# -----------------------------
# 5) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(ONCL_6_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- ONCL_6_coords$site
colnames(geo_dist_km) <- ONCL_6_coords$site

# -----------------------------
# 6) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(ONCL_6_fst)[row(ONCL_6_fst)[upper.tri(ONCL_6_fst)]],
  site2   = colnames(ONCL_6_fst)[col(ONCL_6_fst)[upper.tri(ONCL_6_fst)]],
  fst     = ONCL_6_fst[upper.tri(ONCL_6_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 7) map
# include US + Canada, zoom to point extent
# plot only in RStudio; do not save
# -----------------------------
usa <- map_data("state")
can <- map_data("world", region = "Canada")

xpad <- 1.5
ypad <- 1.0

xlim_map <- c(
  min(ONCL_6_coords$lon) - xpad,
  max(ONCL_6_coords$lon) + xpad
)

ylim_map <- c(
  min(ONCL_6_coords$lat) - ypad,
  max(ONCL_6_coords$lat) + ypad
)

ggplot() +
  geom_polygon(
    data = usa,
    aes(x = long, y = lat, group = group),
    fill = "grey92", color = "black", linewidth = 0.25
  ) +
  geom_polygon(
    data = can,
    aes(x = long, y = lat, group = group),
    fill = "grey86", color = "black", linewidth = 0.25
  ) +
  geom_point(
    data = ONCL_6_coords,
    aes(x = lon, y = lat),
    size = 2.7
  ) +
  geom_text(
    data = ONCL_6_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.10,
    size = 3.3
  ) +
  coord_fixed(
    xlim = xlim_map,
    ylim = ylim_map
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "ONCL-6 sampling sites"
  )


# -----------------------------
# 8) IBD plot
# plot only in RStudio; do not save
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.8, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "ONCL-6 isolation by distance"
  )


# -----------------------------
# 9) quick checks
# -----------------------------
stopifnot(nrow(ONCL_6_coords) == 9)
stopifnot(identical(rownames(ONCL_6_fst), ONCL_6_coords$site))
stopifnot(identical(colnames(ONCL_6_fst), ONCL_6_coords$site))
stopifnot(isTRUE(all.equal(ONCL_6_fst, t(ONCL_6_fst))))
stopifnot(length(fst_vals) == nrow(ONCL_6_coords) * (nrow(ONCL_6_coords) - 1) / 2)

# -----------------------------
# 10) save RData
# -----------------------------
save(
  ONCL_6_fst,
  ONCL_6_coords,
  file = file.path(data_dir, "ONCL-6.RData")
)