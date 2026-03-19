# ========================================
# Ictalurus punctatus
# ICPU-1
# Wabash + Ohio rivers
# site coordinates + FST matrix + IBD plot
# ========================================

library(ggplot2)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ICPU-1/"

# ------------------------------
# 0) site coordinates
# order matches the published FST matrix
# ------------------------------
ICPU_1_coords <- data.frame(
  site = c("WA_268", "WA_189", "WA_136", "WA_19",
           "OH_1478", "OH_1511", "OH_1548"),
  lat = c(39.12247, 38.64305, 38.29545, 37.90736,
          37.23654, 37.12979, 37.22777),
  lon = c(-87.64498, -87.61680, -87.83276, -88.09838,
          -88.47449, -88.42943, -88.95386),
  stringsAsFactors = FALSE
)

# ------------------------------
# 1) helpers
# fill symmetric matrix from row-wise triangle values
# ------------------------------
fill_sym_from_upper <- function(pops, vals, diag_val = 0) {
  n <- length(pops)
  stopifnot(length(vals) == n * (n - 1) / 2)
  
  mat <- matrix(0, nrow = n, ncol = n, dimnames = list(pops, pops))
  
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

# ------------------------------
# 2) pairwise FST matrix
# values from Table 2, entered row-wise from the upper triangle
# ------------------------------
fst_vals <- c(
  0.009, 0.072, 0.123, 0.173, 0.157, 0.184,
  0.029, 0.080, 0.142, 0.102, 0.123,
  0.018, 0.069, 0.049, 0.074,
  0.066, 0.050, 0.087,
  0.021, 0.172,
  0.113
)

ICPU_1_fst <- fill_sym_from_upper(ICPU_1_coords$site, fst_vals, diag_val = 0)

# ------------------------------
# 3) river distance matrix (km)
# values from Table 2, entered row-wise from the lower triangle
# ------------------------------
river_dist_vals <- c(
  79,
  132, 53,
  249, 170, 117,
  379, 299, 246, 129,
  412, 332, 279, 162, 33,
  449, 369, 316, 199, 70, 37
)

river_dist_km <- fill_sym_from_lower(ICPU_1_coords$site, river_dist_vals, diag_val = 0)

# ------------------------------
# 4) IBD plot dataframe
# ------------------------------
ibd_df <- data.frame(
  site1 = rownames(ICPU_1_fst)[row(ICPU_1_fst)[upper.tri(ICPU_1_fst)]],
  site2 = colnames(ICPU_1_fst)[col(ICPU_1_fst)[upper.tri(ICPU_1_fst)]],
  distance_km = river_dist_km[upper.tri(river_dist_km)],
  fst = ICPU_1_fst[upper.tri(ICPU_1_fst)]
)

# ------------------------------
# 5) map of sampling locations
# ------------------------------
ggplot(ICPU_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.08, size = 3.5) +
  coord_equal() +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "ICPU-1 sampling locations"
  )

# ------------------------------
# 6) IBD plot
# ------------------------------
ggplot(ibd_df, aes(x = distance_km, y = fst)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(
    x = "River distance (km)",
    y = expression(F[ST]),
    title = "ICPU-1 isolation by distance"
  )

# ------------------------------
# 7) quick checks
# ------------------------------
stopifnot(identical(rownames(ICPU_1_fst), ICPU_1_coords$site))
stopifnot(identical(colnames(ICPU_1_fst), ICPU_1_coords$site))
stopifnot(isTRUE(all.equal(ICPU_1_fst, t(ICPU_1_fst))))
stopifnot(isTRUE(all.equal(river_dist_km, t(river_dist_km))))

# ------------------------------
# 8) save RData
# ------------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  ICPU_1_fst,
  ICPU_1_coords,
  file = file.path(out_dir, "ICPU-1.RData")
)