# -----------------------------
# ETRA-1 Yazoo darter
# site coordinates + FST matrix
# -----------------------------

library(dplyr)
library(geosphere)
library(ggplot2)

# directory where everything is
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETRA-1/"

# -----------------------------
# 0) site-level coordinates
# -----------------------------
sites <- data.frame(
  site = 1:17,
  pop  = c(1, 2, 3, 4, 4, 5, 6, 6, 6, 7, 8, 9, 10, 11, 11, 11, 11),
  stream = c(
    "Yellow Rabbit Creek",
    "South Chili Creek",
    "Tippah River Tributary",
    "Chewalla Creek Tributary upstream of dam",
    "Chewalla Creek Tributary downstream of dam",
    "Big Spring Creek Tributary",
    "Puskus Creek upstream of dam",
    "Puskus Creek upstream of dam",
    "Puskus Creek upstream of dam",
    "Puskus Creek downstream of dam",
    "Cypress Creek",
    "Taylor Creek Tributary",
    "Morris Creek",
    "Johnson Creek",
    "Otoucalofa Creek Tributary",
    "Gordon Branch",
    "Mill Creek"
  ),
  lat = c(
    34.819, 34.682, 34.708, 34.760, 34.725,
    34.663, 34.395, 34.428, 34.431, 34.445,
    34.393, 34.123, 34.282, 34.123, 34.125,
    34.140, 34.166
  ),
  lon = -c(
    89.105, 89.172, 89.255, 89.332, 89.305,
    89.412, 89.372, 89.394, 89.375, 89.336,
    89.286, 89.641, 89.543, 89.641, 89.610,
    89.549, 89.520
  )
)

# -----------------------------
# 1) population coordinates
# use population IDs as site names
# order matches the FST matrix
# -----------------------------
ETRA_1_coords <- sites %>%
  group_by(pop) %>%
  summarise(
    lat = mean(lat),
    lon = mean(lon),
    n_sites = n(),
    .groups = "drop"
  ) %>%
  mutate(site = as.character(pop)) %>%
  select(site, lat, lon, n_sites)

fst_order <- c("11","10","9","5","4","1","3","2","8","7","6")

ETRA_1_coords <- ETRA_1_coords[match(fst_order, ETRA_1_coords$site), ]

# -----------------------------
# 2) helper: fill symmetric matrix
# from row-wise lower triangle
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
# values entered row-wise from lower triangle
# -----------------------------
fst_vals <- c(
  0.219,
  0.199, 0.170,
  0.280, 0.267, 0.225,
  0.248, 0.229, 0.172, 0.104,
  0.240, 0.201, 0.165, 0.091, 0.040,
  0.254, 0.228, 0.187, 0.102, 0.063, 0.034,
  0.287, 0.276, 0.221, 0.063, 0.086, 0.080, 0.109,
  0.293, 0.285, 0.262, 0.103, 0.095, 0.097, 0.124, 0.114,
  0.271, 0.260, 0.211, 0.069, 0.079, 0.078, 0.085, 0.070, 0.051,
  0.281, 0.265, 0.223, 0.094, 0.010, 0.083, 0.117, 0.088, 0.064, 0.047
)

ETRA_1_fst <- fill_sym_from_lower(fst_order, fst_vals, diag_val = 0)

# -----------------------------
# IBD plot: straight-line distance vs FST
# -----------------------------
dist_km <- geosphere::distm(
  as.matrix(ETRA_1_coords[, c("lon", "lat")]),
  fun = geosphere::distGeo
) / 1000

ibd_df <- data.frame(
  dist = dist_km[upper.tri(dist_km)],
  fst  = ETRA_1_fst[upper.tri(ETRA_1_fst)]
)

ggplot(ibd_df, aes(x = dist, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Geographic distance (km)",
    y = expression(F[ST]),
    title = "IBD: straight-line distance vs FST"
  )

# -----------------------------
# 4) quick checks
# -----------------------------
stopifnot(identical(rownames(ETRA_1_fst), ETRA_1_coords$site))
stopifnot(identical(colnames(ETRA_1_fst), ETRA_1_coords$site))
stopifnot(isTRUE(all.equal(ETRA_1_fst, t(ETRA_1_fst))))

# -----------------------------
# 5) save RData
# -----------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  ETRA_1_fst,
  ETRA_1_coords,
  file = file.path(out_dir, "ETRA-1.RData")
)