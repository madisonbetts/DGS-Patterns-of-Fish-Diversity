# -----------------------------
# GIOR-1 Arroyo chub
# site coordinates + FST matrix
# -----------------------------

library(dplyr)
library(geosphere)
library(ggplot2)

# directory where everything is
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/GIOR-1/"

# -----------------------------
# 0) site-level coordinates
# -----------------------------
sites <- data.frame(
  watershed = c(
    rep("MC", 4),
    rep("LA", 3),
    rep("SG", 2),
    "SA",
    rep("SJ", 2),
    "SM"
  ),
  site = c(
    "Las Virgenes Creek",
    "Above Serra Road Bridge",
    "Above Rindge Dam",
    "Near Cross Creek Road Bridge",
    "Pacoima Canyon",
    "Big Tujunga Creek",
    "unnamed LA site",
    "West Fork",
    "Walnut Creek",
    "Santa Ana River",
    "Bell Canyon–Starr Ranch",
    "Hot Springs Creek",
    "Temecula Creek"
  ),
  pop = c(
    rep("MC", 4),
    "LA/PC",
    "LA/BTC",
    NA,
    "SG/WF",
    "SG/WC",
    "SA",
    "SJ",
    "SJ",
    "SM"
  ),
  lat = c(
    34.09680, 34.04722, 34.07640, 34.04539,
    34.34541, 34.29451, 34.30181,
    34.24319, 34.08722,
    34.03594,
    33.63169, 33.60814,
    33.43408
  ),
  lon = c(
    -118.72845, -118.68972, -118.70230, -118.68703,
    -118.35827, -118.24322, -118.25575,
    -117.87497, -117.84511,
    -117.35670,
    -117.55531, -117.51082,
    -116.85529
  )
)

# -----------------------------
# 1) population coordinates
# only keep populations used in Table 2
# -----------------------------
GIOR_1_coords <- sites %>%
  filter(!is.na(pop)) %>%
  group_by(pop) %>%
  summarise(
    lat = mean(lat),
    lon = mean(lon),
    n_sites = n(),
    .groups = "drop"
  ) %>%
  rename(site = pop)

fst_order <- c("MC", "LA/PC", "LA/BTC", "SG/WF", "SG/WC", "SA", "SJ", "SM")

GIOR_1_coords <- GIOR_1_coords[match(fst_order, GIOR_1_coords$site), ]

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
  0.354,
  0.167, 0.188,
  0.217, 0.246, 0.068,
  0.217, 0.235, 0.064, 0.070,
  0.163, 0.223, 0.048, 0.071, 0.086,
  0.189, 0.246, 0.058, 0.101, 0.073, 0.073,
  0.302, 0.400, 0.208, 0.238, 0.238, 0.215, 0.199
)

GIOR_1_fst <- fill_sym_from_lower(fst_order, fst_vals, diag_val = 0)

# -----------------------------
# IBD plot: straight-line distance vs FST
# -----------------------------
dist_km <- geosphere::distm(
  as.matrix(GIOR_1_coords[, c("lon", "lat")]),
  fun = geosphere::distGeo
) / 1000

ibd_df <- data.frame(
  dist = dist_km[upper.tri(dist_km)],
  fst  = GIOR_1_fst[upper.tri(GIOR_1_fst)]
)

ggplot(ibd_df, aes(x = dist, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(
    x = "Straight-line distance (km)",
    y = expression(F[ST]),
    title = "IBD: straight-line distance vs FST"
  )

# -----------------------------
# 4) quick checks
# -----------------------------
stopifnot(identical(rownames(GIOR_1_fst), GIOR_1_coords$site))
stopifnot(identical(colnames(GIOR_1_fst), GIOR_1_coords$site))
stopifnot(isTRUE(all.equal(GIOR_1_fst, t(GIOR_1_fst))))

# -----------------------------
# 5) save RData
# -----------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  GIOR_1_fst,
  GIOR_1_coords,
  file = file.path(out_dir, "GIOR-1.RData")
)