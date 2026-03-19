# -----------------------------
# ETNI-1 Johnny darter
# site coordinates + FST matrix
# -----------------------------

# -------------------------
# 0) site metadata
# order matches the FST matrix
# -------------------------
ETNI_1_coords <- data.frame(
  site = c(
    "TOM", "BIG", "PM", "WSL", "SOL", "LCO", "MC", "BEA",
    "LYN", "PRC", "BN", "BUT", "NOQ", "NEB", "MCK", "HMO"
  ),
  #lake = c(
  #  "Tomahawk Lake", "Big Lake", "Plum Lake", "White Sand Lake",
  #  "Solberg Lake", "Lac Courte Oreilles", "Moon Chain", "Bearskin Lake",
  #  "Lynx Lake", "Pike/Round Chain", "Butternut Lake", "Butternut Lake",
  #  "Lake Noquebay", "Lake Nebagamon", "McKenzie Lake", "Half Moon Lake"
  #),
  lat = c(
    45.830, 46.154, 46.004, 46.008,
    45.749, 45.892, 45.660, 45.731,
    46.195, 45.928, 45.983, 45.923,
    45.249, 46.513, 45.936, 45.498
  ),
  lon = c(
    -89.661, -89.767, -89.514, -89.827,
    -90.369, -91.438, -89.304, -89.685,
    -89.666, -90.067, -90.515, -88.986,
    -87.924, -91.703, -92.046, -92.438
  ),
  stringsAsFactors = FALSE
)

# -------------------------
# 1) helper: fill symmetric matrix
#    from row-wise lower triangle
# -------------------------
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

# -------------------------
# 2) pairwise FST matrix
# order follows Table 12
# -------------------------
fst_lower <- c(
  0.24,
  0.06, 0.19,
  0.20, 0.08, 0.16,
  0.24, 0.11, 0.20, 0.15,
  0.16, 0.11, 0.11, 0.09, 0.13,
  0.12, 0.31, 0.11, 0.28, 0.31, 0.21,
  0.08, 0.21, 0.07, 0.18, 0.23, 0.15, 0.11,
  0.36, 0.29, 0.32, 0.30, 0.33, 0.27, 0.40, 0.33,
  0.14, 0.10, 0.11, 0.13, 0.12, 0.11, 0.22, 0.14, 0.33,
  0.22, 0.03, 0.17, 0.06, 0.10, 0.09, 0.28, 0.18, 0.26, 0.10,
  0.17, 0.33, 0.18, 0.29, 0.35, 0.22, 0.12, 0.16, 0.39, 0.27, 0.30,
  0.18, 0.24, 0.16, 0.22, 0.23, 0.18, 0.18, 0.13, 0.40, 0.16, 0.21, 0.24,
  0.24, 0.24, 0.18, 0.21, 0.27, 0.16, 0.26, 0.17, 0.41, 0.21, 0.21, 0.30, 0.20,
  0.22, 0.27, 0.19, 0.25, 0.28, 0.17, 0.24, 0.17, 0.43, 0.21, 0.25, 0.28, 0.19, 0.11,
  0.23, 0.18, 0.20, 0.13, 0.25, 0.11, 0.25, 0.21, 0.34, 0.21, 0.17, 0.26, 0.22, 0.18, 0.19
)

ETNI_1_fst <- fill_sym_from_lower(ETNI_1_coords$site, fst_lower, diag_val = 0)

# -----------------------------
# IBD plot: straight-line distance vs FST
# -----------------------------
coords <- ETNI_1_coords[, c("lon", "lat")]

dist_km <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
dist_km <- as.matrix(dist_km)

ibd_df <- data.frame(
  dist = dist_km[upper.tri(dist_km)],
  fst  = ETNI_1_fst[upper.tri(ETNI_1_fst)]
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
# 3) quick checks
# -----------------------------
stopifnot(identical(rownames(ETNI_1_fst), ETNI_1_coords$site))
stopifnot(identical(colnames(ETNI_1_fst), ETNI_1_coords$site))
stopifnot(isTRUE(all.equal(ETNI_1_fst, t(ETNI_1_fst))))
stopifnot(all.equal(ETNI_1_fst["BIG", "TOM"], 0.24))
stopifnot(all.equal(ETNI_1_fst["PM", "BIG"], 0.19))

# -----------------------------
# 4) save RData
# -----------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETNI-1/data"

if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

save(
  ETNI_1_fst,
  ETNI_1_coords,
  file = file.path(save_dir, "ETNI-1.RData")
)