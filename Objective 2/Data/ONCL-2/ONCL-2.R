# -----------------------------
# ONCL-2 cutthroat trout
# site coordinates + FST matrix
# -----------------------------

library(dplyr)
library(geosphere)
library(ggplot2)

# directory where everything is
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ONCL-2/"

# -----------------------------
# 0) site-level coordinates
# read final coordinates from csv
# expected columns: site, lat, lon
# -----------------------------
ONCL_2_coords <- read.csv(
  file.path(save_dir, "ONCL-2_xys.csv"),
  stringsAsFactors = FALSE
) %>%
  select(-WKT)

# order matches the FST matrix in the paper
fst_order <- c(
  "WMR", "EMR1", "EMR2", "MRBC1", "MRBC2", "QCK", "BC", "MS",
  "CC1", "CC2", "TC1", "TC2", "DC1", "DC2", "WC1", "WC2"
)

ONCL_2_coords <- ONCL_2_coords[match(fst_order, ONCL_2_coords$site), ]

# -----------------------------
# 1) helper: fill symmetric matrix
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
# 2) pairwise FST matrix
# values entered row-wise from lower triangle
# transcribed from Table 2
# -----------------------------
fst_vals <- c(
  # EMR1
  0.05,
  
  # EMR2
  0.06, 0.01,
  
  # MRBC1
  0.15, 0.17, 0.15,
  
  # MRBC2
  0.06, 0.02, 0.01, 0.14,
  
  # QCK
  0.08, 0.06, 0.06, 0.16, 0.06,
  
  # BC
  0.05, 0.03, 0.03, 0.14, 0.03, 0.07,
  
  # MS
  0.04, 0.01, 0.01, 0.13, 0.01, 0.05, 0.02,
  
  # CC1
  0.11, 0.10, 0.10, 0.18, 0.11, 0.12, 0.10, 0.08,
  
  # CC2
  0.08, 0.06, 0.06, 0.18, 0.06, 0.09, 0.06, 0.05, 0.03,
  
  # TC1
  0.13, 0.09, 0.07, 0.22, 0.19, 0.14, 0.09, 0.08, 0.15, 0.12,
  
  # TC2
  0.15, 0.12, 0.10, 0.24, 0.19, 0.16, 0.10, 0.10, 0.18, 0.14, 0.04,
  
  # DC1
  0.23, 0.20, 0.18, 0.32, 0.13, 0.23, 0.18, 0.18, 0.27, 0.22, 0.09, 0.05,
  
  # DC2
  0.22, 0.20, 0.18, 0.31, 0.13, 0.22, 0.17, 0.18, 0.25, 0.21, 0.10, 0.05, 0.01,
  
  # WC1
  0.15, 0.15, 0.14, 0.22, 0.13, 0.13, 0.15, 0.12, 0.20, 0.15, 0.14, 0.14, 0.16, 0.16,
  
  # WC2
  0.16, 0.15, 0.13, 0.23, 0.13, 0.14, 0.15, 0.12, 0.20, 0.16, 0.12, 0.14, 0.14, 0.15, 0.03
)

ONCL_2_fst <- fill_sym_from_lower(fst_order, fst_vals, diag_val = 0)

# -----------------------------
# IBD plot: straight-line distance vs FST
# -----------------------------
dist_km <- geosphere::distm(
  as.matrix(ONCL_2_coords[, c("lon", "lat")]),
  fun = geosphere::distGeo
) / 1000

ibd_df <- data.frame(
  dist = dist_km[upper.tri(dist_km)],
  fst  = ONCL_2_fst[upper.tri(ONCL_2_fst)]
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
stopifnot(identical(rownames(ONCL_2_fst), ONCL_2_coords$site))
stopifnot(identical(colnames(ONCL_2_fst), ONCL_2_coords$site))
stopifnot(isTRUE(all.equal(ONCL_2_fst, t(ONCL_2_fst))))

# -----------------------------
# 4) save RData
# -----------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  ONCL_2_fst,
  ONCL_2_coords,
  file = file.path(out_dir, "ONCL-2.RData")
)