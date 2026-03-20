# -----------------------------
# GAAC-1 threespine stickleback
# Oregon sites only
# site coordinates + FST matrix
# -----------------------------

library(dplyr)
library(geosphere)
library(ggplot2)

# directory where everything is
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/GAAC-1/"

# -----------------------------
# site key
# Oregon populations only
# Alaska sites removed: Ra, BP, MiF, MiO
# -----------------------------
# C     = Columbia River
# Sz    = Millport Slough
# Si    = South Jetty
# SiC   = Cushman Slough
# U     = Dean Creek
# WC    = Winchester Creek
# EC    = Eel Creek
# LL    = Lily Lake
# PC    = Pony Creek Reservoir
# U196  = Page Road
# Mo14  = Milk Creek
# Sa39  = Buell-Miller Slough
# L34   = Jont Creek
# Ma32  = Finley Grey Creek Swamp
# W296  = Science Factory
# W320  = Dougren Slough
# CF27  = Lynx Hollow Slough
# W278  = Green Island
# Mc18  = Riverbend
# Mc42  = Walterville Slough
# Mc61  = Leaburg Fish Hatchery
# CR    = Crooked River
# ST    = South Twin Lake
# PL    = Paulina Lake

# -----------------------------
# 0) site-level coordinates
# Table 1 coordinates converted to decimal degrees
# order matches Supplementary Table 6
# -----------------------------
GAAC_1_coords <- data.frame(
  site = c(
    "C", "Sz", "Si", "SiC", "U", "WC", "EC", "LL", "PC", "U196",
    "Mo14", "Sa39", "L34", "Ma32", "W296", "W320", "CF27", "W278",
    "Mc18", "Mc42", "Mc61", "CR", "ST", "PL"
  ),
  lat = c(
    46.243344, 44.887411, 44.002156, 43.989556, 43.692878, 43.277053,
    43.586961, 44.092692, 43.370008, 43.284083, 45.236894, 44.770122,
    44.774533, 44.395581, 44.057083, 43.966997, 43.859742, 44.145006,
    44.078047, 44.070456, 44.134400, 44.290436, 43.710661, 43.713517
  ),
  lon = c(
    -123.902531, -123.996167, -124.133128, -124.045261, -124.000461, -124.319047,
    -124.186100, -124.116336, -124.262089, -123.329628, -122.631772, -122.846378,
    -123.318061, -123.343172, -123.075258, -122.869003, -123.023606, -123.118022,
    -123.026367, -122.798522, -122.608850, -120.846408, -121.767925, -121.272633
  ),
  stringsAsFactors = FALSE
)

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
# transcribed from Supplementary Table 6
# order:
# C, Sz, Si, SiC, U, WC, EC, LL, PC, U196,
# Mo14, Sa39, L34, Ma32, W296, W320, CF27, W278,
# Mc18, Mc42, Mc61, CR, ST, PL
# -----------------------------
fst_order <- GAAC_1_coords$site

fst_vals <- c(
  # Sz
  0.016,
  
  # Si
  0.012, 0.014,
  
  # SiC
  0.011, 0.012, 0.005,
  
  # U
  0.026, 0.029, 0.016, 0.017,
  
  # WC
  0.022, 0.024, 0.009, 0.011, 0.030,
  
  # EC
  0.090, 0.097, 0.058, 0.062, 0.115, 0.097,
  
  # LL
  0.056, 0.059, 0.036, 0.039, 0.057, 0.049, 0.069,
  
  # PC
  0.045, 0.051, 0.028, 0.033, 0.043, 0.031, 0.085, 0.045,
  
  # U196
  0.098, 0.102, 0.077, 0.074, 0.164, 0.148, 0.330, 0.161, 0.162,
  
  # Mo14
  0.109, 0.123, 0.090, 0.091, 0.159, 0.148, 0.299, 0.146, 0.165, 0.424,
  
  # Sa39
  0.044, 0.053, 0.037, 0.037, 0.062, 0.058, 0.151, 0.080, 0.077, 0.212, 0.107,
  
  # L34
  0.152, 0.168, 0.126, 0.126, 0.234, 0.213, 0.406, 0.193, 0.222, 0.580, 0.194, 0.186,
  
  # Ma32
  0.142, 0.157, 0.118, 0.118, 0.212, 0.193, 0.366, 0.178, 0.206, 0.509, 0.160, 0.162, 0.163,
  
  # W296
  0.142, 0.157, 0.118, 0.118, 0.210, 0.192, 0.364, 0.177, 0.206, 0.505, 0.159, 0.161, 0.175, 0.088,
  
  # W320
  0.158, 0.174, 0.128, 0.128, 0.246, 0.221, 0.431, 0.200, 0.228, 0.609, 0.237, 0.201, 0.341, 0.180, 0.095,
  
  # CF27
  0.155, 0.172, 0.127, 0.127, 0.242, 0.218, 0.424, 0.197, 0.225, 0.604, 0.230, 0.199, 0.329, 0.192, 0.135, 0.371,
  
  # W278
  0.147, 0.163, 0.123, 0.123, 0.218, 0.201, 0.375, 0.184, 0.213, 0.515, 0.174, 0.172, 0.171, 0.086, 0.089, 0.160, 0.157,
  
  # Mc18
  0.081, 0.095, 0.079, 0.081, 0.082, 0.085, 0.171, 0.096, 0.131, 0.230, 0.061, 0.056, 0.070, 0.048, 0.045, 0.065, 0.067, 0.045,
  
  # Mc42
  0.052, 0.061, 0.041, 0.042, 0.076, 0.067, 0.175, 0.087, 0.086, 0.265, 0.158, 0.059, 0.248, 0.215, 0.214, 0.268, 0.263, 0.230, 0.047,
  
  # Mc61
  0.057, 0.058, 0.040, 0.039, 0.084, 0.073, 0.202, 0.096, 0.095, 0.281, 0.301, 0.128, 0.416, 0.370, 0.367, 0.430, 0.429, 0.382, 0.172, 0.106,
  
  # CR
  0.128, 0.143, 0.104, 0.104, 0.210, 0.185, 0.362, 0.170, 0.186, 0.528, 0.146, 0.148, 0.218, 0.143, 0.118, 0.299, 0.294, 0.161, 0.032, 0.186, 0.358,
  
  # ST
  0.143, 0.158, 0.120, 0.120, 0.212, 0.195, 0.367, 0.178, 0.208, 0.508, 0.149, 0.161, 0.194, 0.138, 0.114, 0.248, 0.251, 0.152, 0.048, 0.210, 0.376, 0.023,
  
  # PL
  0.128, 0.143, 0.104, 0.104, 0.217, 0.188, 0.377, 0.172, 0.187, 0.561, 0.164, 0.154, 0.267, 0.165, 0.141, 0.348, 0.356, 0.186, 0.035, 0.193, 0.373, 0.085, 0.058
)

GAAC_1_fst <- fill_sym_from_lower(fst_order, fst_vals, diag_val = 0)

# -----------------------------
# IBD plot: straight-line distance vs FST
# -----------------------------
dist_km <- geosphere::distm(
  as.matrix(GAAC_1_coords[, c("lon", "lat")]),
  fun = geosphere::distGeo
) / 1000

ibd_df <- data.frame(
  dist = dist_km[upper.tri(dist_km)],
  fst  = GAAC_1_fst[upper.tri(GAAC_1_fst)]
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
stopifnot(identical(rownames(GAAC_1_fst), GAAC_1_coords$site))
stopifnot(identical(colnames(GAAC_1_fst), GAAC_1_coords$site))
stopifnot(isTRUE(all.equal(GAAC_1_fst, t(GAAC_1_fst))))

# -----------------------------
# 4) save RData
# -----------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  GAAC_1_fst,
  GAAC_1_coords,
  file = file.path(out_dir, "GAAC-1.RData")
)