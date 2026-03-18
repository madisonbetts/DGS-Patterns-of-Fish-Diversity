############################
## RKBS_1 workflow
# WI NORTHERN ROCKBASS 
## 1) site coordinates / metadata
## 2) pairwise FST matrix
## 3) pairwise p-value matrix
## 4) pairwise Dest matrix
## 5) reorder to match site table
## 6) save outputs
############################

# -------------------------
# 1) Table 5 site metadata
# alphabetical order as shown in site table
# -------------------------
AMRU_1_coords <- data.frame(
  lake = c(
    "Bearskin Lake", "Big Lake", "Butternut Lake", "Cranberry Lake",
    "Deerskin Lake", "Half Moon Lake", "Lac Courte Oreilles",
    "Lake Galilee", "Lake Nebagamon", "Lake Noquebay",
    "Long/Herde Lake Chain", "Lynx Lake", "McKenzie Lake",
    "Pelican Lake", "Pike/Round Chain", "Plum Lake",
    "Sawyer Lake", "Solberg Lake", "Somo Lake",
    "Spread Eagles Chain", "Tomahawk Lake", "White Sand Lake"
  ),
  pop = c(
    "BEA", "BIG", "BN", "CRA", "DEE", "HMO", "LCO", "GAL", "NEB", "NOQ",
    "LHC", "LYN", "MCK", "PEL", "PRC", "PM", "SAW", "SOL", "SOM", "SPE",
    "TOM", "WSL"
  ),
  county = c(
    "Oneida", "Vilas", "Price/Ashland", "Price", "Vilas", "Polk", "Sawyer",
    "Ashland", "Bayfield", "Marinette", "Chippewa", "Vilas", "Washburn",
    "Oneida", "Price", "Vilas", "Langlade", "Price", "Lincoln", "Florence",
    "Oneida", "Vilas"
  ),
  lat = c(
    45.731, 46.154, 45.983, 45.621, 45.975, 45.498, 45.892, 46.287, 46.513,
    45.249, 45.241, 46.195, 45.936, 45.503, 45.928, 46.004, 45.247, 45.749,
    45.516, 45.901, 45.830, 46.008
  ),
  lon = c(
    -89.685, -89.767, -90.515, -90.350, -89.171, -92.438, -91.438, -90.602,
    -91.703, -87.924, -91.416, -89.666, -92.046, -89.202, -90.067, -89.514,
    -88.757, -90.369, -89.868, -88.145, -89.661, -89.827
  ),
  mgmt_unit = c(
    "Wisconsin", "Chippewa", "Chippewa", "Chippewa", "Wisconsin", "St. Croix",
    "Chippewa", "Lake Superior", "Lake Superior", "Lake Michigan", "Chippewa",
    "Lake Superior", "St. Croix", "Wisconsin", "Chippewa", "Wisconsin",
    "Lake Michigan", "Chippewa", "Wisconsin", "Lake Michigan", "Wisconsin",
    "Chippewa"
  ),
  stringsAsFactors = FALSE
)

AMRU_1_coords

# -------------------------
# Population order in Tables 11 and 13
# -------------------------
pops <- c(
  "TOM", "PM", "BEA", "SOL", "WSL", "PEL", "LYN", "LCO", "BIG", "SAW",
  "SOM", "CRA", "DEE", "NOQ", "SPE", "BN", "GAL", "PRC", "LHC", "MCK",
  "NEB", "HMO"
)

n <- length(pops)
stopifnot(n == 22)

# -------------------------
# Helper: fill symmetric matrix from row-wise lower triangle
# -------------------------
fill_sym_from_lower <- function(pops, vals, diag_val = 0) {
  n <- length(pops)
  stopifnot(length(vals) == n * (n - 1) / 2)
  
  mat <- matrix(
    0,
    nrow = n,
    ncol = n,
    dimnames = list(pops, pops)
  )
  
  lt_idx <- which(lower.tri(mat), arr.ind = TRUE)
  lt_idx <- lt_idx[order(lt_idx[, 1], lt_idx[, 2]), ]
  
  mat[lt_idx] <- vals
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  diag(mat) <- diag_val
  
  mat
}

# -------------------------
# 2) Table 11: pairwise FST
# below diagonal, row-wise
# -------------------------
fst_lower <- c(
  0.04,
  0.07, 0.10,
  0.15, 0.10, 0.20,
  0.24, 0.22, 0.25, 0.13,
  0.15, 0.11, 0.17, 0.13, 0.24,
  0.14, 0.15, 0.17, 0.13, 0.11, 0.21,
  0.18, 0.15, 0.18, 0.08, 0.04, 0.20, 0.09,
  0.22, 0.18, 0.22, 0.09, 0.03, 0.22, 0.13, 0.04,
  0.28, 0.23, 0.33, 0.18, 0.33, 0.29, 0.26, 0.24, 0.28,
  0.17, 0.21, 0.17, 0.16, 0.30, 0.16, 0.26, 0.22, 0.25, 0.36,
  0.24, 0.18, 0.24, 0.17, 0.12, 0.24, 0.15, 0.10, 0.14, 0.26, 0.37,
  0.28, 0.20, 0.27, 0.30, 0.37, 0.23, 0.26, 0.29, 0.35, 0.37, 0.42, 0.32,
  0.14, 0.09, 0.19, 0.05, 0.19, 0.13, 0.15, 0.12, 0.15, 0.11, 0.18, 0.18, 0.26,
  0.23, 0.18, 0.25, 0.10, 0.04, 0.22, 0.15, 0.05, 0.04, 0.25, 0.32, 0.07, 0.34, 0.14,
  0.18, 0.14, 0.20, 0.11, 0.14, 0.19, 0.13, 0.08, 0.08, 0.21, 0.24, 0.16, 0.26, 0.12, 0.10,
  0.15, 0.11, 0.20, 0.06, 0.14, 0.20, 0.12, 0.05, 0.12, 0.19, 0.25, 0.13, 0.27, 0.09, 0.10, 0.09,
  0.23, 0.18, 0.24, 0.10, 0.07, 0.19, 0.10, 0.03, 0.08, 0.23, 0.28, 0.09, 0.29, 0.13, 0.07, 0.11, 0.09,
  0.24, 0.20, 0.25, 0.11, 0.14, 0.25, 0.17, 0.06, 0.09, 0.21, 0.26, 0.19, 0.34, 0.14, 0.12, 0.09, 0.11, 0.09,
  0.22, 0.19, 0.33, 0.16, 0.41, 0.30, 0.37, 0.30, 0.34, 0.29, 0.28, 0.42, 0.51, 0.17, 0.37, 0.28, 0.21, 0.38, 0.30,
  0.22, 0.21, 0.30, 0.15, 0.33, 0.32, 0.31, 0.21, 0.25, 0.25, 0.29, 0.36, 0.45, 0.16, 0.28, 0.21, 0.15, 0.30, 0.19, 0.12,
  0.38, 0.31, 0.45, 0.33, 0.52, 0.36, 0.47, 0.43, 0.45, 0.32, 0.49, 0.45, 0.49, 0.31, 0.44, 0.34, 0.37, 0.45, 0.40, 0.37, 0.41
)

AMRU_1_fst <- fill_sym_from_lower(pops, fst_lower, diag_val = 0)




# -------------------------
# 5) Reorder matrices to match AMRU_1_coords
# alphabetical order from Table 5
# -------------------------
site_order <- AMRU_1_coords$pop

stopifnot(length(site_order) == n)
stopifnot(all(site_order %in% pops))

AMRU_1_fst     <- AMRU_1_fst[site_order, site_order]
#RKBS_pvalues <- RKBS_pvalues[site_order, site_order]
#RKBS_dest    <- RKBS_dest[site_order, site_order]

# -------------------------
# 6) Save
# -------------------------
out_file <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/AMRU-1/data/AMRU_1.RData"

save(
  AMRU_1_fst,
  AMRU_1_coords,
  file = out_file
)
