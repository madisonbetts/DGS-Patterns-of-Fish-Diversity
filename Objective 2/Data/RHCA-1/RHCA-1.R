# -----------------------------
# RHCA-1 Rhinichthys cataractae complex
# Beaver River, ON removed
# coordinates taken directly from Table 1
# pairwise FST collapsed by stream from Table S8
# -----------------------------

library(ggplot2)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/RHCA-1"

# -----------------------------
# 0) site coordinates
# exact coordinates from Table 1
# site 12 (Beaver River, ON) removed
# -----------------------------
RHCA_1_coords <- data.frame(
  site = as.character(1:11),
  stream = c(
    "Porter Creek, WA",
    "West Fork Satsop River, WA",
    "Bertrand Creek, BC",
    "Brunette River, BC",
    "Coquitlam River, BC",
    "Alouette River, BC",
    "Kanaka Creek, BC",
    "Norrish Creek, BC",
    "Fraser River, BC",
    "Coquihalla River, BC",
    "Beaver Creek, BC"
  ),
  lat = c(
    46.9465,
    47.0598,
    49.0043,
    49.2414,
    49.0416,
    49.2392,
    49.2022,
    49.2349,
    49.3831,
    49.3885,
    49.1009
  ),
  lon = c(
    -123.2954,
    -123.5410,
    -122.5322,
    -122.8959,
    -122.7711,
    -122.5793,
    -122.5413,
    -122.1336,
    -122.4522,
    -121.4334,
    -117.5572
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map of sampling locations
# -----------------------------
ggplot(RHCA_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.20, size = 3.5) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "RHCA-1 sampling locations"
  )

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(RHCA_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- RHCA_1_coords$site
colnames(geo_dist_km) <- RHCA_1_coords$site
diag(geo_dist_km) <- 0

# -----------------------------
# 3) helper function
# takes a vector of lower-triangle values (row-wise order)
# and reconstructs a full symmetric matrix
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
# collapsed from Table S8
#
# original sympatric subsets:
# Coquitlam River = Coquitlam LND / Coquitlam NSD
# Alouette River  = Alouette LND / Alouette NSD
# Kanaka Creek    = Kanaka LND / Kanaka NSD
#
# collapsing rule:
# - allopatric vs sympatric = mean(LND, NSD)
# - sympatric vs sympatric  = mean(all 4 cross-comparisons)
#
# site 12 (Beaver River, ON) removed entirely
# -----------------------------
fst_vals <- c(
  
  # 2 Satsop vs 1 Porter
  0.0227,
  
  # 3 Bertrand vs 1-2
  0.1005, 0.0799,
  
  # 4 Brunette vs 1-3
  0.3463, 0.2721, 0.2292,
  
  # 5 Coquitlam vs 1-4
  mean(c(0.0853, 0.0906)),   # Porter
  mean(c(0.0784, 0.0867)),   # Satsop
  mean(c(0.0475, 0.0554)),   # Bertrand
  mean(c(0.1840, 0.2059)),   # Brunette
  
  # 6 Alouette vs 1-5
  mean(c(0.0755, 0.0728)),   # Porter
  mean(c(0.0868, 0.0840)),   # Satsop
  mean(c(0.0509, 0.0555)),   # Bertrand
  mean(c(0.1874, 0.2136)),   # Brunette
  mean(c(0.0153, 0.0234, 0.0208, 0.0390)), # Coquitlam x Alouette
  
  # 7 Kanaka vs 1-6
  mean(c(0.0747, 0.0833)),   # Porter
  mean(c(0.0792, 0.0715)),   # Satsop
  mean(c(0.0375, 0.0379)),   # Bertrand
  mean(c(0.2640, 0.2237)),   # Brunette
  mean(c(0.0390, 0.0473, 0.0450, 0.0439)), # Coquitlam x Kanaka
  mean(c(0.0374, 0.0454, 0.0403, 0.0531)), # Alouette x Kanaka
  
  # 8 Norrish vs 1-7
  0.1134,
  0.1280,
  0.1144,
  0.3051,
  mean(c(0.0972, 0.1084)),   # Coquitlam
  mean(c(0.0874, 0.0764)),   # Alouette
  mean(c(0.0943, 0.1150)),   # Kanaka
  
  # 9 Fraser vs 1-8
  0.0682,
  0.0926,
  0.0887,
  0.3120,
  mean(c(0.0527, 0.0509)),   # Coquitlam
  mean(c(0.0262, 0.0363)),   # Alouette
  mean(c(0.0565, 0.0736)),   # Kanaka
  0.0808,
  
  # 10 Coquihalla vs 1-9
  0.0883,
  0.1039,
  0.1097,
  0.2660,
  mean(c(0.0596, 0.0587)),   # Coquitlam
  mean(c(0.0484, 0.0476)),   # Alouette
  mean(c(0.0929, 0.1098)),   # Kanaka
  0.0859,
  0.0338,
  
  # 11 Beaver Creek BC vs 1-10
  0.1704,
  0.1549,
  0.1664,
  0.3343,
  mean(c(0.1465, 0.1512)),   # Coquitlam
  mean(c(0.1411, 0.1422)),   # Alouette
  mean(c(0.1492, 0.1620)),   # Kanaka
  0.1833,
  0.1260,
  0.1192
)

RHCA_1_fst <- fill_sym_from_lower(
  pops = RHCA_1_coords$site,
  vals = fst_vals,
  diag_val = 0
)

# set negative FST to 0 and force diagonal = 0
RHCA_1_fst <- pmax(RHCA_1_fst, 0)
diag(RHCA_1_fst) <- 0

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(RHCA_1_fst)[row(RHCA_1_fst)[upper.tri(RHCA_1_fst)]],
  site2   = colnames(RHCA_1_fst)[col(RHCA_1_fst)[upper.tri(RHCA_1_fst)]],
  fst     = RHCA_1_fst[upper.tri(RHCA_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 6) IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "RHCA-1 isolation by distance"
  )

# -----------------------------
# 7) quick checks
# -----------------------------
stopifnot(identical(rownames(RHCA_1_fst), RHCA_1_coords$site))
stopifnot(identical(colnames(RHCA_1_fst), RHCA_1_coords$site))
stopifnot(isTRUE(all.equal(RHCA_1_fst, t(RHCA_1_fst))))
stopifnot(all(diag(RHCA_1_fst) == 0))
stopifnot(all(diag(geo_dist_km) == 0))

# -----------------------------
# 8) save RData
# -----------------------------
out_dir <- file.path(base_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  RHCA_1_fst,
  RHCA_1_coords,
  file = file.path(out_dir, "RHCA-1.RData")
)