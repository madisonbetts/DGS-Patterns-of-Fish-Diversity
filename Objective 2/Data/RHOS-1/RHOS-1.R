# -----------------------------
# RHOS-1 speckled dace
# site coordinates + FST matrix + IBD plot
# -----------------------------

library(ggplot2)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/RHOS-1"

# -----------------------------
# 0) site coordinates
# order matches Table 1 and Appendix B exactly
# -----------------------------
RHOS_1_coords <- data.frame(
  site = c(
    "SEV","SPR","LNK","SPE","JEN","COP","WIL","SR","SK","MOF",
    "LSR","CN","BB","BLU","NFS","SFS","TT","DL","FG","CAN",
    "GVY","SFK","SWFT","TRL","EFT"
  ),
  lat = c(
    42.69421, 42.56339, 42.22160, 42.15257, 42.11791,
    41.99150, 41.86591, 41.59100, 41.81309, 41.63370,
    41.76529, 41.61363, 41.25204, 41.44417, 41.29301,
    41.18689, 41.02586, 40.77490, 40.37700, 40.73800,
    40.68979, 40.85580, 40.98642, 41.05328, 41.00850
  ),
  lon = c(
    -122.0727, -121.8624, -121.7934, -122.0278, -122.3670,
    -122.1902, -122.4648, -122.4380, -122.5923, -122.7490,
    -123.0213, -123.4958, -123.6348, -123.9069, -123.2301,
    -123.2139, -123.6403, -123.3250, -123.3256, -123.0495,
    -122.8576, -122.8851, -122.7089, -122.6965, -122.6201
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map of sampling locations
# -----------------------------
ggplot(RHOS_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.03, size = 3.5) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "RHOS-1 sampling locations"
  )

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(RHOS_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- RHOS_1_coords$site
colnames(geo_dist_km) <- RHOS_1_coords$site

# -----------------------------
# 3) pairwise FST matrix
# Appendix B: FST values are ABOVE the diagonal
# -----------------------------
RHOS_1_fst <- matrix(
  0,
  nrow = nrow(RHOS_1_coords),
  ncol = nrow(RHOS_1_coords),
  dimnames = list(RHOS_1_coords$site, RHOS_1_coords$site)
)

RHOS_1_fst["SEV", c("SPR","LNK","SPE","JEN","COP","WIL","SR","SK","MOF","LSR","CN","BB","BLU","NFS","SFS","TT","DL","FG","CAN","GVY","SFK","SWFT","TRL","EFT")] <-
  c(0.010,0.024,0.016,0.241,0.019,0.026,0.040,0.021,0.032,0.028,0.019,0.039,0.031,0.150,0.182,0.127,0.221,0.204,0.230,0.207,0.222,0.235,0.229,0.190)

RHOS_1_fst["SPR", c("LNK","SPE","JEN","COP","WIL","SR","SK","MOF","LSR","CN","BB","BLU","NFS","SFS","TT","DL","FG","CAN","GVY","SFK","SWFT","TRL","EFT")] <-
  c(0.017,0.010,0.245,0.013,0.017,0.031,0.024,0.032,0.020,0.020,0.023,0.032,0.146,0.181,0.122,0.221,0.205,0.233,0.209,0.224,0.234,0.228,0.194)

RHOS_1_fst["LNK", c("SPE","JEN","COP","WIL","SR","SK","MOF","LSR","CN","BB","BLU","NFS","SFS","TT","DL","FG","CAN","GVY","SFK","SWFT","TRL","EFT")] <-
  c(0.006,0.259,0.014,0.020,0.037,0.024,0.037,0.020,0.016,0.020,0.020,0.147,0.180,0.125,0.211,0.211,0.224,0.202,0.214,0.229,0.220,0.188)

RHOS_1_fst["SPE", c("JEN","COP","WIL","SR","SK","MOF","LSR","CN","BB","BLU","NFS","SFS","TT","DL","FG","CAN","GVY","SFK","SWFT","TRL","EFT")] <-
  c(0.246,0.002,0.005,0.026,0.011,0.024,0.007,0.002,0.015,0.011,0.137,0.171,0.114,0.206,0.197,0.220,0.196,0.209,0.223,0.214,0.178)

RHOS_1_fst["JEN", c("COP","WIL","SR","SK","MOF","LSR","CN","BB","BLU","NFS","SFS","TT","DL","FG","CAN","GVY","SFK","SWFT","TRL","EFT")] <-
  c(0.248,0.215,0.264,0.218,0.238,0.233,0.227,0.249,0.231,0.285,0.320,0.268,0.329,0.316,0.343,0.312,0.337,0.332,0.345,0.298)

RHOS_1_fst["COP", c("WIL","SR","SK","MOF","LSR","CN","BB","BLU","NFS","SFS","TT","DL","FG","CAN","GVY","SFK","SWFT","TRL","EFT")] <-
  c(0.008,0.023,0.008,0.016,0.002,0.002,0.007,0.005,0.121,0.156,0.103,0.196,0.182,0.207,0.186,0.198,0.209,0.201,0.171)

RHOS_1_fst["WIL", c("SR","SK","MOF","LSR","CN","BB","BLU","NFS","SFS","TT","DL","FG","CAN","GVY","SFK","SWFT","TRL","EFT")] <-
  c(0.005,0.017,0.021,0.008,0.009,0.012,0.011,0.123,0.159,0.106,0.193,0.184,0.208,0.183,0.195,0.209,0.202,0.170)

RHOS_1_fst["SR", c("SK","MOF","LSR","CN","BB","BLU","NFS","SFS","TT","DL","FG","CAN","GVY","SFK","SWFT","TRL","EFT")] <-
  c(0.033,0.021,0.032,0.030,0.037,0.025,0.185,0.224,0.166,0.270,0.253,0.283,0.256,0.275,0.285,0.281,0.238)

RHOS_1_fst["SK", c("MOF","LSR","CN","BB","BLU","NFS","SFS","TT","DL","FG","CAN","GVY","SFK","SWFT","TRL","EFT")] <-
  c(0.016,0.004,0.005,0.015,0.009,0.102,0.136,0.087,0.167,0.158,0.178,0.157,0.171,0.179,0.175,0.137)

RHOS_1_fst["MOF", c("LSR","CN","BB","BLU","NFS","SFS","TT","DL","FG","CAN","GVY","SFK","SWFT","TRL","EFT")] <-
  c(0.012,0.024,0.031,0.027,0.137,0.173,0.116,0.207,0.190,0.220,0.198,0.213,0.226,0.216,0.179)

RHOS_1_fst["LSR", c("CN","BB","BLU","NFS","SFS","TT","DL","FG","CAN","GVY","SFK","SWFT","TRL","EFT")] <-
  c(0.005,0.004,0.007,0.100,0.137,0.085,0.169,0.160,0.180,0.162,0.173,0.184,0.176,0.145)

RHOS_1_fst["CN", c("BB","BLU","NFS","SFS","TT","DL","FG","CAN","GVY","SFK","SWFT","TRL","EFT")] <-
  c(0.012,0.003,0.106,0.141,0.090,0.177,0.162,0.189,0.167,0.180,0.189,0.183,0.152)

RHOS_1_fst["BB", c("BLU","NFS","SFS","TT","DL","FG","CAN","GVY","SFK","SWFT","TRL","EFT")] <-
  c(0.003,0.099,0.131,0.087,0.176,0.170,0.188,0.167,0.179,0.189,0.183,0.156)

RHOS_1_fst["BLU", c("NFS","SFS","TT","DL","FG","CAN","GVY","SFK","SWFT","TRL","EFT")] <-
  c(0.109,0.143,0.095,0.188,0.172,0.202,0.178,0.192,0.200,0.196,0.165)

RHOS_1_fst["NFS", c("SFS","TT","DL","FG","CAN","GVY","SFK","SWFT","TRL","EFT")] <-
  c(0.010,0.016,0.054,0.043,0.066,0.052,0.056,0.045,0.048,0.045)

RHOS_1_fst["SFS", c("TT","DL","FG","CAN","GVY","SFK","SWFT","TRL","EFT")] <-
  c(0.043,0.081,0.075,0.094,0.079,0.080,0.067,0.073,0.069)

RHOS_1_fst["TT", c("DL","FG","CAN","GVY","SFK","SWFT","TRL","EFT")] <-
  c(0.056,0.027,0.059,0.052,0.059,0.054,0.047,0.037)

RHOS_1_fst["DL", c("FG","CAN","GVY","SFK","SWFT","TRL","EFT")] <-
  c(0.002,0.010,0.002,0.002,0.009,0.000,0.013)

RHOS_1_fst["FG", c("CAN","GVY","SFK","SWFT","TRL","EFT")] <-
  c(0.080,0.058,0.069,0.064,0.054,0.054)

RHOS_1_fst["CAN", c("GVY","SFK","SWFT","TRL","EFT")] <-
  c(0.008,0.012,0.006,0.019,0.015)

RHOS_1_fst["GVY", c("SFK","SWFT","TRL","EFT")] <-
  c(0.003,0.010,0.006,0.015)

RHOS_1_fst["SFK", c("SWFT","TRL","EFT")] <-
  c(0.008,0.000,0.018)

RHOS_1_fst["SWFT", c("TRL","EFT")] <-
  c(0.000,0.019)

RHOS_1_fst["TRL", "EFT"] <- 0.014

# reflect upper triangle to lower triangle
RHOS_1_fst[lower.tri(RHOS_1_fst)] <- t(RHOS_1_fst)[lower.tri(RHOS_1_fst)]
diag(RHOS_1_fst) <- 0
RHOS_1_fst[RHOS_1_fst < 0] <- 0

# -----------------------------
# 4) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(RHOS_1_fst)[row(RHOS_1_fst)[upper.tri(RHOS_1_fst)]],
  site2   = colnames(RHOS_1_fst)[col(RHOS_1_fst)[upper.tri(RHOS_1_fst)]],
  fst     = RHOS_1_fst[upper.tri(RHOS_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 5) IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "RHOS-1 isolation by distance"
  )

# -----------------------------
# 6) quick checks
# -----------------------------
stopifnot(identical(rownames(RHOS_1_fst), RHOS_1_coords$site))
stopifnot(identical(colnames(RHOS_1_fst), RHOS_1_coords$site))
stopifnot(isTRUE(all.equal(RHOS_1_fst, t(RHOS_1_fst))))

# -----------------------------
# 7) save RData
# -----------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  RHOS_1_fst,
  RHOS_1_coords,
  geo_dist_km,
  ibd_df,
  file = file.path(out_dir, "RHOS-1.RData")
)