# -----------------------------
# ONCL-1 cutthroat trout
# site coordinates + FST matrix + IBD plot
# -----------------------------

library(ggplot2)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ONCL-1/"

# ------------------------------
# 0) site coordinates
# ------------------------------
ONCL_1_coords <- data.frame(
  site = c("4","5","6","7","8","9","10","11","12","13","14","15","17","18","20","21","22","23","24"),
  lat  = c(48.81624,48.82404,48.83410,48.84562,48.85330,48.86019,48.86138,
           48.86132,48.85884,48.85991,48.86753,48.88227,48.86632,48.87812,
           48.89360,48.86600,48.87413,48.88316,48.89112),
  lon  = c(-114.29494,-114.29529,-114.29529,-114.29082,-114.28576,-114.27529,
           -114.25675,-114.23821,-114.22825,-114.21701,-114.20422,-114.19916,
           -114.26359,-114.26480,-114.26995,-114.24179,-114.23364,-114.22832,
           -114.22016),
  stringsAsFactors = FALSE
)

# ------------------------------
# 1) pairwise FST matrix
# ------------------------------
ONCL_1_fst <- matrix(
  0,
  nrow = nrow(ONCL_1_coords),
  ncol = nrow(ONCL_1_coords),
  dimnames = list(ONCL_1_coords$site, ONCL_1_coords$site)
)

ONCL_1_fst["5",  c("4")] <- c(0.000)

ONCL_1_fst["6",  c("4","5")] <- c(0.000, 0.000)

ONCL_1_fst["7",  c("4","5","6")] <- c(0.000, 0.000, 0.000)

ONCL_1_fst["8",  c("4","5","6","7")] <- c(0.000, 0.007, 0.000, 0.000)

ONCL_1_fst["9",  c("4","5","6","7","8")] <- c(0.000, 0.000, 0.000, 0.000, 0.000)

ONCL_1_fst["10", c("4","5","6","7","8","9")] <- c(0.011, 0.000, 0.000, 0.025, 0.037, 0.000)

ONCL_1_fst["11", c("4","5","6","7","8","9","10")] <- c(0.000, 0.000, 0.010, 0.003, 0.008, 0.000, 0.000)

ONCL_1_fst["12", c("4","5","6","7","8","9","10","11")] <- c(0.024, 0.016, 0.049, 0.031, 0.030, 0.015, 0.014, 0.022)

ONCL_1_fst["13", c("4","5","6","7","8","9","10","11","12")] <- c(0.000, 0.000, 0.037, 0.012, 0.028, 0.010, 0.003, 0.000, 0.006)

ONCL_1_fst["14", c("4","5","6","7","8","9","10","11","12","13")] <- c(0.005, 0.000, 0.025, 0.010, 0.014, 0.000, 0.000, 0.000, 0.014, 0.000)

ONCL_1_fst["15", c("4","5","6","7","8","9","10","11","12","13","14")] <- c(0.001, 0.008, 0.011, 0.003, 0.018, 0.011, 0.035, 0.006, 0.021, 0.000, 0.000)

ONCL_1_fst["17", c("4","5","6","7","8","9","10","11","12","13","14","15")] <- c(0.023, 0.000, 0.022, 0.006, 0.021, 0.016, 0.008, 0.020, 0.019, 0.000, 0.005, 0.014)

ONCL_1_fst["18", c("4","5","6","7","8","9","10","11","12","13","14","15","17")] <- c(0.000, 0.000, 0.010, 0.000, 0.000, 0.000, 0.007, 0.010, 0.002, 0.004, 0.000, 0.008, 0.010)

ONCL_1_fst["20", c("4","5","6","7","8","9","10","11","12","13","14","15","17","18")] <- c(0.009, 0.000, 0.038, 0.021, 0.037, 0.015, 0.000, 0.000, 0.007, 0.000, 0.016, 0.014, 0.006, 0.004)

ONCL_1_fst["21", c("4","5","6","7","8","9","10","11","12","13","14","15","17","18","20")] <- c(0.027, 0.049, 0.038, 0.025, 0.054, 0.024, 0.016, 0.025, 0.009, 0.000, 0.038, 0.008, 0.010, 0.037, 0.027)

ONCL_1_fst["22", c("4","5","6","7","8","9","10","11","12","13","14","15","17","18","20","21")] <- c(0.048, 0.059, 0.046, 0.008, 0.055, 0.027, 0.069, 0.042, 0.051, 0.054, 0.069, 0.034, 0.041, 0.049, 0.018, 0.050)

ONCL_1_fst["23", c("4","5","6","7","8","9","10","11","12","13","14","15","17","18","20","21","22")] <- c(0.043, 0.015, 0.074, 0.050, 0.061, 0.039, 0.033, 0.038, 0.088, 0.036, 0.054, 0.060, 0.030, 0.021, 0.060, 0.063, 0.094)

ONCL_1_fst["24", c("4","5","6","7","8","9","10","11","12","13","14","15","17","18","20","21","22","23")] <- c(0.051, 0.081, 0.086, 0.068, 0.098, 0.061, 0.071, 0.066, 0.135, 0.098, 0.118, 0.089, 0.085, 0.085, 0.079, 0.104, 0.096, 0.047)

ONCL_1_fst[upper.tri(ONCL_1_fst)] <- t(ONCL_1_fst)[upper.tri(ONCL_1_fst)]
diag(ONCL_1_fst) <- 0
ONCL_1_fst[ONCL_1_fst < 0] <- 0

# ------------------------------
# 2) map of sampling locations
# ------------------------------
ggplot(ONCL_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), position = position_jitter(width = 0.003, height = 0.003)) +
  coord_equal() +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "ONCL-1 sampling locations"
  )

# ------------------------------
# 3) geographic distance matrix
# straight-line distance in km
# ------------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(ONCL_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- ONCL_1_coords$site
colnames(geo_dist_km) <- ONCL_1_coords$site

# ------------------------------
# 4) pairwise dataframe for IBD plot
# ------------------------------
ibd_df <- data.frame(
  site1 = rownames(ONCL_1_fst)[row(ONCL_1_fst)[upper.tri(ONCL_1_fst)]],
  site2 = colnames(ONCL_1_fst)[col(ONCL_1_fst)[upper.tri(ONCL_1_fst)]],
  fst = ONCL_1_fst[upper.tri(ONCL_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# ------------------------------
# 5) IBD plot
# ------------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(
    x = "Straight-line distance (km)",
    y = expression(F[ST]),
    title = "ONCL-1 isolation by distance"
  )

# ------------------------------
# 6) quick checks
# ------------------------------
stopifnot(identical(rownames(ONCL_1_fst), ONCL_1_coords$site))
stopifnot(identical(colnames(ONCL_1_fst), ONCL_1_coords$site))
stopifnot(isTRUE(all.equal(ONCL_1_fst, t(ONCL_1_fst))))

# ------------------------------
# 7) save RData
# ------------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  ONCL_1_fst,
  ONCL_1_coords,
  file = file.path(out_dir, "ONCL-1.RData")
)