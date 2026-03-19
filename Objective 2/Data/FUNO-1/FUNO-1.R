# -----------------------------
# FUNO-1 Fundulus notatus
# site coordinates + FST matrix
# -----------------------------

library(ggplot2)
library(geosphere)

# directory where everything is
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/FUNO-1/"

# -----------------------------
# 0) site coordinates
# only sites included in the FST matrix
# -----------------------------
site_data <- data.frame(
  site_id = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17),
  fork = c(
    "Middle", "Middle", "Middle", "Middle", "Middle",
    "Eagle", "Eagle", "Eagle",
    "Main", "Main", "Main",
    "North", "North", "South", "North", "South", "South"
  ),
  latitude = c(
    37.81, 37.86, 37.81, 37.79, 37.77,
    37.66, 37.65, 37.65,
    37.57, 37.73, 37.70,
    37.89, 37.94, 37.64, 37.83, 37.66, 37.68
  ),
  longitude = c(
    -88.76, -88.70, -88.61, -88.71, -88.54,
    -88.32, -88.30, -88.26,
    -88.13, -88.35, -88.29,
    -88.33, -88.33, -88.63, -88.33, -88.63, -88.80
  ),
  fundulus_olivaceus_n = c(
    27, 19, 0, 30, 0,
    14, 6, 6,
    1, 1, 10,
    42, 14, 41, 1, 28, 18
  ),
  f_notatus_n = c(
    0, 65, 19, 0, 6,
    2, 11, 14,
    29, 8, 12,
    2, 0, 0, 22, 38, 0
  ),
  hybrid_ancestry_n = c(
    0, 9, 5, 0, 0,
    2, 2, 3,
    3, 0, 4,
    1, 6, 0, 2, 0, 0
  ),
  cumulative_drainage_area_km2 = c(
    3.8, 18.9, 181.1, 15.0, 215.3,
    65.8, 58.7, 148.6,
    3041.9, 1437.0, 2415.1,
    132.9, 29.2, 84.7, 998.1, 657.5, 0.9
  )
)

FUNO_1_coords <- subset(
  site_data,
  site_id %in% c(2, 3, 5, 7, 8, 9, 10, 11, 15, 16),
  select = c(site_id, latitude, longitude)
)

names(FUNO_1_coords) <- c("site", "lat", "lon")
FUNO_1_coords$site <- as.character(FUNO_1_coords$site)

FUNO_1_coords <- FUNO_1_coords[
  match(c("2", "3", "5", "7", "8", "9", "10", "11", "15", "16"), FUNO_1_coords$site),
]

# -----------------------------
# 1) pairwise FST matrix
# -----------------------------
FUNO_1_fst <- matrix(
  0,
  nrow = nrow(FUNO_1_coords),
  ncol = nrow(FUNO_1_coords),
  dimnames = list(FUNO_1_coords$site, FUNO_1_coords$site)
)

FUNO_1_fst["3",  "2"] <- 0.026
FUNO_1_fst["5",  c("2","3")] <- c(0.066, 0.039)
FUNO_1_fst["7",  c("2","3","5")] <- c(0.050, 0.022, 0.060)
FUNO_1_fst["8",  c("2","3","5","7")] <- c(0.040, 0.010, 0.061, 0.012)
FUNO_1_fst["9",  c("2","3","5","7","8")] <- c(0.034, 0.009, 0.026, 0.010, 0.028)
FUNO_1_fst["10", c("2","3","5","7","8","9")] <- c(0.060, 0.035, 0.072, 0.059, 0.027, 0.028)
FUNO_1_fst["11", c("2","3","5","7","8","9","10")] <- c(0.033, 0.001, 0.037, 0.020, 0.009, 0.000, 0.050)
FUNO_1_fst["15", c("2","3","5","7","8","9","10","11")] <- c(0.100, 0.149, 0.176, 0.141, 0.119, 0.130, 0.165, 0.144)
FUNO_1_fst["16", c("2","3","5","7","8","9","10","11","15")] <- c(0.033, 0.034, 0.061, 0.039, 0.015, 0.024, 0.052, 0.024, 0.056)

FUNO_1_fst[upper.tri(FUNO_1_fst)] <- t(FUNO_1_fst)[upper.tri(FUNO_1_fst)]
diag(FUNO_1_fst) <- 0

# -----------------------------
# IBD plot: straight-line distance vs FST
# -----------------------------
dist_km <- geosphere::distm(
  as.matrix(FUNO_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist = dist_km[upper.tri(dist_km)],
  fst  = FUNO_1_fst[upper.tri(FUNO_1_fst)]
)

ggplot(ibd_df, aes(dist, fst)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(
    x = "Straight-line distance (km)",
    y = expression(F[ST]),
    title = "IBD: straight-line distance vs FST"
  )

# -----------------------------
# 2) quick checks
# -----------------------------
stopifnot(identical(rownames(FUNO_1_fst), FUNO_1_coords$site))
stopifnot(identical(colnames(FUNO_1_fst), FUNO_1_coords$site))
stopifnot(isTRUE(all.equal(FUNO_1_fst, t(FUNO_1_fst))))

# -----------------------------
# 3) save RData
# -----------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  FUNO_1_fst,
  FUNO_1_coords,
  file = file.path(out_dir, "FUNO-1.RData")
)