# -----------------------------
# ETSA-1 Kentucky arrow darter
# site coordinates + FST matrix
# -----------------------------

library(ggplot2)
library(geosphere)

# directory where everything is
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETSA-1/"

# -----------------------------
# 0) site coordinates
# -----------------------------
ETSA_1_coords <- data.frame(
  site = as.character(1:9),
  lat = c(
    37.4758623,
    37.0143963,
    37.1159372,
    37.34798,
    37.2199002,
    37.6323547,
    37.6303108,
    37.466745,
    37.7647083
  ),
  lon = c(
    -83.8512715,
    -83.5452998,
    -83.5838696,
    -83.56304,
    -83.4495912,
    -83.3544548,
    -83.6661473,
    -83.1503212,
    -83.570909
  ),
  stringsAsFactors = FALSE
)

ggplot(ETSA_1_coords, aes(lon, lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.03) +
  coord_equal() +
  theme_bw()

# -----------------------------
# 1) pairwise FST matrix
# -----------------------------
ETSA_1_fst <- matrix(
  0,
  9,
  9,
  dimnames = list(ETSA_1_coords$site, ETSA_1_coords$site)
)

ETSA_1_fst["2","1"] <- 0.19

ETSA_1_fst["3", c("1","2")] <- c(0.38, 0.34)

ETSA_1_fst["4", c("1","2","3")] <- c(0.14, 0.13, 0.40)

ETSA_1_fst["5", c("1","2","3","4")] <- c(0.34, 0.41, 0.63, 0.39)

ETSA_1_fst["6", c("1","2","3","4","5")] <- c(0.29, 0.27, 0.47, 0.28, 0.47)

ETSA_1_fst["7", c("1","2","3","4","5","6")] <- c(0.19, 0.21, 0.46, 0.20, 0.44, 0.24)

ETSA_1_fst["8", c("1","2","3","4","5","6","7")] <- c(0.23, 0.29, 0.55, 0.22, 0.40, 0.27, 0.29)

ETSA_1_fst["9", c("1","2","3","4","5","6","7","8")] <- c(0.49, 0.49, 0.74, 0.47, 0.70, 0.52, 0.58, 0.49)

ETSA_1_fst[upper.tri(ETSA_1_fst)] <- t(ETSA_1_fst)[upper.tri(ETSA_1_fst)]
diag(ETSA_1_fst) <- 0

# -----------------------------
# IBD plot: straight-line distance vs FST
# -----------------------------
dist_km <- geosphere::distm(
  as.matrix(ETSA_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist = dist_km[upper.tri(dist_km)],
  fst  = ETSA_1_fst[upper.tri(ETSA_1_fst)]
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
stopifnot(identical(rownames(ETSA_1_fst), ETSA_1_coords$site))
stopifnot(identical(colnames(ETSA_1_fst), ETSA_1_coords$site))
stopifnot(isTRUE(all.equal(ETSA_1_fst, t(ETSA_1_fst))))

# -----------------------------
# 3) save RData
# -----------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  ETSA_1_fst,
  ETSA_1_coords,
  file = file.path(out_dir, "ETSA-1.RData")
)