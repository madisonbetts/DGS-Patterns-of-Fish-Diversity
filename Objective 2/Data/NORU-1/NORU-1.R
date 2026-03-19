#--------------------------
# NORU darter sites and fst 
#--------------------------
library(ggplot2)
library(ggrepel)
library(geosphere)

# sites and coordinates
bp_sites <- data.frame(
  site = 1:18,
  lat = c(
    32.002947,
    32.010479,
    32.000078,
    32.012273,
    31.952733,
    31.976405,
    32.023391,
    32.032717,
    32.023033,
    32.036662,
    31.980709,
    31.981068,
    31.949146,
    31.930854,
    31.897856,
    31.908257,
    31.900725,
    31.901442
  ),
  lon = c(
    -90.902614,
    -90.807137,
    -90.789008,
    -90.767254,
    -90.720120,
    -90.720120,
    -90.720120,
    -90.683863,
    -90.670568,
    -90.653648,
    -90.651231,
    -90.608931,
    -90.584760,
    -90.564214,
    -90.549711,
    -90.543668,
    -90.587177,
    -90.497743
  ),
  stringsAsFactors = FALSE
)

bp_sites

# -----------------------------
# map of estimated sampling sites
# -----------------------------
ggplot(bp_sites, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = site),
    size = 4,
    box.padding = 0.3,
    point.padding = 0.3
  ) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "Estimated sampling locations (Bayou Pierre system)"
  )

# -----------------------------
# FST matrix for Nothonotus rubris
# row/col names match bp_sites$site
# -----------------------------
sites <- as.character(bp_sites$site)

fst <- matrix(0, nrow = length(sites), ncol = length(sites),
              dimnames = list(sites, sites))

fst["2",  "1"] <- 0.058

fst["3",  c("1","2")] <- c(0.053, 0.041)

fst["4",  c("1","2","3")] <- c(0.094, 0.087, 0.077)

fst["5",  c("1","2","3","4")] <- c(0.064, 0.050, 0.050, 0.089)

fst["6",  c("1","2","3","4","5")] <- c(0.089, 0.088, 0.073, 0.115, 0.080)

fst["7",  c("1","2","3","4","5","6")] <- c(0.080, 0.053, 0.055, 0.098, 0.064, 0.102)

fst["8",  c("1","2","3","4","5","6","7")] <- c(0.075, 0.066, 0.057, 0.083, 0.069, 0.092, 0.073)

fst["9",  c("1","2","3","4","5","6","7","8")] <- c(0.075, 0.067, 0.061, 0.084, 0.070, 0.097, 0.071, 0.058)

fst["10", c("1","2","3","4","5","6","7","8","9")] <- c(0.122, 0.116, 0.103, 0.131, 0.112, 0.134, 0.115, 0.098, 0.109)

fst["11", c("1","2","3","4","5","6","7","8","9","10")] <- c(0.066, 0.060, 0.049, 0.075, 0.055, 0.084, 0.074, 0.059, 0.061, 0.104)

fst["12", c("1","2","3","4","5","6","7","8","9","10","11")] <- c(0.070, 0.059, 0.055, 0.086, 0.063, 0.088, 0.075, 0.071, 0.073, 0.112, 0.061)

fst["13", c("1","2","3","4","5","6","7","8","9","10","11","12")] <- c(0.056, 0.050, 0.041, 0.078, 0.047, 0.074, 0.064, 0.057, 0.057, 0.104, 0.045, 0.057)

fst["14", c("1","2","3","4","5","6","7","8","9","10","11","12","13")] <- c(0.074, 0.067, 0.063, 0.104, 0.068, 0.098, 0.081, 0.083, 0.080, 0.131, 0.070, 0.078, 0.062)

fst["15", c("1","2","3","4","5","6","7","8","9","10","11","12","13","14")] <- c(0.084, 0.069, 0.068, 0.101, 0.070, 0.109, 0.087, 0.083, 0.087, 0.128, 0.072, 0.079, 0.066, 0.084)

fst["16", c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15")] <- c(0.092, 0.086, 0.077, 0.102, 0.080, 0.108, 0.098, 0.087, 0.096, 0.133, 0.077, 0.087, 0.074, 0.101, 0.074)

fst["17", c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16")] <- c(0.106, 0.101, 0.092, 0.123, 0.097, 0.118, 0.119, 0.106, 0.110, 0.144, 0.088, 0.106, 0.088, 0.109, 0.081, 0.089)

fst["18", c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")] <- c(0.103, 0.095, 0.089, 0.122, 0.095, 0.122, 0.115, 0.097, 0.098, 0.136, 0.089, 0.096, 0.082, 0.103, 0.079, 0.092, 0.102)

fst[upper.tri(fst)] <- t(fst)[upper.tri(fst)]

fst

# -----------------------------
# pairwise straight-line distance (km)
# great-circle distance from coordinates
# -----------------------------
coords <- bp_sites[, c("lon", "lat")]

geo_dist <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
geo_dist <- as.matrix(geo_dist)

rownames(geo_dist) <- as.character(bp_sites$site)
colnames(geo_dist) <- as.character(bp_sites$site)

round(geo_dist, 2)

# -----------------------------
# build pairwise dataframe
# lower triangle only, no duplicates
# -----------------------------
idx <- lower.tri(fst)

ibd_df <- data.frame(
  site1 = rownames(fst)[row(fst)[idx]],
  site2 = colnames(fst)[col(fst)[idx]],
  fst = fst[idx],
  dist_km = geo_dist[idx]
)

ibd_df

# optional pair labels
ibd_df$pair <- paste(ibd_df$site1, ibd_df$site2, sep = "-")

# -----------------------------
# straight-line distance vs FST
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Straight-line distance (km)",
    y = expression(F[ST]),
    title = "Isolation by distance for Nothonotus rubris"
  )