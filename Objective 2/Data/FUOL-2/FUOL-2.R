# -----------------------------
# FUOL-2 blackspotted topminnow
# site coordinates + FST matrix + IBD plot
# -----------------------------

library(ggplot2)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/FUOL-2"

# -----------------------------
# site key
# -----------------------------
# 1 = Brushy Creek
# 2 = Baker Creek
# 3 = Gasconade Mid-Stream
# 4 = Gasconade Large River
# 5 = Big River
# 6 = Dry Fork
# 7 = Meramec Mid-Stream
# 8 = Meramec Large River

# -----------------------------
# 0) site coordinates
# -----------------------------
FUOL_2_coords <- data.frame(
  site = as.character(1:8),
  lat = c(
    37.333000,  # 1
    37.168683,  # 2
    37.938267,  # 3
    38.563500,  # 4
    37.753833,  # 5
    37.853267,  # 6
    38.481333,  # 7
    38.437533   # 8
  ),
  lon = c(
    -91.957567, # 1
    -92.614433, # 2
    -91.979650, # 3
    -91.598167, # 4
    -90.884617, # 5
    -91.667350, # 6
    -90.630433, # 7
    -90.346350  # 8
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map of sampling locations
# -----------------------------
ggplot(FUOL_2_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.04, size = 3.5) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "FUOL-2 sampling locations"
  )

# -----------------------------
# 2) geographic distance matrix
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(FUOL_2_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- FUOL_2_coords$site
colnames(geo_dist_km) <- FUOL_2_coords$site

# -----------------------------
# 3) pairwise FST matrix
# -----------------------------
FUOL_2_fst <- matrix(
  0,
  nrow = nrow(FUOL_2_coords),
  ncol = nrow(FUOL_2_coords),
  dimnames = list(FUOL_2_coords$site, FUOL_2_coords$site)
)

FUOL_2_fst["2", "1"] <- 0.2424

FUOL_2_fst["3", c("1","2")] <- c(0.0796, 0.2084)

FUOL_2_fst["4", c("1","2","3")] <- c(0.1120, 0.1798, 0.0574)

FUOL_2_fst["5", c("1","2","3","4")] <- c(0.3628, 0.3511, 0.3432, 0.2131)

FUOL_2_fst["6", c("1","2","3","4","5")] <- c(0.3349, 0.3221, 0.3097, 0.1726, 0.1096)

FUOL_2_fst["7", c("1","2","3","4","5","6")] <- c(0.3299, 0.3155, 0.3062, 0.1669, 0.1017, 0.0557)

FUOL_2_fst["8", c("1","2","3","4","5","6","7")] <- c(0.3148, 0.2977, 0.2895, 0.1541, 0.0903, 0.0440, 0.0144)

# reflect to upper triangle
FUOL_2_fst[upper.tri(FUOL_2_fst)] <- t(FUOL_2_fst)[upper.tri(FUOL_2_fst)]
diag(FUOL_2_fst) <- 0

# -----------------------------
# 4) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(FUOL_2_fst)[row(FUOL_2_fst)[upper.tri(FUOL_2_fst)]],
  site2   = colnames(FUOL_2_fst)[col(FUOL_2_fst)[upper.tri(FUOL_2_fst)]],
  fst     = FUOL_2_fst[upper.tri(FUOL_2_fst)],
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
    title = "FUOL-2 isolation by distance"
  )

# -----------------------------
# 6) quick checks
# -----------------------------
stopifnot(identical(rownames(FUOL_2_fst), FUOL_2_coords$site))
stopifnot(identical(colnames(FUOL_2_fst), FUOL_2_coords$site))
stopifnot(isTRUE(all.equal(FUOL_2_fst, t(FUOL_2_fst))))

# -----------------------------
# 7) save RData (ONLY fst + coords)
# -----------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  FUOL_2_fst,
  FUOL_2_coords,
  file = file.path(out_dir, "FUOL-2.RData")
)