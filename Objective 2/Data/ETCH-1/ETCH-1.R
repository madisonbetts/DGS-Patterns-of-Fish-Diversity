# ETCH-1.R
# Etheostoma chienense
# Source: Kattawar & Piller 2020, Comparative population genetics...

library(ggplot2)
library(maps)
library(geosphere)

base_path <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"
study_code <- "ETCH-1"
study_dir <- file.path(base_path, study_code)
data_dir <- file.path(study_dir, "data")

dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# site coordinates (Table 2)
# -----------------------------
ETCH_1_coords <- data.frame(
  site = 1:9,
  #site_name = c(
  #  "BDC at Pea Ridge Road",
  #  "Jackson Creek at Lawrence Road",
  #  "BDC at Route 1283 Bridge",
  #  "BDC at Route 307 Bridge",
  #  "BDC at Howell Road Bridge",
  #  "BDC at US 51 Bridge",
  #  "BDC 1 km upstream of KY 239",
  #  "Little BDC at KY 1125",
  #  "Little BDC at KY 1706"
  #),
  lat = c(
    36.56878, 36.58436, 36.60468, 36.61221, 36.62812,
    36.62875, 36.62227, 36.53374, 36.52709
  ),
  lon = c(
    -88.74555, -88.79022, -88.81670, -88.87130, -88.92690,
    -88.96368, -89.01769, -88.96138, -88.94261
  )
)

# -----------------------------
# pairwise FST matrix (Table 3)
# negative values set to 0
# -----------------------------
ETCH_1_fst <- matrix(0, nrow = 9, ncol = 9)
rownames(ETCH_1_fst) <- colnames(ETCH_1_fst) <- as.character(1:9)

ETCH_1_fst[2,1] <- 0.0271
ETCH_1_fst[3,1] <- 0.0159
ETCH_1_fst[3,2] <- 0.0132
ETCH_1_fst[4,1] <- 0.0108
ETCH_1_fst[4,2] <- 0.0002
ETCH_1_fst[4,3] <- -0.0010
ETCH_1_fst[5,1] <- 0.0063
ETCH_1_fst[5,2] <- 0.0072
ETCH_1_fst[5,3] <- 0.0005
ETCH_1_fst[5,4] <- -0.0085
ETCH_1_fst[6,1] <- 0.0360
ETCH_1_fst[6,2] <- 0.0296
ETCH_1_fst[6,3] <- 0.0237
ETCH_1_fst[6,4] <- 0.0293
ETCH_1_fst[6,5] <- -0.0076
ETCH_1_fst[7,1] <- 0.0651
ETCH_1_fst[7,2] <- 0.0819
ETCH_1_fst[7,3] <- 0.0588
ETCH_1_fst[7,4] <- 0.0609
ETCH_1_fst[7,5] <- 0.0490
ETCH_1_fst[7,6] <- -0.0539
ETCH_1_fst[8,1] <- 0.0064
ETCH_1_fst[8,2] <- -0.0009
ETCH_1_fst[8,3] <- -0.0007
ETCH_1_fst[8,4] <- -0.0033
ETCH_1_fst[8,5] <- 0.0002
ETCH_1_fst[8,6] <- 0.0023
ETCH_1_fst[8,7] <- 0.0543
ETCH_1_fst[9,1] <- 0.0133
ETCH_1_fst[9,2] <- -0.0041
ETCH_1_fst[9,3] <- -0.0003
ETCH_1_fst[9,4] <- 0.0012
ETCH_1_fst[9,5] <- -0.0083
ETCH_1_fst[9,6] <- 0.0253
ETCH_1_fst[9,7] <- 0.0716
ETCH_1_fst[9,8] <- -0.0241

ETCH_1_fst <- ETCH_1_fst + t(ETCH_1_fst)
diag(ETCH_1_fst) <- 0
ETCH_1_fst[ETCH_1_fst < 0] <- 0

# -----------------------------
# straight-line distance matrix (km)
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = ETCH_1_coords[, c("lon", "lat")],
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- colnames(geo_dist_km) <- as.character(ETCH_1_coords$site_id)

# -----------------------------
# pairwise dataframe for IBD plot
# -----------------------------
ETCH_1_ibd <- data.frame(
  site1   = rownames(ETCH_1_fst)[row(ETCH_1_fst)[upper.tri(ETCH_1_fst)]],
  site2   = colnames(ETCH_1_fst)[col(ETCH_1_fst)[upper.tri(ETCH_1_fst)]],
  fst     = ETCH_1_fst[upper.tri(ETCH_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# map of sites (US + Canada context)
# -----------------------------
# -----------------------------
# map of sites (KY + TN only)
# -----------------------------
states_map <- map_data("state")

ky_tn <- subset(states_map, region %in% c("kentucky", "tennessee"))

map_xlim <- c(-89, -88.6)
map_ylim <- c(36.4, 36.7)

ggplot() +
  geom_polygon(
    data = ky_tn,
    aes(x = long, y = lat, group = group),
    fill = "grey92", color = "grey40", linewidth = 0.3
  ) +
  geom_point(
    data = ETCH_1_coords,
    aes(x = lon, y = lat),
    size = 2.2
  ) +
  geom_text(
    data = ETCH_1_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.01,
    size = 3
  ) +
  coord_fixed(xlim = map_xlim, ylim = map_ylim) +
  theme_classic() +
  labs(
    title = "ETCH-1 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )
# -----------------------------
# IBD plot
# -----------------------------
ggplot(ETCH_1_ibd, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.3, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    title = "ETCH-1 IBD QC plot",
    x = "Euclidean distance (km)",
    y = expression(F[ST])
  )



save(
  ETCH_1_fst,
  ETCH_1_coords,
  #geo_dist_km,
  #ETCH_1_ibd,
  file = file.path(data_dir, "ETCH-1.RData")
)
