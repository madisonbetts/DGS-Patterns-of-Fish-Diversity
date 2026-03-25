# ETOO-1.R
# Etheostoma cf. oophylax (East Fork only)
# Source: Kattawar & Piller 2020, Comparative population genetics...

library(ggplot2)
library(maps)
library(geosphere)

base_path <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"
study_code <- "ETOO-1"
study_dir <- file.path(base_path, study_code)
data_dir <- file.path(study_dir, "data")

dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# site coordinates (Table 2; East Fork only)
# -----------------------------
ETOO_1_coords <- data.frame(
  site = 1:8,
  #site_name = c(
  #  "White Oak Creek at Route 1497",
  #  "EFC Trib at Highway 1497",
  #  "Farley Branch at Tom Taylor Road",
  #  "Middle Fork at Martin Chapel Road",
  #  "Clayton Creek",
  #  "EFC Route 94 Bridge",
  #  "EFC Route 1536 Bridge",
  #  "Bee Creek at Highway 2075"
  #),
  lat = c(
    36.54566, 36.57431, 36.56494, 36.57780,
    36.59217, 36.61257, 36.62041, 36.62622
  ),
  lon = c(
    -88.30559, -88.29654, -88.34426, -88.32805,
    -88.25590, -88.28836, -88.25498, -88.30142
  )
)

# -----------------------------
# pairwise FST matrix (Table 4)
# negative values set to 0
# -----------------------------
ETOO_1_fst <- matrix(0, nrow = 8, ncol = 8)
rownames(ETOO_1_fst) <- colnames(ETOO_1_fst) <- as.character(1:8)

ETOO_1_fst[2,1] <- -0.0217
ETOO_1_fst[3,1] <- -0.1707
ETOO_1_fst[3,2] <- -0.1767
ETOO_1_fst[4,1] <- -0.0057
ETOO_1_fst[4,2] <- 0.0181
ETOO_1_fst[4,3] <- -0.1044
ETOO_1_fst[5,1] <- 0.0190
ETOO_1_fst[5,2] <- 0.0434
ETOO_1_fst[5,3] <- -0.0637
ETOO_1_fst[5,4] <- 0.0128
ETOO_1_fst[6,1] <- 0.0033
ETOO_1_fst[6,2] <- 0.0276
ETOO_1_fst[6,3] <- -0.1904
ETOO_1_fst[6,4] <- -0.0010
ETOO_1_fst[6,5] <- 0.0073
ETOO_1_fst[7,1] <- 0.0123
ETOO_1_fst[7,2] <- 0.0361
ETOO_1_fst[7,3] <- -0.0579
ETOO_1_fst[7,4] <- 0.0082
ETOO_1_fst[7,5] <- 0.0033
ETOO_1_fst[7,6] <- 0.0039
ETOO_1_fst[8,1] <- 0.0061
ETOO_1_fst[8,2] <- 0.0295
ETOO_1_fst[8,3] <- -0.1095
ETOO_1_fst[8,4] <- -0.0018
ETOO_1_fst[8,5] <- 0.0108
ETOO_1_fst[8,6] <- 0.0054
ETOO_1_fst[8,7] <- 0.0066

ETOO_1_fst <- ETOO_1_fst + t(ETOO_1_fst)
diag(ETOO_1_fst) <- 0
ETOO_1_fst[ETOO_1_fst < 0] <- 0

# -----------------------------
# straight-line distance matrix (km)
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = ETOO_1_coords[, c("lon", "lat")],
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- colnames(geo_dist_km) <- as.character(ETOO_1_coords$site)

# -----------------------------
# pairwise dataframe for IBD plot
# -----------------------------
ETOO_1_ibd <- data.frame(
  site1   = rownames(ETOO_1_fst)[row(ETOO_1_fst)[upper.tri(ETOO_1_fst)]],
  site2   = colnames(ETOO_1_fst)[col(ETOO_1_fst)[upper.tri(ETOO_1_fst)]],
  fst     = ETOO_1_fst[upper.tri(ETOO_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# map of sites (KY + TN only)
# -----------------------------
states_map <- map_data("state")

ky_tn <- subset(states_map, region %in% c("kentucky", "tennessee"))

map_xlim <- c(-88.5, -88.2)
map_ylim <- c(36.5, 36.7)

ggplot() +
  geom_polygon(
    data = ky_tn,
    aes(x = long, y = lat, group = group),
    fill = "grey92", color = "grey40", linewidth = 0.3
  ) +
  geom_point(
    data = ETOO_1_coords,
    aes(x = lon, y = lat),
    size = 2.2
  ) +
  geom_text(
    data = ETOO_1_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.01,
    size = 3
  ) +
  coord_fixed(xlim = map_xlim, ylim = map_ylim) +
  theme_classic() +
  labs(
    title = "ETOO-1 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )

# -----------------------------
# IBD plot
# -----------------------------
ggplot(ETOO_1_ibd, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.3, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    title = "ETOO-1 IBD QC plot",
    x = "Euclidean distance (km)",
    y = expression(F[ST])
  )

save(
  ETOO_1_fst,
  ETOO_1_coords,
  #geo_dist_km,
  #ETOO_1_ibd,
  file = file.path(data_dir, "ETOO-1.RData")
)