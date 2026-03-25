# POSP-2.R
# Polyodon spathula
# Source: Kratina et al. 2023, Transactions of the American Fisheries Society
# Pairwise FST values transcribed from Table 5
# Representative reach coordinates inferred from Figure 1 and dam/section locations

library(ggplot2)
library(maps)
library(geosphere)

base_path <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"
study_code <- "POSP-2"
study_dir <- file.path(base_path, study_code)
data_dir <- file.path(study_dir, "data")

dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# representative section coordinates
# section order follows Table 5:
# 1 = CSA, 2 = JBR, 3 = MFR, 4 = CL, 5 = LAR, 6 = TOM
# coords are representative reach centroids for grouped river sections
# -----------------------------
POSP_2_coords <- data.frame(
  site = 1:6,
  #section = c("CSA", "JBR", "MFR", "CL", "LAR", "TOM"),
  lat = c(
    32.58,   # CSA
    32.52,   # JBR
    32.18,   # MFR
    31.82,   # CL
    31.35,   # LAR
    31.50    # TOM
  ),
  lon = c(
    -86.33,  # CSA
    -86.64,  # JBR
    -87.10,  # MFR
    -87.53,  # CL
    -87.78,  # LAR
    -88.13   # TOM
  )
)

# -----------------------------
# pairwise FST matrix (Table 5)
# negative values set to 0
# -----------------------------
POSP_2_fst <- matrix(c(
  0.00000, 0.00322, 0.00013, 0.00279, 0.00397, -0.00017,
  0.00322, 0.00000, 0.00056, 0.00127, -0.00066, 0.00050,
  0.00013, 0.00056, 0.00000, 0.00064, 0.00057, -0.00043,
  0.00279, 0.00127, 0.00064, 0.00000, 0.00109, 0.00125,
  0.00397, -0.00066, 0.00057, 0.00109, 0.00000, -0.00007,
  -0.00017, 0.00050, -0.00043, 0.00125, -0.00007, 0.00000
), nrow = 6, byrow = TRUE)

rownames(POSP_2_fst) <- colnames(POSP_2_fst) <- as.character(1:6)
POSP_2_fst[POSP_2_fst < 0] <- 0
diag(POSP_2_fst) <- 0

# -----------------------------
# straight-line distance matrix (km)
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = POSP_2_coords[, c("lon", "lat")],
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- colnames(geo_dist_km) <- as.character(POSP_2_coords$site)

# -----------------------------
# pairwise dataframe for IBD plot
# -----------------------------
POSP_2_ibd <- data.frame(
  site1   = rownames(POSP_2_fst)[row(POSP_2_fst)[upper.tri(POSP_2_fst)]],
  site2   = colnames(POSP_2_fst)[col(POSP_2_fst)[upper.tri(POSP_2_fst)]],
  fst     = POSP_2_fst[upper.tri(POSP_2_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# map of sites (AL only)
# -----------------------------
states_map <- map_data("state")
al_map <- subset(states_map, region %in% c("alabama"))

map_xlim <- c(-88.6, -85.8)
map_ylim <- c(30.8, 33.4)

ggplot() +
  geom_polygon(
    data = al_map,
    aes(x = long, y = lat, group = group),
    fill = "grey92", color = "grey40", linewidth = 0.3
  ) +
  geom_point(
    data = POSP_2_coords,
    aes(x = lon, y = lat),
    size = 2.2
  ) +
  geom_text(
    data = POSP_2_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.05,
    size = 3
  ) +
  coord_fixed(xlim = map_xlim, ylim = map_ylim) +
  theme_classic() +
  labs(
    title = "POSP-2 representative reach points",
    x = "Longitude",
    y = "Latitude"
  )

# -----------------------------
# IBD plot
# -----------------------------
ggplot(POSP_2_ibd, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.3, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    title = "POSP-2 IBD QC plot",
    x = "Euclidean distance (km)",
    y = expression(F[ST])
  )

save(
  POSP_2_fst,
  POSP_2_coords,
  #geo_dist_km,
  #POSP_2_ibd,
  file = file.path(data_dir, "POSP-2.RData")
)
