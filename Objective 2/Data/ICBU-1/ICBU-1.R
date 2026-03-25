# ICBU-1.R
# Ictiobus bubalus
# Source: Kratina et al. 2023, Transactions of the American Fisheries Society
# Pairwise FST values transcribed from Table 5
# Representative reach coordinates inferred from Figure 1 and dam/section locations

library(ggplot2)
library(maps)
library(geosphere)

base_path <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"
study_code <- "ICBU-1"
study_dir <- file.path(base_path, study_code)
data_dir <- file.path(study_dir, "data")

dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# representative section coordinates
# section order follows Table 5 / Table 3:
# 1 = JBR, 2 = MFR, 3 = CL, 4 = LAR
# coords are representative reach centroids for grouped river sections
# -----------------------------
ICBU_1_coords <- data.frame(
  site = 1:4,
  #section = c("JBR", "MFR", "CL", "LAR"),
  lat = c(
    32.52,   # JBR
    32.18,   # MFR
    31.82,   # CL
    31.35    # LAR
  ),
  lon = c(
    -86.64,  # JBR
    -87.10,  # MFR
    -87.53,  # CL
    -87.78   # LAR
  )
)

# -----------------------------
# pairwise FST matrix (Table 5, upper diagonal)
# negative values set to 0
# -----------------------------
ICBU_1_fst <- matrix(c(
  0.00000, 0.00232, 0.00051, 0.00149,
  0.00232, 0.00000, 0.00203, 0.00206,
  0.00051, 0.00203, 0.00000, 0.00045,
  0.00149, 0.00206, 0.00045, 0.00000
), nrow = 4, byrow = TRUE)

rownames(ICBU_1_fst) <- colnames(ICBU_1_fst) <- as.character(1:4)
ICBU_1_fst[ICBU_1_fst < 0] <- 0
diag(ICBU_1_fst) <- 0

# -----------------------------
# straight-line distance matrix (km)
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = ICBU_1_coords[, c("lon", "lat")],
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- colnames(geo_dist_km) <- as.character(ICBU_1_coords$site)

# -----------------------------
# pairwise dataframe for IBD plot
# -----------------------------
ICBU_1_ibd <- data.frame(
  site1   = rownames(ICBU_1_fst)[row(ICBU_1_fst)[upper.tri(ICBU_1_fst)]],
  site2   = colnames(ICBU_1_fst)[col(ICBU_1_fst)[upper.tri(ICBU_1_fst)]],
  fst     = ICBU_1_fst[upper.tri(ICBU_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# map of sites (AL only)
# -----------------------------
states_map <- map_data("state")
al_map <- subset(states_map, region %in% c("alabama"))

map_xlim <- c(-88.2, -86.0)
map_ylim <- c(31.0, 33.0)

ggplot() +
  geom_polygon(
    data = al_map,
    aes(x = long, y = lat, group = group),
    fill = "grey92", color = "grey40", linewidth = 0.3
  ) +
  geom_point(
    data = ICBU_1_coords,
    aes(x = lon, y = lat),
    size = 2.2
  ) +
  geom_text(
    data = ICBU_1_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.05,
    size = 3
  ) +
  coord_fixed(xlim = map_xlim, ylim = map_ylim) +
  theme_classic() +
  labs(
    title = "ICBU-1 representative reach points",
    x = "Longitude",
    y = "Latitude"
  )

# -----------------------------
# IBD plot
# -----------------------------
ggplot(ICBU_1_ibd, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.3, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    title = "ICBU-1 IBD QC plot",
    x = "Euclidean distance (km)",
    y = expression(F[ST])
  )

save(
  ICBU_1_fst,
  ICBU_1_coords,
  #geo_dist_km,
  #ICBU_1_ibd,
  file = file.path(data_dir, "ICBU-1.RData")
)
