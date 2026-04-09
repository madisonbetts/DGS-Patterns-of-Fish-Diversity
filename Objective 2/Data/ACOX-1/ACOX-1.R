# ============================================================
# ACOX-1 | Atlantic sturgeon
# updated with user-provided coordinates (hardcoded)
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(geosphere)
})

study_code <- "ACOX-1"
study_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ACOX-1"
out_dir <- file.path(study_dir, "data")

if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

# -----------------------------
# coordinates (USER PROVIDED)
# -----------------------------
ACOX_1_coords_full <- data.frame(
  site = 1:5,
  site_name = c("Edisto River", "Savannah River", "Ogeechee River", "Altamaha River", "Satilla River"),
  lat = c(33.0581517, 32.4278842, 32.1050282, 31.66553, 31.087228),
  lon = c(-80.5867839, -81.2131894, -81.3983092, -81.8555354, -81.9258256),
  stringsAsFactors = FALSE
)

ACOX_1_coords <- ACOX_1_coords_full %>%
  select(site, lat, lon)

# -----------------------------
# FST matrix (real, from paper)
# -----------------------------
ACOX_1_fst <- matrix(0, nrow = 5, ncol = 5)
rownames(ACOX_1_fst) <- as.character(1:5)
colnames(ACOX_1_fst) <- as.character(1:5)

ACOX_1_fst[2, 1] <- 0.024
ACOX_1_fst[3, 1] <- 0.033
ACOX_1_fst[3, 2] <- 0.028
ACOX_1_fst[4, 1] <- 0.024
ACOX_1_fst[4, 2] <- 0.013
ACOX_1_fst[4, 3] <- 0.039
ACOX_1_fst[5, 1] <- 0.028
ACOX_1_fst[5, 2] <- 0.022
ACOX_1_fst[5, 3] <- 0.020
ACOX_1_fst[5, 4] <- 0.019

ACOX_1_fst <- ACOX_1_fst + t(ACOX_1_fst)
diag(ACOX_1_fst) <- 0
ACOX_1_fst[ACOX_1_fst < 0] <- 0

# -----------------------------
# map
# -----------------------------
world_df <- ggplot2::map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

ggplot() +
  geom_polygon(data = world_df, aes(x = long, y = lat, group = group),
               fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_point(data = ACOX_1_coords_full, aes(x = lon, y = lat), size = 2) +
  geom_text(data = ACOX_1_coords_full,
            aes(x = lon, y = lat, label = paste0(site, " ", site_name)),
            nudge_y = 0.12, size = 3) +
  coord_quickmap() +
  theme_bw()

# -----------------------------
# IBD
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(ACOX_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  fst = ACOX_1_fst[upper.tri(ACOX_1_fst)]
)

ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw()

# -----------------------------
# save
# -----------------------------
save(
  ACOX_1_fst,
  ACOX_1_coords,
  file = file.path(out_dir, "ACOX-1.RData")
)
