# ============================================================
# COBA-4 | Rocky Mountain sculpin (Cottus sp.)
# Ruppert et al. 2017, Conservation Genetics
# ============================================================

library(ggplot2)
library(maps)
library(mapdata)
library(geosphere)
library(dplyr)

# -----------------------------
# 1) robust working directory
# -----------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"

study_code <- "COBA-4"
study_dir <- file.path(base_dir, study_code)

if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE)
setwd(study_dir)

if (!dir.exists("data")) dir.create("data")

# -----------------------------
# 2) site coordinates
# -----------------------------
COBA_4_coords <- data.frame(
  site = 1:37,
  lat = c(49.027, 49.021, 49.025, 49.028, 49.031, 49.034, 49.037, 49.018, 49.047, 49.073, 49.089, 49.11, 49.128, 49.143, 49.161, 49.177, 49.191, 49.201, 49.209, 49.108, 49.127, 49.143, 49.156, 49.175, 49.197, 49.224, 49.248, 49.271, 49.04, 49.05, 49.064, 49.077, 49.093, 49.107, 49.122, 49.137, 49.149),
  lon = c(-114.565, -114.479, -114.476, -114.473, -114.47, -114.468, -114.466, -113.813, -113.781, -113.725, -113.673, -113.629, -113.598, -113.574, -113.548, -113.519, -113.496, -113.476, -113.458, -113.509, -113.47, -113.435, -113.399, -113.369, -113.358, -113.366, -113.388, -113.414, -112.91, -112.892, -112.874, -112.858, -112.828, -112.792, -112.747, -112.692, -112.565)
)

# -----------------------------
# 3) basin metadata
# -----------------------------
COBA_4_site_meta <- data.frame(
  site = 1:37,
  basin = c('Flathead', 'Flathead', 'Flathead', 'Flathead', 'Flathead', 'Flathead', 'Flathead', 'Lee Creek', 'Lee Creek', 'Lee Creek', 'Lee Creek', 'Lee Creek', 'Lee Creek', 'Lee Creek', 'Lee Creek', 'Lee Creek', 'Lee Creek', 'Lee Creek', 'Lee Creek', 'St. Mary River', 'St. Mary River', 'St. Mary River', 'St. Mary River', 'St. Mary River', 'St. Mary River', 'St. Mary River', 'St. Mary River', 'St. Mary River', 'North Milk River', 'North Milk River', 'North Milk River', 'North Milk River', 'North Milk River', 'North Milk River', 'North Milk River', 'North Milk River', 'North Milk River'),
  stringsAsFactors = FALSE
)

# -----------------------------
# 4) pairwise FST matrix
# values transcribed row-by-row from COBA_FST.xlsx
# lower diagonal only, then mirrored to upper triangle
# diagonal filled with zero
# -----------------------------
fill_row <- function(mat, row, vals) {
  stopifnot(length(vals) == (row - 1))
  mat[row, seq_len(row - 1)] <- vals
  mat
}

COBA_4_fst <- matrix(0, nrow = 37, ncol = 37)

COBA_4_fst <- fill_row(COBA_4_fst,  2, c(0.02))
COBA_4_fst <- fill_row(COBA_4_fst,  3, c(0.01, 0.02))
COBA_4_fst <- fill_row(COBA_4_fst,  4, c(0.01, 0.00, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst,  5, c(0.02, 0.00, 0.00, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst,  6, c(0.05, 0.01, 0.01, 0.00, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst,  7, c(0.05, 0.00, 0.02, 0.00, 0.00, 0.01))
COBA_4_fst <- fill_row(COBA_4_fst,  8, c(0.44, 0.42, 0.43, 0.42, 0.43, 0.41, 0.44))
COBA_4_fst <- fill_row(COBA_4_fst,  9, c(0.41, 0.40, 0.41, 0.40, 0.41, 0.40, 0.42, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 10, c(0.43, 0.41, 0.42, 0.41, 0.43, 0.41, 0.44, 0.01, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 11, c(0.41, 0.39, 0.40, 0.39, 0.41, 0.40, 0.42, 0.02, 0.01, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 12, c(0.40, 0.39, 0.40, 0.39, 0.40, 0.39, 0.41, 0.01, 0.00, 0.00, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 13, c(0.40, 0.39, 0.40, 0.39, 0.40, 0.39, 0.41, 0.01, 0.00, 0.00, 0.00, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 14, c(0.40, 0.38, 0.39, 0.38, 0.40, 0.38, 0.40, 0.02, 0.01, 0.00, 0.01, 0.00, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 15, c(0.43, 0.41, 0.42, 0.41, 0.43, 0.42, 0.44, 0.03, 0.02, 0.01, 0.01, 0.01, 0.00, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 16, c(0.40, 0.38, 0.39, 0.38, 0.40, 0.38, 0.40, 0.02, 0.01, 0.00, 0.01, 0.00, 0.00, 0.00, 0.01))
COBA_4_fst <- fill_row(COBA_4_fst, 17, c(0.39, 0.37, 0.39, 0.37, 0.39, 0.37, 0.40, 0.02, 0.01, 0.00, 0.01, 0.00, 0.00, 0.00, 0.01, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 18, c(0.42, 0.40, 0.42, 0.40, 0.43, 0.41, 0.44, 0.05, 0.04, 0.00, 0.01, 0.01, 0.01, 0.01, 0.01, 0.00, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 19, c(0.37, 0.35, 0.36, 0.35, 0.37, 0.35, 0.37, 0.07, 0.07, 0.05, 0.04, 0.04, 0.05, 0.04, 0.04, 0.04, 0.02, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 20, c(0.38, 0.37, 0.38, 0.37, 0.39, 0.37, 0.39, 0.07, 0.06, 0.05, 0.04, 0.03, 0.04, 0.04, 0.03, 0.03, 0.02, 0.00, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 21, c(0.37, 0.35, 0.36, 0.35, 0.37, 0.35, 0.37, 0.06, 0.05, 0.04, 0.04, 0.03, 0.04, 0.02, 0.03, 0.03, 0.02, 0.00, 0.00, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 22, c(0.37, 0.35, 0.36, 0.35, 0.37, 0.36, 0.38, 0.06, 0.05, 0.04, 0.03, 0.03, 0.04, 0.02, 0.03, 0.02, 0.01, 0.00, 0.00, 0.00, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 23, c(0.36, 0.35, 0.36, 0.35, 0.37, 0.36, 0.37, 0.06, 0.04, 0.04, 0.03, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.00, 0.00, 0.00, 0.00, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 24, c(0.37, 0.35, 0.36, 0.35, 0.37, 0.36, 0.38, 0.05, 0.05, 0.03, 0.03, 0.02, 0.03, 0.03, 0.03, 0.02, 0.01, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 25, c(0.37, 0.36, 0.37, 0.35, 0.38, 0.36, 0.38, 0.05, 0.04, 0.03, 0.03, 0.02, 0.03, 0.01, 0.02, 0.01, 0.02, 0.00, 0.01, 0.00, 0.00, 0.00, 0.00, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 26, c(0.37, 0.35, 0.37, 0.35, 0.38, 0.36, 0.38, 0.05, 0.04, 0.03, 0.02, 0.02, 0.02, 0.02, 0.01, 0.02, 0.01, 0.00, 0.00, 0.00, 0.01, 0.00, 0.00, 0.00, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 27, c(0.37, 0.35, 0.36, 0.35, 0.37, 0.35, 0.38, 0.06, 0.06, 0.03, 0.04, 0.03, 0.04, 0.04, 0.04, 0.03, 0.02, 0.00, 0.01, 0.01, 0.01, 0.01, 0.01, 0.00, 0.01, 0.01))
COBA_4_fst <- fill_row(COBA_4_fst, 28, c(0.39, 0.37, 0.38, 0.37, 0.39, 0.37, 0.39, 0.06, 0.06, 0.04, 0.04, 0.03, 0.04, 0.04, 0.03, 0.02, 0.01, 0.00, 0.00, 0.00, 0.01, 0.00, 0.01, 0.00, 0.01, 0.01, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 29, c(0.47, 0.45, 0.45, 0.45, 0.46, 0.45, 0.47, 0.18, 0.19, 0.16, 0.16, 0.15, 0.16, 0.15, 0.17, 0.14, 0.13, 0.12, 0.15, 0.13, 0.13, 0.13, 0.15, 0.13, 0.13, 0.16, 0.11, 0.12))
COBA_4_fst <- fill_row(COBA_4_fst, 30, c(0.46, 0.45, 0.45, 0.45, 0.46, 0.45, 0.47, 0.19, 0.19, 0.16, 0.16, 0.15, 0.16, 0.15, 0.17, 0.14, 0.13, 0.12, 0.16, 0.14, 0.13, 0.13, 0.15, 0.13, 0.13, 0.16, 0.12, 0.12, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 31, c(0.45, 0.43, 0.44, 0.43, 0.45, 0.44, 0.45, 0.19, 0.19, 0.16, 0.16, 0.15, 0.16, 0.14, 0.17, 0.14, 0.14, 0.12, 0.15, 0.13, 0.12, 0.12, 0.14, 0.12, 0.12, 0.15, 0.12, 0.12, 0.01, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 32, c(0.46, 0.45, 0.45, 0.45, 0.46, 0.45, 0.47, 0.18, 0.18, 0.15, 0.16, 0.14, 0.15, 0.14, 0.16, 0.14, 0.12, 0.11, 0.15, 0.13, 0.13, 0.13, 0.15, 0.12, 0.13, 0.15, 0.11, 0.11, 0.00, 0.00, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 33, c(0.43, 0.42, 0.42, 0.42, 0.44, 0.43, 0.44, 0.18, 0.17, 0.14, 0.14, 0.14, 0.14, 0.14, 0.15, 0.13, 0.12, 0.10, 0.14, 0.12, 0.12, 0.12, 0.13, 0.11, 0.11, 0.13, 0.11, 0.10, 0.01, 0.00, 0.00, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 34, c(0.49, 0.47, 0.47, 0.47, 0.49, 0.48, 0.49, 0.21, 0.21, 0.18, 0.18, 0.18, 0.18, 0.16, 0.19, 0.16, 0.16, 0.13, 0.18, 0.16, 0.15, 0.15, 0.17, 0.15, 0.15, 0.18, 0.14, 0.14, 0.02, 0.00, 0.00, 0.00, 0.01))
COBA_4_fst <- fill_row(COBA_4_fst, 35, c(0.47, 0.45, 0.46, 0.46, 0.47, 0.46, 0.47, 0.22, 0.22, 0.18, 0.19, 0.18, 0.19, 0.17, 0.19, 0.17, 0.16, 0.15, 0.19, 0.16, 0.16, 0.16, 0.17, 0.15, 0.15, 0.18, 0.15, 0.15, 0.01, 0.00, 0.00, 0.01, 0.00, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 36, c(0.47, 0.46, 0.46, 0.46, 0.48, 0.47, 0.48, 0.23, 0.23, 0.19, 0.20, 0.19, 0.20, 0.17, 0.20, 0.17, 0.17, 0.14, 0.18, 0.16, 0.15, 0.15, 0.17, 0.15, 0.15, 0.18, 0.14, 0.14, 0.01, 0.01, 0.00, 0.01, 0.01, 0.00, 0.00))
COBA_4_fst <- fill_row(COBA_4_fst, 37, c(0.47, 0.46, 0.46, 0.46, 0.47, 0.46, 0.48, 0.22, 0.22, 0.19, 0.20, 0.19, 0.20, 0.17, 0.20, 0.17, 0.17, 0.15, 0.19, 0.17, 0.16, 0.16, 0.18, 0.16, 0.16, 0.19, 0.15, 0.15, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00))

COBA_4_fst[upper.tri(COBA_4_fst)] <- t(COBA_4_fst)[upper.tri(COBA_4_fst)]
diag(COBA_4_fst) <- 0
COBA_4_fst[COBA_4_fst < 0] <- 0

rownames(COBA_4_fst) <- as.character(1:37)
colnames(COBA_4_fst) <- as.character(1:37)

stopifnot(isTRUE(all.equal(COBA_4_fst, t(COBA_4_fst))))

# -----------------------------
# 5) geographic distance
# -----------------------------
coord_mat <- as.matrix(COBA_4_coords[, c("lon", "lat")])
COBA_4_dist_km <- geosphere::distm(coord_mat, fun = geosphere::distHaversine) / 1000
rownames(COBA_4_dist_km) <- as.character(COBA_4_coords$site)
colnames(COBA_4_dist_km) <- as.character(COBA_4_coords$site)

# -----------------------------
# 6) pairwise dataframe
# -----------------------------
fst_to_pairs <- function(fst_mat, dist_mat, meta_df) {
  fst_df <- as.data.frame(as.table(fst_mat))
  names(fst_df) <- c("site1", "site2", "fst")
  
  dist_df <- as.data.frame(as.table(dist_mat))
  names(dist_df) <- c("site1", "site2", "dist_km")
  
  out <- fst_df %>%
    left_join(dist_df, by = c("site1", "site2")) %>%
    mutate(site1 = as.integer(site1), site2 = as.integer(site2)) %>%
    filter(site1 < site2) %>%
    left_join(meta_df, by = c("site1" = "site")) %>%
    rename(basin1 = basin) %>%
    left_join(meta_df, by = c("site2" = "site")) %>%
    rename(basin2 = basin)
  
  out
}

COBA_4_pairs <- fst_to_pairs(COBA_4_fst, COBA_4_dist_km, COBA_4_site_meta)

# -----------------------------
# 7) map plot
# -----------------------------
world_df <- map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

x_pad <- max(0.5, diff(range(COBA_4_coords$lon)) * 0.15)
y_pad <- max(0.5, diff(range(COBA_4_coords$lat)) * 0.15)

ggplot() +
  geom_polygon(
    data = world_df,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = COBA_4_coords,
    aes(x = lon, y = lat),
    size = 2
  ) +
  geom_text(
    data = COBA_4_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = y_pad * 0.05,
    size = 2.7
  ) +
  coord_quickmap(
    xlim = c(min(COBA_4_coords$lon) - x_pad, max(COBA_4_coords$lon) + x_pad),
    ylim = c(min(COBA_4_coords$lat) - y_pad, max(COBA_4_coords$lat) + y_pad)
  ) +
  labs(
    title = "COBA-4 sampling sites",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw()

# -----------------------------
# 8) global IBD
# -----------------------------
ggplot(COBA_4_pairs, aes(x = dist_km, y = fst)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(title = "COBA-4 Global IBD",
       x = "Geographic distance (km)",
       y = "Pairwise FST")

# -----------------------------
# 9) basin-specific IBD
# -----------------------------
ggplot(
  COBA_4_pairs %>% filter(basin1 == basin2),
  aes(x = dist_km, y = fst)
) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ basin1) +
  theme_bw() +
  labs(title = "COBA-4 Basin-specific IBD",
       x = "Geographic distance (km)",
       y = "Pairwise FST")

# -----------------------------
# 10) save
# -----------------------------
save(
  COBA_4_fst,
  COBA_4_coords,
  file = "data/COBA-4.RData"
)
