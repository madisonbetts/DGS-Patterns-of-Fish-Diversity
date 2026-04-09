# ============================================================
# SAFO-10 | Brook Trout
# Salvelinus fontinalis
# Torterotot et al. 2014, Transactions of the American Fisheries Society
# DOI: 10.1080/00028487.2014.952449
# ============================================================
# Workflow notes
# - site coordinates are reported directly in Table 1
# - pairwise FST values are transcribed directly from Supplementary Table S.3
# - pairwise geographic distances reported in Supplementary Table S.3 are not
#   used here because Objective 2 workflows use site coordinates and rebuild
#   the distance matrix in R
# - negative FST values are set to 0
# - plots are shown in RStudio and are not saved
# ============================================================

library(ggplot2)
library(dplyr)
library(geosphere)
library(maps)

study_code <- "SAFO-10"
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"
study_dir <- file.path(base_dir, study_code)
out_dir <- file.path(study_dir, "data")

if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

site_key <- data.frame(
  site = 1:25,
  site_code = c(
    "MS1","T1","MS2","T2","T3","T4","T5","T6","MS3","T7",
    "T8","MS4","T9","T10","T11","MS5","T12","T13","MS6","T14",
    "T15","T16","T17","MS7","T18"
  ),
  stringsAsFactors = FALSE
)

site_lookup <- data.frame(
  site = 1:25,
  site_code = site_key$site_code,
  lat = c(
    48 + 36/60 + 27.40/3600,
    48 + 38/60 + 58.50/3600,
    48 + 37/60 + 21.50/3600,
    48 + 38/60 + 58.00/3600,
    48 + 39/60 + 11.90/3600,
    48 + 39/60 + 24.20/3600,
    48 + 39/60 + 32.20/3600,
    48 + 39/60 + 29.40/3600,
    48 + 40/60 + 18.70/3600,
    48 + 40/60 + 28.20/3600,
    48 + 42/60 + 59.20/3600,
    48 + 43/60 +  3.30/3600,
    48 + 42/60 + 55.40/3600,
    48 + 43/60 + 15.80/3600,
    48 + 43/60 + 24.00/3600,
    48 + 45/60 + 40.10/3600,
    48 + 46/60 + 31.30/3600,
    48 + 46/60 + 24.10/3600,
    48 + 46/60 + 32.50/3600,
    48 + 46/60 + 27.80/3600,
    48 + 46/60 + 35.50/3600,
    48 + 47/60 + 33.40/3600,
    48 + 47/60 + 38.10/3600,
    48 + 47/60 + 36.40/3600,
    48 + 47/60 + 30.50/3600
  ),
  lon = -c(
    70 + 56/60 + 12.60/3600,
    70 + 55/60 + 21.10/3600,
    70 + 55/60 + 43.20/3600,
    70 + 55/60 + 18.90/3600,
    70 + 54/60 + 58.10/3600,
    70 + 54/60 + 28.70/3600,
    70 + 53/60 + 36.50/3600,
    70 + 53/60 + 19.80/3600,
    70 + 56/60 +  6.00/3600,
    70 + 56/60 + 17.20/3600,
    70 + 55/60 + 50.50/3600,
    70 + 55/60 + 37.40/3600,
    70 + 54/60 + 55.80/3600,
    70 + 55/60 + 20.80/3600,
    70 + 55/60 + 11.50/3600,
    70 + 53/60 + 29.40/3600,
    70 + 51/60 + 18.80/3600,
    70 + 50/60 + 56.50/3600,
    70 + 50/60 + 28.40/3600,
    70 + 50/60 + 19.10/3600,
    70 + 49/60 + 17.50/3600,
    70 + 48/60 + 30.90/3600,
    70 + 47/60 + 59.80/3600,
    70 + 48/60 +  7.60/3600,
    70 + 48/60 + 26.30/3600
  ),
  stringsAsFactors = FALSE
)

SAFO_10_coords <- site_lookup %>%
  select(site, lat, lon)

fill_row_upper <- function(mat, row, vals) {
  stopifnot(length(vals) == ncol(mat) - row)
  mat[row, (row + 1):ncol(mat)] <- vals
  mat
}

SAFO_10_fst <- matrix(0, nrow = 25, ncol = 25)

SAFO_10_fst <- fill_row_upper(SAFO_10_fst,  1, c(0.00,0.01,0.01,0.07,0.07,0.06,0.07,0.01,0.03,0.02,0.01,0.02,0.02,0.02,0.01,0.02,0.02,0.03,0.05,0.04,0.04,0.08,0.03,0.03))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst,  2, c(0.03,0.02,0.07,0.07,0.07,0.09,0.02,0.04,0.04,0.03,0.03,0.03,0.03,0.03,0.04,0.03,0.05,0.06,0.05,0.05,0.09,0.04,0.05))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst,  3, c(0.01,0.10,0.09,0.07,0.09,0.00,0.03,0.01,0.00,0.02,0.02,0.02,0.00,0.01,0.01,0.02,0.04,0.02,0.04,0.09,0.01,0.02))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst,  4, c(0.05,0.04,0.04,0.05,0.01,0.03,0.03,0.02,0.02,0.03,0.03,0.01,0.01,0.01,0.02,0.05,0.04,0.04,0.07,0.02,0.03))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst,  5, c(0.00,0.02,0.03,0.10,0.10,0.12,0.10,0.09,0.13,0.12,0.09,0.09,0.10,0.09,0.14,0.13,0.11,0.13,0.12,0.11))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst,  6, c(0.02,0.03,0.09,0.10,0.11,0.10,0.09,0.12,0.11,0.09,0.08,0.09,0.08,0.13,0.12,0.10,0.12,0.11,0.11))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst,  7, c(0.02,0.08,0.09,0.10,0.08,0.08,0.11,0.10,0.08,0.07,0.08,0.08,0.13,0.11,0.10,0.13,0.10,0.10))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst,  8, c(0.09,0.10,0.11,0.10,0.09,0.12,0.11,0.09,0.08,0.08,0.09,0.14,0.13,0.11,0.14,0.11,0.11))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst,  9, c(0.03,0.02,0.00,0.02,0.01,0.01,0.00,0.01,0.01,0.03,0.04,0.03,0.04,0.09,0.01,0.02))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst, 10, c(0.04,0.03,0.05,0.05,0.05,0.03,0.04,0.04,0.05,0.06,0.06,0.07,0.10,0.04,0.05))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst, 11, c(0.01,0.04,0.04,0.04,0.02,0.03,0.03,0.03,0.06,0.04,0.05,0.10,0.02,0.03))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst, 12, c(0.01,0.02,0.03,0.01,0.02,0.02,0.03,0.05,0.04,0.05,0.10,0.02,0.03))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst, 13, c(0.03,0.03,0.02,0.03,0.02,0.05,0.07,0.05,0.06,0.10,0.04,0.05))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst, 14, c(0.00,0.02,0.04,0.02,0.07,0.06,0.04,0.08,0.13,0.04,0.06))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst, 15, c(0.02,0.03,0.03,0.08,0.05,0.05,0.08,0.13,0.04,0.06))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst, 16, c(0.01,0.01,0.02,0.03,0.02,0.04,0.09,0.00,0.02))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst, 17, c(0.01,0.02,0.04,0.02,0.02,0.06,0.01,0.02))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst, 18, c(0.02,0.05,0.04,0.05,0.09,0.02,0.03))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst, 19, c(0.05,0.04,0.02,0.06,0.02,0.00))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst, 20, c(0.03,0.05,0.10,0.02,0.04))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst, 21, c(0.05,0.10,0.01,0.03))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst, 22, c(0.01,0.03,0.01))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst, 23, c(0.08,0.06))
SAFO_10_fst <- fill_row_upper(SAFO_10_fst, 24, c(0.01))

SAFO_10_fst[lower.tri(SAFO_10_fst)] <- t(SAFO_10_fst)[lower.tri(SAFO_10_fst)]
diag(SAFO_10_fst) <- 0
SAFO_10_fst[SAFO_10_fst < 0] <- 0

rownames(SAFO_10_fst) <- as.character(site_key$site)
colnames(SAFO_10_fst) <- as.character(site_key$site)

stopifnot(isTRUE(all.equal(SAFO_10_fst, t(SAFO_10_fst))))

world_df <- map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

x_pad <- max(0.15, diff(range(site_lookup$lon)) * 0.10)
y_pad <- max(0.15, diff(range(site_lookup$lat)) * 0.10)

p_map <- ggplot() +
  geom_polygon(
    data = world_df,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = site_lookup,
    aes(x = lon, y = lat),
    size = 2
  ) +
  geom_text(
    data = site_lookup,
    aes(x = lon, y = lat, label = site),
    nudge_y = y_pad * 0.03,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(min(site_lookup$lon) - x_pad, max(site_lookup$lon) + x_pad),
    ylim = c(min(site_lookup$lat) - y_pad, max(site_lookup$lat) + y_pad)
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "SAFO-10 sampling sites",
    subtitle = "Coordinates reported in Table 1"
  ) +
  theme_bw()

print(p_map)

geo_dist <- geosphere::distm(
  x = as.matrix(site_lookup[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist[upper.tri(geo_dist)],
  fst = SAFO_10_fst[upper.tri(SAFO_10_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Geographic distance (km)",
    y = "Pairwise FST",
    title = "SAFO-10 IBD plot"
  ) +
  theme_bw()

print(p_ibd)

save(
  SAFO_10_fst,
  SAFO_10_coords,
  file = file.path(out_dir, "SAFO-10.RData")
)
