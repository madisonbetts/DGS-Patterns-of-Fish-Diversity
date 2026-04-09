# -----------------------------
# ACFU-1 | Lake sturgeon
# Acipenser fulvescens
# DeHaan et al. 2006, Transactions of the American Fisheries Society
# DOI: 10.1577/T05-213.1
# -----------------------------
# Workflow notes
# - pairwise FST values are transcribed from Table 2 in the paper
# - site coordinates below are hard-coded from the user-provided
#   ACFU-1_coords.csv file
# - negative FST values are set to 0
# - plots are shown in RStudio and are not saved
# -----------------------------

library(ggplot2)
library(dplyr)
library(geosphere)
library(maps)

study_code <- "ACFU-1"
out_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ACFU-1/data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) site key
# paper order from Figure 1 / Table 2:
# 1 Manistee
# 2 Muskegon
# 3 Peshtigo
# 4 Fox
# 5 Oconto
# 6 Menominee
# 7 Wolf
# 8 Black Lake
# 9 St. Clair River / Lake St. Clair
# 10 Bad River
# 11 Sturgeon River
# -----------------------------
site_key <- data.frame(
  site = 1:11,
  site_name = c(
    "Manistee River",
    "Muskegon River",
    "Peshtigo River",
    "Fox River",
    "Oconto River",
    "Menominee River",
    "Wolf River",
    "Black Lake",
    "St. Clair River / Lake St. Clair",
    "Bad River",
    "Sturgeon River"
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) pairwise FST matrix
# transcribed from Table 2 (microsatellite FST values above diagonal)
# order is the paper order listed above
# -----------------------------
fill_row <- function(mat, row, vals) {
  stopifnot(length(vals) == row - 1)
  mat[row, seq_along(vals)] <- vals
  mat
}

fst_mat <- matrix(0, nrow = 11, ncol = 11)

fst_mat <- fill_row(fst_mat,  2, c(0.021))
fst_mat <- fill_row(fst_mat,  3, c(0.055, 0.030))
fst_mat <- fill_row(fst_mat,  4, c(0.069, 0.025, 0.024))
fst_mat <- fill_row(fst_mat,  5, c(0.066, 0.030, 0.000, 0.007))
fst_mat <- fill_row(fst_mat,  6, c(0.070, 0.051, 0.036, 0.048, 0.024))
fst_mat <- fill_row(fst_mat,  7, c(0.074, 0.034, 0.023, 0.006, 0.007, 0.044))
fst_mat <- fill_row(fst_mat,  8, c(0.047, 0.025, 0.041, 0.056, 0.046, 0.044, 0.048))
fst_mat <- fill_row(fst_mat,  9, c(0.040, 0.011, 0.040, 0.038, 0.033, 0.059, 0.037, 0.020))
fst_mat <- fill_row(fst_mat, 10, c(0.147, 0.086, 0.121, 0.127, 0.123, 0.115, 0.130, 0.080, 0.074))
fst_mat <- fill_row(fst_mat, 11, c(0.109, 0.069, 0.076, 0.084, 0.070, 0.065, 0.087, 0.052, 0.056, 0.081))

fst_mat[upper.tri(fst_mat)] <- t(fst_mat)[upper.tri(fst_mat)]
diag(fst_mat) <- 0
fst_mat[fst_mat < 0] <- 0

ACFU_1_fst <- fst_mat
rownames(ACFU_1_fst) <- as.character(site_key$site)
colnames(ACFU_1_fst) <- as.character(site_key$site)

stopifnot(isTRUE(all.equal(ACFU_1_fst, t(ACFU_1_fst))))

# -----------------------------
# 3) site coordinates
# UPDATED: hard-coded from user-provided ACFU-1_coords.csv
# -----------------------------
site_lookup <- data.frame(
  site = 1:11,
  site_name = site_key$site_name,
  lat = c(
    44.25,  # Manistee
    43.23,  # Muskegon
    45.05,  # Peshtigo
    44.25,  # Fox
    44.90,  # Oconto
    45.10,  # Menominee
    44.75,  # Wolf
    45.47,  # Black Lake
    42.82,  # St. Clair
    46.60,  # Bad
    45.70   # Sturgeon
  ),
  lon = c(
    -86.32,
    -86.25,
    -87.75,
    -88.50,
    -87.87,
    -87.60,
    -88.70,
    -84.30,
    -82.48,
    -90.65,
    -87.00
  ),
  stringsAsFactors = FALSE
)

ACFU_1_coords <- site_lookup %>%
  select(site, lat, lon)

stopifnot(identical(as.character(ACFU_1_coords$site), rownames(ACFU_1_fst)))

# -----------------------------
# 4) map plot
# include US and Canada and zoom to point extent
# -----------------------------
world_df <- map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

x_pad <- max(1.0, diff(range(site_lookup$lon)) * 0.15)
y_pad <- max(1.0, diff(range(site_lookup$lat)) * 0.15)

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
    title = "ACFU-1 sampling sites",
    subtitle = "User-updated approximate coordinates"
  ) +
  theme_bw()

print(p_map)

# -----------------------------
# 5) IBD plot
# broad-scale straight-line geographic distance among approximate reach midpoints
# -----------------------------
geo_dist <- geosphere::distm(
  x = as.matrix(site_lookup[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist[upper.tri(geo_dist)],
  fst = ACFU_1_fst[upper.tri(ACFU_1_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Geographic distance (km)",
    y = "Pairwise FST",
    title = "ACFU-1 IBD plot"
  ) +
  theme_bw()

print(p_ibd)

# -----------------------------
# 6) save objects
# -----------------------------
save(
  ACFU_1_fst,
  ACFU_1_coords,
  file = file.path(out_dir, "ACFU-1.RData")
)

