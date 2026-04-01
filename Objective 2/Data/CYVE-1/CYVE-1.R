
# ========================================
# CYVE-1
# Cyprinella venusta
# Blacktail Shiner
#
# Source:
# Husemann et al. 2012. Comparative biogeography reveals differences
# in population genetic structure of five species of stream fishes.
#
# This workflow uses:
#   - site coordinates from Table 1 in the main paper
#   - pairwise FST and stream distances from Appendix 3
#
# Important note on coordinates:
# The genetic analyses are reported at six catchment-level populations
# (Trinity, Bosq1, Bosq2, Little1, Little2, Little3), but Table 1 lists
# nine site coordinates because some catchments had more than one sampled
# site. For catchments represented by multiple Table 1 sites, this script
# uses the mean latitude/longitude of the listed sites as the best-available
# catchment coordinate:
#   - Bosq1 = mean(NBOS3, NBOS5)
#   - Little2 = mean(COWH, SALA)
#   - Little3 = mean(ROCK, LAMP1)
#
# Builds:
#   CYVE_1_coords
#   CYVE_1_fst
#   CYVE_1_rivdists
#
# Saves to:
#   /Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/CYVE-1/data/CYVE-1.RData
#
# Plots in RStudio:
#   - site map
#   - isolation by distance plot using stream distance
# ========================================

library(dplyr)
library(ggplot2)
library(maps)

study_code <- "CYVE-1"

base_root <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"
base_dir  <- file.path(base_root, study_code)
data_dir  <- file.path(base_dir, "data")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) catchment coordinates
# -----------------------------
# Table 1 coordinates from the paper
site_tbl <- data.frame(
  map_code = c("BEAR", "ROCK", "MBOS", "NBOS3", "NBOS5", "COWH", "SALA", "LAMP2", "LAMP1"),
  pop      = c("Trinity", "Little3", "Bosq2", "Bosq1", "Bosq1", "Little2", "Little2", "Little1", "Little3"),
  lat      = c(32.59442, 30.94494, 31.50748, 31.97692, 31.63760, 31.28327, 30.91275, 31.37802, 31.11558),
  lon      = c(-97.51018, -97.99117, -97.35624, -98.03974, -97.36640, -97.88241, -97.60105, -98.18063, -98.05432),
  stringsAsFactors = FALSE
)

# catchment-level coordinates used in the genetic matrices
coords_lookup <- site_tbl %>%
  group_by(pop) %>%
  summarise(
    lat = mean(lat),
    lon = mean(lon),
    .groups = "drop"
  )

pop_order <- c("Trinity", "Bosq1", "Bosq2", "Little1", "Little2", "Little3")

coords_lookup <- coords_lookup %>%
  slice(match(pop_order, pop))

CYVE_1_coords <- data.frame(
  site = as.character(seq_along(pop_order)),
  lat = coords_lookup$lat,
  lon = coords_lookup$lon,
  stringsAsFactors = FALSE
)

site_lookup <- data.frame(
  site = CYVE_1_coords$site,
  pop = pop_order,
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) pairwise FST matrix
# transcribed from Appendix 3
# any FST < 0 set to 0 per workflow standard
# -----------------------------
CYVE_1_fst <- matrix(
  c(
    0.0, 0.0426, 0.0206, 0.0751, 0.0847, 0.1176,
    0.0426, 0.0, 0.0021, 0.1255, 0.0824, 0.0295,
    0.0206, 0.0021, 0.0, 0.0284, 0.0308, 0.0639,
    0.0751, 0.1255, 0.0284, 0.0, 0.0039, 0.2053,
    0.0847, 0.0824, 0.0308, 0.0039, 0.0, 0.1682,
    0.1176, 0.0295, 0.0639, 0.2053, 0.1682, 0.0
  ),
  nrow = 6,
  byrow = TRUE,
  dimnames = list(CYVE_1_coords$site, CYVE_1_coords$site)
)

CYVE_1_fst[CYVE_1_fst < 0] <- 0
diag(CYVE_1_fst) <- 0

# -----------------------------
# 3) pairwise stream distance matrix (km)
# transcribed from Appendix 3
# -----------------------------
CYVE_1_rivdists <- matrix(
  c(
    0.0, 1350.4000, 1356.0, 1466.4000, 1413.9000, 1456.0,
    1350.4000, 0.0, 54.4000, 503.0, 450.5000, 492.6000,
    1356.0, 54.4000, 0.0, 508.6000, 456.1000, 498.2000,
    1466.4000, 503.0, 508.6000, 0.0, 184.6000, 33.2000,
    1413.9000, 450.5000, 456.1000, 184.6000, 0.0, 174.2000,
    1456.0, 492.6000, 498.2000, 33.2000, 174.2000, 0.0
  ),
  nrow = 6,
  byrow = TRUE,
  dimnames = list(CYVE_1_coords$site, CYVE_1_coords$site)
)

diag(CYVE_1_rivdists) <- 0

# -----------------------------
# 4) checks
# -----------------------------
stopifnot(identical(rownames(CYVE_1_fst), CYVE_1_coords$site))
stopifnot(identical(colnames(CYVE_1_fst), CYVE_1_coords$site))
stopifnot(identical(rownames(CYVE_1_rivdists), CYVE_1_coords$site))
stopifnot(identical(colnames(CYVE_1_rivdists), CYVE_1_coords$site))
stopifnot(isTRUE(all.equal(CYVE_1_fst, t(CYVE_1_fst))))
stopifnot(isTRUE(all.equal(CYVE_1_rivdists, t(CYVE_1_rivdists))))

# -----------------------------
# 5) pairwise dataframe for IBD
# -----------------------------
ibd_df <- data.frame(
  site1 = rownames(CYVE_1_fst)[row(CYVE_1_fst)[upper.tri(CYVE_1_fst)]],
  site2 = colnames(CYVE_1_fst)[col(CYVE_1_fst)[upper.tri(CYVE_1_fst)]],
  fst = CYVE_1_fst[upper.tri(CYVE_1_fst)],
  rivdist_km = CYVE_1_rivdists[upper.tri(CYVE_1_rivdists)],
  stringsAsFactors = FALSE
)

ibd_df$site1_pop <- site_lookup$pop[match(ibd_df$site1, site_lookup$site)]
ibd_df$site2_pop <- site_lookup$pop[match(ibd_df$site2, site_lookup$site)]

# -----------------------------
# 6) map
# US + Canada, zoomed to extent of points
# -----------------------------
world_map <- map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

plot_sites <- left_join(CYVE_1_coords, site_lookup, by = "site")

lon_rng <- range(plot_sites$lon, na.rm = TRUE)
lat_rng <- range(plot_sites$lat, na.rm = TRUE)

x_pad <- max(1.5, diff(lon_rng) * 0.20)
y_pad <- max(1.0, diff(lat_rng) * 0.20)

xlim_use <- c(lon_rng[1] - x_pad, lon_rng[2] + x_pad)
ylim_use <- c(lat_rng[1] - y_pad, lat_rng[2] + y_pad)

map_plot <- ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "grey96",
    color = "grey65",
    linewidth = 0.25
  ) +
  geom_point(
    data = plot_sites,
    aes(x = lon, y = lat),
    size = 3
  ) +
  geom_text(
    data = plot_sites,
    aes(x = lon, y = lat, label = paste0(site, ". ", pop)),
    nudge_y = 0.08,
    size = 3.4
  ) +
  coord_fixed(xlim = xlim_use, ylim = ylim_use) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = paste0(study_code, " sampling locations")
  )

print(map_plot)

# -----------------------------
# 7) isolation by distance plot
# -----------------------------
ibd_plot <- ggplot(ibd_df, aes(x = rivdist_km, y = fst)) +
  geom_point(size = 2.8, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Stream distance (km)",
    y = expression(F[ST]),
    title = paste0(study_code, " isolation by distance")
  )

print(ibd_plot)

# -----------------------------
# 8) save objects
# -----------------------------
CYVE_1_rivdists = CYVE_1_rivdists / 1000

save(
  CYVE_1_fst,
  CYVE_1_coords,
  CYVE_1_rivdists,
  file = file.path(data_dir, paste0(study_code, ".RData"))
)
