
# ============================================================
# LECO-2 | Northern leatherside chub
# Lepidomeda copei
# Blakney et al. 2014, Conservation Genetics
# Bear River watershed populations only
# FST from Supplemental Table 3b
# Coordinates from Supplemental Table 1 (UTM NAD83)
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
  library(dplyr)
  library(geosphere)
})

study_code <- "LECO-2"
study_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/LECO-2"
out_dir <- file.path(study_dir, "data")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

# -----------------------------
# 1) site metadata from Supplemental Table 1
# Bear River watershed only
# UTM = NAD83, zone 12 as reported in the supplement
# -----------------------------
site_meta_raw <- data.frame(
  site = 1:9,
  site_code = c("DMN", "DRY", "LCH", "MIL", "MLC", "MUD", "SLR", "TWN", "YLW"),
  site_name = c(
    "Deadman Creek",
    "Dry Fork Smiths Fork",
    "LaChapelle Creek",
    "Mill Creek (MIL)",
    "Mill Creek (MLC)",
    "Muddy Creek",
    "Sulphur Creek",
    "Twin Creek",
    "Yellow Creek"
  ),
  state = c("UT", "WY", "WY", "WY", "WY", "WY", "WY", "WY", "UT"),
  watershed = rep("Bear River", 9),
  subbasin = c(
    "upper Bear River", "middle Bear River", "upper Bear River", "middle Bear River",
    "upper Bear River", "middle Bear River", "upper Bear River", "middle Bear River",
    "upper Bear River"
  ),
  utm_zone = rep(12, 9),
  easting = c(517232, 510500, 518443, 505203, 508823, 508535, 515185, 518484, 502980),
  northing = c(4531725, 4690882, 4548485, 4672296, 4545855, 4670572, 4545586, 4629805, 4538488),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) convert UTM NAD83 to lat/lon
# -----------------------------
site_sf <- site_meta_raw %>%
  st_as_sf(coords = c("easting", "northing"), crs = 26912, remove = FALSE) %>%
  st_transform(4326) %>%
  arrange(site)

coords_mat <- st_coordinates(site_sf)

LECO_2_coords_full <- site_sf %>%
  st_drop_geometry() %>%
  mutate(
    lon = coords_mat[, "X"],
    lat = coords_mat[, "Y"]
  ) %>%
  select(site, site_code, site_name, state, watershed, subbasin, utm_zone, easting, northing, lat, lon)

LECO_2_coords <- LECO_2_coords_full %>%
  select(site, lat, lon)

stopifnot(nrow(LECO_2_coords) == 9)

# -----------------------------
# 3) pairwise FST from Supplemental Table 3b
# below diagonal in the paper
# make symmetric matrix with 0 diagonal
# -----------------------------
LECO_2_fst <- matrix(0, nrow = 9, ncol = 9)

LECO_2_fst[2, 1] <- 0.125
LECO_2_fst[3, 1] <- 0.101
LECO_2_fst[3, 2] <- 0.120
LECO_2_fst[4, 1] <- 0.073
LECO_2_fst[4, 2] <- 0.070
LECO_2_fst[4, 3] <- 0.058
LECO_2_fst[5, 1] <- 0.058
LECO_2_fst[5, 2] <- 0.100
LECO_2_fst[5, 3] <- 0.021
LECO_2_fst[5, 4] <- 0.047
LECO_2_fst[6, 1] <- 0.080
LECO_2_fst[6, 2] <- 0.071
LECO_2_fst[6, 3] <- 0.044
LECO_2_fst[6, 4] <- 0.012
LECO_2_fst[6, 5] <- 0.025
LECO_2_fst[7, 1] <- 0.072
LECO_2_fst[7, 2] <- 0.095
LECO_2_fst[7, 3] <- 0.023
LECO_2_fst[7, 4] <- 0.032
LECO_2_fst[7, 5] <- 0.007
LECO_2_fst[7, 6] <- 0.021
LECO_2_fst[8, 1] <- 0.114
LECO_2_fst[8, 2] <- 0.118
LECO_2_fst[8, 3] <- 0.053
LECO_2_fst[8, 4] <- 0.082
LECO_2_fst[8, 5] <- 0.041
LECO_2_fst[8, 6] <- 0.052
LECO_2_fst[8, 7] <- 0.031
LECO_2_fst[9, 1] <- 0.071
LECO_2_fst[9, 2] <- 0.109
LECO_2_fst[9, 3] <- 0.023
LECO_2_fst[9, 4] <- 0.043
LECO_2_fst[9, 5] <- 0.008
LECO_2_fst[9, 6] <- 0.029
LECO_2_fst[9, 7] <- 0.009
LECO_2_fst[9, 8] <- 0.033

LECO_2_fst <- LECO_2_fst + t(LECO_2_fst)
diag(LECO_2_fst) <- 0
LECO_2_fst[is.na(LECO_2_fst)] <- 0
LECO_2_fst[LECO_2_fst < 0] <- 0

rownames(LECO_2_fst) <- as.character(LECO_2_coords$site)
colnames(LECO_2_fst) <- as.character(LECO_2_coords$site)

stopifnot(isTRUE(all.equal(LECO_2_fst, t(LECO_2_fst))))
stopifnot(identical(rownames(LECO_2_fst), as.character(LECO_2_coords$site)))

# -----------------------------
# 4) map plot
# show US and Canada, zoom to point extent
# -----------------------------
world_df <- ggplot2::map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

plot_df <- LECO_2_coords_full %>%
  mutate(label = paste0(site, " ", site_code))

x_pad <- max(0.5, diff(range(plot_df$lon)) * 0.10)
y_pad <- max(0.5, diff(range(plot_df$lat)) * 0.10)

p_map <- ggplot() +
  geom_polygon(
    data = world_df,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = plot_df,
    aes(x = lon, y = lat),
    size = 2
  ) +
  geom_text(
    data = plot_df,
    aes(x = lon, y = lat, label = label),
    nudge_y = 0.05,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(min(plot_df$lon) - x_pad, max(plot_df$lon) + x_pad),
    ylim = c(min(plot_df$lat) - y_pad, max(plot_df$lat) + y_pad)
  ) +
  labs(
    title = "LECO-2 sampling sites",
    subtitle = "Bear River watershed populations",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw()

print(p_map)

# -----------------------------
# 5) IBD plot
# use straight-line distance for QC plotting only
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(LECO_2_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  fst = LECO_2_fst[upper.tri(LECO_2_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "LECO-2 IBD plot",
    x = "Geographic distance (km)",
    y = "Pairwise FST"
  ) +
  theme_bw()

print(p_ibd)

# -----------------------------
# 6) save outputs
# -----------------------------
save(
  LECO_2_fst,
  LECO_2_coords,
  file = file.path(out_dir, "LECO-2.RData")
)

# optional debug export
# write.csv(LECO_2_coords_full, file.path(out_dir, "LECO-2_site_lookup_debug.csv"), row.names = FALSE)
