
# ============================================================
# LECO-1 | Northern leatherside chub
# Lepidomeda copei
# Blakney et al. 2014, Conservation Genetics
# Snake River watershed populations only
# FST from Supplemental Table 3a
# Coordinates from Supplemental Table 1 (UTM NAD83)
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
  library(dplyr)
  library(geosphere)
})

study_code <- "LECO-1"
study_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/LECO-1"
out_dir <- file.path(study_dir, "data")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

# -----------------------------
# 1) site metadata from Supplemental Table 1
# Snake River watershed only
# UTM = NAD83, zones 11 or 12 as reported in the supplement
# -----------------------------
site_meta_raw <- data.frame(
  site = 1:12,
  site_code = c("BVD", "GSE", "JKN", "PAC", "POL", "SQW", "STP", "TIN", "TRL", "TRA", "TRT", "UNT"),
  site_name = c(
    "Beaverdam Creek",
    "Goose Creek",
    "Jackknife Creek",
    "Pacific Creek",
    "Pole Creek",
    "Squaw Creek",
    "Stump Creek",
    "Tincup Creek",
    "Trail Creek",
    "Trapper Creek",
    "Trout Creek",
    "Unnamed tributary"
  ),
  state = c("ID", "ID", "ID", "WY", "ID/UT", "ID", "ID", "ID", "ID", "ID", "NV", "ID"),
  watershed = rep("Snake River", 12),
  subbasin = c(
    "middle Snake River", "middle Snake River", "upper Snake River", "upper Snake River",
    "middle Snake River", "upper Snake River", "upper Snake River", "upper Snake River",
    "upper Snake River", "middle Snake River", "middle Snake River", "upper Snake River"
  ),
  utm_zone = c(11, 12, 12, 12, 12, 12, 12, 12, 12, 11, 11, 12),
  easting = c(744460, 257390, 486930, 540767, 252910, 490680, 493772, 477070, 485473, 747799, 738986, 476378),
  northing = c(4655851, 4661722, 4765184, 4857677, 4653139, 4766506, 4737976, 4758657, 4765645, 4671263, 4644966, 4759522),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) convert UTM NAD83 to lat/lon
# -----------------------------
utm11 <- site_meta_raw %>%
  filter(utm_zone == 11) %>%
  st_as_sf(coords = c("easting", "northing"), crs = 26911, remove = FALSE) %>%
  st_transform(4326)

utm12 <- site_meta_raw %>%
  filter(utm_zone == 12) %>%
  st_as_sf(coords = c("easting", "northing"), crs = 26912, remove = FALSE) %>%
  st_transform(4326)

site_sf <- bind_rows(utm11, utm12) %>%
  arrange(site)

coords_mat <- st_coordinates(site_sf)

LECO_1_coords_full <- site_sf %>%
  st_drop_geometry() %>%
  mutate(
    lon = coords_mat[, "X"],
    lat = coords_mat[, "Y"]
  ) %>%
  select(site, site_code, site_name, state, watershed, subbasin, utm_zone, easting, northing, lat, lon)

LECO_1_coords <- LECO_1_coords_full %>%
  select(site, lat, lon)

stopifnot(nrow(LECO_1_coords) == 12)

# -----------------------------
# 3) pairwise FST from Supplemental Table 3a
# below diagonal in the paper
# make symmetric matrix with 0 diagonal
# -----------------------------
LECO_1_fst <- matrix(0, nrow = 12, ncol = 12)

LECO_1_fst[2, 1]  <- 0.036
LECO_1_fst[3, 1]  <- 0.126
LECO_1_fst[3, 2]  <- 0.116
LECO_1_fst[4, 1]  <- 0.127
LECO_1_fst[4, 2]  <- 0.110
LECO_1_fst[4, 3]  <- 0.107
LECO_1_fst[5, 1]  <- 0.065
LECO_1_fst[5, 2]  <- 0.042
LECO_1_fst[5, 3]  <- 0.101
LECO_1_fst[5, 4]  <- 0.112
LECO_1_fst[6, 1]  <- 0.122
LECO_1_fst[6, 2]  <- 0.121
LECO_1_fst[6, 3]  <- 0.026
LECO_1_fst[6, 4]  <- 0.116
LECO_1_fst[6, 5]  <- 0.112
LECO_1_fst[7, 1]  <- 0.106
LECO_1_fst[7, 2]  <- 0.107
LECO_1_fst[7, 3]  <- 0.020
LECO_1_fst[7, 4]  <- 0.107
LECO_1_fst[7, 5]  <- 0.082
LECO_1_fst[7, 6]  <- 0.051
LECO_1_fst[8, 1]  <- 0.136
LECO_1_fst[8, 2]  <- 0.147
LECO_1_fst[8, 3]  <- 0.055
LECO_1_fst[8, 4]  <- 0.087
LECO_1_fst[8, 5]  <- 0.137
LECO_1_fst[8, 6]  <- 0.071
LECO_1_fst[8, 7]  <- 0.039
LECO_1_fst[9, 1]  <- 0.119
LECO_1_fst[9, 2]  <- 0.112
LECO_1_fst[9, 3]  <- 0.002
LECO_1_fst[9, 4]  <- 0.094
LECO_1_fst[9, 5]  <- 0.089
LECO_1_fst[9, 6]  <- 0.034
LECO_1_fst[9, 7]  <- 0.046
LECO_1_fst[9, 8]  <- 0.071
LECO_1_fst[10, 1] <- 0.045
LECO_1_fst[10, 2] <- 0.027
LECO_1_fst[10, 3] <- 0.093
LECO_1_fst[10, 4] <- 0.080
LECO_1_fst[10, 5] <- 0.047
LECO_1_fst[10, 6] <- 0.090
LECO_1_fst[10, 7] <- 0.067
LECO_1_fst[10, 8] <- 0.095
LECO_1_fst[10, 9] <- 0.082
LECO_1_fst[11, 1] <- 0.052
LECO_1_fst[11, 2] <- 0.033
LECO_1_fst[11, 3] <- 0.085
LECO_1_fst[11, 4] <- 0.085
LECO_1_fst[11, 5] <- 0.031
LECO_1_fst[11, 6] <- 0.089
LECO_1_fst[11, 7] <- 0.070
LECO_1_fst[11, 8] <- 0.099
LECO_1_fst[11, 9] <- 0.082
LECO_1_fst[11,10] <- 0.029
LECO_1_fst[12, 1] <- 0.166
LECO_1_fst[12, 2] <- 0.162
LECO_1_fst[12, 3] <- 0.074
LECO_1_fst[12, 4] <- 0.092
LECO_1_fst[12, 5] <- 0.164
LECO_1_fst[12, 6] <- 0.096
LECO_1_fst[12, 7] <- 0.056
LECO_1_fst[12, 8] <- 0.013
LECO_1_fst[12, 9] <- 0.081
LECO_1_fst[12,10] <- 0.122
LECO_1_fst[12,11] <- 0.123

LECO_1_fst <- LECO_1_fst + t(LECO_1_fst)
diag(LECO_1_fst) <- 0
LECO_1_fst[is.na(LECO_1_fst)] <- 0
LECO_1_fst[LECO_1_fst < 0] <- 0

rownames(LECO_1_fst) <- as.character(LECO_1_coords$site)
colnames(LECO_1_fst) <- as.character(LECO_1_coords$site)

stopifnot(isTRUE(all.equal(LECO_1_fst, t(LECO_1_fst))))
stopifnot(identical(rownames(LECO_1_fst), as.character(LECO_1_coords$site)))

# -----------------------------
# 4) map plot
# show US and Canada, zoom to point extent
# -----------------------------
world_df <- ggplot2::map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

plot_df <- LECO_1_coords_full %>%
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
    nudge_y = 0.08,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(min(plot_df$lon) - x_pad, max(plot_df$lon) + x_pad),
    ylim = c(min(plot_df$lat) - y_pad, max(plot_df$lat) + y_pad)
  ) +
  labs(
    title = "LECO-1 sampling sites",
    subtitle = "Snake River watershed populations",
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
  x = as.matrix(LECO_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  fst = LECO_1_fst[upper.tri(LECO_1_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "LECO-1 IBD plot",
    x = "Geographic distance (km)",
    y = "Pairwise FST"
  ) +
  theme_bw()

print(p_ibd)

# -----------------------------
# 6) save outputs
# -----------------------------
save(
  LECO_1_fst,
  LECO_1_coords,
  file = file.path(out_dir, "LECO-1.RData")
)

# optional debug export
# write.csv(LECO_1_coords_full, file.path(out_dir, "LECO-1_site_lookup_debug.csv"), row.names = FALSE)
