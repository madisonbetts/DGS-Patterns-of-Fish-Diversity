# ============================================================
# CHSA-1 — Clinch Dace
# Chrosomus sp. cf. saylori
# Bourquin et al. 2023, Fishes 8:365
#
# Objective 2 extraction workflow
# - marker type: microsatellites
# - FST matrix transcribed from Table 5 in the paper
# - site coordinates estimated from named stream centroids
# - output objects:
#     CHSA_1_fst
#     CHSA_1_coords
# - save to:
#     .../Objective 2/Data/CHSA-1/data/CHSA-1.RData
# ============================================================

# -----------------------------
# 0) libraries
# -----------------------------
library(dplyr)
library(ggplot2)
library(maps)
library(geosphere)

# -----------------------------
# 1) setup
# -----------------------------
study_code <- "CHSA-1"

base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"
study_dir <- file.path(base_dir, study_code)
out_dir   <- file.path(study_dir, "data")

if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------
# 2) stream key
# using stream names from Table 5 / results
# note: paper text earlier says "Big Creek" and "Hurricane Creek",
# but Table 5 and the STRUCTURE figure use Big Lick Creek and
# Hurricane Fork; those are used here for consistency.
# -----------------------------
stream_key <- data.frame(
  site = 1:7,
  stream = c(
    "Big Lick Creek",
    "Greasy Creek",
    "Hart Creek",
    "Hurricane Fork",
    "Lewis Creek",
    "Middle Creek",
    "Pine Creek"
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 3) estimated site coordinates
# updated from uploaded CSV of stream centroids
# -----------------------------
CHSA_1_coords <- data.frame(
  site = 1:7,
  lat = c(
    37.069,
    37.1753088,
    37.0188432,
    36.9930858,
    37.0591661,
    37.138276,
    37.081836
  ),
  lon = c(
    -81.92485,
    -81.7313567,
    -82.1000614,
    -82.162691,
    -81.988355,
    -81.7497582,
    -81.9117191
  )
)

# keep only standard columns
CHSA_1_coords <- CHSA_1_coords[, c("site", "lat", "lon")]

# -----------------------------
# 4) pairwise FST matrix
# transcribed from Table 5
# order:
# 1 Big Lick Creek
# 2 Greasy Creek
# 3 Hart Creek
# 4 Hurricane Fork
# 5 Lewis Creek
# 6 Middle Creek
# 7 Pine Creek
# -----------------------------
CHSA_1_fst <- matrix(
  c(
    0.000, 0.186, 0.305, 0.105, 0.243, 0.157, 0.191,
    0.186, 0.000, 0.310, 0.299, 0.257, 0.207, 0.175,
    0.305, 0.310, 0.000, 0.272, 0.416, 0.313, 0.337,
    0.105, 0.299, 0.272, 0.000, 0.272, 0.335, 0.371,
    0.243, 0.257, 0.416, 0.272, 0.000, 0.526, 0.462,
    0.157, 0.207, 0.313, 0.335, 0.526, 0.000, 0.036,
    0.191, 0.175, 0.337, 0.371, 0.462, 0.036, 0.000
  ),
  nrow = 7,
  byrow = TRUE
)

rownames(CHSA_1_fst) <- as.character(CHSA_1_coords$site)
colnames(CHSA_1_fst) <- as.character(CHSA_1_coords$site)

# workflow rules
diag(CHSA_1_fst) <- 0
CHSA_1_fst[CHSA_1_fst < 0] <- 0

# -----------------------------
# 5) sanity checks
# -----------------------------
stopifnot(isSymmetric(CHSA_1_fst))
stopifnot(identical(rownames(CHSA_1_fst), as.character(CHSA_1_coords$site)))
stopifnot(identical(colnames(CHSA_1_fst), as.character(CHSA_1_coords$site)))

print(stream_key)
print(CHSA_1_coords)
print(round(CHSA_1_fst, 3))

# -----------------------------
# 6) map
# include US and Canada
# zoom to point extent
# plot only; do not save
# -----------------------------
us_map <- map_data("state")
can_map <- map_data("world") %>%
  filter(region == "Canada")

map_df <- CHSA_1_coords %>%
  left_join(stream_key, by = "site")

xpad <- 0.5
ypad <- 0.5

xlim_map <- range(map_df$lon) + c(-xpad, xpad)
ylim_map <- range(map_df$lat) + c(-ypad, ypad)

ggplot() +
  geom_polygon(
    data = can_map,
    aes(x = long, y = lat, group = group),
    fill = "grey95", color = "grey60", linewidth = 0.3
  ) +
  geom_polygon(
    data = us_map,
    aes(x = long, y = lat, group = group),
    fill = "grey90", color = "grey50", linewidth = 0.3
  ) +
  geom_point(
    data = map_df,
    aes(x = lon, y = lat),
    size = 2.8
  ) +
  geom_text(
    data = map_df,
    aes(x = lon, y = lat, label = paste0(site, ": ", stream)),
    nudge_y = 0.10,
    size = 3.0
  ) +
  coord_fixed(
    xlim = xlim_map,
    ylim = ylim_map
  ) +
  theme_bw() +
  labs(
    title = "CHSA-1 sampling streams",
    x = "Longitude",
    y = "Latitude"
  )

# -----------------------------
# 7) IBD plot
# geodesic distance between stream centroids
# plot only; do not save
# -----------------------------
coord_mat <- as.matrix(CHSA_1_coords[, c("lon", "lat")])
geo_dist_km <- geosphere::distm(coord_mat, fun = geosphere::distHaversine) / 1000

ibd_df <- data.frame(
  site1 = character(),
  site2 = character(),
  distance_km = numeric(),
  fst = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:(nrow(CHSA_1_fst) - 1)) {
  for (j in (i + 1):ncol(CHSA_1_fst)) {
    ibd_df <- rbind(
      ibd_df,
      data.frame(
        site1 = rownames(CHSA_1_fst)[i],
        site2 = colnames(CHSA_1_fst)[j],
        distance_km = geo_dist_km[i, j],
        fst = CHSA_1_fst[i, j],
        stringsAsFactors = FALSE
      )
    )
  }
}

print(ibd_df)

ggplot(ibd_df, aes(x = distance_km, y = fst)) +
  geom_point(size = 2.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  geom_text(
    aes(label = paste0(site1, "-", site2)),
    nudge_y = 0.012,
    size = 3
  ) +
  theme_bw() +
  labs(
    title = "CHSA-1 IBD plot",
    x = "Geodesic distance (km)",
    y = "Pairwise FST"
  )

# -----------------------------
# 8) assign and save
# -----------------------------
obj_prefix <- gsub("-", "_", study_code)
fst_name <- paste0(obj_prefix, "_fst")
coords_name <- paste0(obj_prefix, "_coords")

assign(fst_name, CHSA_1_fst, envir = .GlobalEnv)
assign(coords_name, CHSA_1_coords, envir = .GlobalEnv)

save(
  list = c(fst_name, coords_name),
  file = file.path(out_dir, paste0(study_code, ".RData"))
)

# -----------------------------
# 9) optional reload check
# -----------------------------
# load(file.path(out_dir, paste0(study_code, ".RData")))
# ls(pattern = "CHSA_1")