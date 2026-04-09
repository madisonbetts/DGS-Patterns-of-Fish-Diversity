# ============================================================
# ETSO-1 — Etheostoma osburni (Candy Darter)
# McBaine et al. 2023
# Objective 2 extraction workflow
#
# Sources used:
# - main paper Table 3 microsatellite FST values
# - supplement Table S1 site coordinates
#
# Notes:
# - this study used microsatellites, not SNPs
# - using the 4 stream populations as the units:
#   1 = Cripple Creek
#   2 = Dismal Creek
#   3 = Laurel Creek
#   4 = Stony Creek
# - FST values are the microsatellite values BELOW the diagonal in Table 3
# - coordinates are stream centroids derived from Table S1 site coordinates
# - no negative FST values were reported, but rule is enforced anyway
# - plots are shown in RStudio and not saved
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
study_code <- "ETSO-1"

base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"
study_dir <- file.path(base_dir, study_code)
out_dir   <- file.path(study_dir, "data")

if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------
# 2) stream/pop metadata
# -----------------------------
pop_key <- data.frame(
  site = 1:4,
  stream = c("Cripple Creek", "Dismal Creek", "Laurel Creek", "Stony Creek"),
  stringsAsFactors = FALSE
)

# -----------------------------
# 3) site coordinates
# from supplement Table S1
# coordinates are decimal degrees and only reported to 0.1 degree,
# so stream centroids are necessarily coarse but are based directly
# on the provided within-stream sampling locations
# -----------------------------

# Table S1 values entered explicitly
site_tbl <- data.frame(
  stream = c(
    rep("Cripple Creek", 3),
    rep("Dismal Creek", 3),
    rep("Stony Creek", 5),
    rep("Laurel Creek", 5)
  ),
  location_within_stream = c(
    "Lower", "Center", "Upper",
    "Below Dismal Falls", "Above Dismal Falls 1", "Above Dismal Falls 2",
    "Lower 1", "Lower 2", "Center 1", "Center 2", "Upper",
    "Lower 1", "Lower 2", "Center", "Upper 1", "Upper 2"
  ),
  lat = c(
    36.9, 36.9, 36.9,
    37.2, 37.2, 37.2,
    37.4, 37.4, 37.4, 37.4, 37.4,
    37.2, 37.3, 37.3, 37.3, 37.3
  ),
  lon = c(
    -80.9, -81.0, -81.0,
    -80.9, -80.9, -80.9,
    -80.7, -80.7, -80.6, -80.6, -80.6,
    -81.1, -81.1, -81.1, -81.1, -81.1
  ),
  stringsAsFactors = FALSE
)

# stream centroids used for Objective 2 coords df
coords_centroids <- site_tbl %>%
  group_by(stream) %>%
  summarise(
    lat = mean(lat),
    lon = mean(lon),
    .groups = "drop"
  )

ETSO_1_coords <- pop_key %>%
  left_join(coords_centroids, by = "stream") %>%
  select(site, lat, lon)

# -----------------------------
# 4) pairwise FST matrix
# microsatellite values from Table 3, below diagonal
# -----------------------------
ETSO_1_fst <- matrix(
  c(
    0.00, 0.28, 0.19, 0.17,
    0.28, 0.00, 0.33, 0.31,
    0.19, 0.33, 0.00, 0.25,
    0.17, 0.31, 0.25, 0.00
  ),
  nrow = 4,
  byrow = TRUE
)

rownames(ETSO_1_fst) <- ETSO_1_coords$site
colnames(ETSO_1_fst) <- ETSO_1_coords$site

# enforce workflow rules
diag(ETSO_1_fst) <- 0
ETSO_1_fst[ETSO_1_fst < 0] <- 0

# -----------------------------
# 5) sanity checks
# -----------------------------
stopifnot(isSymmetric(ETSO_1_fst))
stopifnot(identical(rownames(ETSO_1_fst), as.character(ETSO_1_coords$site)))
stopifnot(identical(colnames(ETSO_1_fst), as.character(ETSO_1_coords$site)))

print(pop_key)
print(ETSO_1_coords)
print(round(ETSO_1_fst, 3))

# -----------------------------
# 6) map plot
# show US + Canada and zoom to point extent
# do not save, just plot in RStudio
# -----------------------------
us_map <- map_data("state")
can_map <- map_data("world") %>%
  filter(region == "Canada")

xpad <- 2
ypad <- 2

xlim_map <- range(ETSO_1_coords$lon) + c(-xpad, xpad)
ylim_map <- range(ETSO_1_coords$lat) + c(-ypad, ypad)

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
    data = ETSO_1_coords,
    aes(x = lon, y = lat),
    size = 2.8
  ) +
  geom_text(
    data = ETSO_1_coords %>% left_join(pop_key, by = "site"),
    aes(x = lon, y = lat, label = paste0(site, ": ", stream)),
    nudge_y = 0.18,
    size = 3.2
  ) +
  coord_fixed(
    xlim = xlim_map,
    ylim = ylim_map
  ) +
  theme_bw() +
  labs(
    title = "ETSO-1 sampling streams",
    x = "Longitude",
    y = "Latitude"
  )

# -----------------------------
# 7) IBD plot
# using geodesic distance between stream centroids
# do not save, just plot in RStudio
# -----------------------------
coord_mat <- as.matrix(ETSO_1_coords[, c("lon", "lat")])
geo_dist_m <- geosphere::distm(coord_mat, fun = geosphere::distHaversine)
geo_dist_km <- geo_dist_m / 1000

ibd_df <- data.frame(
  site1 = character(),
  site2 = character(),
  distance_km = numeric(),
  fst = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:(nrow(ETSO_1_fst) - 1)) {
  for (j in (i + 1):ncol(ETSO_1_fst)) {
    ibd_df <- rbind(
      ibd_df,
      data.frame(
        site1 = rownames(ETSO_1_fst)[i],
        site2 = colnames(ETSO_1_fst)[j],
        distance_km = geo_dist_km[i, j],
        fst = ETSO_1_fst[i, j],
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
    nudge_y = 0.01,
    size = 3
  ) +
  theme_bw() +
  labs(
    title = "ETSO-1 IBD plot",
    x = "Geodesic distance (km)",
    y = "Pairwise FST"
  )

# -----------------------------
# 8) assign and save
# -----------------------------
obj_prefix <- gsub("-", "_", study_code)
fst_name <- paste0(obj_prefix, "_fst")
coords_name <- paste0(obj_prefix, "_coords")

assign(fst_name, ETSO_1_fst, envir = .GlobalEnv)
assign(coords_name, ETSO_1_coords, envir = .GlobalEnv)

save(
  list = c(fst_name, coords_name),
  file = file.path(out_dir, paste0(study_code, ".RData"))
)

# -----------------------------
# 9) optional quick reload check
# -----------------------------
# load(file.path(out_dir, paste0(study_code, ".RData")))
# ls(pattern = "ETSO_1")