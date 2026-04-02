# ============================================================
# CYEL-1 — Cycleptus elongatus (Blue sucker)
# Bessert & Orti 2008, Conservation Genetics 9:821–832
#
# Objective 2 extraction workflow
# - pairwise FST from Table 3 (Missouri River sites only)
# - site coordinates from Table 1
# - Missouri site B excluded because Table 3 explicitly excludes it
#   due to small sample size
# - sites retained: A, C, D, E, F, G, H, I
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
study_code <- "CYEL-1"

base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"
study_dir <- file.path(base_dir, study_code)
out_dir   <- file.path(study_dir, "data")

if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------
# 2) site key
# Missouri River sites in Table 3
# Table 3 excludes site B due to small sample size
# -----------------------------
site_key <- data.frame(
  site = 1:8,
  site_code = c("A", "C", "D", "E", "F", "G", "H", "I"),
  site_name = c(
    "Missouri River, Missouri",
    "Missouri River, Nebraska",
    "Platte River, Nebraska",
    "Missouri River, South Dakota",
    "Missouri River, North Dakota",
    "Yellowstone River, Montana",
    "Missouri River, Montana (H)",
    "Missouri River, Montana (I)"
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 3) site coordinates
# transcribed from Table 1
# central locality (CL) converted from DMS to decimal degrees
# order follows Table 3 Missouri FST matrix: A, C, D, E, F, G, H, I
# -----------------------------
CYEL_1_coords <- data.frame(
  site = 1:8,
  lat = c(
    38.6800000,   # A  Missouri/Missouri
    40.8136111,   # C  Missouri/Nebraska
    41.0869444,   # D  Platte/Nebraska
    42.7805556,   # E  Missouri/South Dakota
    47.4283333,   # F  Missouri/North Dakota
    46.7266667,   # G  Yellowstone/Montana
    48.0616667,   # H  Missouri/Montana
    47.7308333    # I  Missouri/Montana
  ),
  lon = c(
    -90.6688889,  # A
    -95.8452778,  # C
    -96.1680556,  # D
    -97.9808333,  # E
    -101.4094444, # F
    -96.5175000,  # G
    -106.9283333, # H
    -109.9908333  # I
  )
)

# standard Objective 2 format
CYEL_1_coords <- CYEL_1_coords[, c("site", "lat", "lon")]

# -----------------------------
# 4) pairwise FST matrix
# from Table 3, Missouri River
# below diagonal in the published table
# sites in order: A, C, D, E, F, G, H, I
# -----------------------------
CYEL_1_fst <- matrix(
  c(
    0.0000, 0.0007, 0.0039, 0.0005, 0.0113, 0.0150, 0.0109, 0.0371,
    0.0007, 0.0000, 0.0030, 0.0009, 0.0108, 0.0093, 0.0105, 0.0291,
    0.0039, 0.0030, 0.0000, 0.0029, 0.0068, 0.0095, 0.0074, 0.0363,
    0.0005, 0.0009, 0.0029, 0.0000, 0.0035, 0.0013, 0.0006, 0.0108,
    0.0113, 0.0108, 0.0068, 0.0035, 0.0000, 0.0056, 0.0001, 0.0094,
    0.0150, 0.0093, 0.0095, 0.0013, 0.0056, 0.0000, 0.0011, 0.0102,
    0.0109, 0.0105, 0.0074, 0.0006, 0.0001, 0.0011, 0.0000, 0.0214,
    0.0371, 0.0291, 0.0363, 0.0108, 0.0094, 0.0102, 0.0214, 0.0000
  ),
  nrow = 8,
  byrow = TRUE
)

rownames(CYEL_1_fst) <- CYEL_1_coords$site
colnames(CYEL_1_fst) <- CYEL_1_coords$site

# workflow rules
diag(CYEL_1_fst) <- 0
CYEL_1_fst[CYEL_1_fst < 0] <- 0

# -----------------------------
# 5) sanity checks
# -----------------------------
stopifnot(isSymmetric(CYEL_1_fst))
stopifnot(identical(rownames(CYEL_1_fst), as.character(CYEL_1_coords$site)))
stopifnot(identical(colnames(CYEL_1_fst), as.character(CYEL_1_coords$site)))

print(site_key)
print(CYEL_1_coords)
print(round(CYEL_1_fst, 4))

# -----------------------------
# 6) map
# include US + Canada
# zoom to site extent
# plot only; do not save
# -----------------------------
us_map <- map_data("state")
can_map <- map_data("world") %>%
  filter(region == "Canada")

map_df <- CYEL_1_coords %>%
  left_join(site_key, by = "site")

xpad <- 5
ypad <- 4

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
    aes(x = lon, y = lat, label = paste0(site, ": ", site_code)),
    nudge_y = 0.35,
    size = 3.2
  ) +
  coord_fixed(
    xlim = xlim_map,
    ylim = ylim_map
  ) +
  theme_bw() +
  labs(
    title = "CYEL-1 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )

# -----------------------------
# 7) IBD plot
# using geodesic distance among site central localities
# plot only; do not save
# -----------------------------
coord_mat <- as.matrix(CYEL_1_coords[, c("lon", "lat")])
geo_dist_km <- geosphere::distm(coord_mat, fun = geosphere::distHaversine) / 1000

ibd_df <- data.frame(
  site1 = character(),
  site2 = character(),
  distance_km = numeric(),
  fst = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:(nrow(CYEL_1_fst) - 1)) {
  for (j in (i + 1):ncol(CYEL_1_fst)) {
    ibd_df <- rbind(
      ibd_df,
      data.frame(
        site1 = rownames(CYEL_1_fst)[i],
        site2 = colnames(CYEL_1_fst)[j],
        distance_km = geo_dist_km[i, j],
        fst = CYEL_1_fst[i, j],
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
    nudge_y = 0.0012,
    size = 2.8
  ) +
  theme_bw() +
  labs(
    title = "CYEL-1 IBD plot",
    x = "Geodesic distance (km)",
    y = "Pairwise FST"
  )

# -----------------------------
# 8) assign and save
# -----------------------------
obj_prefix <- gsub("-", "_", study_code)
fst_name <- paste0(obj_prefix, "_fst")
coords_name <- paste0(obj_prefix, "_coords")

assign(fst_name, CYEL_1_fst, envir = .GlobalEnv)
assign(coords_name, CYEL_1_coords, envir = .GlobalEnv)

save(
  list = c(fst_name, coords_name),
  file = file.path(out_dir, paste0(study_code, ".RData"))
)

# -----------------------------
# 9) optional reload check
# -----------------------------
# load(file.path(out_dir, paste0(study_code, ".RData")))
# ls(pattern = "CYEL_1")
