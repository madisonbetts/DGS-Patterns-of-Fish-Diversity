# ============================================================
# COBA-3 — Cottus bairdi (Mottled sculpin)
# Lamphere & Blum 2012, Ecology of Freshwater Fish 21:75–86
#
# Objective 2 extraction workflow
# - pairwise FST from Table 2 (top matrix; sites A–H)
# - site coordinates from Table 1
# - 8 sites spanning 5.6 km in the Nantahala River, NC
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
study_code <- "COBA-3"

base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"
study_dir <- file.path(base_dir, study_code)
out_dir   <- file.path(study_dir, "data")

if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------
# 2) site key
# Table 1 / Fig. 1 order
# -----------------------------
site_key <- data.frame(
  site = 1:8,
  site_code = LETTERS[1:8],
  site_name = paste("Site", LETTERS[1:8]),
  stringsAsFactors = FALSE
)

# -----------------------------
# 3) site coordinates
# transcribed from Table 1
# original format in paper: latitude (W), longitude (N)
# interpreted conventionally here as:
#   lat = N
#   lon = W
# DMS converted to decimal degrees
# -----------------------------
COBA_3_coords <- data.frame(
  site = 1:8,
  lat = c(
    35.0428333,  # A 35 02 34.20 N
    35.0445972,  # B 35 02 40.55 N
    35.0519500,  # C 35 03 07.22 N
    35.0569000,  # D 35 03 24.84 N
    35.0596056,  # E 35 03 34.58 N
    35.0676778,  # F 35 04 03.64 N
    35.0727472,  # G 35 04 21.89 N
    35.0752417   # H 35 04 30.87 N
  ),
  lon = c(
    -83.5084083, # A 83 30 30.27 W
    -83.5101083, # B 83 30 36.39 W
    -83.5134500, # C 83 30 48.42 W
    -83.5122472, # D 83 30 44.09 W
    -83.5165028, # E 83 30 59.41 W
    -83.5209056, # F 83 31 15.26 W
    -83.5285028, # G 83 31 42.61 W
    -83.5302778  # H 83 31 49.00 W
  )
)

# standard Objective 2 format
COBA_3_coords <- COBA_3_coords[, c("site", "lat", "lon")]

# -----------------------------
# 4) pairwise FST matrix
# from Table 2 (top matrix), sites A–H
# values reported below diagonal in the paper
# -----------------------------
COBA_3_fst <- matrix(
  c(
    0.000, 0.003, 0.008, 0.018, 0.021, 0.039, 0.071, 0.059,
    0.003, 0.000, 0.006, 0.013, 0.004, 0.025, 0.055, 0.051,
    0.008, 0.006, 0.000, 0.000, 0.000, 0.027, 0.056, 0.053,
    0.018, 0.013, 0.000, 0.000, 0.003, 0.018, 0.045, 0.043,
    0.021, 0.004, 0.000, 0.003, 0.000, 0.009, 0.039, 0.041,
    0.039, 0.025, 0.027, 0.018, 0.009, 0.000, 0.007, 0.014,
    0.071, 0.055, 0.056, 0.045, 0.039, 0.007, 0.000, 0.005,
    0.059, 0.051, 0.053, 0.043, 0.041, 0.014, 0.005, 0.000
  ),
  nrow = 8,
  byrow = TRUE
)

rownames(COBA_3_fst) <- COBA_3_coords$site
colnames(COBA_3_fst) <- COBA_3_coords$site

# workflow rules
diag(COBA_3_fst) <- 0
COBA_3_fst[COBA_3_fst < 0] <- 0

# -----------------------------
# 5) sanity checks
# -----------------------------
stopifnot(isSymmetric(COBA_3_fst))
stopifnot(identical(rownames(COBA_3_fst), as.character(COBA_3_coords$site)))
stopifnot(identical(colnames(COBA_3_fst), as.character(COBA_3_coords$site)))

print(site_key)
print(COBA_3_coords)
print(round(COBA_3_fst, 3))

# -----------------------------
# 6) map
# include US + Canada
# zoom to site extent
# plot only; do not save
# -----------------------------
us_map <- map_data("state")
can_map <- map_data("world") %>%
  filter(region == "Canada")

map_df <- COBA_3_coords %>%
  left_join(site_key, by = "site")

xpad <- 1.5
ypad <- 1.0

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
    nudge_y = 0.03,
    size = 3.2
  ) +
  coord_fixed(
    xlim = xlim_map,
    ylim = ylim_map
  ) +
  theme_bw() +
  labs(
    title = "COBA-3 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )

# -----------------------------
# 7) IBD plot
# using geodesic distance among site central localities
# plot only; do not save
# -----------------------------
coord_mat <- as.matrix(COBA_3_coords[, c("lon", "lat")])
geo_dist_km <- geosphere::distm(coord_mat, fun = geosphere::distHaversine) / 1000

ibd_df <- data.frame(
  site1 = character(),
  site2 = character(),
  distance_km = numeric(),
  fst = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:(nrow(COBA_3_fst) - 1)) {
  for (j in (i + 1):ncol(COBA_3_fst)) {
    ibd_df <- rbind(
      ibd_df,
      data.frame(
        site1 = rownames(COBA_3_fst)[i],
        site2 = colnames(COBA_3_fst)[j],
        distance_km = geo_dist_km[i, j],
        fst = COBA_3_fst[i, j],
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
    nudge_y = 0.0015,
    size = 2.8
  ) +
  theme_bw() +
  labs(
    title = "COBA-3 IBD plot",
    x = "Geodesic distance (km)",
    y = "Pairwise FST"
  )

# -----------------------------
# 8) assign and save
# -----------------------------
obj_prefix <- gsub("-", "_", study_code)
fst_name <- paste0(obj_prefix, "_fst")
coords_name <- paste0(obj_prefix, "_coords")

assign(fst_name, COBA_3_fst, envir = .GlobalEnv)
assign(coords_name, COBA_3_coords, envir = .GlobalEnv)

save(
  list = c(fst_name, coords_name),
  file = file.path(out_dir, paste0(study_code, ".RData"))
)

# -----------------------------
# 9) optional reload check
# -----------------------------
# load(file.path(out_dir, paste0(study_code, ".RData")))
# ls(pattern = "COBA_3")
