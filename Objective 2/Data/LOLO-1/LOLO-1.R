############################
## LOLO-1 burbot
## Lota lota maculosa
## US sites only (Wnn removed)
##
## Sites retained:
## 1 = Id4 (Kootenai River, ID)
## 2 = Sak (Missouri River, ND)
## 3 = Mis (Leech Lake, MN)
## 4 = Nei (Oneida Lake, NY)
##
## Saved objects:
## - LOLO_1_fst
## - LOLO_1_coords
############################

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(geosphere)
  library(dplyr)
  library(tigris)
  library(sf)
})

options(tigris_use_cache = TRUE)

# -----------------------------
# STEP 0: paths
# -----------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"
dataset_code <- "LOLO-1"
dataset_dir <- file.path(base_dir, dataset_code)
data_dir <- file.path(dataset_dir, "data")

if (!dir.exists(dataset_dir)) dir.create(dataset_dir, recursive = TRUE)
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# -----------------------------
# STEP 1: coordinates
# 4 US sites (Wnn removed)
# -----------------------------
LOLO_1_coords <- data.frame(
  site_id = 1:4,
  lat = c(
    48 + 58/60,  # 1 = Id4
    47 + 37/60,  # 2 = Sak
    47 +  9/60,  # 3 = Mis
    43 + 10/60   # 4 = Nei
  ),
  lon = c(
    -(116 + 32/60),
    -(102 + 14/60),
    -( 94 + 25/60),
    -( 75 + 55/60)
  ),
  stringsAsFactors = FALSE
)

# site key
# 1 = Id4 = Kootenai River, Idaho
# 2 = Sak = Missouri River, North Dakota
# 3 = Mis = Leech Lake, Minnesota
# 4 = Nei = Oneida Lake, New York

print(LOLO_1_coords)

# -----------------------------
# STEP 2: pairwise FST matrix
# from Table 5 (Id2/Id4 treated as Id4)
# -----------------------------
site_ids <- as.character(LOLO_1_coords$site_id)

LOLO_1_fst <- matrix(
  0,
  nrow = 4,
  ncol = 4,
  dimnames = list(site_ids, site_ids)
)

LOLO_1_fst["1", "2"] <- 0.198
LOLO_1_fst["1", "3"] <- 0.246
LOLO_1_fst["1", "4"] <- 0.282
LOLO_1_fst["2", "3"] <- 0.111
LOLO_1_fst["2", "4"] <- 0.128
LOLO_1_fst["3", "4"] <- 0.192

LOLO_1_fst[lower.tri(LOLO_1_fst)] <- t(LOLO_1_fst)[lower.tri(LOLO_1_fst)]
LOLO_1_fst[LOLO_1_fst < 0] <- 0
diag(LOLO_1_fst) <- 0

print(round(LOLO_1_fst, 4))

# -----------------------------
# STEP 3: Euclidean distances
# -----------------------------
LOLO_1_geo_dist_km <- geosphere::distm(
  x = LOLO_1_coords[, c("lon", "lat")],
  fun = geosphere::distGeo
) / 1000

rownames(LOLO_1_geo_dist_km) <- site_ids
colnames(LOLO_1_geo_dist_km) <- site_ids

# -----------------------------
# STEP 4: IBD dataframe
# -----------------------------
pair_idx <- lower.tri(LOLO_1_fst)

LOLO_1_ibd_df <- data.frame(
  site1 = as.integer(rownames(LOLO_1_fst)[row(LOLO_1_fst)[pair_idx]]),
  site2 = as.integer(colnames(LOLO_1_fst)[col(LOLO_1_fst)[pair_idx]]),
  fst = LOLO_1_fst[pair_idx],
  dist_km = LOLO_1_geo_dist_km[pair_idx],
  stringsAsFactors = FALSE
)

# -----------------------------
# STEP 5: map sites
# -----------------------------
state_names <- c(
  "Idaho",
  "North Dakota",
  "Minnesota",
  "New York"
)

us_states <- tigris::states(cb = TRUE, year = 2022, class = "sf") %>%
  dplyr::filter(NAME %in% state_names)

xlim <- range(LOLO_1_coords$lon) + c(-2.5, 2.5)
ylim <- range(LOLO_1_coords$lat) + c(-1.5, 1.5)

site_map <- ggplot() +
  geom_sf(data = us_states, fill = "grey95", color = "black", linewidth = 0.3) +
  geom_point(data = LOLO_1_coords, aes(x = lon, y = lat), size = 2.8) +
  geom_text_repel(
    data = LOLO_1_coords,
    aes(x = lon, y = lat, label = site_id),
    size = 3.3
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_bw() +
  labs(title = "LOLO-1 sampling sites")

print(site_map)

# -----------------------------
# STEP 6: IBD plot
# -----------------------------
ibd_plot <- ggplot(LOLO_1_ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(
    x = "Distance (km)",
    y = expression(F[ST]),
    title = "IBD: LOLO-1"
  )

print(ibd_plot)

# -----------------------------
# STEP 7: save
# -----------------------------
save(
  LOLO_1_fst,
  LOLO_1_coords,
  file = file.path(data_dir, "LOLO-1.RData")
)