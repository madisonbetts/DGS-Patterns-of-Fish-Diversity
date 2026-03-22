############################
## PEFL-3 yellow perch
## Perca flavescens
############################

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(geosphere)
  library(dplyr)
  library(sf)
  library(maps)
})

# -----------------------------
# STEP 0: paths
# -----------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"
dataset_code <- "PEFL-3"
dataset_dir <- file.path(base_dir, dataset_code)
data_dir <- file.path(dataset_dir, "data")

if (!dir.exists(dataset_dir)) dir.create(dataset_dir, recursive = TRUE)
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# -----------------------------
# STEP 1: coordinates
# approximate locations read from Fig. 1 and site names in text/Table 2
# -----------------------------
PEFL_3_coords <- data.frame(
  site_id = 1:8,
  lat = c(
    42.43,  # 1 = A = Lake St. Clair
    41.85,  # 2 = B = Monroe, Michigan
    42.26,  # 3 = C = Erieau, Ontario
    41.76,  # 4 = D = Fairport, Ohio
    41.81,  # 5 = E = Perry, Ohio
    42.16,  # 6 = F = Erie, Pennsylvania
    42.48,  # 7 = G = Dunkirk, New York
    43.25   # 8 = H = Lake Ontario out-group
  ),
  lon = c(
    -82.80,
    -83.38,
    -81.90,
    -81.27,
    -81.14,
    -80.08,
    -79.33,
    -76.55
  ),
  stringsAsFactors = FALSE
)

# site key
# 1 = A = Lake St. Clair
# 2 = B = Monroe, Michigan
# 3 = C = Erieau, Ontario
# 4 = D = Fairport, Ohio
# 5 = E = Perry, Ohio
# 6 = F = Erie, Pennsylvania
# 7 = G = Dunkirk, New York
# 8 = H = Lake Ontario

print(PEFL_3_coords)

# -----------------------------
# STEP 2: pairwise FST matrix
# 2009 values from Table 5 (above diagonal in paper)
# -----------------------------
site_ids <- as.character(PEFL_3_coords$site_id)

PEFL_3_fst <- matrix(
  0,
  nrow = 8,
  ncol = 8,
  dimnames = list(site_ids, site_ids)
)

# row A
PEFL_3_fst["1", "2"] <- 0.030
PEFL_3_fst["1", "3"] <- 0.044
PEFL_3_fst["1", "4"] <- 0.065
PEFL_3_fst["1", "5"] <- 0.052
PEFL_3_fst["1", "6"] <- 0.044
PEFL_3_fst["1", "7"] <- 0.040
PEFL_3_fst["1", "8"] <- 0.033

# row B
PEFL_3_fst["2", "3"] <- 0.042
PEFL_3_fst["2", "4"] <- 0.069
PEFL_3_fst["2", "5"] <- 0.048
PEFL_3_fst["2", "6"] <- 0.055
PEFL_3_fst["2", "7"] <- 0.061
PEFL_3_fst["2", "8"] <- 0.051

# row C
PEFL_3_fst["3", "4"] <- 0.056
PEFL_3_fst["3", "5"] <- 0.047
PEFL_3_fst["3", "6"] <- 0.037
PEFL_3_fst["3", "7"] <- 0.033
PEFL_3_fst["3", "8"] <- 0.045

# row D
PEFL_3_fst["4", "5"] <- 0.002
PEFL_3_fst["4", "6"] <- 0.030
PEFL_3_fst["4", "7"] <- 0.038
PEFL_3_fst["4", "8"] <- 0.056

# row E
PEFL_3_fst["5", "6"] <- 0.028
PEFL_3_fst["5", "7"] <- 0.037
PEFL_3_fst["5", "8"] <- 0.049

# row F
PEFL_3_fst["6", "7"] <- 0.014
PEFL_3_fst["6", "8"] <- 0.034

# row G
PEFL_3_fst["7", "8"] <- 0.034

PEFL_3_fst[lower.tri(PEFL_3_fst)] <- t(PEFL_3_fst)[lower.tri(PEFL_3_fst)]
PEFL_3_fst[PEFL_3_fst < 0] <- 0
diag(PEFL_3_fst) <- 0

print(round(PEFL_3_fst, 4))

# -----------------------------
# STEP 3: Euclidean distance matrix
# for IBD plot only
# -----------------------------
PEFL_3_geo_dist_km <- geosphere::distm(
  x = PEFL_3_coords[, c("lon", "lat")],
  fun = geosphere::distGeo
) / 1000

rownames(PEFL_3_geo_dist_km) <- site_ids
colnames(PEFL_3_geo_dist_km) <- site_ids

# -----------------------------
# STEP 4: long dataframe for IBD plot
# -----------------------------
pair_idx <- lower.tri(PEFL_3_fst)

PEFL_3_ibd_df <- data.frame(
  site1 = as.integer(rownames(PEFL_3_fst)[row(PEFL_3_fst)[pair_idx]]),
  site2 = as.integer(colnames(PEFL_3_fst)[col(PEFL_3_fst)[pair_idx]]),
  fst = PEFL_3_fst[pair_idx],
  dist_km = PEFL_3_geo_dist_km[pair_idx],
  stringsAsFactors = FALSE
)

print(PEFL_3_ibd_df)

# -----------------------------
# STEP 5: basemap and site plot
# Canada + USA polygons around lower Great Lakes
# -----------------------------
base_map <- sf::st_as_sf(
  maps::map("world", regions = c("Canada", "USA"), fill = TRUE, plot = FALSE)
)

xlim <- range(PEFL_3_coords$lon) + c(-0.8, 0.8)
ylim <- range(PEFL_3_coords$lat) + c(-0.5, 0.5)

site_map <- ggplot() +
  geom_sf(
    data = base_map,
    fill = "grey95",
    color = "black",
    linewidth = 0.25
  ) +
  geom_point(
    data = PEFL_3_coords,
    aes(x = lon, y = lat),
    size = 2.8
  ) +
  geom_text_repel(
    data = PEFL_3_coords,
    aes(x = lon, y = lat, label = site_id),
    size = 3.3,
    min.segment.length = 0,
    box.padding = 0.25,
    point.padding = 0.2
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "PEFL-2 sampling sites"
  )

print(site_map)

# -----------------------------
# STEP 6: IBD plot
# -----------------------------
ibd_plot <- ggplot(PEFL_3_ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(pairwise ~ F[ST]),
    title = "Isolation by distance for PEFL-2"
  )

print(ibd_plot)

# -----------------------------
# STEP 7: save only requested objects
# -----------------------------
save(
  PEFL_3_fst,
  PEFL_3_coords,
  file = file.path(data_dir, "PEFL-3.RData")
)

#message("Saved RData to: ", file.path(data_dir, "PEFL-2.RData"))
