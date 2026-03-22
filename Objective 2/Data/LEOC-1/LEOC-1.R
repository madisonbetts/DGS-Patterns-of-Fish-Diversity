############################
## LEOC-1 spotted gar
## Lepisosteus oculatus
## Glass et al. 2015 Conserv Genet
##
## Inputs extracted from:
## - Table 1: site names and coordinates
## - Table 4a: pairwise FST among 8 sample sites
##
## Saved objects:
## - LEOC_1_fst
## - LEOC_1_coords
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
# STEP 0: define output paths
# -----------------------------
# make dataset folder and data subfolder if they do not already exist
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"
dataset_code <- "LEOC-1"
dataset_dir <- file.path(base_dir, dataset_code)
data_dir <- file.path(dataset_dir, "data")

if (!dir.exists(dataset_dir)) dir.create(dataset_dir, recursive = TRUE)
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# -----------------------------
# STEP 1: build coordinate dataframe
# -----------------------------
# numeric site IDs match the FST matrix order exactly
LEOC_1_coords <- data.frame(
  site_id = 1:8,
  lat = c(
    29.9151,   # 1
    30.7749,   # 2
    33.1439,   # 3
    36.97966,  # 4
    41.8689,   # 5
    41.9655,   # 6
    42.2873,   # 7
    42.6145    # 8
  ),
  lon = c(
    -90.5728,  # 1
    -88.6888,  # 2
    -90.4110,  # 3
    -90.2047,  # 4
    -84.9427,  # 5
    -82.5094,  # 6
    -81.8978,  # 7
    -80.4503   # 8
  ),
  stringsAsFactors = FALSE
)

# site key
# 1 = Louisiana
# 2 = Mississippi 1
# 3 = Mississippi 2
# 4 = Missouri
# 5 = Michigan
# 6 = Point Pelee
# 7 = Rondeau Bay
# 8 = Long Point Bay

print(LEOC_1_coords)

# -----------------------------
# STEP 2: build pairwise FST matrix
# -----------------------------
# values are from Table 4a in lower-triangle order
# negative FST values are set to 0 per workflow
site_ids <- as.character(LEOC_1_coords$site_id)

LEOC_1_fst <- matrix(
  0,
  nrow = nrow(LEOC_1_coords),
  ncol = nrow(LEOC_1_coords),
  dimnames = list(site_ids, site_ids)
)

fst_vals <- c(
  0.0813,
  0.0541, 0.0774,
  0.0437, 0.0706, 0.0043,
  0.2254, 0.1326, 0.2166, 0.1923,
  0.2425, 0.1994, 0.2004, 0.1894, 0.2774,
  0.1668, 0.0899, 0.1455, 0.1238, 0.0785, 0.1202,
  0.1801, 0.0732, 0.1123, 0.0968, 0.2451, 0.2369, 0.1035
)

lt_idx <- which(lower.tri(LEOC_1_fst), arr.ind = TRUE)
lt_idx <- lt_idx[order(lt_idx[, 1], lt_idx[, 2]), , drop = FALSE]

LEOC_1_fst[lt_idx] <- fst_vals
LEOC_1_fst[upper.tri(LEOC_1_fst)] <- t(LEOC_1_fst)[upper.tri(LEOC_1_fst)]
LEOC_1_fst[LEOC_1_fst < 0] <- 0
diag(LEOC_1_fst) <- 0

print(round(LEOC_1_fst, 4))

# -----------------------------
# STEP 3: calculate Euclidean distance matrix
# -----------------------------
# this is only used to make the IBD plot
# it is not saved into the RData object
LEOC_1_geo_dist_km <- geosphere::distm(
  x = LEOC_1_coords[, c("lon", "lat")],
  fun = geosphere::distGeo
) / 1000

rownames(LEOC_1_geo_dist_km) <- site_ids
colnames(LEOC_1_geo_dist_km) <- site_ids

# -----------------------------
# STEP 4: convert pairwise data into long format for IBD plotting
# -----------------------------
pair_idx <- lower.tri(LEOC_1_fst)

LEOC_1_ibd_df <- data.frame(
  site1 = as.integer(rownames(LEOC_1_fst)[row(LEOC_1_fst)[pair_idx]]),
  site2 = as.integer(colnames(LEOC_1_fst)[col(LEOC_1_fst)[pair_idx]]),
  fst = LEOC_1_fst[pair_idx],
  dist_km = LEOC_1_geo_dist_km[pair_idx],
  stringsAsFactors = FALSE
)

print(LEOC_1_ibd_df)

# -----------------------------
# STEP 5: plot sampling sites
# -----------------------------
# use tigris state polygons only
# labels are numeric site IDs
state_names <- c(
  "Louisiana",
  "Mississippi",
  "Missouri",
  "Michigan",
  "Ohio",
  "Pennsylvania",
  "New York",
  "Indiana",
  "Illinois"
)

us_states <- tigris::states(cb = TRUE, year = 2022, class = "sf") %>%
  dplyr::filter(NAME %in% state_names)

xlim <- range(LEOC_1_coords$lon) + c(-2.5, 2.5)
ylim <- range(LEOC_1_coords$lat) + c(-1.5, 1.5)

site_map <- ggplot() +
  geom_sf(
    data = us_states,
    fill = "grey95",
    color = "black",
    linewidth = 0.3
  ) +
  geom_point(
    data = LEOC_1_coords,
    aes(x = lon, y = lat),
    size = 2.8
  ) +
  geom_text_repel(
    data = LEOC_1_coords,
    aes(x = lon, y = lat, label = site_id),
    size = 3.3,
    min.segment.length = 0,
    box.padding = 0.25,
    point.padding = 0.2
  ) +
  coord_sf(
    xlim = xlim,
    ylim = ylim,
    expand = FALSE
  ) +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "LEOC-1 sampling sites"
  )

print(site_map)

# -----------------------------
# STEP 6: make IBD plot
# -----------------------------
# plot pairwise FST against Euclidean geographic distance
ibd_plot <- ggplot(LEOC_1_ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(pairwise ~ F[ST]),
    title = "Isolation by distance for LEOC-1"
  )

print(ibd_plot)

# -----------------------------
# STEP 7: save only the requested objects
# -----------------------------
save(
  LEOC_1_fst,
  LEOC_1_coords,
  file = file.path(data_dir, "LEOC-1.RData")
)

message("Saved RData to: ", file.path(data_dir, "LEOC-1.RData"))