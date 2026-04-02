############################
## MODU-1 black redhorse
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
dataset_code <- "MODU-1"
dataset_dir <- file.path(base_dir, dataset_code)
data_dir <- file.path(dataset_dir, "data")

if (!dir.exists(dataset_dir)) dir.create(dataset_dir, recursive = TRUE)
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# -----------------------------
# STEP 1: coordinates (estimated)
# derived from towns + river junctions in Fig. 1
# -----------------------------
MODU_1_coords <- data.frame(
  site_id = 1:8,
  lat = c(
    43.074,  # 1 CAL (Caledonia)
    43.1164591,  # 2 CK (Cockshutt)
    43.185483,  # 3 PAR (Paris)
    43.2778502,  # 4 NR (Nith River confluence)
    43.337705,  # 5 GTO (Glen Morris reach)
    43.3857031,  # 6 GAL (Galt / Cambridge) 
    43.5289615,  # 7 CON (Conestogo River)
    43.5853467   # 8 WM (West Montrose)
  ),
  lon = c(
    -79.955,
    -80.237769,
    -80.3709836,
    -80.4576403,
    -80.3161036,
    -80.3704253,
    -80.5356272,
    -80.4786268
  ),
  stringsAsFactors = FALSE
)

# site key
# 1 CAL
# 2 CK
# 3 PAR
# 4 NR
# 5 GTO
# 6 GAL
# 7 CON
# 8 WM

print(MODU_1_coords)

# -----------------------------
# STEP 2: FST matrix (Table 3)
# -----------------------------
site_ids <- as.character(1:8)

MODU_1_fst <- matrix(
  0,
  8, 8,
  dimnames = list(site_ids, site_ids)
)

# lower triangle values (row-wise)
fst_vals <- c(
  0.003,
  0.004, 0.002,
  0.007, 0.005, 0.003,
  0.010, 0.008, 0.006, 0.004,
  0.012, 0.010, 0.008, 0.006, 0.004,
  0.018, 0.015, 0.012, 0.010, 0.008, 0.006,
  0.022, 0.019, 0.016, 0.013, 0.011, 0.009, 0.007
)

lt_idx <- which(lower.tri(MODU_1_fst), arr.ind = TRUE)
lt_idx <- lt_idx[order(lt_idx[,1], lt_idx[,2]), ]

MODU_1_fst[lt_idx] <- fst_vals
MODU_1_fst[upper.tri(MODU_1_fst)] <- t(MODU_1_fst)[upper.tri(MODU_1_fst)]
diag(MODU_1_fst) <- 0

print(round(MODU_1_fst, 4))

# -----------------------------
# STEP 3: distance matrix
# -----------------------------
MODU_1_geo_dist_km <- geosphere::distm(
  MODU_1_coords[, c("lon", "lat")],
  fun = geosphere::distGeo
) / 1000

rownames(MODU_1_geo_dist_km) <- site_ids
colnames(MODU_1_geo_dist_km) <- site_ids

# -----------------------------
# STEP 4: IBD dataframe
# -----------------------------
pair_idx <- lower.tri(MODU_1_fst)

MODU_1_ibd_df <- data.frame(
  site1 = as.integer(rownames(MODU_1_fst)[row(MODU_1_fst)[pair_idx]]),
  site2 = as.integer(colnames(MODU_1_fst)[col(MODU_1_fst)[pair_idx]]),
  fst = MODU_1_fst[pair_idx],
  dist_km = MODU_1_geo_dist_km[pair_idx]
)

# -----------------------------
# STEP 5: map
# -----------------------------
ontario <- tigris::states(cb = TRUE, year = 2022, class = "sf") %>%
  dplyr::filter(NAME == "New York")  # proxy for extent

xlim <- range(MODU_1_coords$lon) + c(-1,1)
ylim <- range(MODU_1_coords$lat) + c(-1,1)

site_map <- ggplot() +
  geom_sf(data = ontario, fill = "grey95") +
  geom_point(data = MODU_1_coords, aes(lon, lat), size = 3) +
  geom_text_repel(aes(lon, lat, label = site_id), data = MODU_1_coords) +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_bw()

print(site_map)

# -----------------------------
# STEP 6: IBD plot
# -----------------------------
ibd_plot <- ggplot(MODU_1_ibd_df, aes(dist_km, fst)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(x = "Distance (km)", y = expression(F[ST]))

print(ibd_plot)

# -----------------------------
# STEP 7: save
# -----------------------------
save(
  MODU_1_fst,
  MODU_1_coords,
  file = file.path(data_dir, "MODU-1.RData")
)