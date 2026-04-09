# ============================================================
# SAFO-12 | Brook Trout
# Salvelinus fontinalis
# Wood et al. 2018. North American Journal of Fisheries Management
#
# Updated Objective 2 extraction workflow
# - uses the user-mapped coordinates
# - averages FST values across repeated annual samples for
#   Lamothe Hollow and Beaver Creek 1
# - plots shown in RStudio only; not saved
# - saves SAFO_12_fst and SAFO_12_coords to data/SAFO-12.RData
#
# Notes:
# 1) Beaver Creek 2 before restoration was excluded from FST analysis
#    in the paper due to small sample size.
# 2) Lamothe Hollow and Beaver Creek 1 were sampled before and after
#    culvert removal. Here they are collapsed to single sites by averaging
#    the pairwise FST estimates across years, per user request.
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(geosphere)
})

study_code <- "SAFO-12"
study_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SAFO-12"
out_dir <- file.path(study_dir, "data")

if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

# -----------------------------
# 1) site metadata
# user-mapped coordinates, with repeated sites collapsed
# by averaging pre- and post-restoration coordinates
# -----------------------------
SAFO_12_coords <- data.frame(
  site_id = 1:14,
  site_name = c(
    "Buck Run Upper",
    "Buck Run Middle",
    "Black Run 1 Upper",
    "Black Run 1 Middle",
    "Upper Rocky",
    "Lower Rocky",
    "Oates Run",
    "Spruce Lower",
    "Spruce Upper",
    "Black Run 2 Lower",
    "Black Run 2 Upper",
    "Lamothe Hollow",
    "Beaver Creek 1",
    "Beaver Creek 2 after restoration"
  ),
  lat = c(
    38.5418000,
    38.5404392,
    38.5610040,
    38.5561007,
    38.4877602,
    38.4853000,
    38.4539926,
    38.4594474,
    38.4594233,
    38.4435186,
    38.4449671,
    38.4778684,
    38.5185054,
    38.5220759
  ),
  lon = c(
    -79.9455944,
    -79.9383749,
    -79.9273320,
    -79.9265441,
    -79.9811807,
    -79.9621474,
    -79.9550531,
    -79.9597438,
    -79.9693078,
    -79.9655373,
    -79.9755932,
    -79.9511240,
    -79.9443357,
    -79.9402474
  ),
  locality_label = c(
    "User-mapped coordinate",
    "User-mapped coordinate",
    "User-mapped coordinate",
    "User-mapped coordinate",
    "User-mapped coordinate",
    "User-mapped coordinate",
    "User-mapped coordinate",
    "User-mapped coordinate",
    "User-mapped coordinate",
    "User-mapped coordinate",
    "User-mapped coordinate",
    "Mean of pre- and post-restoration user-mapped coordinates",
    "Mean of pre- and post-restoration user-mapped coordinates",
    "User-mapped coordinate"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  SAFO_12_coords,
  file = file.path(study_dir, "SAFO-12_collapsed_coords.csv"),
  row.names = FALSE
)

SAFO_12_coords_out <- SAFO_12_coords %>%
  select(site_id, lat, lon)

# -----------------------------
# 2) pairwise FST matrix
# values from Table 3 in the paper, collapsed across repeated
# annual samples by averaging:
# - Lamothe Hollow pre and post
# - Beaver Creek 1 pre and post
# final order:
# 1 Buck Run Upper
# 2 Buck Run Middle
# 3 Black Run 1 Upper
# 4 Black Run 1 Middle
# 5 Upper Rocky
# 6 Lower Rocky
# 7 Oates Run
# 8 Spruce Lower
# 9 Spruce Upper
# 10 Black Run 2 Lower
# 11 Black Run 2 Upper
# 12 Lamothe Hollow (mean of pre/post)
# 13 Beaver Creek 1 (mean of pre/post)
# 14 Beaver Creek 2 after restoration
# -----------------------------
SAFO_12_fst <- matrix(0, nrow = 14, ncol = 14)
rownames(SAFO_12_fst) <- as.character(1:14)
colnames(SAFO_12_fst) <- as.character(1:14)

# row 1
SAFO_12_fst[1,2]  <- 0.008
SAFO_12_fst[1,3]  <- 0.038
SAFO_12_fst[1,4]  <- 0.019
SAFO_12_fst[1,5]  <- 0.040
SAFO_12_fst[1,6]  <- 0.021
SAFO_12_fst[1,7]  <- 0.018
SAFO_12_fst[1,8]  <- 0.040
SAFO_12_fst[1,9]  <- 0.037
SAFO_12_fst[1,10] <- 0.020
SAFO_12_fst[1,11] <- 0.014
SAFO_12_fst[1,12] <- 0.063
SAFO_12_fst[1,13] <- 0.049
SAFO_12_fst[1,14] <- 0.054

# row 2
SAFO_12_fst[2,3]  <- 0.044
SAFO_12_fst[2,4]  <- 0.013
SAFO_12_fst[2,5]  <- 0.019
SAFO_12_fst[2,6]  <- 0.016
SAFO_12_fst[2,7]  <- 0.025
SAFO_12_fst[2,8]  <- 0.022
SAFO_12_fst[2,9]  <- 0.024
SAFO_12_fst[2,10] <- 0.017
SAFO_12_fst[2,11] <- 0.030
SAFO_12_fst[2,12] <- 0.043
SAFO_12_fst[2,13] <- 0.0605
SAFO_12_fst[2,14] <- 0.078

# row 3
SAFO_12_fst[3,4]  <- 0.020
SAFO_12_fst[3,5]  <- 0.059
SAFO_12_fst[3,6]  <- 0.015
SAFO_12_fst[3,7]  <- 0.047
SAFO_12_fst[3,8]  <- 0.031
SAFO_12_fst[3,9]  <- 0.031
SAFO_12_fst[3,10] <- 0.030
SAFO_12_fst[3,11] <- 0.028
SAFO_12_fst[3,12] <- 0.0495
SAFO_12_fst[3,13] <- 0.029
SAFO_12_fst[3,14] <- 0.059

# row 4
SAFO_12_fst[4,5]  <- 0.029
SAFO_12_fst[4,6]  <- 0.000
SAFO_12_fst[4,7]  <- 0.027
SAFO_12_fst[4,8]  <- 0.013
SAFO_12_fst[4,9]  <- 0.006
SAFO_12_fst[4,10] <- 0.009
SAFO_12_fst[4,11] <- 0.022
SAFO_12_fst[4,12] <- 0.041
SAFO_12_fst[4,13] <- 0.0235
SAFO_12_fst[4,14] <- 0.041

# row 5
SAFO_12_fst[5,6]  <- 0.001
SAFO_12_fst[5,7]  <- 0.029
SAFO_12_fst[5,8]  <- 0.010
SAFO_12_fst[5,9]  <- 0.028
SAFO_12_fst[5,10] <- 0.015
SAFO_12_fst[5,11] <- 0.034
SAFO_12_fst[5,12] <- 0.069
SAFO_12_fst[5,13] <- 0.066
SAFO_12_fst[5,14] <- 0.095

# row 6
SAFO_12_fst[6,7]  <- 0.006
SAFO_12_fst[6,8]  <- 0.000
SAFO_12_fst[6,9]  <- 0.000
SAFO_12_fst[6,10] <- 0.006
SAFO_12_fst[6,11] <- 0.018
SAFO_12_fst[6,12] <- 0.0285
SAFO_12_fst[6,13] <- 0.0345
SAFO_12_fst[6,14] <- 0.059

# row 7
SAFO_12_fst[7,8]  <- 0.020
SAFO_12_fst[7,9]  <- 0.017
SAFO_12_fst[7,10] <- 0.010
SAFO_12_fst[7,11] <- 0.013
SAFO_12_fst[7,12] <- 0.068
SAFO_12_fst[7,13] <- 0.0735
SAFO_12_fst[7,14] <- 0.095

# row 8
SAFO_12_fst[8,9]  <- 0.000
SAFO_12_fst[8,10] <- 0.000
SAFO_12_fst[8,11] <- 0.020
SAFO_12_fst[8,12] <- 0.049
SAFO_12_fst[8,13] <- 0.060
SAFO_12_fst[8,14] <- 0.089

# row 9
SAFO_12_fst[9,10] <- 0.003
SAFO_12_fst[9,11] <- 0.024
SAFO_12_fst[9,12] <- 0.0575
SAFO_12_fst[9,13] <- 0.070
SAFO_12_fst[9,14] <- 0.087

# row 10
SAFO_12_fst[10,11] <- 0.000
SAFO_12_fst[10,12] <- 0.059
SAFO_12_fst[10,13] <- 0.058
SAFO_12_fst[10,14] <- 0.090

# row 11
SAFO_12_fst[11,12] <- 0.0705
SAFO_12_fst[11,13] <- 0.0575
SAFO_12_fst[11,14] <- 0.091

# row 12
SAFO_12_fst[12,13] <- 0.053375
SAFO_12_fst[12,14] <- 0.099

# row 13
SAFO_12_fst[13,14] <- 0.0285

SAFO_12_fst <- SAFO_12_fst + t(SAFO_12_fst)
diag(SAFO_12_fst) <- 0
SAFO_12_fst[SAFO_12_fst < 0] <- 0

stopifnot(isTRUE(all.equal(SAFO_12_fst, t(SAFO_12_fst))))
stopifnot(identical(rownames(SAFO_12_fst), as.character(SAFO_12_coords_out$site_id)))

# -----------------------------
# 3) map plot
# include USA and Canada and zoom to point extent
# -----------------------------
world_df <- ggplot2::map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

plot_df <- SAFO_12_coords %>%
  mutate(label = paste0(site_id, " ", site_name))

# tighter bounding box around the points
x_pad <- max(0.01, diff(range(plot_df$lon)) * 0.03)
y_pad <- max(0.01, diff(range(plot_df$lat)) * 0.03)

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
    size = 2.2
  ) +
  geom_text(
    data = plot_df,
    aes(x = lon, y = lat, label = label),
    nudge_y = 0.003,
    size = 2.8,
    check_overlap = TRUE
  ) +
  coord_quickmap(
    xlim = c(min(plot_df$lon) - x_pad, max(plot_df$lon) + x_pad),
    ylim = c(min(plot_df$lat) - y_pad, max(plot_df$lat) + y_pad),
    expand = FALSE
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  labs(
    title = "SAFO-12 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )

print(p_map)

# -----------------------------
# 4) IBD plot
# straight-line distance for QC plotting only
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(SAFO_12_coords_out[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  fst = SAFO_12_fst[upper.tri(SAFO_12_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "SAFO-12 IBD plot",
    x = "Geographic distance (km)",
    y = "Pairwise FST"
  ) +
  theme_bw()

print(p_ibd)

# -----------------------------
# 5) save outputs
# -----------------------------
SAFO_12_coords <- SAFO_12_coords[ , c("site_id", "lat", "lon")]
names(SAFO_12_coords)[names(SAFO_12_coords) == "site_id"] <- "site"

# export 
save(
  SAFO_12_fst,
  SAFO_12_coords,
  file = file.path(out_dir, "SAFO-12.RData")
)
