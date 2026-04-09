# ============================================================
# SACO-6 | Bull Trout
# Salvelinus confluentus
# Nyce et al. 2013. North American Journal of Fisheries Management
#
# Updated Objective 2 extraction workflow
# - pairwise FST matrix hard-coded from Table 3
# - coordinates updated from the cleaned user-provided CSV
# - plots shown in RStudio only; not saved
# - saves SACO_6_fst and SACO_6_coords to data/SACO-6.RData
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(geosphere)
})

study_code <- "SACO-6"
study_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SACO-6"
out_dir <- file.path(study_dir, "data")

if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

# -----------------------------
# 1) site metadata
# site numbering / names from Table 2
# coordinates updated from the user-cleaned SACO-6 csv
# -----------------------------
SACO_6_coords <- data.frame(
  site = 1:17,
  site_name = c(
    "Clifford Creek",
    "Martin Creek",
    "Meadow Creek",
    "Moose Creek",
    "Orphan Creek",
    "Star Creek",
    "Swift Creek",
    "Tolan Creek",
    "Warm Springs Creek",
    "East Fork main stem",
    "North Fork Sheephead Creek",
    "Sheephead Creek",
    "Slate Creek",
    "Daly Creek",
    "Skalkaho Creek",
    "Burnt Fork Creek",
    "Willow Creek"
  ),
  lat = c(
    45.9354903,
    45.9730667,
    45.8583098,
    45.9466000,
    45.8790238,
    45.9255063,
    45.8542137,
    45.8115087,
    45.7989586,
    45.8422629,
    45.7563000,
    45.7335192,
    45.6932871,
    46.1906065,
    46.1635144,
    46.4713981,
    46.2937087
  ),
  lon = c(
    -113.6394692,
    -113.8009370,
    -113.8087730,
    -113.7181000,
    -113.6952877,
    -113.5985506,
    -113.7639758,
    -113.8643046,
    -114.0665407,
    -113.9571169,
    -114.4781000,
    -114.4746091,
    -114.2136926,
    -113.8830046,
    -114.0277859,
    -113.9552402,
    -113.9948846
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) pairwise FST matrix
# hard-coded from Table 3; lower triangle mirrored to symmetric matrix
# -----------------------------
SACO_6_fst <- matrix(0, nrow = 17, ncol = 17)
rownames(SACO_6_fst) <- as.character(1:17)
colnames(SACO_6_fst) <- as.character(1:17)

lower_vals <- list(
  c(2,1,0.071),
  c(3,1,0.045), c(3,2,0.082),
  c(4,1,0.041), c(4,2,0.059), c(4,3,0.032),
  c(5,1,0.103), c(5,2,0.120), c(5,3,0.073), c(5,4,0.073),
  c(6,1,0.043), c(6,2,0.067), c(6,3,0.047), c(6,4,0.048), c(6,5,0.095),
  c(7,1,0.068), c(7,2,0.098), c(7,3,0.043), c(7,4,0.060), c(7,5,0.085), c(7,6,0.065),
  c(8,1,0.124), c(8,2,0.165), c(8,3,0.105), c(8,4,0.116), c(8,5,0.188), c(8,6,0.147), c(8,7,0.148),
  c(9,1,0.045), c(9,2,0.090), c(9,3,0.041), c(9,4,0.038), c(9,5,0.106), c(9,6,0.041), c(9,7,0.066), c(9,8,0.113),
  c(10,1,0.016), c(10,2,0.056), c(10,3,0.019), c(10,4,0.018), c(10,5,0.071), c(10,6,0.023), c(10,7,0.030), c(10,8,0.106), c(10,9,0.020),
  c(11,1,0.083), c(11,2,0.140), c(11,3,0.076), c(11,4,0.098), c(11,5,0.156), c(11,6,0.089), c(11,7,0.106), c(11,8,0.131), c(11,9,0.074), c(11,10,0.070),
  c(12,1,0.087), c(12,2,0.122), c(12,3,0.074), c(12,4,0.102), c(12,5,0.161), c(12,6,0.097), c(12,7,0.091), c(12,8,0.126), c(12,9,0.083), c(12,10,0.071), c(12,11,0.041),
  c(13,1,0.088), c(13,2,0.114), c(13,3,0.078), c(13,4,0.089), c(13,5,0.139), c(13,6,0.089), c(13,7,0.100), c(13,8,0.122), c(13,9,0.071), c(13,10,0.068), c(13,11,0.046), c(13,12,0.050),
  c(14,1,0.080), c(14,2,0.111), c(14,3,0.066), c(14,4,0.072), c(14,5,0.119), c(14,6,0.086), c(14,7,0.092), c(14,8,0.105), c(14,9,0.063), c(14,10,0.060), c(14,11,0.042), c(14,12,0.066), c(14,13,0.050),
  c(15,1,0.073), c(15,2,0.109), c(15,3,0.059), c(15,4,0.075), c(15,5,0.119), c(15,6,0.082), c(15,7,0.078), c(15,8,0.093), c(15,9,0.061), c(15,10,0.054), c(15,11,0.042), c(15,12,0.054), c(15,13,0.046), c(15,14,0.009),
  c(16,1,0.097), c(16,2,0.133), c(16,3,0.095), c(16,4,0.096), c(16,5,0.136), c(16,6,0.102), c(16,7,0.130), c(16,8,0.143), c(16,9,0.093), c(16,10,0.081), c(16,11,0.063), c(16,12,0.098), c(16,13,0.079), c(16,14,0.030), c(16,15,0.057),
  c(17,1,0.240), c(17,2,0.284), c(17,3,0.212), c(17,4,0.244), c(17,5,0.300), c(17,6,0.227), c(17,7,0.262), c(17,8,0.212), c(17,9,0.239), c(17,10,0.224), c(17,11,0.215), c(17,12,0.208), c(17,13,0.221), c(17,14,0.195), c(17,15,0.179), c(17,16,0.209)
)

for (x in lower_vals) {
  i <- x[1]; j <- x[2]; v <- x[3]
  SACO_6_fst[i, j] <- v
  SACO_6_fst[j, i] <- v
}

diag(SACO_6_fst) <- 0
SACO_6_fst[SACO_6_fst < 0] <- 0

# -----------------------------
# 3) map plot
# -----------------------------
world_df <- ggplot2::map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

plot_df <- SACO_6_coords %>%
  mutate(label = paste0(site, " ", site_name))

x_pad <- max(0.6, diff(range(plot_df$lon)) * 0.10)
y_pad <- max(0.4, diff(range(plot_df$lat)) * 0.10)

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
    nudge_y = 0.04,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(min(plot_df$lon) - x_pad, max(plot_df$lon) + x_pad),
    ylim = c(min(plot_df$lat) - y_pad, max(plot_df$lat) + y_pad)
  ) +
  labs(
    title = "SACO-6 sampling sites",
    subtitle = "Updated with cleaned user-provided coordinates",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw()

print(p_map)

# -----------------------------
# 4) IBD plot
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(SACO_6_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  fst = SACO_6_fst[upper.tri(SACO_6_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "SACO-6 IBD plot",
    x = "Geographic distance (km)",
    y = "Pairwise FST"
  ) +
  theme_bw()

print(p_ibd)

# -----------------------------
# 5) save outputs
# -----------------------------
save(
  SACO_6_fst,
  SACO_6_coords,
  file = file.path(out_dir, "SACO-6.RData")
)
