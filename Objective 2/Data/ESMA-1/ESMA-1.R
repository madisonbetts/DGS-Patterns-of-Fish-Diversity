# -----------------------------
# ESMA-1 — Muskellunge
# Johnson et al. 2025
# updated with user-mapped coordinates
# -----------------------------

setwd("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ESMA-1")

library(geosphere)
library(ggplot2)
library(dplyr)

study_code <- "ESMA-1"
out_dir <- file.path(getwd(), "data")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) coordinates
# updated from user-mapped CSV
# order matches FST matrix below
# -----------------------------
ESMA_1_coords <- data.frame(
  site_id = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
  site_name = c("Buckhannon River", "East Lynn Lake", "Kimsey Run Lake", "Little Kanawha River", "Middle Island Creek", "Monongahela River", "North Bend Lake", "Sandy Creek", "Stonewall Jackson Lake"),
  lat = c(38.8195000, 38.1052347, 38.9581000, 39.0237416, 39.4385570, 39.5850975, 39.2211037, 38.9231905, 39.0035682),
  lon = c(-80.2139000, -82.3583653, -78.8183000, -81.3736106, -81.0081708, -80.0151775, -81.0829904, -81.7239845, -80.4727139),
  locality_label = c("Buckhannon River near Alton / broodstock reach proxy", "East Lynn Lake reservoir centroid / gauge proxy", "Kimsey Run Lake boat ramp / lake proxy", "Little Kanawha River at Glenville / broodstock reach proxy", "Middle Island Creek lower-middle reach centroid proxy", "Monongahela River at Fairmont / WV confluence reach proxy", "North Bend Lake centroid / state park proxy", "Sandy Creek near Thornton / broodstock reach proxy", "Stonewall Jackson Lake centroid"),
  stringsAsFactors = FALSE
)

write.csv(
  ESMA_1_coords,
  file = file.path(getwd(), "ESMA-1_coords.csv"),
  row.names = FALSE
)

# -----------------------------
# 2) pairwise FST matrix
# transcribed from Figure 6A
# order:
# Buckhannon, East Lynn, Kimsey Run, Little Kanawha,
# Middle Island Creek, Monongahela, North Bend,
# Sandy Creek, Stonewall Jackson
# -----------------------------
ESMA_1_fst <- matrix(0, nrow = 9, ncol = 9)
rownames(ESMA_1_fst) <- colnames(ESMA_1_fst) <- as.character(ESMA_1_coords$site_id)

ESMA_1_fst[1,2] <- 0.055
ESMA_1_fst[1,3] <- 0.084
ESMA_1_fst[1,4] <- 0.107
ESMA_1_fst[1,5] <- 0.088
ESMA_1_fst[1,6] <- 0.042
ESMA_1_fst[1,7] <- 0.083
ESMA_1_fst[1,8] <- 0.111
ESMA_1_fst[1,9] <- 0.045

ESMA_1_fst[2,3] <- 0.036
ESMA_1_fst[2,4] <- 0.066
ESMA_1_fst[2,5] <- 0.022
ESMA_1_fst[2,6] <- 0.024
ESMA_1_fst[2,7] <- 0.018
ESMA_1_fst[2,8] <- 0.026
ESMA_1_fst[2,9] <- 0.052

ESMA_1_fst[3,4] <- 0.070
ESMA_1_fst[3,5] <- 0.083
ESMA_1_fst[3,6] <- 0.007
ESMA_1_fst[3,7] <- 0.089
ESMA_1_fst[3,8] <- 0.040
ESMA_1_fst[3,9] <- 0.010

ESMA_1_fst[4,5] <- 0.077
ESMA_1_fst[4,6] <- 0.081
ESMA_1_fst[4,7] <- 0.076
ESMA_1_fst[4,8] <- 0.061
ESMA_1_fst[4,9] <- 0.098

ESMA_1_fst[5,6] <- 0.065
ESMA_1_fst[5,7] <- 0.013
ESMA_1_fst[5,8] <- 0.029
ESMA_1_fst[5,9] <- 0.103

ESMA_1_fst[6,7] <- 0.062
ESMA_1_fst[6,8] <- 0.048
ESMA_1_fst[6,9] <- 0.013

ESMA_1_fst[7,8] <- 0.039
ESMA_1_fst[7,9] <- 0.096

ESMA_1_fst[8,9] <- 0.083

ESMA_1_fst <- ESMA_1_fst + t(ESMA_1_fst)
diag(ESMA_1_fst) <- 0
ESMA_1_fst[is.na(ESMA_1_fst)] <- 0
ESMA_1_fst[ESMA_1_fst < 0] <- 0

# -----------------------------
# 3) straight-line distance matrix (km)
# -----------------------------
coords_mat <- as.matrix(ESMA_1_coords[, c("lon", "lat")])

ESMA_1_dist <- geosphere::distm(coords_mat, fun = geosphere::distHaversine) / 1000
rownames(ESMA_1_dist) <- colnames(ESMA_1_dist) <- as.character(ESMA_1_coords$site_id)

# -----------------------------
# 4) IBD dataframe
# -----------------------------
ibd_df <- data.frame(
  fst = ESMA_1_fst[upper.tri(ESMA_1_fst)],
  dist_km = ESMA_1_dist[upper.tri(ESMA_1_dist)]
)

# -----------------------------
# 5) IBD plot
# -----------------------------
p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(
    x = "Straight-line distance (km)",
    y = "Pairwise FST",
    title = "ESMA-1 IBD (straight-line distance)"
  )

print(p_ibd)

# -----------------------------
# 6) site map
# -----------------------------
world_df <- ggplot2::map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

plot_df <- ESMA_1_coords %>%
  mutate(label = paste0(site_id, " ", site_name))

x_pad <- max(0.8, diff(range(plot_df$lon)) * 0.18)
y_pad <- max(0.8, diff(range(plot_df$lat)) * 0.18)

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
    size = 2.5
  ) +
  geom_text(
    data = plot_df,
    aes(x = lon, y = lat, label = label),
    nudge_y = 0.10,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(min(plot_df$lon) - x_pad, max(plot_df$lon) + x_pad),
    ylim = c(min(plot_df$lat) - y_pad, max(plot_df$lat) + y_pad)
  ) +
  theme_bw() +
  labs(
    title = "ESMA-1 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )

print(p_map)

# -----------------------------
# 7) save outputs
# -----------------------------
save(
  ESMA_1_fst,
  ESMA_1_coords,
  file = file.path(out_dir, "ESMA-1.RData")
)
