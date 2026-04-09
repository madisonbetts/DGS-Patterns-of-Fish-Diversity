# ============================================================
# ACDE-1 | Gulf sturgeon
# Acipenser desotoi
# Dugo et al. 2004 J. Appl. Ichthyol.
# updated with user-provided coordinates
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(geosphere)
})

study_code <- "ACDE-1"
study_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ACDE-1"
out_dir <- file.path(study_dir, "data")

if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

# -----------------------------
# 1) site metadata
# updated from user-provided coordinates
# order must match the FST matrix below:
# Pearl, Pascagoula, Escambia, Yellow, Choctawhatchee, Apalachicola
# -----------------------------
ACDE_1_coords_full <- data.frame(
  site_id = c(1, 2, 3, 4, 5, 6),
  site_name = c("Pearl River", "Pascagoula River", "Escambia River", "Yellow River", "Choctawhatchee River", "Apalachicola River"),
  lat = c(30.3870780, 30.5163794, 30.7763723, 30.7787562, 30.7765890, 30.3181756),
  lon = c(-89.7410808, -88.6027757, -87.3089417, -86.6256816, -85.8258226, -85.0428197),
  locality_label = c("West Pearl River / Rigolets lower holding-area proxy", "West Pascagoula River at Hwy 90 at Gautier holding-area proxy", "Lower Escambia River near mouth / Hwy 90 proxy", "Lower Yellow River near Blackwater Bay / Milton proxy", "Lower Choctawhatchee River / bay mouth proxy", "Lower Apalachicola River near Apalachicola proxy"),
  stringsAsFactors = FALSE
)

ACDE_1_coords <- ACDE_1_coords_full %>%
  select(site_id, lat, lon)

write.csv(
  ACDE_1_coords_full,
  file = file.path(study_dir, "ACDE-1_coords.csv"),
  row.names = FALSE
)

# -----------------------------
# 2) pairwise FST matrix
# Table 3 in the paper
# values reported below the diagonal in this order:
# Pearl, Pascagoula, Escambia, Yellow, Choctawhatchee, Apalachicola
# -----------------------------
ACDE_1_fst <- matrix(0, nrow = 6, ncol = 6)
rownames(ACDE_1_fst) <- as.character(1:6)
colnames(ACDE_1_fst) <- as.character(1:6)

ACDE_1_fst[2, 1] <- 0.033
ACDE_1_fst[3, 1] <- 0.147
ACDE_1_fst[3, 2] <- 0.133
ACDE_1_fst[4, 1] <- 0.176
ACDE_1_fst[4, 2] <- 0.149
ACDE_1_fst[4, 3] <- 0.070
ACDE_1_fst[5, 1] <- 0.224
ACDE_1_fst[5, 2] <- 0.186
ACDE_1_fst[5, 3] <- 0.063
ACDE_1_fst[5, 4] <- 0.076
ACDE_1_fst[6, 1] <- 0.156
ACDE_1_fst[6, 2] <- 0.131
ACDE_1_fst[6, 3] <- 0.074
ACDE_1_fst[6, 4] <- 0.101
ACDE_1_fst[6, 5] <- 0.075

ACDE_1_fst <- ACDE_1_fst + t(ACDE_1_fst)
diag(ACDE_1_fst) <- 0
ACDE_1_fst[is.na(ACDE_1_fst)] <- 0
ACDE_1_fst[ACDE_1_fst < 0] <- 0

stopifnot(isTRUE(all.equal(ACDE_1_fst, t(ACDE_1_fst))))
stopifnot(identical(rownames(ACDE_1_fst), as.character(ACDE_1_coords$site_id)))

# -----------------------------
# 3) map plot
# -----------------------------
world_df <- ggplot2::map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

plot_df <- ACDE_1_coords_full %>%
  mutate(label = paste0(site_id, " ", site_name))

x_pad <- max(1.0, diff(range(plot_df$lon)) * 0.20)
y_pad <- max(0.8, diff(range(plot_df$lat)) * 0.25)

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
    nudge_y = 0.10,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(min(plot_df$lon) - x_pad, max(plot_df$lon) + x_pad),
    ylim = c(min(plot_df$lat) - y_pad, max(plot_df$lat) + y_pad)
  ) +
  theme_bw() +
  labs(
    title = "ACDE-1 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )

print(p_map)

# -----------------------------
# 4) IBD plot
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(ACDE_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  fst = ACDE_1_fst[upper.tri(ACDE_1_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "ACDE-1 IBD plot",
    x = "Geographic distance (km)",
    y = "Pairwise FST"
  ) +
  theme_bw()

print(p_ibd)

# -----------------------------
# 5) save outputs
# -----------------------------
save(
  ACDE_1_fst,
  ACDE_1_coords,
  file = file.path(out_dir, "ACDE-1.RData")
)
