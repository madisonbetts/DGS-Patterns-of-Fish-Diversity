# ============================================================
# MAST-1 | Silver Chub
# Macrhybopsis storeriana
# Elbassiouny et al. 2023. North American Journal of Fisheries Management
#
# Objective 2 extraction workflow
# - representative site coordinates hard-coded from Table 1
# - pairwise FST matrix hard-coded from Figure 4B
# - plots shown in RStudio only; not saved
# - saves MAST_1_fst and MAST_1_coords
#   to data/MAST-1.RData
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(geosphere)
})

study_code <- "MAST-1"
study_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/MAST-1"
out_dir <- file.path(study_dir, "data")

if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

# -----------------------------
# 1) site metadata
# representative GPS coordinates from Table 1
# -----------------------------
MAST_1_coords <- data.frame(
  site_id = 1:7,
  site_name = c(
    "Assiniboine River",
    "Lake Erie",
    "Lower Mississippi River",
    "Upper Mississippi River",
    "Missouri River",
    "Ohio River",
    "Wabash River"
  ),
  lat = c(
    49.8674,
    41.8852,
    31.0271,
    44.5836,
    40.3649,
    38.7953,
    38.9348
  ),
  lon = c(
    -97.4375,
    -82.8812,
    -91.5785,
    -92.5102,
    -95.6392,
    -84.1742,
    -87.5163
  ),
  locality_label = rep("Representative GPS from Table 1", 7),
  stringsAsFactors = FALSE
)

write.csv(
  MAST_1_coords,
  file = file.path(study_dir, "MAST-1_coords.csv"),
  row.names = FALSE
)

# object saved in RData should just have numeric site IDs with matching order
MAST_1_coords_out <- MAST_1_coords %>%
  select(site_id, lat, lon)

# -----------------------------
# 2) pairwise FST matrix
# values from Figure 4B as shown in the paper
# order:
# Assiniboine, Lake Erie, Lower Mississippi, Upper Mississippi,
# Missouri, Ohio, Wabash
# -----------------------------
MAST_1_fst <- matrix(0, nrow = 7, ncol = 7)
rownames(MAST_1_fst) <- as.character(1:7)
colnames(MAST_1_fst) <- as.character(1:7)

MAST_1_fst[1,2] <- 0.05
MAST_1_fst[1,3] <- 0.03
MAST_1_fst[1,4] <- 0.00
MAST_1_fst[1,5] <- 0.01
MAST_1_fst[1,6] <- 0.01
MAST_1_fst[1,7] <- 0.02

MAST_1_fst[2,3] <- 0.06
MAST_1_fst[2,4] <- 0.04
MAST_1_fst[2,5] <- 0.05
MAST_1_fst[2,6] <- 0.04
MAST_1_fst[2,7] <- 0.05

MAST_1_fst[3,4] <- 0.00
MAST_1_fst[3,5] <- 0.02
MAST_1_fst[3,6] <- 0.01
MAST_1_fst[3,7] <- 0.00

MAST_1_fst[4,5] <- 0.00
MAST_1_fst[4,6] <- 0.00
MAST_1_fst[4,7] <- 0.00

MAST_1_fst[5,6] <- 0.00
MAST_1_fst[5,7] <- 0.01

MAST_1_fst[6,7] <- 0.00

MAST_1_fst <- MAST_1_fst + t(MAST_1_fst)
diag(MAST_1_fst) <- 0
MAST_1_fst[MAST_1_fst < 0] <- 0

stopifnot(isTRUE(all.equal(MAST_1_fst, t(MAST_1_fst))))
stopifnot(identical(rownames(MAST_1_fst), as.character(MAST_1_coords_out$site_id)))

# -----------------------------
# 3) map plot
# include USA and Canada and zoom to point extent
# -----------------------------
world_df <- ggplot2::map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

plot_df <- MAST_1_coords %>%
  mutate(label = paste0(site_id, " ", site_name))

x_pad <- max(3, diff(range(plot_df$lon)) * 0.12)
y_pad <- max(2, diff(range(plot_df$lat)) * 0.12)

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
    nudge_y = 0.6,
    size = 3,
    check_overlap = TRUE
  ) +
  coord_quickmap(
    xlim = c(min(plot_df$lon) - x_pad, max(plot_df$lon) + x_pad),
    ylim = c(min(plot_df$lat) - y_pad, max(plot_df$lat) + y_pad)
  ) +
  theme_bw() +
  labs(
    title = "MAST-1 sampling locations",
    x = "Longitude",
    y = "Latitude"
  )

print(p_map)

# -----------------------------
# 4) IBD plot
# straight-line distance for QC plotting only
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(MAST_1_coords_out[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  fst = MAST_1_fst[upper.tri(MAST_1_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(
    title = "MAST-1 IBD plot",
    x = "Geographic distance (km)",
    y = "Pairwise FST"
  )

print(p_ibd)

# -----------------------------
# 5) save outputs
# -----------------------------
save(
  MAST_1_fst,
  MAST_1_coords_out,
  file = file.path(out_dir, "MAST-1.RData")
)
