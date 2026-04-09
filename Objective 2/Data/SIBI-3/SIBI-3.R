# ============================================================
# SIBI-3 | Mohave tui chub
# Siphateles bicolor mohavensis
# Chen, Parmenter & May 2013, Conservation Genetics
#
# Objective 2 workflow
# - pairwise Weir and Cockerham's FST values transcribed from Table 2
# - coordinates hardcoded from sample collection coordinates in the paper
# - plots shown in RStudio only; not saved
# - saves SIBI_3_fst and SIBI_3_coords to data/SIBI-3.RData
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(geosphere)
})

study_code <- "SIBI-3"
study_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SIBI-3"
out_dir <- file.path(study_dir, "data")

if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

# -----------------------------
# 1) site metadata
# order follows Table 2:
# Camp Cady, China Lake, Lake Tuendae, MC Spring
# coordinates taken directly from sample collection text
# -----------------------------
SIBI_3_coords_check <- data.frame(
  site = 1:4,
  site_code = c("CC", "CL", "LT", "MC"),
  site_name = c(
    "Camp Cady",
    "China Lake",
    "Lake Tuendae",
    "MC Spring"
  ),
  lat = c(
    34 + 56/60 + 12/3600,
    35 + 42/60 +  0/3600,
    35 +  8/60 + 36/3600,
    35 +  8/60 + 27/3600
  ),
  lon = c(
    -(116 + 36/60 + 42/3600),
    -(117 + 37/60 + 48/3600),
    -(116 +  6/60 + 15/3600),
    -(116 +  6/60 + 15/3600)
  ),
  stringsAsFactors = FALSE
)

SIBI_3_coords <- SIBI_3_coords_check %>%
  transmute(site = site, lat = lat, lon = lon)

# -----------------------------
# 2) pairwise FST matrix
# values transcribed from Table 2
# order:
# Camp Cady, China Lake, Lake Tuendae, MC Spring
# -----------------------------
site_ids <- as.character(SIBI_3_coords$site)

SIBI_3_fst <- matrix(0, nrow = 4, ncol = 4)
rownames(SIBI_3_fst) <- site_ids
colnames(SIBI_3_fst) <- site_ids

upper_vals <- list(
  c(1, 2, 0.07),
  c(1, 3, 0.10),
  c(1, 4, 0.26),
  c(2, 3, 0.02),
  c(2, 4, 0.17),
  c(3, 4, 0.17)
)

for (x in upper_vals) {
  i <- x[1]; j <- x[2]; v <- x[3]
  SIBI_3_fst[i, j] <- v
  SIBI_3_fst[j, i] <- v
}

diag(SIBI_3_fst) <- 0
SIBI_3_fst[is.na(SIBI_3_fst)] <- 0
SIBI_3_fst[SIBI_3_fst < 0] <- 0

stopifnot(identical(as.character(SIBI_3_coords$site), rownames(SIBI_3_fst)))

# -----------------------------
# 3) map plot
# include USA and Canada, then zoom to point extent
# -----------------------------
world_df <- ggplot2::map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

plot_df <- SIBI_3_coords_check %>%
  mutate(label = paste0(site, " ", site_code))

x_pad <- max(1.5, diff(range(plot_df$lon)) * 0.12)
y_pad <- max(1.0, diff(range(plot_df$lat)) * 0.12)

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
    nudge_y = 0.08,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(min(plot_df$lon) - x_pad, max(plot_df$lon) + x_pad),
    ylim = c(min(plot_df$lat) - y_pad, max(plot_df$lat) + y_pad)
  ) +
  labs(
    title = "SIBI-3 sampling sites",
    subtitle = "Coordinates from sample collection section",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw()

print(p_map)

# -----------------------------
# 4) IBD plot
# straight-line distance for QC only
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(SIBI_3_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  fst = SIBI_3_fst[upper.tri(SIBI_3_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "SIBI-3 IBD plot",
    x = "Geographic distance (km)",
    y = "Pairwise FST"
  ) +
  theme_bw()

print(p_ibd)

# -----------------------------
# 5) save outputs
# -----------------------------
save(
  SIBI_3_fst,
  SIBI_3_coords,
  file = file.path(out_dir, "SIBI-3.RData")
)

write.csv(
  SIBI_3_coords_check,
  file = file.path(out_dir, "SIBI-3_coords_check.csv"),
  row.names = FALSE
)
