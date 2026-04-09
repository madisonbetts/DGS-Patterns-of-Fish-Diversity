# ============================================================
# GICY-1 | Humpback Chub / Gila spp. comparison dataset
# Dzul et al. 2025. North American Journal of Fisheries Management
#
# Expanded Objective 2 workflow
# - includes all groups reported in Table 3:
#   BTC, WGC, HAV, LCR-MIG, LCR-RES, DGP, RTC
# - triangular FST table flipped into a full symmetric matrix
# - coordinates taken directly from Table 1
# - plots shown in RStudio only; not saved
# - saves GICY_1_fst and GICY_1_coords to data/GICY-1.RData
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(geosphere)
})

study_code <- "GICY-1"
study_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/GICY-1"
out_dir <- file.path(study_dir, "data")

if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

# -----------------------------
# 1) site metadata
# coordinates taken directly from Table 1
# order follows Table 3:
# BTC, WGC, HAV, LCR-MIG, LCR-RES, DGP, RTC
# -----------------------------
GICY_1_coords <- data.frame(
  site = 1:7,
  site_abbr = c("BTC", "WGC", "HAV", "LCR-MIG", "LCR-RES", "DGP", "RTC"),
  lat = c(
    33.194344,  # BTC
    35.916233,  # WGC
    36.296067,  # HAV
    36.169014,  # LCR-MIG
    36.173867,  # LCR-RES
    39.296017,  # DGP
    34.474674   # RTC
  ),
  lon = c(
    -104.350733, # BTC
    -113.334649, # WGC
    -112.736870, # HAV
    -111.811668, # LCR-MIG
    -111.708984, # LCR-RES
    -110.062061, # DGP
    -111.800951  # RTC
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) pairwise FST matrix
# Table 3 is triangular, so fill the upper triangle and mirror it
# across the diagonal.
# order: BTC, WGC, HAV, LCR-MIG, LCR-RES, DGP, RTC
# -----------------------------
site_ids <- as.character(GICY_1_coords$site)

GICY_1_fst <- matrix(0, nrow = 7, ncol = 7)
rownames(GICY_1_fst) <- site_ids
colnames(GICY_1_fst) <- site_ids

upper_vals <- list(
  c(1,2,0.212),
  c(1,3,0.251),
  c(1,4,0.197),
  c(1,5,0.158),
  c(1,6,0.273),
  c(1,7,0.335),

  c(2,3,0.014),
  c(2,4,0.005),
  c(2,5,0.005),
  c(2,6,0.037),
  c(2,7,0.116),

  c(3,4,0.011),
  c(3,5,0.009),
  c(3,6,0.050),
  c(3,7,0.149),

  c(4,5,0.002),
  c(4,6,0.034),
  c(4,7,0.106),

  c(5,6,0.031),
  c(5,7,0.082),

  c(6,7,0.168)
)

for (x in upper_vals) {
  i <- x[1]; j <- x[2]; v <- x[3]
  GICY_1_fst[i, j] <- v
  GICY_1_fst[j, i] <- v
}

diag(GICY_1_fst) <- 0
GICY_1_fst[GICY_1_fst < 0] <- 0

# -----------------------------
# 3) map plot
# -----------------------------
world_df <- ggplot2::map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

plot_df <- GICY_1_coords %>%
  mutate(label = paste0(site, " ", site_abbr))

x_pad <- max(0.6, diff(range(plot_df$lon)) * 0.12)
y_pad <- max(0.6, diff(range(plot_df$lat)) * 0.12)

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
    nudge_y = 0.15,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(min(plot_df$lon) - x_pad, max(plot_df$lon) + x_pad),
    ylim = c(min(plot_df$lat) - y_pad, max(plot_df$lat) + y_pad)
  ) +
  labs(
    title = "GICY-1 sampling sites",
    subtitle = "All groups from Tables 1 and 3",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw()

print(p_map)

# -----------------------------
# 4) IBD plot
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(GICY_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  fst = GICY_1_fst[upper.tri(GICY_1_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "GICY-1 IBD plot",
    x = "Geographic distance (km)",
    y = "Pairwise FST"
  ) +
  theme_bw()

print(p_ibd)

# -----------------------------
# 5) save outputs
# -----------------------------
# drop site_abbrv
GICY_1_coords= select(GICY_1_coords, c('site','lat','lon'))

# export
save(
  GICY_1_fst,
  GICY_1_coords,
  file = file.path(out_dir, "GICY-1.RData")
)

write.csv(GICY_1_coords, file.path(out_dir, "GICY-1_coords.csv"), row.names = FALSE)
