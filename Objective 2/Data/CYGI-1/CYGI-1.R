# -----------------------------
# CYGI-1 Tallapoosa shiner
# Cyprinella gibbsi
# reservoir tributaries only from Fluker et al. 2014
# pairwise microsatellite FST from Table S1
# -----------------------------

library(ggplot2)
library(ggrepel)
library(dplyr)
library(maps)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/CYGI-1"
data_dir <- file.path(save_dir, "data")

# -----------------------------
# 1) site coordinates
# most downstream sample locality per tributary
# -----------------------------
CYGI_1_coords <- data.frame(
  site = as.character(1:4),
  code = c("OK", "SA", "MA", "EL"),
  site_name = c(
    "Oakachoy Creek at AL Hwy 259",
    "Sandy Creek at Co. Rd. 89",
    "Manoy Creek at US Hwy 280",
    "Elkahatchee Creek at AL Hwy 22"
  ),
  lat = c(
    32.83389,
    32.78028,
    32.87333,
    32.90500
  ),
  lon = c(
    -86.04028,
    -85.65306,
    -85.80917,
    -86.01028
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) FST matrix
# Table S1, Cyprinella gibbsi values are ABOVE diagonal
# order: OK, SA, MA, EL
# -----------------------------
CYGI_1_fst <- matrix(
  c(
    0.000, 0.026, 0.059, 0.053,
    0.026, 0.000, 0.047, 0.036,
    0.059, 0.047, 0.000, 0.082,
    0.053, 0.036, 0.082, 0.000
  ),
  nrow = 4,
  byrow = TRUE
)

rownames(CYGI_1_fst) <- CYGI_1_coords$site
colnames(CYGI_1_fst) <- CYGI_1_coords$site

CYGI_1_fst[CYGI_1_fst < 0] <- 0
diag(CYGI_1_fst) <- 0

# -----------------------------
# 3) inspect
# -----------------------------
print(CYGI_1_coords)
print(round(CYGI_1_fst, 3))

cat("\nSite key:\n")
cat("1 = OK = Oakachoy Creek\n")
cat("2 = SA = Sandy Creek\n")
cat("3 = MA = Manoy Creek\n")
cat("4 = EL = Elkahatchee Creek\n")

# -----------------------------
# 4) map with US + Canada
# -----------------------------
world_map <- map_data("world") |>
  dplyr::filter(region %in% c("USA", "Canada"))

states_map <- map_data("state")

x_pad <- 8
y_pad <- 6

xlim_use <- range(CYGI_1_coords$lon) + c(-x_pad, x_pad)
ylim_use <- range(CYGI_1_coords$lat) + c(-y_pad, y_pad)

p_map <- ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "gray98",
    color = "gray80",
    linewidth = 0.2
  ) +
  geom_polygon(
    data = states_map,
    aes(x = long, y = lat, group = group),
    fill = "gray95",
    color = "gray70",
    linewidth = 0.2
  ) +
  geom_point(
    data = CYGI_1_coords,
    aes(x = lon, y = lat),
    size = 2.8
  ) +
  geom_text_repel(
    data = CYGI_1_coords,
    aes(x = lon, y = lat, label = paste0(site, ": ", code)),
    size = 3.2,
    max.overlaps = 100
  ) +
  coord_quickmap(
    xlim = xlim_use,
    ylim = ylim_use,
    expand = FALSE
  ) +
  theme_classic() +
  labs(
    title = "CYGI-1 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )

print(p_map)

# -----------------------------
# 5) IBD plot
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(CYGI_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- CYGI_1_coords$site
colnames(geo_dist_km) <- CYGI_1_coords$site

ibd_df <- data.frame(
  site1   = rownames(CYGI_1_fst)[row(CYGI_1_fst)[upper.tri(CYGI_1_fst)]],
  site2   = colnames(CYGI_1_fst)[col(CYGI_1_fst)[upper.tri(CYGI_1_fst)]],
  fst     = CYGI_1_fst[upper.tri(CYGI_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    title = "CYGI-1 isolation by distance",
    x = "Euclidean distance among site centroids (km)",
    y = expression(F[ST])
  )

print(p_ibd)

# -----------------------------
# 6) save
# -----------------------------
save(
  CYGI_1_fst,
  CYGI_1_coords,
  file = file.path(data_dir, "CYGI-1.RData")
)