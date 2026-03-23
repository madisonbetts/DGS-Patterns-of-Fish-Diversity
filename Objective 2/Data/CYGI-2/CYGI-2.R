# -----------------------------
# CYGI-2 Tallapoosa shiner
# Cyprinella gibbsi
# river tributaries only from Fluker et al. 2014
# pairwise microsatellite FST from Table S2
# -----------------------------

library(ggplot2)
library(ggrepel)
library(dplyr)
library(maps)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/CYGI-2"
data_dir <- file.path(save_dir, "data")

# -----------------------------
# 1) site coordinates
# microsatellite sample localities from Table 1
# -----------------------------
CYGI_2_coords <- data.frame(
  site = as.character(1:4),
  code = c("EM", "SW", "EG", "CH"),
  site_name = c(
    "Emuckfaw Creek at Co. Rd. 98",
    "Sweetwater Creek at Sweetwater Rd.",
    "Eagle Creek at Russell Ida Rd.",
    "Little Chatahospee Creek at Co. Rd. 66"
  ),
  lat = c(
    33.07889,
    32.99184,
    32.95336,
    32.92917
  ),
  lon = c(
    -85.69472,
    -85.69430,
    -85.71649,
    -85.53056
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) FST matrix
# Table S2, Cyprinella gibbsi values are ABOVE diagonal
# order: EM, SW, EG, CH
# -----------------------------
CYGI_2_fst <- matrix(
  c(
    0.000, 0.015, 0.072, 0.025,
    0.015, 0.000, 0.026, -0.002,
    0.072, 0.026, 0.000, 0.038,
    0.025, -0.002, 0.038, 0.000
  ),
  nrow = 4,
  byrow = TRUE
)

rownames(CYGI_2_fst) <- CYGI_2_coords$site
colnames(CYGI_2_fst) <- CYGI_2_coords$site

CYGI_2_fst[CYGI_2_fst < 0] <- 0
diag(CYGI_2_fst) <- 0

# -----------------------------
# 3) inspect
# -----------------------------
print(CYGI_2_coords)
print(round(CYGI_2_fst, 3))

cat("\nSite key:\n")
cat("1 = EM = Emuckfaw Creek\n")
cat("2 = SW = Sweetwater Creek\n")
cat("3 = EG = Eagle Creek\n")
cat("4 = CH = Little Chatahospee Creek\n")

# -----------------------------
# 4) map with US + Canada
# -----------------------------
world_map <- map_data("world") |>
  dplyr::filter(region %in% c("USA", "Canada"))

states_map <- map_data("state")

x_pad <- 8
y_pad <- 6

xlim_use <- range(CYGI_2_coords$lon) + c(-x_pad, x_pad)
ylim_use <- range(CYGI_2_coords$lat) + c(-y_pad, y_pad)

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
    data = CYGI_2_coords,
    aes(x = lon, y = lat),
    size = 2.8
  ) +
  geom_text_repel(
    data = CYGI_2_coords,
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
    title = "CYGI-2 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )

print(p_map)

# -----------------------------
# 5) IBD plot
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(CYGI_2_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- CYGI_2_coords$site
colnames(geo_dist_km) <- CYGI_2_coords$site

ibd_df <- data.frame(
  site1   = rownames(CYGI_2_fst)[row(CYGI_2_fst)[upper.tri(CYGI_2_fst)]],
  site2   = colnames(CYGI_2_fst)[col(CYGI_2_fst)[upper.tri(CYGI_2_fst)]],
  fst     = CYGI_2_fst[upper.tri(CYGI_2_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    title = "CYGI-2 isolation by distance",
    x = "Euclidean distance among site centroids (km)",
    y = expression(F[ST])
  )

print(p_ibd)

# -----------------------------
# 6) save
# -----------------------------
save(
  CYGI_2_fst,
  CYGI_2_coords,
  file = file.path(data_dir, "CYGI-2.RData")
)