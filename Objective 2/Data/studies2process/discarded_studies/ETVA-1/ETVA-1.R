# -----------------------------
# ETVA-1 | Variegate darter
# Etheostoma variatum
# Argentina et al. 2018 | Freshwater Biology
# Objective 2 workflow
# -----------------------------

# Notes:
# 1) Published pairwise FST are reported among STRUCTURE population clusters,
#    not among all 34 individual sites.
# 2) Therefore, this workflow treats each published cluster as the sampling unit
#    for Objective 2 IBD.
# 3) Site coordinates below come directly from Appendix S1.
# 4) Cluster coordinates used for IBD are simple centroids of the site coordinates
#    belonging to each published cluster included in Table 3.
# 5) Monongahela is excluded because Table 3 reports pairwise FST excluding that cluster.

library(dplyr)
library(geosphere)
library(ggplot2)
library(maps)

# -----------------------------
# 1) site coordinates from Appendix S1
# -----------------------------
site_coords <- tibble::tribble(
  ~site_id, ~site_name, ~cluster,       ~lat,      ~lon,
  1,  "Bell Run",                     "Allegheny", 41.97970, -78.23660,
  2,  "Allegheny River",              "Allegheny", 42.07052, -78.42996,
  3,  "Stillwater Creek",             "Allegheny", 42.04372, -79.23126,
  4,  "Brokenstraw Creek",            "Allegheny", 41.83528, -79.25975,
  5,  "French Creek",                 "Allegheny", 42.01589, -79.50438,
  6,  "South Branch French Creek 1",  "Allegheny", 41.90269, -79.85466,
  7,  "South Branch French Creek 2",  "Allegheny", 41.89770, -79.83500,
  8,  "Mahoning Creek",               "Allegheny", 40.94880, -79.31980,
  9,  "Cheat River",                  "Monongahela", 39.12250, -79.67510,
  10, "Middle Fork Beaver Creek",     "Allegheny", 40.77480, -80.77980,
  11, "Back Fork Elk River",          "Kanawha",   38.47905, -80.38910,
  12, "Birch River",                  "Kanawha",   38.59257, -80.87021,
  13, "Dry Fork",                     "Tug",       37.41460, -81.78312,
  14, "Tug Fork",                     "Tug",       37.46841, -81.86341,
  15, "Dismal River 1",               "Levisa",    37.24127, -81.87148,
  16, "Dismal River 2",               "Levisa",    37.24624, -82.02314,
  17, "Slate Creek",                  "Levisa",    37.27739, -82.09370,
  18, "Levisa Fork 1",                "Levisa",    37.23333, -82.10251,
  19, "Levisa Fork 2",                "Levisa",    37.23682, -82.06055,
  20, "Levisa Fork 3",                "Levisa",    37.30709, -82.15599,
  21, "Levisa Fork 4",                "Levisa",    37.39688, -82.25520,
  22, "Levisa Fork 22-24 centroid",   "Levisa",    37.41162, -82.42752,
  25, "Levisa Fork 5",                "Levisa",    37.46396, -82.52574,
  26, "Russell Fork",                 "Levisa",    37.36608, -82.41172,
  27, "Tygarts Creek",                "Tygarts",   38.46640, -83.05080,
  28, "Big Darby Creek",              "Lower Ohio",39.80716, -83.15501,
  29, "Salt Creek OH",                "Lower Ohio",39.34880, -82.67760,
  30, "Little Miami River",           "Lower Ohio",39.37889, -84.25159,
  31, "Salt Creek IN",                "Lower Ohio",39.40330, -85.20640,
  32, "Whitewater River",             "Lower Ohio",39.31440, -84.90500,
  33, "Middle Fork Red River",        "Kentucky",  37.81510, -83.71840,
  34, "Redbird River",                "Kentucky",  37.18310, -83.60030
)

# -----------------------------
# 2) derive cluster centroids for the published FST units
# -----------------------------
cluster_order <- c("Allegheny", "Kanawha", "Tug", "Levisa", "Tygarts", "Lower Ohio", "Kentucky")

ETVA_1_coords <- site_coords %>%
  filter(cluster %in% cluster_order) %>%
  group_by(cluster) %>%
  summarise(
    lat = mean(lat),
    lon = mean(lon),
    .groups = "drop"
  ) %>%
  mutate(site = seq_len(n())) %>%
  select(site, lat, lon, cluster)

# keep only the columns requested in the saved object
ETVA_1_coords <- ETVA_1_coords %>% select(site, lat, lon)

# separate lookup for plotting labels only
cluster_lookup <- tibble(
  site = seq_along(cluster_order),
  cluster = cluster_order
)

# -----------------------------
# 3) pairwise FST among published clusters (Table 3)
# -----------------------------
fst_names <- cluster_order
ETVA_1_fst <- matrix(0, nrow = length(fst_names), ncol = length(fst_names),
                     dimnames = list(as.character(seq_along(fst_names)), as.character(seq_along(fst_names))))

# row 2: Kanawha vs Allegheny
ETVA_1_fst[2, 1] <- 0.041

# row 3: Tug
ETVA_1_fst[3, 1] <- 0.046
ETVA_1_fst[3, 2] <- 0.056

# row 4: Levisa
ETVA_1_fst[4, 1] <- 0.030
ETVA_1_fst[4, 2] <- 0.041
ETVA_1_fst[4, 3] <- 0.023

# row 5: Tygarts
ETVA_1_fst[5, 1] <- 0.091
ETVA_1_fst[5, 2] <- 0.121
ETVA_1_fst[5, 3] <- 0.114
ETVA_1_fst[5, 4] <- 0.101

# row 6: Lower Ohio
ETVA_1_fst[6, 1] <- 0.028
ETVA_1_fst[6, 2] <- 0.047
ETVA_1_fst[6, 3] <- 0.038
ETVA_1_fst[6, 4] <- 0.023
ETVA_1_fst[6, 5] <- 0.095

# row 7: Kentucky
ETVA_1_fst[7, 1] <- 0.067
ETVA_1_fst[7, 2] <- 0.058
ETVA_1_fst[7, 3] <- 0.065
ETVA_1_fst[7, 4] <- 0.047
ETVA_1_fst[7, 5] <- 0.139
ETVA_1_fst[7, 6] <- 0.042

# mirror and clean
ETVA_1_fst <- ETVA_1_fst + t(ETVA_1_fst)
diag(ETVA_1_fst) <- 0
ETVA_1_fst[ETVA_1_fst < 0] <- 0

# -----------------------------
# 4) quick checks
# -----------------------------
stopifnot(all(dim(ETVA_1_fst) == c(7, 7)))
stopifnot(nrow(ETVA_1_coords) == 7)
stopifnot(all(rownames(ETVA_1_fst) == as.character(ETVA_1_coords$site)))
stopifnot(all(colnames(ETVA_1_fst) == as.character(ETVA_1_coords$site)))

# -----------------------------
# 5) IBD data frame
# -----------------------------
coord_mat <- as.matrix(ETVA_1_coords[, c("lon", "lat")])
geo_dist_km <- geosphere::distm(coord_mat, fun = geosphere::distHaversine) / 1000

ibd_df <- data.frame(
  site1 = row(ETVA_1_fst)[lower.tri(ETVA_1_fst)],
  site2 = col(ETVA_1_fst)[lower.tri(ETVA_1_fst)],
  fst = ETVA_1_fst[lower.tri(ETVA_1_fst)],
  dist_km = geo_dist_km[lower.tri(geo_dist_km)]
) %>%
  mutate(
    fst_linear = fst / (1 - fst),
    label1 = cluster_lookup$cluster[match(site1, cluster_lookup$site)],
    label2 = cluster_lookup$cluster[match(site2, cluster_lookup$site)]
  )

# -----------------------------
# 6) map plot (US + Canada context)
# -----------------------------
world_map <- map_data("world")
plot_coords <- ETVA_1_coords %>%
  left_join(cluster_lookup, by = "site")

xpad <- 2
ypad <- 1.5

map_plot <- ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "grey95", color = "grey70", linewidth = 0.2
  ) +
  geom_point(
    data = plot_coords,
    aes(x = lon, y = lat),
    size = 2.8
  ) +
  geom_text(
    data = plot_coords,
    aes(x = lon, y = lat, label = paste(site, cluster, sep = ": ")),
    size = 3,
    nudge_y = 0.15,
    check_overlap = TRUE
  ) +
  coord_quickmap(
    xlim = range(plot_coords$lon) + c(-xpad, xpad),
    ylim = range(plot_coords$lat) + c(-ypad, ypad)
  ) +
  labs(
    title = "ETVA-1 sampling units",
    subtitle = "Cluster centroids derived from Appendix S1 site coordinates",
    x = "Longitude", y = "Latitude"
  ) +
  theme_bw()

print(map_plot)

# -----------------------------
# 7) IBD plot
# -----------------------------
ibd_plot <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  labs(
    title = "ETVA-1 IBD",
    subtitle = "Cluster-level pairwise FST from published Table 3",
    x = "Geographic distance (km)",
    y = "Pairwise FST"
  ) +
  theme_bw()

print(ibd_plot)

# -----------------------------
# 8) save objects
# -----------------------------
out_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETVA-1/data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save(
  ETVA_1_fst,
  ETVA_1_coords,
  file = file.path(out_dir, "ETVA-1.RData")
)
