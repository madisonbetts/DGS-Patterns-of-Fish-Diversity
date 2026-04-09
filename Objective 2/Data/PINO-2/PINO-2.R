# ========================================
# PINO-2
# Pimephales notatus
# bluntnose minnow
#
# Collapsed from fine-scale sampling reaches in:
# Smith 2016 thesis on low-head dams in the
# Vermilion / North Fork Vermilion rivers.
#
# To make the dataset usable for Objective 2 IBD workflows,
# sites are collapsed into 3 geographic groups:
#   1 = Vermilion River mainstem
#   2 = North Fork Vermilion below Ellsworth Park Dam
#   3 = North Fork Vermilion above Ellsworth Park Dam
#
# Pairwise FST for collapsed groups is the mean of all
# pairwise site-by-site FST values between groups
# after setting negative values to 0.
# ========================================

library(ggplot2)
library(geosphere)
library(maps)

study_code <- "PINO-2"

base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/PINO-2"
data_dir <- file.path(base_dir, "data")

dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) original site coordinates
# Appendix Table I
# -----------------------------
site_coords <- data.frame(
  site = c(
    "V_BD1", "V_BD2", "V_P1", "V_P2", "V_R1", "V_R2",
    "NF_BD1", "NF_BD2", "NF_P1", "NF_P2", "NF_R1", "NF_R2"
  ),
  lat = c(
    40.117908, 40.121214, 40.122522, 40.117217, 40.120531, 40.122811,
    40.120542, 40.122756, 40.124717, 40.128669, 40.132097, 40.129844
  ),
  lon = c(
    -87.629033, -87.631983, -87.633308, -87.650892, -87.653156, -87.660347,
    -87.641531, -87.639333, -87.639153, -87.642339, -87.647511, -87.650286
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) collapsed groups
# -----------------------------
group_lookup <- data.frame(
  site = site_coords$site,
  group_id = c(
    "1", "1", "1", "1", "1", "1",
    "2", "2",
    "3", "3", "3", "3"
  ),
  group_name = c(
    rep("Vermilion River", 6),
    rep("North Fork below dam", 2),
    rep("North Fork above dam", 4)
  ),
  stringsAsFactors = FALSE
)

PINO_2_coords <- aggregate(
  cbind(lat, lon) ~ group_id + group_name,
  data = merge(site_coords, group_lookup, by = "site", sort = FALSE),
  FUN = mean
)

PINO_2_coords <- PINO_2_coords[order(as.numeric(PINO_2_coords$group_id)), ]
PINO_2_coords <- data.frame(
  site_id = PINO_2_coords$group_id,
  lat = PINO_2_coords$lat,
  lon = PINO_2_coords$lon,
  stringsAsFactors = FALSE
)

site_lookup <- data.frame(
  site_id = c("1", "2", "3"),
  site_name = c(
    "Vermilion River",
    "North Fork below dam",
    "North Fork above dam"
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 3) original BLS pairwise FST
# Table 2.6
# negatives set to 0 per workflow
# -----------------------------
orig_sites <- c(
  "V_BD1", "V_BD2", "V_P1", "V_P2", "V_R1", "V_R2",
  "NF_BD1", "NF_BD2", "NF_P1", "NF_P2", "NF_R1", "NF_R2"
)

PINO_raw_fst <- matrix(0, nrow = length(orig_sites), ncol = length(orig_sites),
                       dimnames = list(orig_sites, orig_sites))

fill_pair <- function(a, b, val) {
  PINO_raw_fst[a, b] <<- val
  PINO_raw_fst[b, a] <<- val
}

fill_pair("V_BD2","V_BD1", 0.0001)
fill_pair("V_P1","V_BD1",-0.0005)
fill_pair("V_P1","V_BD2", 0.0055)
fill_pair("V_P2","V_BD1", 0.0017)
fill_pair("V_P2","V_BD2", 0.0021)
fill_pair("V_P2","V_P1", 0.0063)
fill_pair("V_R1","V_BD1",-0.0018)
fill_pair("V_R1","V_BD2",-0.0025)
fill_pair("V_R1","V_P1", 0.0166)
fill_pair("V_R1","V_P2",-0.0024)
fill_pair("V_R2","V_BD1", 0.0013)
fill_pair("V_R2","V_BD2",-0.0003)
fill_pair("V_R2","V_P1", 0.0038)
fill_pair("V_R2","V_P2",-0.0020)
fill_pair("V_R2","V_R1",-0.0023)
fill_pair("NF_BD1","V_BD1",-0.0029)
fill_pair("NF_BD1","V_BD2",-0.0028)
fill_pair("NF_BD1","V_P1", 0.0030)
fill_pair("NF_BD1","V_P2",-0.0002)
fill_pair("NF_BD1","V_R1",-0.0052)
fill_pair("NF_BD1","V_R2", 0.0014)
fill_pair("NF_BD2","V_BD1",-0.0038)
fill_pair("NF_BD2","V_BD2",-0.0015)
fill_pair("NF_BD2","V_P1", 0.0034)
fill_pair("NF_BD2","V_P2", 0.0001)
fill_pair("NF_BD2","V_R1",-0.0021)
fill_pair("NF_BD2","V_R2", 0.0013)
fill_pair("NF_BD2","NF_BD1",-0.0019)
fill_pair("NF_P1","V_BD1", 0.0227)
fill_pair("NF_P1","V_BD2", 0.0157)
fill_pair("NF_P1","V_P1", 0.0372)
fill_pair("NF_P1","V_P2", 0.0215)
fill_pair("NF_P1","V_R1", 0.0234)
fill_pair("NF_P1","V_R2", 0.0180)
fill_pair("NF_P1","NF_BD1", 0.0146)
fill_pair("NF_P1","NF_BD2", 0.0149)
fill_pair("NF_P2","V_BD1",-0.0118)
fill_pair("NF_P2","V_BD2",-0.0153)
fill_pair("NF_P2","V_P1",-0.0145)
fill_pair("NF_P2","V_P2", 0.0017)
fill_pair("NF_P2","V_R1",-0.0026)
fill_pair("NF_P2","V_R2",-0.0059)
fill_pair("NF_P2","NF_BD1", 0.0007)
fill_pair("NF_P2","NF_BD2",-0.0186)
fill_pair("NF_P2","NF_P1",-0.0132)
fill_pair("NF_R1","V_BD1", 0.0136)
fill_pair("NF_R1","V_BD2", 0.0098)
fill_pair("NF_R1","V_P1", 0.0185)
fill_pair("NF_R1","V_P2", 0.0142)
fill_pair("NF_R1","V_R1", 0.0075)
fill_pair("NF_R1","V_R2", 0.0112)
fill_pair("NF_R1","NF_BD1", 0.0077)
fill_pair("NF_R1","NF_BD2", 0.0115)
fill_pair("NF_R1","NF_P1",-0.0264)
fill_pair("NF_R1","NF_P2", 0.0068)
fill_pair("NF_R2","V_BD1", 0.0185)
fill_pair("NF_R2","V_BD2", 0.0186)
fill_pair("NF_R2","V_P1", 0.0289)
fill_pair("NF_R2","V_P2", 0.0166)
fill_pair("NF_R2","V_R1", 0.0151)
fill_pair("NF_R2","V_R2", 0.0154)
fill_pair("NF_R2","NF_BD1", 0.0182)
fill_pair("NF_R2","NF_BD2", 0.0171)
fill_pair("NF_R2","NF_P1",-0.0035)
fill_pair("NF_R2","NF_P2", 0.0111)
fill_pair("NF_R2","NF_R1", 0.0019)

PINO_raw_fst[PINO_raw_fst < 0] <- 0
diag(PINO_raw_fst) <- 0

# -----------------------------
# 4) collapse to 3 groups
# mean pairwise FST between groups
# -----------------------------
groups <- split(group_lookup$site, group_lookup$group_id)

group_ids <- c("1", "2", "3")
PINO_2_fst <- matrix(0, nrow = 3, ncol = 3, dimnames = list(group_ids, group_ids))

for (i in 1:3) {
  for (j in 1:3) {
    if (i == j) {
      PINO_2_fst[i, j] <- 0
    } else {
      vals <- as.vector(PINO_raw_fst[groups[[group_ids[i]]], groups[[group_ids[j]]]])
      PINO_2_fst[i, j] <- mean(vals, na.rm = TRUE)
    }
  }
}

diag(PINO_2_fst) <- 0

# -----------------------------
# 5) geographic distance matrix
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(PINO_2_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- PINO_2_coords$site_id
colnames(geo_dist_km) <- PINO_2_coords$site_id

# -----------------------------
# 6) IBD dataframe
# -----------------------------
ibd_df <- data.frame(
  site1 = rownames(PINO_2_fst)[row(PINO_2_fst)[upper.tri(PINO_2_fst)]],
  site2 = colnames(PINO_2_fst)[col(PINO_2_fst)[upper.tri(PINO_2_fst)]],
  fst = PINO_2_fst[upper.tri(PINO_2_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  stringsAsFactors = FALSE
)

# -----------------------------
# 7) map
# -----------------------------
state_map <- map_data("state")

plot_sites <- merge(PINO_2_coords, site_lookup, by = "site_id", sort = FALSE)

lon_rng <- range(plot_sites$lon, na.rm = TRUE)
lat_rng <- range(plot_sites$lat, na.rm = TRUE)

x_pad <- max(0.08, diff(lon_rng) * 0.35)
y_pad <- max(0.06, diff(lat_rng) * 0.35)

xlim_use <- c(lon_rng[1] - x_pad, lon_rng[2] + x_pad)
ylim_use <- c(lat_rng[1] - y_pad, lat_rng[2] + y_pad)

map_plot <- ggplot() +
  geom_polygon(
    data = state_map,
    aes(x = long, y = lat, group = group),
    fill = "grey96",
    color = "grey55",
    linewidth = 0.25
  ) +
  geom_point(
    data = plot_sites,
    aes(x = lon, y = lat),
    size = 3
  ) +
  geom_text(
    data = plot_sites,
    aes(x = lon, y = lat, label = paste0(site_id, ". ", site_name)),
    nudge_y = 0.003,
    size = 3.3
  ) +
  coord_fixed(xlim = xlim_use, ylim = ylim_use) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "PINO-2 sampling locations"
  )

print(map_plot)

ggsave(
  filename = file.path(base_dir, "PINO-2_map.png"),
  plot = map_plot,
  width = 7,
  height = 5.2,
  dpi = 300
)

# -----------------------------
# 8) IBD plot
# -----------------------------
ibd_plot <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "PINO-2 isolation by distance"
  )

print(ibd_plot)

ggsave(
  filename = file.path(base_dir, "PINO-2_IBD.png"),
  plot = ibd_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# -----------------------------
# 9) checks + save
# -----------------------------
stopifnot(identical(rownames(PINO_2_fst), PINO_2_coords$site_id))
stopifnot(identical(colnames(PINO_2_fst), PINO_2_coords$site_id))
stopifnot(isTRUE(all.equal(PINO_2_fst, t(PINO_2_fst))))
stopifnot(isTRUE(all.equal(geo_dist_km, t(geo_dist_km))))

save(
  PINO_2_fst,
  PINO_2_coords,
  file = file.path(data_dir, "PINO-2.RData")
)
