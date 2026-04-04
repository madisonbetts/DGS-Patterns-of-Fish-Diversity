# ============================================================
# MIDO-5 | Smallmouth Bass
# Micropterus dolomieu
# Schall et al. 2017, North American Journal of Fisheries Management
# DOI: 10.1080/02755947.2017.1327902
# ============================================================
# Workflow notes
# - pairwise FST values are transcribed directly from Table 3
# - three low-N sites were excluded by the paper and are not included here:
#   WBW (n = 3), WSC (n = 16), and KC (n = 17)
# - coordinates below are best-available approximate site coordinates
#   inferred from Figure 1, Table 2 site names, and named towns / creek mouths
# - these are approximate georeferences for Objective 2 use and should not be
#   treated as exact electrofishing coordinates
# - negative FST values are set to 0
# - plots are shown in RStudio and are not saved
# ============================================================

library(ggplot2)
library(dplyr)
library(geosphere)
library(maps)

study_code <- "MIDO-5"
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"
study_dir <- file.path(base_dir, study_code)
out_dir <- file.path(study_dir, "data")

if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

# -----------------------------
# 1) site key
# order follows Table 2 / Table 3 in the paper
# -----------------------------
site_key <- data.frame(
  site = 1:22,
  abbrev = c(
    "ARF","BEC","BW","CHIL","DAN","FA","HAL","JRH","LOY","PCM","PCR",
    "PEC","SMA","SRH","SRM","SSN","SW","TC","UBE","WBM","WSM","WYA"
  ),
  site_name = c(
    "Allegheny River, Franklin",
    "Bald Eagle Creek, Castanea",
    "North Branch Susquehanna River, Berwick",
    "Chillisquaque Creek",
    "North Branch Susquehanna River, Danville",
    "North Branch Susquehanna River, Falls",
    "Susquehanna River, Halifax",
    "Juniata River, Howe Township",
    "Loyalsock Creek",
    "Pine Creek Mouth",
    "Pine Creek, Ramsey",
    "Penns Creek",
    "Susquehanna River, Mahantango",
    "Susquehanna River, Harrisburg",
    "Susquehanna River, Marietta",
    "Susquehanna River, Shady Nook",
    "Swatara Creek",
    "Tuscarora Creek",
    "Upper Bald Eagle Creek",
    "West Branch Mahantango Creek",
    "West Branch Susquehanna River, McElhattan",
    "Wyalusing Creek"
  ),
  basin = c(
    "Out of basin",
    "West Branch",
    "Middle",
    "West Branch",
    "Middle",
    "Middle",
    "Lower",
    "Juniata",
    "West Branch",
    "West Branch",
    "West Branch",
    "Lower",
    "Lower",
    "Lower",
    "Lower",
    "Lower",
    "Lower",
    "Juniata",
    "West Branch",
    "Lower",
    "West Branch",
    "Middle"
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) approximate site coordinates
# best-available estimates from Figure 1 + named localities / creek mouths
# -----------------------------
site_lookup <- data.frame(
  site = site_key$site,
  abbrev = site_key$abbrev,
  site_name = site_key$site_name,
  lat = c(
    41.3895000,
    41.0990995,
    41.0540000,
    40.9677589,
    40.9422223,
    41.4542578,
    40.4690000,
    40.4687949,
    41.2470000,
    41.1684111,
    41.2601184,
    40.7727400,
    40.6600000,
    40.3178920,
    40.0508468,
    40.8400268,
    40.2400414,
    40.4853123,
    40.8186924,
    40.6389700,
    41.1400000,
    41.6970200
  ),
  lon = c(
    -79.8203300,
    -77.5035821,
    -76.2340000,
    -76.8127626,
    -76.6068199,
    -75.8846376,
    -76.9310000,
    -77.1133561,
    -76.9860000,
    -77.2835185,
    -77.3334812,
    -76.8642670,
    -76.9160000,
    -76.9077303,
    -76.5402070,
    -76.8157738,
    -76.7278286,
    -77.4889675,
    -78.0306650,
    -76.9391400,
    -77.3600000,
    -76.2307700
  ),
  stringsAsFactors = FALSE
)

MIDO_5_coords <- site_lookup %>%
  select(site, lat, lon)

# -----------------------------
# 3) pairwise FST matrix
# lower diagonal transcribed from Table 3
# then mirrored to upper triangle
# -----------------------------
fill_row <- function(mat, row, vals) {
  stopifnot(length(vals) == row - 1)
  mat[row, seq_along(vals)] <- vals
  mat
}

MIDO_5_fst <- matrix(0, nrow = 22, ncol = 22)

MIDO_5_fst <- fill_row(MIDO_5_fst,  2, c(0.047))
MIDO_5_fst <- fill_row(MIDO_5_fst,  3, c(0.077, 0.023))
MIDO_5_fst <- fill_row(MIDO_5_fst,  4, c(0.076, 0.025, 0.000))
MIDO_5_fst <- fill_row(MIDO_5_fst,  5, c(0.077, 0.020, 0.002, 0.004))
MIDO_5_fst <- fill_row(MIDO_5_fst,  6, c(0.071, 0.018, 0.001, 0.005, 0.001))
MIDO_5_fst <- fill_row(MIDO_5_fst,  7, c(0.064, 0.009, 0.000, 0.002, 0.000, 0.002))
MIDO_5_fst <- fill_row(MIDO_5_fst,  8, c(0.077, 0.023, 0.007, 0.011, 0.009, 0.008, 0.009))
MIDO_5_fst <- fill_row(MIDO_5_fst,  9, c(0.056, 0.007, 0.000, 0.000, 0.000, 0.001, 0.000, 0.005))
MIDO_5_fst <- fill_row(MIDO_5_fst, 10, c(0.062, 0.000, 0.012, 0.012, 0.004, 0.010, 0.003, 0.010, 0.000))
MIDO_5_fst <- fill_row(MIDO_5_fst, 11, c(0.077, 0.011, 0.024, 0.024, 0.016, 0.022, 0.016, 0.019, 0.003, 0.000))
MIDO_5_fst <- fill_row(MIDO_5_fst, 12, c(0.064, 0.020, 0.000, 0.004, 0.000, 0.002, 0.000, 0.000, 0.000, 0.007, 0.019))
MIDO_5_fst <- fill_row(MIDO_5_fst, 13, c(0.079, 0.015, 0.001, 0.010, 0.003, 0.000, 0.004, 0.000, 0.003, 0.003, 0.010, 0.001))
MIDO_5_fst <- fill_row(MIDO_5_fst, 14, c(0.070, 0.016, 0.000, 0.001, 0.004, 0.005, 0.000, 0.003, 0.000, 0.008, 0.021, 0.000, 0.003))
MIDO_5_fst <- fill_row(MIDO_5_fst, 15, c(0.074, 0.018, 0.000, 0.001, 0.004, 0.007, 0.002, 0.003, 0.000, 0.002, 0.009, 0.000, 0.004, 0.000))
MIDO_5_fst <- fill_row(MIDO_5_fst, 16, c(0.071, 0.025, 0.000, 0.003, 0.000, 0.004, 0.001, 0.000, 0.000, 0.010, 0.020, 0.000, 0.002, 0.000, 0.000))
MIDO_5_fst <- fill_row(MIDO_5_fst, 17, c(0.070, 0.016, 0.020, 0.026, 0.032, 0.026, 0.020, 0.020, 0.010, 0.010, 0.017, 0.018, 0.018, 0.012, 0.011, 0.028))
MIDO_5_fst <- fill_row(MIDO_5_fst, 18, c(0.079, 0.026, 0.006, 0.000, 0.010, 0.015, 0.008, 0.012, 0.004, 0.018, 0.028, 0.006, 0.018, 0.001, 0.004, 0.005, 0.023))
MIDO_5_fst <- fill_row(MIDO_5_fst, 19, c(0.115, 0.063, 0.053, 0.054, 0.058, 0.052, 0.049, 0.048, 0.044, 0.051, 0.063, 0.038, 0.059, 0.042, 0.036, 0.053, 0.029, 0.037))
MIDO_5_fst <- fill_row(MIDO_5_fst, 20, c(0.071, 0.014, 0.001, 0.005, 0.003, 0.007, 0.003, 0.000, 0.000, 0.001, 0.012, 0.000, 0.003, 0.002, 0.000, 0.000, 0.021, 0.014, 0.064))
MIDO_5_fst <- fill_row(MIDO_5_fst, 21, c(0.064, 0.001, 0.006, 0.015, 0.009, 0.010, 0.003, 0.012, 0.000, 0.002, 0.012, 0.010, 0.008, 0.006, 0.006, 0.011, 0.010, 0.014, 0.044, 0.007))
MIDO_5_fst <- fill_row(MIDO_5_fst, 22, c(0.069, 0.028, 0.005, 0.007, 0.008, 0.006, 0.005, 0.019, 0.009, 0.027, 0.044, 0.009, 0.019, 0.012, 0.018, 0.013, 0.038, 0.012, 0.058, 0.015, 0.016))

MIDO_5_fst[upper.tri(MIDO_5_fst)] <- t(MIDO_5_fst)[upper.tri(MIDO_5_fst)]
diag(MIDO_5_fst) <- 0
MIDO_5_fst[MIDO_5_fst < 0] <- 0

rownames(MIDO_5_fst) <- as.character(site_key$site)
colnames(MIDO_5_fst) <- as.character(site_key$site)

stopifnot(isTRUE(all.equal(MIDO_5_fst, t(MIDO_5_fst))))

# -----------------------------
# 4) map plot
# include US and Canada and zoom to point extent
# -----------------------------
world_df <- map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

x_pad <- max(0.5, diff(range(site_lookup$lon)) * 0.12)
y_pad <- max(0.5, diff(range(site_lookup$lat)) * 0.12)

p_map <- ggplot() +
  geom_polygon(
    data = world_df,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = site_lookup,
    aes(x = lon, y = lat),
    size = 2
  ) +
  geom_text(
    data = site_lookup,
    aes(x = lon, y = lat, label = site),
    nudge_y = y_pad * 0.03,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(min(site_lookup$lon) - x_pad, max(site_lookup$lon) + x_pad),
    ylim = c(min(site_lookup$lat) - y_pad, max(site_lookup$lat) + y_pad)
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "MIDO-5 sampling sites",
    subtitle = "Best-available approximate site coordinates inferred from map + site names"
  ) +
  theme_bw()

print(p_map)

# -----------------------------
# 5) IBD plot
# straight-line geographic distances among approximate site coordinates
# -----------------------------
geo_dist <- geosphere::distm(
  x = as.matrix(site_lookup[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist[upper.tri(geo_dist)],
  fst = MIDO_5_fst[upper.tri(MIDO_5_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Geographic distance (km)",
    y = "Pairwise FST",
    title = "MIDO-5 IBD plot"
  ) +
  theme_bw()

print(p_ibd)

# -----------------------------
# 6) save objects
# -----------------------------
save(
  MIDO_5_fst,
  MIDO_5_coords,
  file = file.path(out_dir, "MIDO-5.RData")
)
