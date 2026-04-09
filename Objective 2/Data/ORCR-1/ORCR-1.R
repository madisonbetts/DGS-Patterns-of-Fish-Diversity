# -----------------------------
# ORCR-1 Oregon chub
# Oregonichthys crameri
# DeHaan et al. 2012
#
# Collapsed across temporal replicate samples from the same
# physical site (1/1A, 13/13A, 14/14A, 15/15A).
# Pairwise FST values among collapsed sites are the mean of all
# corresponding raw pairwise values from Table 3.
# -----------------------------

library(ggplot2)
library(ggrepel)
library(geosphere)
library(maps)

# -----------------------------
# paths
# -----------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ORCR-1"
data_dir <- file.path(save_dir, "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 0) physical sites / coordinates
#
# The paper provides the site map and names, but not a coordinate
# table for each population. Coordinates below are best-available
# approximate site anchors inferred from Figure 1 and the named
# localities in the Willamette basin.
#
# numeric site IDs follow the collapsed physical-site order:
# 1  = Geren Island North Channel
# 2  = Stayton Public Works Pond
# 3  = Warren Gray Slough
# 4  = Santiam I-5 Channel Pond
# 5  = Finley NWR Gray Creek Swamp
# 6  = Shetzline South Pond
# 7  = Big Island
# 8  = CF Willamette Side Channel
# 9  = Elijah Bristow Berry Slough
# 10 = Elijah Bristow North Slough
# 11 = Dexter Reservoir RV Alcove
# 12 = Dexter Reservoir Alcove - The Pit
# 13 = EF Minnow Creek Pond
# 14 = Hospital Pond
# 15 = Shady Dell Pond
# 16 = Buckhead Creek
# -----------------------------
ORCR_1_coords <- data.frame(
  site = as.character(1:16),
  site_name = c(
    "Geren Island North Channel",
    "Stayton Public Works Pond",
    "Warren Gray Slough",
    "Santiam I-5 Channel Pond",
    "Finley NWR Gray Creek Swamp",
    "Shetzline South Pond",
    "Big Island",
    "CF Willamette Side Channel",
    "Elijah Bristow Berry Slough",
    "Elijah Bristow North Slough",
    "Dexter Reservoir RV Alcove",
    "Dexter Reservoir Alcove - The Pit",
    "EF Minnow Creek Pond",
    "Hospital Pond",
    "Shady Dell Pond",
    "Buckhead Creek"
  ),
  lat = c(
    44.7810, 44.8000, 44.7460, 44.7050,
    44.4070, 44.1090, 44.0280, 43.7800,
    43.9220, 43.9300, 43.9050, 43.8970,
    43.8850, 43.8460, 43.8000, 43.7860
  ),
  lon = c(
    -122.7480, -122.7830, -122.8560, -123.0080,
    -123.3010, -122.9050, -123.0200, -123.0850,
    -122.7810, -122.7700, -122.8120, -122.7970,
    -122.6900, -122.6400, -122.6000, -122.5550
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) plot sites
# US + Canada context, zoomed to site extent
# -----------------------------
usa <- map_data("state")
can <- map_data("world", region = "Canada")

xpad <- 2
ypad <- 1.5

xlim_map <- c(min(ORCR_1_coords$lon) - xpad, max(ORCR_1_coords$lon) + xpad)
ylim_map <- c(min(ORCR_1_coords$lat) - ypad, max(ORCR_1_coords$lat) + ypad)

ggplot() +
  geom_polygon(
    data = usa,
    aes(x = long, y = lat, group = group),
    fill = "grey92", color = "black", linewidth = 0.25
  ) +
  geom_polygon(
    data = can,
    aes(x = long, y = lat, group = group),
    fill = "grey86", color = "black", linewidth = 0.25
  ) +
  geom_point(
    data = ORCR_1_coords,
    aes(x = lon, y = lat),
    size = 2.7
  ) +
  ggrepel::geom_text_repel(
    data = ORCR_1_coords,
    aes(x = lon, y = lat, label = site),
    size = 3.2,
    seed = 1,
    min.segment.length = 0
  ) +
  coord_fixed(xlim = xlim_map, ylim = ylim_map) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "ORCR-1 sampling sites"
  )

# -----------------------------
# 2) raw sample labels from Table 3
# temporal replicate samples retained at this stage
# -----------------------------
raw_labels <- c(
  "1", "1A", "2", "3", "4", "5", "6", "7", "8",
  "9", "10", "11", "12", "13", "13A", "14", "14A", "15", "15A", "16"
)

# map raw labels to collapsed physical sites
collapse_map <- c(
  "1" = "1",
  "1A" = "1",
  "2" = "2",
  "3" = "3",
  "4" = "4",
  "5" = "5",
  "6" = "6",
  "7" = "7",
  "8" = "8",
  "9" = "9",
  "10" = "10",
  "11" = "11",
  "12" = "12",
  "13" = "13",
  "13A" = "13",
  "14" = "14",
  "14A" = "14",
  "15" = "15",
  "15A" = "15",
  "16" = "16"
)

# -----------------------------
# 3) raw pairwise FST matrix from Table 3
# -----------------------------
raw_fst <- matrix(
  0,
  nrow = length(raw_labels),
  ncol = length(raw_labels),
  dimnames = list(raw_labels, raw_labels)
)

set_pair <- function(a, b, val) {
  raw_fst[a, b] <<- val
  raw_fst[b, a] <<- val
}

set_pair("2","1",0.040)
set_pair("2","1A",0.047)

set_pair("3","1",0.019);  set_pair("3","1A",0.025); set_pair("3","2",0.029)

set_pair("4","1",0.044);  set_pair("4","1A",0.050); set_pair("4","2",0.057); set_pair("4","3",0.036)

set_pair("5","1",0.148);  set_pair("5","1A",0.139); set_pair("5","2",0.138); set_pair("5","3",0.138); set_pair("5","4",0.141)

set_pair("6","1",0.177);  set_pair("6","1A",0.169); set_pair("6","2",0.215); set_pair("6","3",0.186); set_pair("6","4",0.181); set_pair("6","5",0.250)

set_pair("7","1",0.069);  set_pair("7","1A",0.067); set_pair("7","2",0.096); set_pair("7","3",0.073); set_pair("7","4",0.064); set_pair("7","5",0.135); set_pair("7","6",0.130)

set_pair("8","1",0.079);  set_pair("8","1A",0.071); set_pair("8","2",0.099); set_pair("8","3",0.067); set_pair("8","4",0.079); set_pair("8","5",0.150); set_pair("8","6",0.156); set_pair("8","7",0.079)

set_pair("9","1",0.069);  set_pair("9","1A",0.062); set_pair("9","2",0.072); set_pair("9","3",0.049); set_pair("9","4",0.054); set_pair("9","5",0.089); set_pair("9","6",0.154); set_pair("9","7",0.058); set_pair("9","8",0.059)

set_pair("10","1",0.059); set_pair("10","1A",0.052); set_pair("10","2",0.063); set_pair("10","3",0.040); set_pair("10","4",0.051); set_pair("10","5",0.103); set_pair("10","6",0.164); set_pair("10","7",0.059); set_pair("10","8",0.055); set_pair("10","9",0.010)

set_pair("11","1",0.072); set_pair("11","1A",0.064); set_pair("11","2",0.083); set_pair("11","3",0.051); set_pair("11","4",0.068); set_pair("11","5",0.122); set_pair("11","6",0.169); set_pair("11","7",0.085); set_pair("11","8",0.053); set_pair("11","9",0.028); set_pair("11","10",0.017)

set_pair("12","1",0.060); set_pair("12","1A",0.055); set_pair("12","2",0.063); set_pair("12","3",0.047); set_pair("12","4",0.054); set_pair("12","5",0.101); set_pair("12","6",0.157); set_pair("12","7",0.064); set_pair("12","8",0.059); set_pair("12","9",0.016); set_pair("12","10",0.012); set_pair("12","11",0.025)

set_pair("13","1",0.056); set_pair("13","1A",0.053); set_pair("13","2",0.063); set_pair("13","3",0.056); set_pair("13","4",0.047); set_pair("13","5",0.112); set_pair("13","6",0.173); set_pair("13","7",0.056); set_pair("13","8",0.074); set_pair("13","9",0.024); set_pair("13","10",0.031); set_pair("13","11",0.058); set_pair("13","12",0.031)

set_pair("13A","1",0.052); set_pair("13A","1A",0.046); set_pair("13A","2",0.069); set_pair("13A","3",0.054); set_pair("13A","4",0.049); set_pair("13A","5",0.113); set_pair("13A","6",0.149); set_pair("13A","7",0.051); set_pair("13A","8",0.063); set_pair("13A","9",0.017); set_pair("13A","10",0.026); set_pair("13A","11",0.050); set_pair("13A","12",0.023); set_pair("13A","13",0.001)

set_pair("14","1",0.055); set_pair("14","1A",0.044); set_pair("14","2",0.073); set_pair("14","3",0.047); set_pair("14","4",0.049); set_pair("14","5",0.116); set_pair("14","6",0.162); set_pair("14","7",0.054); set_pair("14","8",0.056); set_pair("14","9",0.019); set_pair("14","10",0.014); set_pair("14","11",0.037); set_pair("14","12",0.021); set_pair("14","13",0.020); set_pair("14","13A",0.014)

set_pair("14A","1",0.049); set_pair("14A","1A",0.037); set_pair("14A","2",0.070); set_pair("14A","3",0.044); set_pair("14A","4",0.051); set_pair("14A","5",0.101); set_pair("14A","6",0.152); set_pair("14A","7",0.054); set_pair("14A","8",0.052); set_pair("14A","9",0.024); set_pair("14A","10",0.017); set_pair("14A","11",0.033); set_pair("14A","12",0.023); set_pair("14A","13",0.024); set_pair("14A","13A",0.020); set_pair("14A","14",0.001)

set_pair("15","1",0.071); set_pair("15","1A",0.056); set_pair("15","2",0.086); set_pair("15","3",0.066); set_pair("15","4",0.074); set_pair("15","5",0.112); set_pair("15","6",0.156); set_pair("15","7",0.060); set_pair("15","8",0.069); set_pair("15","9",0.023); set_pair("15","10",0.031); set_pair("15","11",0.047); set_pair("15","12",0.033); set_pair("15","13",0.035); set_pair("15","13A",0.025); set_pair("15","14",0.020); set_pair("15","14A",0.021)

set_pair("15A","1",0.071); set_pair("15A","1A",0.058); set_pair("15A","2",0.085); set_pair("15A","3",0.063); set_pair("15A","4",0.069); set_pair("15A","5",0.120); set_pair("15A","6",0.156); set_pair("15A","7",0.052); set_pair("15A","8",0.056); set_pair("15A","9",0.025); set_pair("15A","10",0.031); set_pair("15A","11",0.050); set_pair("15A","12",0.035); set_pair("15A","13",0.031); set_pair("15A","13A",0.023); set_pair("15A","14",0.021); set_pair("15A","14A",0.022); set_pair("15A","15",0.001)

set_pair("16","1",0.087); set_pair("16","1A",0.072); set_pair("16","2",0.113); set_pair("16","3",0.086); set_pair("16","4",0.094); set_pair("16","5",0.131); set_pair("16","6",0.167); set_pair("16","7",0.070); set_pair("16","8",0.081); set_pair("16","9",0.037); set_pair("16","10",0.050); set_pair("16","11",0.062); set_pair("16","12",0.048); set_pair("16","13",0.040); set_pair("16","13A",0.028); set_pair("16","14",0.025); set_pair("16","14A",0.027); set_pair("16","15",0.014); set_pair("16","15A",0.018)

diag(raw_fst) <- 0

# -----------------------------
# 4) collapse temporal replicates to physical sites
# same-site temporal comparisons become zero by construction
# between-site values are averaged across all raw comparisons
# -----------------------------
collapsed_sites <- as.character(1:16)

ORCR_1_fst <- matrix(
  0,
  nrow = length(collapsed_sites),
  ncol = length(collapsed_sites),
  dimnames = list(collapsed_sites, collapsed_sites)
)

for (i in seq_along(collapsed_sites)) {
  for (j in seq_along(collapsed_sites)) {

    si <- collapsed_sites[i]
    sj <- collapsed_sites[j]

    raw_i <- names(collapse_map)[collapse_map == si]
    raw_j <- names(collapse_map)[collapse_map == sj]

    if (si == sj) {
      ORCR_1_fst[si, sj] <- 0
    } else {
      vals <- as.vector(raw_fst[raw_i, raw_j, drop = FALSE])
      ORCR_1_fst[si, sj] <- mean(vals, na.rm = TRUE)
    }
  }
}

diag(ORCR_1_fst) <- 0
ORCR_1_fst[ORCR_1_fst < 0] <- 0

# -----------------------------
# 5) geographic distance matrix (km)
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(ORCR_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- ORCR_1_coords$site
colnames(geo_dist_km) <- ORCR_1_coords$site

# -----------------------------
# 6) IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(ORCR_1_fst)[row(ORCR_1_fst)[upper.tri(ORCR_1_fst)]],
  site2   = colnames(ORCR_1_fst)[col(ORCR_1_fst)[upper.tri(ORCR_1_fst)]],
  fst     = ORCR_1_fst[upper.tri(ORCR_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  stringsAsFactors = FALSE
)

ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.8, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "ORCR-1 isolation by distance"
  )

# -----------------------------
# 7) checks
# -----------------------------
stopifnot(identical(rownames(ORCR_1_fst), ORCR_1_coords$site))
stopifnot(identical(colnames(ORCR_1_fst), ORCR_1_coords$site))
stopifnot(isTRUE(all.equal(ORCR_1_fst, t(ORCR_1_fst))))
stopifnot(all(diag(ORCR_1_fst) == 0))

# -----------------------------
# 8) save
# final coords object follows your usual format
# -----------------------------
ORCR_1_coords <- ORCR_1_coords[, c("site", "lat", "lon")]

save(
  ORCR_1_fst,
  ORCR_1_coords,
  file = file.path(data_dir, "ORCR-1.RData")
)
