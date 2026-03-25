# ETTU-1.R
# Etheostoma tuscumbia
# Source: Fluker et al. 2014, Evolution
# Microsatellite GenePop workflow using HWxtest::genepop.to.genind()

library(ggplot2)
library(maps)
library(geosphere)
library(HWxtest)
library(adegenet)
library(hierfstat)

base_path <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"
study_code <- "ETTU-1"
study_dir <- file.path(base_path, study_code)
data_dir <- file.path(study_dir, "data")

dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

ETTU_1_gen <- genepop.to.genind(
  name = file.path(study_dir, "E_tuscumbia_GENEPOP.txt"),
  quiet = TRUE,
  ncode = 3
)

ettu_levels <- levels(adegenet::pop(ETTU_1_gen))
if (length(ettu_levels) != 13) {
  stop("Expected 13 populations in E_tuscumbia_GENEPOP.txt.")
}

levels(adegenet::pop(ETTU_1_gen)) <- c(
  "TS", "BF", "WH", "PY", "LM", "TH", "PK",
  "BD", "KL", "BY", "BR", "MV", "WC"
)

ETTU_1_coords <- data.frame(
  site = 1:13,
  #site_name = c(
  #  "Tuscumbia Spring",
  #  "Buffler Spring",
  #  "Wheeler Spring",
  #  "Pryor Spring",
  #  "Limestone Creek",
  #  "Thorsen Spring",
  #  "Pickens Spring",
  #  "Beaverdam Spring",
  #  "Kelly Spring",
  #  "Byrd Spring",
  #  "Braham Spring",
  #  "Meridianville Spring",
  #  "Unnamed Spring"
  #),
  lat = c(
    34.72970,
    34.83330,
    34.65220,
    34.67560,
    34.68420,
    34.64000,
    34.66690,
    34.70280,
    34.81560,
    34.66420,
    34.70670,
    34.84530,
    34.92750
  ),
  lon = c(
    -87.70330,
    -86.94750,
    -87.25220,
    -86.95000,
    -86.87830,
    -86.80920,
    -86.81280,
    -86.82940,
    -86.71250,
    -86.58250,
    -86.60050,
    -86.56830,
    -86.39420
  )
)

ETTU_1_hier <- hierfstat::genind2hierfstat(ETTU_1_gen)
ETTU_1_hier <- hierfstat::genind2hierfstat(ETTU_1_gen)
ETTU_1_fst_named <- hierfstat::pairwise.WCfst(ETTU_1_hier)

ETTU_1_fst <- as.matrix(ETTU_1_fst_named)
diag(ETTU_1_fst) <- 0
ETTU_1_fst[ETTU_1_fst < 0] <- 0

ettu_order <- c("TS", "BF", "WH", "PY", "LM", "TH", "PK", "BD", "KL", "BY", "BR", "MV", "WC")
ETTU_1_fst <- ETTU_1_fst[ettu_order, ettu_order]
rownames(ETTU_1_fst) <- colnames(ETTU_1_fst) <- as.character(1:13)

geo_dist_km <- geosphere::distm(
  x = ETTU_1_coords[, c("lon", "lat")],
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- colnames(geo_dist_km) <- as.character(ETTU_1_coords$site)

ETTU_1_ibd <- data.frame(
  site1   = rownames(ETTU_1_fst)[row(ETTU_1_fst)[upper.tri(ETTU_1_fst)]],
  site2   = colnames(ETTU_1_fst)[col(ETTU_1_fst)[upper.tri(ETTU_1_fst)]],
  fst     = ETTU_1_fst[upper.tri(ETTU_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

states_map <- map_data("state")
al_tn <- subset(states_map, region %in% c("alabama", "tennessee"))

map_xlim <- c(-87.9, -86.2)
map_ylim <- c(34.5, 35.0)

ggplot() +
  geom_polygon(
    data = al_tn,
    aes(x = long, y = lat, group = group),
    fill = "grey92", color = "grey40", linewidth = 0.3
  ) +
  geom_point(
    data = ETTU_1_coords,
    aes(x = lon, y = lat),
    size = 2.2
  ) +
  geom_text(
    data = ETTU_1_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.015,
    size = 3
  ) +
  coord_fixed(xlim = map_xlim, ylim = map_ylim) +
  theme_classic() +
  labs(
    title = "ETTU-1 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )

ggplot(ETTU_1_ibd, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.3, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    title = "ETTU-1 IBD QC plot",
    x = "Euclidean distance (km)",
    y = expression(F[ST])
  )

save(
  ETTU_1_fst,
  ETTU_1_coords,
  #geo_dist_km,
  #ETTU_1_ibd,
  file = file.path(data_dir, "ETTU-1.RData")
)
