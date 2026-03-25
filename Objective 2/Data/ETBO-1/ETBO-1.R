# ETBO-1.R
# Etheostoma boschungi
# Source: Fluker et al. 2014, Evolution
# Microsatellite GenePop workflow using HWxtest::genepop.to.genind()

library(ggplot2)
library(maps)
library(geosphere)
library(HWxtest)
library(adegenet)
library(hierfstat)

base_path <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"
study_code <- "ETBO-1"
study_dir <- file.path(base_path, study_code)
data_dir <- file.path(study_dir, "data")

dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

ETBO_1_gen <- genepop.to.genind(
  name = file.path(study_dir, "E_boschungi_GENEPOP.txt"),
  quiet = TRUE,
  ncode = 3
)

etbo_prefix <- gsub("[0-9]+$", "", adegenet::indNames(ETBO_1_gen))

etbo_pop_map <- c(
  TBC = "CB",
  TMC = "DD",
  BCH = "CH",
  TLS = "SH",
  BNF = "NF",
  BSF = "SF",
  TSC = "SW",
  FBF = "FL"
)

if (!all(etbo_prefix %in% names(etbo_pop_map))) {
  stop("Unexpected ETBO individual prefixes detected in GenePop file.")
}

adegenet::pop(ETBO_1_gen) <- factor(
  etbo_pop_map[etbo_prefix],
  levels = c("CB", "DD", "CH", "SH", "NF", "SF", "SW", "FL")
)

ETBO_1_coords <- data.frame(
  site = 1:8,
  #site_name = c(
  #  "Cooper Branch",
  #  "Dodd Site",
  #  "Chief Creek",
  #  "Shoal Creek",
  #  "North Fork Buffalo River",
  #  "South Fork Buffalo River",
  #  "Swan Creek",
  #  "Brier Fork Flint River"
  #),
  lat = c(
    35.01589,
    35.06060,
    35.35972,
    35.32694,
    35.42472,
    35.36194,
    34.83260,
    35.00708
  ),
  lon = c(
    -87.82322,
    -87.77250,
    -87.41917,
    -87.27278,
    -87.28139,
    -87.25528,
    -86.94953,
    -86.66580
  )
)

# -----------------------------
# pairwise FST from microsatellite data
# -----------------------------
ETBO_1_hier <- hierfstat::genind2hierfstat(ETBO_1_gen)
ETBO_1_fst_named <- hierfstat::pairwise.WCfst(ETBO_1_hier)

ETBO_1_fst <- as.matrix(ETBO_1_fst_named)
diag(ETBO_1_fst) <- 0
ETBO_1_fst[ETBO_1_fst < 0] <- 0

etbo_order <- c("CB", "DD", "CH", "SH", "NF", "SF", "SW", "FL")
ETBO_1_fst <- ETBO_1_fst[etbo_order, etbo_order]

rownames(ETBO_1_fst) <- colnames(ETBO_1_fst) <- as.character(1:8)

geo_dist_km <- geosphere::distm(
  x = ETBO_1_coords[, c("lon", "lat")],
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- colnames(geo_dist_km) <- as.character(ETBO_1_coords$site)

ETBO_1_ibd <- data.frame(
  site1   = rownames(ETBO_1_fst)[row(ETBO_1_fst)[upper.tri(ETBO_1_fst)]],
  site2   = colnames(ETBO_1_fst)[col(ETBO_1_fst)[upper.tri(ETBO_1_fst)]],
  fst     = ETBO_1_fst[upper.tri(ETBO_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

states_map <- map_data("state")
al_tn <- subset(states_map, region %in% c("alabama", "tennessee"))

map_xlim <- c(-88.1, -86.4)
map_ylim <- c(34.7, 35.6)

ggplot() +
  geom_polygon(
    data = al_tn,
    aes(x = long, y = lat, group = group),
    fill = "grey92", color = "grey40", linewidth = 0.3
  ) +
  geom_point(
    data = ETBO_1_coords,
    aes(x = lon, y = lat),
    size = 2.2
  ) +
  geom_text(
    data = ETBO_1_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.03,
    size = 3
  ) +
  coord_fixed(xlim = map_xlim, ylim = map_ylim) +
  theme_classic() +
  labs(
    title = "ETBO-1 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )

ggplot(ETBO_1_ibd, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.3, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    title = "ETBO-1 IBD QC plot",
    x = "Euclidean distance (km)",
    y = expression(F[ST])
  )

save(
  ETBO_1_fst,
  ETBO_1_coords,
  #geo_dist_km,
  #ETBO_1_ibd,
  file = file.path(data_dir, "ETBO-1.RData")
)
