
# -----------------------------
# SAFO-4 brook trout
# site coordinates + FST matrix + IBD plot
# -----------------------------

library(ggplot2)
library(geosphere)
library(adegenet)
library(hierfstat)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SAFO-4"

# -----------------------------
# 0) input data path
# expects the Dryad genotype file in the species folder
# -----------------------------
geno_file <- file.path("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SAFO-4/doi_10_5061_dryad_3n5tb2rsd__v20240704/Stack_et_al_microsat_dataset.csv")

if (!file.exists(geno_file)) {
  stop("Could not find Stack_et_al_microsat_dataset.csv in SAFO-4 folder.")
}

safo4_raw <- read.csv(geno_file, stringsAsFactors = FALSE)

# -----------------------------
# 1) site coordinates
# tightened approximate reach centroids digitized from Figure 1
# explicit reach coordinates were not tabulated in the paper / README
# reach order is downstream (1) to upstream (11)
# -----------------------------
SAFO_4_coords <- data.frame(
  site = as.character(1:11),
  lat = c(
    40.5954, 40.5880, 40.5802, 40.5721, 40.5600,
    40.5329, 40.5148, 40.5100, 40.5042, 40.4921,
    40.4523
  ),
  lon = c(
    -105.7770, -105.7790, -105.7820, -105.7855, -105.7920,
    -105.8025, -105.7950, -105.7885, -105.7798, -105.7900,
    -105.8165
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) map of sampling locations
# -----------------------------
ggplot(SAFO_4_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.01, size = 3.5) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "SAFO-4 sampling reaches"
  )

# -----------------------------
# 3) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(SAFO_4_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- SAFO_4_coords$site
colnames(geo_dist_km) <- SAFO_4_coords$site

# -----------------------------
# 4) build genind from paired microsatellite columns
# -----------------------------
loci <- unique(sub("_[12]$", "", names(safo4_raw)[-(1:2)]))

collapse_alleles <- function(a1, a2) {
  ifelse(
    is.na(a1) | is.na(a2),
    NA_character_,
    sprintf("%03d%03d", as.integer(a1), as.integer(a2))
  )
}

geno_df <- as.data.frame(
  setNames(
    lapply(loci, function(loc) {
      collapse_alleles(
        safo4_raw[[paste0(loc, "_1")]],
        safo4_raw[[paste0(loc, "_2")]]
      )
    }),
    loci
  ),
  stringsAsFactors = FALSE
)

SAFO_4_gen <- adegenet::df2genind(
  X = geno_df,
  ploidy = 2,
  sep = "",
  ncode = 3,
  ind.names = safo4_raw$Sample,
  pop = factor(safo4_raw$Reach, levels = 1:11),
  NA.char = NA,
  type = "codom"
)

# -----------------------------
# 5) pairwise FST matrix
# calculated from Dryad microsatellite genotypes
# -----------------------------
SAFO_4_hier <- hierfstat::genind2hierfstat(SAFO_4_gen)
SAFO_4_fst <- as.matrix(hierfstat::pairwise.WCfst(SAFO_4_hier))

SAFO_4_fst <- SAFO_4_fst[as.character(1:11), as.character(1:11)]
diag(SAFO_4_fst) <- 0
SAFO_4_fst[SAFO_4_fst < 0] <- 0

rownames(SAFO_4_fst) <- colnames(SAFO_4_fst) <- as.character(1:11)

# -----------------------------
# 6) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(SAFO_4_fst)[row(SAFO_4_fst)[upper.tri(SAFO_4_fst)]],
  site2   = colnames(SAFO_4_fst)[col(SAFO_4_fst)[upper.tri(SAFO_4_fst)]],
  fst     = SAFO_4_fst[upper.tri(SAFO_4_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 7) IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "SAFO-4 isolation by distance"
  )

# -----------------------------
# 8) quick checks
# -----------------------------
stopifnot(identical(rownames(SAFO_4_fst), SAFO_4_coords$site))
stopifnot(identical(colnames(SAFO_4_fst), SAFO_4_coords$site))
stopifnot(isTRUE(all.equal(SAFO_4_fst, t(SAFO_4_fst))))

# -----------------------------
# 9) save RData
# -----------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  SAFO_4_fst,
  SAFO_4_coords,
  file = file.path(out_dir, "SAFO-4.RData")
)
