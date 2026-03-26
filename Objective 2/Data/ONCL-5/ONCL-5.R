# ============================================================
# ONCL-5
# Oncorhynchus clarkii (cutthroat trout)
# Harris et al. 2022 Conservation Genetics
#
# Workflow:
#   1) read SNP genotypes
#   2) filter loci / individuals
#   3) calculate pairwise Weir & Cockerham FST
#   4) build best-available site coordinates
#   5) plot sites + IBD
#   6) save ONCL-5.RData in /data subfolder
#
# NOTE:
# Site coordinates below are best-available approximate pooled-site centroids.
# The paper maps the 5 pooled sampling sites but does not tabulate explicit
# lat/lon. These coordinates were anchored to named features and tightened
# against the study figure.
# ============================================================

rm(list = ls())

# -----------------------------
# packages
# -----------------------------
libs <- c(
  "adegenet", "hierfstat", "pegas", "ggplot2",
  "dplyr", "maps", "geosphere"
)

for (p in libs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Package not installed: ", p)
  }
}

library(adegenet)
library(hierfstat)
library(pegas)
library(ggplot2)
library(dplyr)
library(maps)
library(geosphere)

# -----------------------------
# paths
# -----------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"
study_dir <- file.path(base_dir, "ONCL-5")
data_dir  <- file.path(study_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# recursively find genotype file anywhere inside ONCL-5
candidate_files <- list.files(
  study_dir,
  pattern = "COGE-D-21-00171_CUT_SNP_Genotypes\\.csv$|CUT_SNP_Genotypes\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(candidate_files) == 0) {
  candidate_files <- list.files(
    study_dir,
    pattern = "\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
}

if (length(candidate_files) == 0) {
  stop("Could not find a genotype csv anywhere inside ", study_dir)
}

infile <- candidate_files[1]
message("Using genotype file: ", infile)

# -----------------------------
# 1) read data
# -----------------------------
raw_dat <- read.csv(infile, check.names = FALSE, stringsAsFactors = FALSE)

if (!all(c("Population", "Sample.Name") %in% names(raw_dat))) {
  stop("Expected columns 'Population' and 'Sample.Name' were not found in the genotype file.")
}

site_levels <- c("BG", "GD", "LPPA", "NE", "LPPB")
site_lookup <- data.frame(
  site_num  = as.character(1:5),
  site_code = site_levels,
  site_name = c(
    "Baker Gulch",
    "Grand Ditch",
    "La Poudre Pass Creek Above Long Draw Reservoir",
    "Neota Creek",
    "La Poudre Pass Creek Below Long Draw Reservoir"
  ),
  stringsAsFactors = FALSE
)

raw_dat <- raw_dat %>%
  mutate(
    Population = trimws(as.character(Population)),
    Population = factor(Population, levels = site_levels),
    Sample.Name = make.unique(as.character(Sample.Name), sep = "_dup")
  ) %>%
  filter(!is.na(Population))

if (nrow(raw_dat) == 0) {
  stop("No rows remained after filtering to expected populations: ", paste(site_levels, collapse = ", "))
}

loci_all <- setdiff(names(raw_dat), c("Population", "Sample.Name"))
raw_dat[loci_all] <- lapply(raw_dat[loci_all], function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "na", "nan")] <- NA
  x
})

# -----------------------------
# 2) filtering
# -----------------------------

# loci with <=20% missing
locus_missing <- sapply(raw_dat[loci_all], function(x) mean(is.na(x)))
keep_loci <- names(locus_missing)[locus_missing <= 0.20]

if (length(keep_loci) == 0) {
  stop("No loci passed the <=20% missing filter.")
}

geno1 <- raw_dat[, c("Population", "Sample.Name", keep_loci), drop = FALSE]

# individuals with <=20% missing across retained loci
ind_missing <- apply(is.na(geno1[, keep_loci, drop = FALSE]), 1, mean)
geno2 <- geno1[ind_missing <= 0.20, , drop = FALSE]

if (nrow(geno2) == 0) {
  stop("No individuals passed the <=20% missing-data filter.")
}

# remove monomorphic loci
is_monomorphic <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(TRUE)
  alleles <- unlist(strsplit(x, split = ":", fixed = TRUE))
  length(unique(alleles)) <= 1
}

mono_flag <- sapply(geno2[keep_loci], is_monomorphic)
keep_loci2 <- keep_loci[!mono_flag]

if (length(keep_loci2) == 0) {
  stop("No loci remained after removing monomorphic loci.")
}

geno3 <- geno2[, c("Population", "Sample.Name", keep_loci2), drop = FALSE]

# optional HWE filter
extract_hw_p <- function(res) {
  if (is.null(res)) return(NA_real_)
  if (!is.null(res$p.value)) return(as.numeric(res$p.value[1]))
  if (is.data.frame(res) || is.matrix(res)) {
    cn <- tolower(colnames(res))
    if (length(cn) > 0) {
      hit <- which(grepl("p", cn))[1]
      if (!is.na(hit)) return(as.numeric(res[1, hit]))
    }
  }
  if (is.numeric(res)) return(as.numeric(res[1]))
  NA_real_
}

hwe_alpha <- 0.05 / (length(keep_loci2) * length(site_levels))

hwe_fail_counts <- sapply(keep_loci2, function(loc) {
  sum(sapply(site_levels, function(pop_i) {
    x <- geno3[[loc]][geno3$Population == pop_i]
    x <- x[!is.na(x)]

    if (length(x) < 5) return(FALSE)

    alleles <- unlist(strsplit(x, split = ":", fixed = TRUE))
    if (length(unique(alleles)) <= 1) return(FALSE)

    x2 <- gsub(":", "/", x, fixed = TRUE)
    res <- tryCatch(
      pegas::hw.test(pegas::as.loci(data.frame(geno = x2)), B = 1000),
      error = function(e) NULL
    )
    p <- extract_hw_p(res)
    !is.na(p) && p < hwe_alpha
  }), na.rm = TRUE)
})

keep_loci3 <- keep_loci2[hwe_fail_counts < 3]
if (length(keep_loci3) == 0) {
  warning("HWE filter removed all loci; falling back to post-missing / post-monomorphic loci.")
  keep_loci3 <- keep_loci2
}

geno_final <- geno3[, c("Population", "Sample.Name", keep_loci3), drop = FALSE]

message("Individuals retained: ", nrow(geno_final))
message("Loci retained: ", length(keep_loci3))
message("Population counts after filtering:")
print(table(geno_final$Population, useNA = "ifany"))

# -----------------------------
# 3) convert SNP strings to genind + hierfstat
# -----------------------------
allele_map <- c(A = 1L, C = 2L, G = 3L, T = 4L)

split_gt <- function(x) {
  if (is.na(x) || !nzchar(x)) return(c(NA_character_, NA_character_))
  parts <- strsplit(x, split = ":", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(c(NA_character_, NA_character_))
  parts
}

fmt_genind <- function(x) {
  sp <- split_gt(x)
  if (any(is.na(sp)) || any(!sp %in% names(allele_map))) return(NA_character_)
  paste0(sprintf("%03d", allele_map[sp[1]]), "/", sprintf("%03d", allele_map[sp[2]]))
}

fmt_hier <- function(x) {
  sp <- split_gt(x)
  if (any(is.na(sp)) || any(!sp %in% names(allele_map))) return(NA_integer_)
  as.integer(paste0(sprintf("%03d", allele_map[sp[1]]), sprintf("%03d", allele_map[sp[2]])))
}

geno_df <- as.data.frame(
  lapply(geno_final[keep_loci3], function(x) vapply(x, fmt_genind, character(1))),
  stringsAsFactors = FALSE
)

sample_ids <- geno_final$Sample.Name
npop <- factor(geno_final$Population, levels = site_levels)

ONCL_5_gen <- adegenet::df2genind(
  X = geno_df,
  ploidy = 2,
  sep = "/",
  ncode = 3,
  ind.names = sample_ids,
  pop = npop,
  NA.char = NA,
  type = "codom"
)

ONCL_5_hier <- data.frame(
  pop = as.numeric(npop),
  as.data.frame(lapply(geno_final[keep_loci3], function(x) vapply(x, fmt_hier, integer(1)))),
  row.names = sample_ids,
  check.names = FALSE
)

# -----------------------------
# 4) pairwise FST matrix
# -----------------------------
ONCL_5_fst <- hierfstat::pairwise.WCfst(ONCL_5_hier)
ONCL_5_fst <- as.matrix(ONCL_5_fst)
ONCL_5_fst[is.na(ONCL_5_fst)] <- 0
ONCL_5_fst[ONCL_5_fst < 0] <- 0
diag(ONCL_5_fst) <- 0

# hierfstat returns populations as numeric ids; remap to site codes
present_pop_ids <- sort(unique(ONCL_5_hier$pop))
present_site_codes <- site_levels[present_pop_ids]

if (nrow(ONCL_5_fst) != length(present_site_codes)) {
  stop("Mismatch between FST matrix dimensions and population ids.")
}

rownames(ONCL_5_fst) <- present_site_codes
colnames(ONCL_5_fst) <- present_site_codes

# expand to full 5 x 5 matrix in study order
fst_full <- matrix(0, nrow = length(site_levels), ncol = length(site_levels))
rownames(fst_full) <- colnames(fst_full) <- site_levels
fst_full[present_site_codes, present_site_codes] <- ONCL_5_fst

ONCL_5_fst <- fst_full
rownames(ONCL_5_fst) <- colnames(ONCL_5_fst) <- as.character(seq_along(site_levels))

# -----------------------------
# 5) site coordinates
# -----------------------------
ONCL_5_coords <- data.frame(
  site = as.character(1:5),
  lat = c(
    40.3260,  # Baker Gulch
    40.4010,  # Grand Ditch pooled reaches centroid
    40.4635,  # La Poudre Pass Creek above Long Draw Reservoir
    40.4675,  # Neota Creek
    40.5039   # La Poudre Pass Creek below Long Draw Reservoir
  ),
  lon = c(
    -105.8567,
    -105.8430,
    -105.8210,
    -105.8540,
    -105.7694
  ),
  site_code = site_levels,
  site_name = site_lookup$site_name,
  stringsAsFactors = FALSE
)

# -----------------------------
# 6) geographic distance matrix
# -----------------------------
coord_mat <- as.matrix(ONCL_5_coords[, c("lon", "lat")])
geo_dist_m <- geosphere::distm(coord_mat, fun = geosphere::distHaversine)
geo_dist_km <- geo_dist_m / 1000
rownames(geo_dist_km) <- colnames(geo_dist_km) <- ONCL_5_coords$site

# -----------------------------
# 7) pairwise dataframe for IBD
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(ONCL_5_fst)[row(ONCL_5_fst)[upper.tri(ONCL_5_fst)]],
  site2   = colnames(ONCL_5_fst)[col(ONCL_5_fst)[upper.tri(ONCL_5_fst)]],
  fst     = ONCL_5_fst[upper.tri(ONCL_5_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 8) map plot
# -----------------------------
states_map <- map_data("state")
canada_map <- map_data("world", region = "Canada")

xpad <- 0.7
ymin <- min(ONCL_5_coords$lat) - 0.5
ymax <- max(ONCL_5_coords$lat) + 0.5
xmin <- min(ONCL_5_coords$lon) - xpad
xmax <- max(ONCL_5_coords$lon) + xpad

ggplot() +
  geom_polygon(
    data = canada_map,
    aes(x = long, y = lat, group = group),
    fill = "grey95", color = "grey70", linewidth = 0.2
  ) +
  geom_polygon(
    data = states_map,
    aes(x = long, y = lat, group = group),
    fill = "grey95", color = "grey70", linewidth = 0.2
  ) +
  geom_point(
    data = ONCL_5_coords,
    aes(x = lon, y = lat),
    size = 3
  ) +
  geom_text(
    data = ONCL_5_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.03,
    size = 4
  ) +
  coord_quickmap(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "ONCL-5 sampling sites"
  )

# -----------------------------
# 9) IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Geographic distance (km)",
    y = expression(F[ST]),
    title = "Isolation by distance: ONCL-5"
  )

# -----------------------------
# 10) save outputs
# -----------------------------

save(
  ONCL_5_fst,
  ONCL_5_coords,
  file = file.path(data_dir, "ONCL-5.RData")
)

# -----------------------------
# quick print
# -----------------------------
print(ONCL_5_coords)
print(round(ONCL_5_fst, 4))
