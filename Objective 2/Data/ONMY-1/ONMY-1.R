# -----------------------------
# ONMY trout (Andrews et al. 2023)
# 11 sites + pairwise FST from GTseq SNP genotypes
# -----------------------------

library(readxl)
library(dplyr)
library(hierfstat)
library(geosphere)
library(ggplot2)

# ------------------------------
# paths
# ------------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ONMY-1"

xlsx_file <- file.path(base_dir, "mec16810-sup-0001-tables1-s10.xlsx")

# -----------------------------
# 1) read Table S1
# first row is notes, second row is header
# -----------------------------
raw <- read_excel(
  xlsx_file,
  sheet = "TableS1_Samples",
  skip = 1
)

# -----------------------------
# 2) keep metadata + genotype columns
# Location is the population
# genotype columns begin in column K
# -----------------------------
raw <- raw %>%
  mutate(
    Location = as.character(Location)
  ) %>%
  filter(!is.na(Location))

geno_raw <- raw[, 11:ncol(raw), drop = FALSE]

# -----------------------------
# 3) site coordinates
# from Table 1 in the paper
# for populations with two coordinate rows reported under one
# population label, use the mean of those coordinates
# -----------------------------
ONMY_1_coords <- data.frame(
  location = c(
    "Big Jacks Creek",
    "Dry Creek",
    "Duncan Creek",
    "Fawn Creek",
    "Keithley Creek",
    "Little Jacks Creek",
    "Little Weser Creek",
    "Mann Creek",
    "S.F. Callahan Creek",
    "Trail Creek",
    "Williams Creek"
  ),
  lat = c(
    42.56278,
    mean(c(43.68899, 43.71776)),
    42.54793,
    44.38234,
    mean(c(44.55295, 44.55282)),
    42.72889,
    44.51247,
    mean(c(44.52606, 44.54770)),
    48.42030,
    48.56912,
    42.87577
  ),
  lon = c(
    -116.043,
    mean(c(-116.174, -116.135)),
    -116.028,
    -116.059,
    mean(c(-116.885, -116.885)),
    -116.105,
    -116.339,
    mean(c(-116.934, -116.987)),
    -116.031,
    -116.388,
    -116.929
  ),
  stringsAsFactors = FALSE
)

# convert site labels to numeric indices
ONMY_1_coords <- ONMY_1_coords %>%
  arrange(location) %>%
  mutate(site = as.character(seq_len(n()))) %>%
  select(site, location, lat, lon)

# lookup table
site_lookup <- ONMY_1_coords %>%
  select(location, site)

# -----------------------------
# 4) map of sampling locations
# -----------------------------
ggplot(ONMY_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.08, size = 3.5) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "ONMY trout sampling locations"
  )

# -----------------------------
# 5) function to convert SNP genotype strings
# (e.g. AA, AG, TT, 0) into hierfstat diploid codes
# A=1, C=2, G=3, T=4
# genotype code = low*10 + high
# -----------------------------
encode_snp_genotype <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("0", "00", "", "NA")] <- NA
  
  out <- rep(NA_real_, length(x))
  
  allele_map <- c(A = 1, C = 2, G = 3, T = 4)
  
  ok <- !is.na(x) & nchar(x) == 2
  g <- x[ok]
  
  a1 <- substring(g, 1, 1)
  a2 <- substring(g, 2, 2)
  
  n1 <- unname(allele_map[a1])
  n2 <- unname(allele_map[a2])
  
  bad <- is.na(n1) | is.na(n2)
  code <- rep(NA_real_, length(g))
  
  if (any(!bad)) {
    lo <- pmin(n1[!bad], n2[!bad])
    hi <- pmax(n1[!bad], n2[!bad])
    code[!bad] <- lo * 10 + hi
  }
  
  out[ok] <- code
  out
}

# -----------------------------
# 6) build hierfstat input directly
# -----------------------------
hf_loci <- lapply(geno_raw, encode_snp_genotype)
hf_dat  <- as.data.frame(hf_loci, check.names = FALSE)

raw2 <- raw %>%
  left_join(site_lookup, by = c("Location" = "location"))

hf <- data.frame(
  pop = as.numeric(raw2$site),
  hf_dat,
  check.names = FALSE
)

# remove loci that are all NA
hf <- hf[, c(TRUE, colSums(!is.na(hf[, -1, drop = FALSE])) > 0)]

# -----------------------------
# 7) pairwise FST
# Weir & Cockerham
# -----------------------------
ONMY_1_fst <- pairwise.WCfst(hf)

# set negative FST to 0
ONMY_1_fst <- pmax(ONMY_1_fst, 0)

# convert row/col names to numeric indices
n <- nrow(ONMY_1_fst)
rownames(ONMY_1_fst) <- as.character(seq_len(n))
colnames(ONMY_1_fst) <- as.character(seq_len(n))

# force diagonal = 0
diag(ONMY_1_fst) <- 0

# -----------------------------
# 8) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(ONMY_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- ONMY_1_coords$site
colnames(geo_dist_km) <- ONMY_1_coords$site
diag(geo_dist_km) <- 0

# -----------------------------
# 9) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(ONMY_1_fst)[row(ONMY_1_fst)[upper.tri(ONMY_1_fst)]],
  site2   = colnames(ONMY_1_fst)[col(ONMY_1_fst)[upper.tri(ONMY_1_fst)]],
  fst     = ONMY_1_fst[upper.tri(ONMY_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 10) IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "ONMY trout isolation by distance"
  )

# -----------------------------
# 11) quick checks
# -----------------------------
stopifnot(identical(rownames(ONMY_1_fst), ONMY_1_coords$site))
stopifnot(identical(colnames(ONMY_1_fst), ONMY_1_coords$site))
stopifnot(isTRUE(all.equal(ONMY_1_fst, t(ONMY_1_fst))))
stopifnot(all(diag(ONMY_1_fst) == 0))
stopifnot(all(diag(geo_dist_km) == 0))

# -----------------------------
# 12) save RData
# -----------------------------
out_dir <- file.path(base_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  ONMY_1_fst,
  ONMY_1_coords,
  file = file.path(out_dir, "ONMY-1.RData")
)