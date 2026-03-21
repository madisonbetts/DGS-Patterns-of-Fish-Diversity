# -----------------------------
# SAFO-2 Brook trout
# Kanno et al. 2011
# cluster-level microsatellite FST from raw genotypes
# -----------------------------

library(readxl)
library(dplyr)
library(hierfstat)
library(geosphere)
library(ggplot2)

# ------------------------------
# paths
# ------------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SAFO-2"

jeff_file   <- file.path(base_dir, "doi_10_5061_dryad_5f8s2__v20110603/Brook Trout Microsatellite and Spatial Locations/Jefferson Hill-Spruce Brook Microsatellite.txt")
kent_file   <- file.path(base_dir, "doi_10_5061_dryad_5f8s2__v20110603/Brook Trout Microsatellite and Spatial Locations/Kent Falls Brook Microsatellite.txt")
coords_file <- file.path(base_dir, "cluster_xys.xlsx")

# -----------------------------
# 1) read microsatellite files
# -----------------------------
jeff_raw <- read.delim(
  jeff_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

kent_raw <- read.delim(
  kent_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

raw <- bind_rows(jeff_raw, kent_raw)

# -----------------------------
# 2) read cluster coordinates
# preserve spreadsheet order
# -----------------------------
SAFO_2_coords <- read_excel(coords_file) %>%
  rename(
    cluster_ID = 1,
    lon = 2,
    lat = 3
  ) %>%
  mutate(
    cluster_ID = as.character(cluster_ID),
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  )

# assign numeric site IDs in spreadsheet order
SAFO_2_coords <- SAFO_2_coords %>%
  mutate(site = as.character(seq_len(n()))) %>%
  select(site, cluster_ID, lat, lon)

# -----------------------------
# 3) clean genotype table
# -----------------------------
raw <- raw %>%
  mutate(
    Fish_ID    = as.character(Fish_ID),
    Cluster_ID = as.character(Cluster_ID),
    Reach_ID   = as.character(Reach_ID)
  ) %>%
  filter(Cluster_ID %in% SAFO_2_coords$cluster_ID)

# -----------------------------
# 4) extract genotype columns
# first 3 cols are metadata:
# Fish_ID, Cluster_ID, Reach_ID
# -----------------------------
geno_raw <- raw[, 4:ncol(raw), drop = FALSE]

# coerce genotype columns to numeric
geno_raw[] <- lapply(geno_raw, as.numeric)

# remove columns that are entirely NA
geno_raw <- geno_raw[, colSums(!is.na(geno_raw)) > 0, drop = FALSE]

# -----------------------------
# 5) rebuild locus names and create hierfstat input
# each locus coded as allele1*1000 + allele2
# -----------------------------
orig_names <- colnames(geno_raw)
base_names <- sub("[ab]$", "", orig_names)

tab <- table(base_names)
valid_loci <- names(tab)[tab == 2]

hf_loci <- vector("list", length(valid_loci))
names(hf_loci) <- valid_loci

for (loc in valid_loci) {
  tmp <- geno_raw[, base_names == loc, drop = FALSE]
  
  # keep allele a then b in original order
  tmp <- tmp[, order(colnames(tmp)), drop = FALSE]
  
  a1 <- as.numeric(tmp[[1]])
  a2 <- as.numeric(tmp[[2]])
  
  # treat 0 as missing if present
  a1[a1 == 0] <- NA
  a2[a2 == 0] <- NA
  
  # sort alleles within genotype so 132/155 == 155/132
  lo <- pmin(a1, a2, na.rm = FALSE)
  hi <- pmax(a1, a2, na.rm = FALSE)
  
  hf_loci[[loc]] <- ifelse(
    is.na(lo) | is.na(hi),
    NA,
    lo * 1000 + hi
  )
}

hf_dat <- as.data.frame(hf_loci, check.names = FALSE)

# -----------------------------
# 6) link individuals to numeric site IDs
# using cluster order from coordinates spreadsheet
# -----------------------------
site_lookup <- SAFO_2_coords %>%
  select(cluster_ID, site)

raw2 <- raw %>%
  left_join(site_lookup, by = c("Cluster_ID" = "cluster_ID"))

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
SAFO_2_fst <- pairwise.WCfst(hf)

# set negative FST to 0
SAFO_2_fst <- pmax(SAFO_2_fst, 0)

# convert row/col names to numeric indices
n <- nrow(SAFO_2_fst)
rownames(SAFO_2_fst) <- as.character(seq_len(n))
colnames(SAFO_2_fst) <- as.character(seq_len(n))

# force diagonal = 0
diag(SAFO_2_fst) <- 0

# -----------------------------
# 8) map of sampling locations
# -----------------------------
ggplot(SAFO_2_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.002, size = 3.5) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "SAFO-2 sampling locations"
  )

# -----------------------------
# 9) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(SAFO_2_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- SAFO_2_coords$site
colnames(geo_dist_km) <- SAFO_2_coords$site
diag(geo_dist_km) <- 0

# -----------------------------
# 10) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(SAFO_2_fst)[row(SAFO_2_fst)[upper.tri(SAFO_2_fst)]],
  site2   = colnames(SAFO_2_fst)[col(SAFO_2_fst)[upper.tri(SAFO_2_fst)]],
  fst     = SAFO_2_fst[upper.tri(SAFO_2_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 11) IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "SAFO-2 isolation by distance"
  )

# -----------------------------
# 12) quick checks
# -----------------------------
stopifnot(identical(rownames(SAFO_2_fst), SAFO_2_coords$site))
stopifnot(identical(colnames(SAFO_2_fst), SAFO_2_coords$site))
stopifnot(isTRUE(all.equal(SAFO_2_fst, t(SAFO_2_fst))))
stopifnot(all(diag(SAFO_2_fst) == 0))
stopifnot(all(diag(geo_dist_km) == 0))

# -----------------------------
# 13) save RData
# -----------------------------
out_dir <- file.path(base_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  SAFO_2_fst,
  SAFO_2_coords,
  file = file.path(out_dir, "SAFO-2.RData")
)