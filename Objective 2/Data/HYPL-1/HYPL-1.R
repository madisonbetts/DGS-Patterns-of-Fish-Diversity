# -----------------------------
# HYPL-1 Plains Minnow
# from Dryad microsatellite dataset
# full workflow
# -----------------------------

library(readxl)
library(dplyr)
library(adegenet)
library(hierfstat)
library(geosphere)
library(ggplot2)

# ------------------------------
# paths
# ------------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/HYPL-1"

xlsx_file <- file.path(base_dir, "Mol_ecology_Dryad.xlsx")

# -----------------------------
# 1) read sheet
# skip the first weird row
# next row contains the usable headers
# -----------------------------
raw0 <- read_excel(
  xlsx_file,
  sheet = "Final_HYBPLA",
  skip = 1,
  col_names = FALSE
)

true_names <- as.character(unlist(raw0[1, ]))
raw <- raw0[-1, , drop = FALSE]

# -----------------------------
# 2) clean column names
# fill blank locus headers from previous nonblank name
# -----------------------------
true_names[is.na(true_names)] <- ""

for (i in seq_along(true_names)) {
  if (true_names[i] == "" && i > 1) {
    true_names[i] <- true_names[i - 1]
  }
}

true_names <- gsub("\\s+", "", true_names)
true_names <- make.unique(true_names, sep = "_")

colnames(raw) <- true_names
colnames(raw)[1:6] <- c("site", "id", "drainage", "fragment", "lat", "lon")

# drop completely empty rows
raw <- raw[rowSums(!is.na(raw)) > 0, , drop = FALSE]

# -----------------------------
# 3) clean data types
# -----------------------------
raw <- raw %>%
  mutate(
    site     = as.character(site),
    id       = as.character(id),
    drainage = as.character(drainage),
    fragment = as.character(fragment),
    lat      = as.numeric(lat),
    lon      = as.numeric(lon)
  ) %>%
  filter(!is.na(site), !is.na(lat), !is.na(lon))

# -----------------------------
# 4) site coordinates
# mean lat/lon per site
# then convert site labels to numeric indices
# -----------------------------
HYPL_1_coords <- raw %>%
  group_by(site) %>%
  summarise(
    lat = mean(lat, na.rm = TRUE),
    lon = mean(lon, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(site)

HYPL_1_coords$site <- as.character(seq_len(nrow(HYPL_1_coords)))

# -----------------------------
# 5) map of sampling locations
# -----------------------------
ggplot(HYPL_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.5, size = 3.5) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "HYPL-1 sampling locations"
  )

# -----------------------------
# 6) extract microsatellite columns
# after lon column
# each locus is stored as two adjacent allele columns
# -----------------------------
loc_start <- which(colnames(raw) == "lon") + 1
geno_raw <- raw[, loc_start:ncol(raw), drop = FALSE]

# convert to numeric
geno_raw[] <- lapply(geno_raw, as.numeric)

# remove columns that are entirely NA
geno_raw <- geno_raw[, colSums(!is.na(geno_raw)) > 0, drop = FALSE]

# -----------------------------
# 7) rebuild locus names and collapse allele pairs
# into genotype strings "a/b"
# -----------------------------
orig_names <- colnames(geno_raw)
base_names <- sub("_[0-9]+$", "", orig_names)

tab <- table(base_names)
valid_loci <- names(tab)[tab == 2]

geno_list <- vector("list", length(valid_loci))
names(geno_list) <- valid_loci

for (loc in valid_loci) {
  tmp <- geno_raw[, base_names == loc, drop = FALSE]
  
  a1 <- as.numeric(tmp[[1]])
  a2 <- as.numeric(tmp[[2]])
  
  # treat 0 as missing
  a1[a1 == 0] <- NA
  a2[a2 == 0] <- NA
  
  geno_list[[loc]] <- ifelse(
    is.na(a1) | is.na(a2),
    NA,
    paste0(a1, "/", a2)
  )
}

geno <- as.data.frame(geno_list, stringsAsFactors = FALSE)

# -----------------------------
# 8) build genind object
# -----------------------------
genind_obj <- df2genind(
  geno,
  ploidy = 2,
  sep = "/",
  type = "codom",
  pop = as.factor(raw$site),
  NA.char = NA
)

# -----------------------------
# 9) convert to hierfstat format
# -----------------------------
hf <- genind2hierfstat(genind_obj)

# -----------------------------
# 10) pairwise FST
# Weir & Cockerham
# -----------------------------
HYPL_1_fst <- pairwise.WCfst(hf)

# set negative FST to 0
HYPL_1_fst <- pmax(HYPL_1_fst, 0)

# reorder to match original site order used in coords
site_levels <- raw %>%
  distinct(site) %>%
  arrange(site) %>%
  pull(site)

HYPL_1_fst <- HYPL_1_fst[site_levels, site_levels]

# convert row/col names to numeric indices
n <- nrow(HYPL_1_fst)
rownames(HYPL_1_fst) <- as.character(seq_len(n))
colnames(HYPL_1_fst) <- as.character(seq_len(n))

# force diagonal = 0
diag(HYPL_1_fst) <- 0

# -----------------------------
# 11) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(HYPL_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- HYPL_1_coords$site
colnames(geo_dist_km) <- HYPL_1_coords$site

# force diagonal = 0
diag(geo_dist_km) <- 0

# -----------------------------
# 12) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(HYPL_1_fst)[row(HYPL_1_fst)[upper.tri(HYPL_1_fst)]],
  site2   = colnames(HYPL_1_fst)[col(HYPL_1_fst)[upper.tri(HYPL_1_fst)]],
  fst     = HYPL_1_fst[upper.tri(HYPL_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 13) IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "HYPL-1 isolation by distance"
  )

# -----------------------------
# 14) quick checks
# -----------------------------
stopifnot(identical(rownames(HYPL_1_fst), HYPL_1_coords$site))
stopifnot(identical(colnames(HYPL_1_fst), HYPL_1_coords$site))
stopifnot(isTRUE(all.equal(HYPL_1_fst, t(HYPL_1_fst))))
stopifnot(all(diag(HYPL_1_fst) == 0))
stopifnot(all(diag(geo_dist_km) == 0))

# -----------------------------
# 15) save RData
# -----------------------------
out_dir <- file.path(base_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  HYPL_1_fst,
  HYPL_1_coords,
  file = file.path(out_dir, "HYPL-1.RData")
)