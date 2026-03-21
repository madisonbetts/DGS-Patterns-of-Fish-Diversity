# -----------------------------
# NOAT-1
# from Dryad microsatellite dataset
# full workflow
# -----------------------------

library(readxl)
library(dplyr)
library(hierfstat)
library(geosphere)
library(ggplot2)

# ------------------------------
# paths
# ------------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/NOAT-1"

xlsx_file <- file.path(base_dir, "Mol_ecology_Dryad.xlsx")

# -----------------------------
# 1) read sheet
# skip the blank first row
# next row contains the true headers
# -----------------------------
raw0 <- read_excel(
  xlsx_file,
  sheet = "NOTATHER",
  skip = 1,
  col_names = FALSE
)

true_names <- as.character(unlist(raw0[1, ]))
raw <- raw0[-1, , drop = FALSE]

# -----------------------------
# 2) clean column names before dplyr
# -----------------------------
true_names[is.na(true_names)] <- ""

for (i in seq_along(true_names)) {
  if (true_names[i] == "" && i > 1) {
    true_names[i] <- true_names[i - 1]
  }
}

true_names <- gsub("\\s+", "", true_names)
true_names <- ifelse(grepl("^[0-9]", true_names), paste0("L", true_names), true_names)
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
    site = as.character(site),
    id   = as.character(id),
    lat  = as.numeric(lat),
    lon  = as.numeric(lon)
  ) %>%
  filter(!is.na(site), !is.na(lat), !is.na(lon))

# -----------------------------
# 4) site coordinates
# mean lat/lon per site
# then convert site labels to numeric indices
# -----------------------------
site_lookup <- raw %>%
  distinct(site) %>%
  arrange(site) %>%
  mutate(site_num = as.character(seq_len(n())))

NOAT_1_coords <- raw %>%
  group_by(site) %>%
  summarise(
    lat = mean(lat, na.rm = TRUE),
    lon = mean(lon, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(site) %>%
  left_join(site_lookup, by = "site") %>%
  transmute(
    site = site_num,
    lat  = lat,
    lon  = lon
  )

# -----------------------------
# 5) map of sampling locations
# -----------------------------
ggplot(NOAT_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.03, size = 3.5) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "NOAT-1 sampling locations"
  )

# -----------------------------
# 6) extract microsatellite columns
# after lon column
# each locus is stored as two adjacent allele columns
# -----------------------------
loc_start <- which(colnames(raw) == "lon") + 1
geno_raw <- raw[, loc_start:ncol(raw), drop = FALSE]

geno_raw[] <- lapply(geno_raw, as.numeric)
geno_raw <- geno_raw[, colSums(!is.na(geno_raw)) > 0, drop = FALSE]

# -----------------------------
# 7) rebuild locus names and create hierfstat input directly
# each locus coded as allele1*1000 + allele2
# -----------------------------
orig_names <- colnames(geno_raw)
base_names <- sub("_[0-9]+$", "", orig_names)

tab <- table(base_names)
valid_loci <- names(tab)[tab == 2]

hf_loci <- vector("list", length(valid_loci))
names(hf_loci) <- valid_loci

for (loc in valid_loci) {
  tmp <- geno_raw[, base_names == loc, drop = FALSE]
  
  a1 <- as.numeric(tmp[[1]])
  a2 <- as.numeric(tmp[[2]])
  
  # treat 0 as missing
  a1[a1 == 0] <- NA
  a2[a2 == 0] <- NA
  
  # sort alleles within genotype so 118/120 == 120/118
  lo <- pmin(a1, a2, na.rm = FALSE)
  hi <- pmax(a1, a2, na.rm = FALSE)
  
  hf_loci[[loc]] <- ifelse(
    is.na(lo) | is.na(hi),
    NA,
    lo * 1000 + hi
  )
}

hf_dat <- as.data.frame(hf_loci)

# add population column as numeric index
raw2 <- raw %>%
  left_join(site_lookup, by = "site")

hf <- data.frame(
  pop = as.numeric(raw2$site_num),
  hf_dat,
  check.names = FALSE
)

# -----------------------------
# 8) pairwise FST
# Weir & Cockerham
# -----------------------------
NOAT_1_fst <- pairwise.WCfst(hf)

# set negative FST to 0
NOAT_1_fst <- pmax(NOAT_1_fst, 0)

# convert row/col names to numeric indices
n <- nrow(NOAT_1_fst)
rownames(NOAT_1_fst) <- as.character(seq_len(n))
colnames(NOAT_1_fst) <- as.character(seq_len(n))

# force diagonal = 0
diag(NOAT_1_fst) <- 0

# -----------------------------
# 9) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(NOAT_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- NOAT_1_coords$site
colnames(geo_dist_km) <- NOAT_1_coords$site
diag(geo_dist_km) <- 0

# -----------------------------
# 10) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(NOAT_1_fst)[row(NOAT_1_fst)[upper.tri(NOAT_1_fst)]],
  site2   = colnames(NOAT_1_fst)[col(NOAT_1_fst)[upper.tri(NOAT_1_fst)]],
  fst     = NOAT_1_fst[upper.tri(NOAT_1_fst)],
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
    title = "NOAT-1 isolation by distance"
  )

# -----------------------------
# 12) quick checks
# -----------------------------
stopifnot(identical(rownames(NOAT_1_fst), NOAT_1_coords$site))
stopifnot(identical(colnames(NOAT_1_fst), NOAT_1_coords$site))
stopifnot(isTRUE(all.equal(NOAT_1_fst, t(NOAT_1_fst))))
stopifnot(all(diag(NOAT_1_fst) == 0))
stopifnot(all(diag(geo_dist_km) == 0))

# -----------------------------
# 13) save RData
# -----------------------------
out_dir <- file.path(base_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  NOAT_1_fst,
  NOAT_1_coords,
  file = file.path(out_dir, "NOAT-1.RData")
)