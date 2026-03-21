# -----------------------------
# CYLU-1 Red Shiner
# from Dryad microsatellite dataset
# full corrected workflow
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
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/CYLU-1"

xlsx_file <- file.path(base_dir, "Mol_ecology_Dryad.xlsx")

# -----------------------------
# 1) read sheet
# skip the blank first row
# next row contains the true headers
# -----------------------------
raw0 <- read_excel(
  xlsx_file,
  sheet = "Final_CYPLUT",
  skip = 1,
  col_names = FALSE
)

true_names <- as.character(unlist(raw0[1, ]))
raw <- raw0[-1, , drop = FALSE]

# -----------------------------
# 2) clean column names before dplyr
# -----------------------------
true_names[is.na(true_names)] <- ""

# fill blank locus names from previous nonblank header
for (i in seq_along(true_names)) {
  if (true_names[i] == "" && i > 1) {
    true_names[i] <- true_names[i - 1]
  }
}

# make syntactically safe and unique
true_names <- gsub("\\s+", "", true_names)
true_names <- ifelse(grepl("^[0-9]", true_names), paste0("L", true_names), true_names)
true_names <- make.unique(true_names, sep = "_")

colnames(raw) <- true_names
colnames(raw)[1:6] <- c("site", "id", "lat", "lon", "drainage", "fragment")

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
# -----------------------------
CYLU_1_coords <- raw %>%
  group_by(site) %>%
  summarise(
    lat = mean(lat, na.rm = TRUE),
    lon = mean(lon, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(site)

# -----------------------------
# 5) map of sampling locations
# -----------------------------
ggplot(CYLU_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.03, size = 3.5) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "CYLU-1 sampling locations"
  )

# -----------------------------
# 6) extract microsatellite columns
# after fragment column
# each locus is stored as two adjacent allele columns
# -----------------------------
loc_start <- which(colnames(raw) == "fragment") + 1
geno_raw <- raw[, loc_start:ncol(raw), drop = FALSE]

# convert to numeric
geno_raw[] <- lapply(geno_raw, as.numeric)

# remove columns that are entirely NA
geno_raw <- geno_raw[, colSums(!is.na(geno_raw)) > 0, drop = FALSE]

# -----------------------------
# 7) rebuild locus names from original headers
# and collapse allele pairs into genotype strings "a/b"
# -----------------------------
orig_names <- colnames(geno_raw)
base_names <- sub("_[0-9]+$", "", orig_names)

# keep loci with exactly 2 columns
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
CYLU_1_fst <- pairwise.WCfst(hf)
CYLU_1_fst[CYLU_1_fst < 0] <- 0 # make negative fsts = 0 

# reorder to match coords
CYLU_1_fst <- CYLU_1_fst[CYLU_1_coords$site, CYLU_1_coords$site]

# -----------------------------
# 11) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(CYLU_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- CYLU_1_coords$site
colnames(geo_dist_km) <- CYLU_1_coords$site

# -----------------------------
# 12) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(CYLU_1_fst)[row(CYLU_1_fst)[upper.tri(CYLU_1_fst)]],
  site2   = colnames(CYLU_1_fst)[col(CYLU_1_fst)[upper.tri(CYLU_1_fst)]],
  fst     = CYLU_1_fst[upper.tri(CYLU_1_fst)],
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
    title = "CYLU-1 isolation by distance"
  )

# -----------------------------
# 14) quick checks
# -----------------------------
stopifnot(identical(rownames(CYLU_1_fst), CYLU_1_coords$site))
stopifnot(identical(colnames(CYLU_1_fst), CYLU_1_coords$site))
stopifnot(isTRUE(all.equal(CYLU_1_fst, t(CYLU_1_fst))))

# -----------------------------
# 15) save RData
# -----------------------------
out_dir <- file.path(base_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  CYLU_1_fst,
  CYLU_1_coords,
  file = file.path(out_dir, "CYLU-1.RData")
)