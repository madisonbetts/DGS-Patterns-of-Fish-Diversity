# ============================================================
# CAFU-1
# Owens sucker | Catostomus fumeiventris
# Campbell et al. 2023, Transactions of the American Fisheries Society
# DOI: 10.1002/tafs.10407
#
# Objective 2 workflow:
# - read plink-pruned.vcf as genlight
# - rebuild the paper's 108-individual population-genetics subset from Supplemental Table S1
# - retain Owens sucker individuals only
# - assign populations from location metadata
# - hard-code Owens site coordinates inferred from Fig. 1B
# - calculate pairwise FST with dartR::gl.fst.pop
# - set FST < 0 to 0
# - build matching coords df with numeric site IDs
# - plot sites and IBD in RStudio (do not save plots)
# - save CAFU_1_fst and CAFU_1_coords to data/CAFU-1.RData
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(stringr)
  library(tidyr)
  library(vcfR)
  library(dartR)
  library(ggplot2)
  library(maps)
  library(geosphere)
})

# -----------------------------
# 0) paths
# -----------------------------
study_code <- "CAFU-1"
species_code <- "CAFU"

base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/CAFU-1"
base_dir <- path.expand(base_dir)

data_dir <- file.path(base_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

vcf_file <- file.path(base_dir, "plink-pruned.vcf")
s1_file  <- file.path(base_dir, "Supplemental Table S1.xlsx")
s2_file  <- file.path(base_dir, "Supplemental Table S2.xlsx")
pdf_file <- file.path(base_dir, "TAFS_TAFS10407.pdf")

stopifnot(file.exists(vcf_file), file.exists(s1_file))

# -----------------------------
# 1) read VCF as genlight
# user requested: do not use gl.filter.callrate
# -----------------------------
snps <- dartR::gl.filter.monomorphs(
  vcfR::vcfR2genlight(
    vcfR::read.vcfR(vcf_file, verbose = FALSE)
  )
)

# -----------------------------
# 2) read Supplemental Table S1
# note:
# - the sheet contains all 307 sequenced samples
# - the pruned VCF contains the 108-individual population-genetics subset
# - we therefore reconstruct the 108-individual focal dataset first
# -----------------------------
s1_raw <- readxl::read_excel(s1_file)
names(s1_raw) <- make.names(names(s1_raw), unique = TRUE)

sample_col   <- names(s1_raw)[str_detect(names(s1_raw), regex("^Sample.?ID$", ignore_case = TRUE))][1]
common_col   <- names(s1_raw)[str_detect(names(s1_raw), regex("Species.?Common.?Name", ignore_case = TRUE))][1]
species_col  <- names(s1_raw)[str_detect(names(s1_raw), regex("Scientific.?Name", ignore_case = TRUE))][1]
location_col <- names(s1_raw)[str_detect(names(s1_raw), regex("^Location$", ignore_case = TRUE))][1]
more_col     <- names(s1_raw)[str_detect(names(s1_raw), regex("More.?Specific", ignore_case = TRUE))][1]
filt_col     <- names(s1_raw)[str_detect(names(s1_raw), regex("Filtered.?Read.?Count", ignore_case = TRUE))][1]

if (any(is.na(c(sample_col, common_col, species_col, location_col, more_col, filt_col)))) {
  stop("Could not identify all required columns in Supplemental Table S1.xlsx")
}

s1 <- s1_raw %>%
  transmute(
    sample_id = as.character(.data[[sample_col]]),
    common_name = as.character(.data[[common_col]]),
    species = as.character(.data[[species_col]]),
    location = as.character(.data[[location_col]]),
    more_info = as.character(.data[[more_col]]),
    filtered_reads = suppressWarnings(as.numeric(.data[[filt_col]]))
  ) %>%
  filter(!is.na(sample_id), sample_id != "") %>%
  mutate(
    common_name = str_squish(common_name),
    species = str_squish(species),
    location = str_squish(location)
  )

# -----------------------------
# 3) rebuild the focal 108-individual dataset used for the
# paper's population-genetic analyses:
# 27 Wall Canyon + 37 Warner + 44 Owens = 108
# Warner is filtered to >= 90,000 aligned/filtered reads
# -----------------------------
focal_108 <- s1 %>%
  filter(species %in% c("Catostomus new species",
                        "Catostomus warnerensis",
                        "Catostomus fumeiventris")) %>%
  filter(species != "Catostomus warnerensis" | filtered_reads >= 90000)

if (nrow(focal_108) != nInd(snps)) {
  stop(
    "Reconstructed focal metadata set has ", nrow(focal_108),
    " rows but genlight has ", nInd(snps),
    " individuals. Check VCF/S1 consistency."
  )
}

id_map <- tibble(
  ind_name = indNames(snps),
  row_order = seq_len(nInd(snps))
) %>%
  bind_cols(focal_108)

# -----------------------------
# 4) retain Owens sucker only
# -----------------------------
meta_owens <- id_map %>%
  filter(species == "Catostomus fumeiventris")

if (nrow(meta_owens) == 0) {
  stop("No Owens sucker individuals found in mapped metadata.")
}

snps_owens <- snps[meta_owens$row_order]
indNames(snps_owens) <- meta_owens$sample_id

# -----------------------------
# 5) normalize Owens population labels
# the spreadsheet uses Bishop, CA suffixes and "South Fork"
# whereas the paper figure uses:
# Lower Hot Creek, Lower Rock Creek, Lower Horton Creek,
# Horton Creek, South Fork Bishop Creek
# -----------------------------
normalize_owens_location <- function(x) {
  x <- as.character(x)
  x <- str_squish(x)

  case_when(
    str_detect(x, regex("^Lower Hot Creek", ignore_case = TRUE)) ~ "Lower Hot Creek",
    str_detect(x, regex("^Lower Rock Creek", ignore_case = TRUE)) ~ "Lower Rock Creek",
    str_detect(x, regex("^Lower Horton Creek", ignore_case = TRUE)) ~ "Lower Horton Creek",
    str_detect(x, regex("^Horton Creek", ignore_case = TRUE)) ~ "Horton Creek",
    str_detect(x, regex("^South Fork", ignore_case = TRUE)) ~ "South Fork Bishop Creek",
    TRUE ~ x
  )
}

meta_owens <- meta_owens %>%
  mutate(location = normalize_owens_location(location))

# quick sanity check against the expected five Owens localities
expected_owens_locs <- c(
  "Lower Hot Creek",
  "Lower Rock Creek",
  "Lower Horton Creek",
  "Horton Creek",
  "South Fork Bishop Creek"
)

unexpected_locs <- setdiff(unique(meta_owens$location), expected_owens_locs)
if (length(unexpected_locs) > 0) {
  stop("Unexpected Owens locations after normalization: ",
       paste(unexpected_locs, collapse = ", "))
}

# -----------------------------
# 6) hard-code site coordinates
# S1 does not provide numeric Owens coordinates; these are
# approximate site centroids inferred from Fig. 1B in the paper
# -----------------------------
site_key <- tibble(
  location = expected_owens_locs,
  lat = c(37.72, 37.56, 37.55, 37.47, 37.36),
  lon = c(-118.85, -118.73, -118.94, -119.06, -118.55)
) %>%
  left_join(
    meta_owens %>% count(location, name = "n_ind"),
    by = "location"
  ) %>%
  mutate(n_ind = ifelse(is.na(n_ind), 0L, n_ind)) %>%
  filter(n_ind > 0) %>%
  mutate(site_id = seq_len(n())) %>%
  select(site_id, location, n_ind, lat, lon)

missing_locs <- setdiff(unique(meta_owens$location), site_key$location)
if (length(missing_locs) > 0) {
  stop("Some Owens locations are not in the hard-coded site key: ",
       paste(missing_locs, collapse = ", "))
}

# attach coords to individual metadata
coord_df <- meta_owens %>%
  left_join(site_key %>% select(location, lat, lon), by = "location")

if (any(is.na(coord_df$lat) | is.na(coord_df$lon))) {
  bad <- coord_df %>%
    filter(is.na(lat) | is.na(lon)) %>%
    pull(sample_id)
  stop("Some Owens individuals did not receive coordinates: ",
       paste(bad, collapse = ", "))
}

# -----------------------------
# 7) assign populations
# -----------------------------
pop_lookup <- setNames(site_key$site_id, site_key$location)
pop(snps_owens) <- factor(pop_lookup[coord_df$location], levels = site_key$site_id)

snps_owens@other$latlon <- as.matrix(coord_df[, c("lat", "lon")])
rownames(snps_owens@other$latlon) <- coord_df$sample_id
colnames(snps_owens@other$latlon) <- c("lat", "lon")
snps_owens@other$sample_meta <- coord_df

# -----------------------------
# 8) pairwise FST
# -----------------------------
CAFU_1_fst <- as.matrix(dartR::gl.fst.pop(snps_owens))
storage.mode(CAFU_1_fst) <- "numeric"
CAFU_1_fst[is.na(CAFU_1_fst)] <- 0
CAFU_1_fst <- (CAFU_1_fst + t(CAFU_1_fst)) / 2
CAFU_1_fst[CAFU_1_fst < 0] <- 0
diag(CAFU_1_fst) <- 0

current_levels <- levels(pop(snps_owens))
CAFU_1_fst <- CAFU_1_fst[current_levels, current_levels, drop = FALSE]
rownames(CAFU_1_fst) <- colnames(CAFU_1_fst) <- as.character(site_key$site_id)

CAFU_1_coords <- site_key %>%
  select(site_id, lat, lon)

stopifnot(identical(as.character(CAFU_1_coords$site_id), rownames(CAFU_1_fst)))

# -----------------------------
# 9) quick PCA for QC
# the paper shows Owens separating strongly from the other
# lineages, and internally Lower Hot Creek vs. other Owens
# localities are also distinct
# -----------------------------
cafu_pca <- glPca(snps_owens, nf = min(20, nInd(snps_owens) - 1))
pca_df <- data.frame(
  PC1 = cafu_pca$scores[, 1],
  PC2 = cafu_pca$scores[, 2],
  location = coord_df$location
)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 2) +
  labs(
    x = "PC1",
    y = "PC2",
    title = "CAFU-1 PCA"
  ) +
  theme_bw()

print(p_pca)

# -----------------------------
# 10) map plot
# include US and Canada and zoom to point extent
# -----------------------------
world_df <- ggplot2::map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

xpad <- max(0.5, diff(range(site_key$lon)) * 0.25)
ypad <- max(0.5, diff(range(site_key$lat)) * 0.25)

p_map <- ggplot() +
  geom_polygon(
    data = world_df,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = site_key,
    aes(x = lon, y = lat),
    size = 2.8
  ) +
  geom_text(
    data = site_key,
    aes(x = lon, y = lat, label = paste0(site_id, ". ", location)),
    nudge_y = ypad * 0.04,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(min(site_key$lon) - xpad, max(site_key$lon) + xpad),
    ylim = c(min(site_key$lat) - ypad, max(site_key$lat) + ypad),
    expand = FALSE
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "CAFU-1 sampling sites"
  ) +
  theme_bw()

print(p_map)

# -----------------------------
# 11) IBD plot
# -----------------------------
geo_dist <- geosphere::distm(
  x = as.matrix(CAFU_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist[upper.tri(geo_dist)],
  fst = CAFU_1_fst[upper.tri(CAFU_1_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Geographic distance (km)",
    y = "Pairwise FST",
    title = "CAFU-1 IBD plot"
  ) +
  theme_bw()

print(p_ibd)

# -----------------------------
# 12) save objects
# -----------------------------

# rename "site_id" as "site"
colnames(CAFU_1_coords)[colnames(CAFU_1_coords) == "site_id"] <- "site"

# save to path
save(
  CAFU_1_fst,
  CAFU_1_coords,
  file = file.path(data_dir, "CAFU-1.RData")
)

cat("Saved: ", file.path(data_dir, "CAFU-1.RData"), "\n", sep = "")
