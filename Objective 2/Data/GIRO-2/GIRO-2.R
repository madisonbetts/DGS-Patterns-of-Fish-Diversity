# ============================================================
# GIRO-2 | Roundtail chub
# Gila robusta
# Pilger et al. 2022, The American Naturalist
# Rebuilt from the actual Gn10_dat.gen BIMr input file
# ============================================================

suppressPackageStartupMessages({
  library(hierfstat)
  library(ggplot2)
  library(dplyr)
  library(geosphere)
})

study_code <- "GIRO-2"
study_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/GIRO-2"
out_dir <- file.path(study_dir, "data")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

# -----------------------------
# 1) read full 16-site metadata
# -----------------------------
metadata_file <- "16site_metadata.csv"

if (!file.exists(metadata_file)) {
  stop("Metadata file not found in GIRO-2 directory.")
}

site_metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

required_cols <- c("study_site_no", "site_code", "site_name", "stream", "lat", "lon", "confidence")
stopifnot(all(required_cols %in% names(site_metadata)))

site_metadata$study_site_no <- as.integer(site_metadata$study_site_no)
site_metadata <- site_metadata[, required_cols]

# -----------------------------
# 1b) standardize site metadata to the Figure 1 / README numbering
# keep uploaded XYs but standardize names by study_site_no
# -----------------------------
site_ref <- data.frame(
  study_site_no = 1:16,
  site_code = c("SSC", "BLUC", "MB", "BA", "RS", "GF", "TURC", "SAP",
                "GV", "LEF", "LC", "HB", "WF", "MF", "BCA", "UEF"),
  site_name = c("Sunset Canal", "Blue Creek", "Middle Box", "Bird Area",
                "Riverside", "Gila Farm", "Turkey Creek", "Sapillo Creek",
                "Grapevine", "Lower East Fork", "Little Creek",
                "Removal Reach", "West Fork", "Middle Fork",
                "Black Canyon", "Upper East Fork"),
  stream = c("Gila River", "Blue Creek", "Gila River", "Gila River",
             "Gila River", "Gila River", "Turkey Creek", "Sapillo Creek",
             "East Fork Gila River", "East Fork Gila River", "Little Creek",
             "Gila River", "West Fork Gila River", "Middle Fork Gila River",
             "Black Canyon", "East Fork Gila River"),
  stringsAsFactors = FALSE
)

site_metadata <- site_metadata %>%
  select(study_site_no, lat, lon, confidence) %>%
  right_join(site_ref, by = "study_site_no") %>%
  select(study_site_no, site_code, site_name, stream, lat, lon, confidence) %>%
  arrange(study_site_no)

print(site_metadata)

# -----------------------------
# 2) read the actual genotype file
# IMPORTANT:
# this is not standard GenePop format
# it is a BIMr tab-delimited file with a short header:
# [individuals]=
# [populations]=
# [loci]=
# [alleles]=
# [genotypes]=
# then rows:
# individual_id_comma <tab> pop_id <tab> locus1 ... locus10
#
# For G. robusta, the BIMr k4 group count matches the four
# actual sampled sites, but we still reconstruct populations
# from the individual ID prefixes so the workflow is explicit.
# -----------------------------
geno_file <- "Gn10_dat.gen"
stopifnot(file.exists(geno_file))

raw_lines <- readLines(geno_file, warn = FALSE)
raw_lines <- gsub("\r", "", raw_lines)

header_vals <- function(tag) {
  i <- grep(paste0("^\\[", tag, "\\]="), raw_lines)
  stopifnot(length(i) == 1)
  trimws(sub(paste0("^\\[", tag, "\\]="), "", raw_lines[i]))
}

n_ind_expected  <- as.integer(header_vals("individuals"))
n_pop_expected  <- as.integer(header_vals("populations"))
n_loci_expected <- as.integer(header_vals("loci"))

geno_start <- grep("^\\[genotypes\\]=", raw_lines)
stopifnot(length(geno_start) == 1)

geno_lines <- raw_lines[(geno_start + 1):length(raw_lines)]
geno_lines <- geno_lines[nzchar(trimws(geno_lines))]

split_tab <- strsplit(geno_lines, "\t", fixed = TRUE)
stopifnot(all(lengths(split_tab) == (2 + n_loci_expected)))

geno_df <- data.frame(
  individual_raw = vapply(split_tab, `[`, character(1), 1),
  bimr_cluster = as.integer(vapply(split_tab, `[`, character(1), 2)),
  stringsAsFactors = FALSE
)

locus_mat <- do.call(rbind, lapply(split_tab, function(x) x[3:(2 + n_loci_expected)]))
colnames(locus_mat) <- paste0("L", seq_len(n_loci_expected))
geno_df <- cbind(geno_df, as.data.frame(locus_mat, stringsAsFactors = FALSE))

geno_df$individual <- trimws(sub(",$", "", geno_df$individual_raw))
geno_df$raw_prefix <- sub("[0-9].*$", "", geno_df$individual)

stopifnot(nrow(geno_df) == n_ind_expected)
stopifnot(length(unique(geno_df$bimr_cluster)) == n_pop_expected)

# -----------------------------
# 3) recover true site identities from individual ID prefixes
#
# Observed prefixes in Gn10_dat.gen:
# WFGn   -> WF   (West Fork; site 13)
# MFGn   -> MF   (Middle Fork; site 14)
# UEFGn  -> UEF  (Upper East Fork; site 16)
# GVGn   -> GV   (Grapevine; site 9)
#
# This matches the four G. robusta sample locations reported
# in the supplement / figure labels: 9, 13, 14, 16.
# -----------------------------
prefix_to_site_code <- c(
  "WFGn"  = "WF",
  "MFGn"  = "MF",
  "UEFGn" = "UEF",
  "GVGn"  = "GV"
)

geno_df$site_code <- unname(prefix_to_site_code[geno_df$raw_prefix])
if (any(is.na(geno_df$site_code))) {
  stop("At least one individual prefix could not be translated to a site_code.")
}

pop_lookup <- geno_df %>%
  group_by(site_code) %>%
  summarise(
    raw_prefix = first(raw_prefix),
    n_ind = n(),
    bimr_clusters = paste(sort(unique(bimr_cluster)), collapse = ","),
    first_individual = first(individual),
    .groups = "drop"
  )

expected_codes <- c("GV", "WF", "MF", "UEF")
stopifnot(setequal(pop_lookup$site_code, expected_codes))

site_code_to_study_no <- c(
  "GV"  = 9,
  "WF"  = 13,
  "MF"  = 14,
  "UEF" = 16
)

pop_lookup <- pop_lookup %>%
  mutate(
    study_site_no = unname(site_code_to_study_no[site_code]),
    site = match(site_code, expected_codes)
  ) %>%
  left_join(
    site_metadata %>%
      select(study_site_no, site_code_meta = site_code, site_name, stream, lat, lon, confidence),
    by = "study_site_no"
  ) %>%
  arrange(site)

stopifnot(!any(is.na(pop_lookup$study_site_no)))
stopifnot(!any(is.na(pop_lookup$lat)))
stopifnot(!any(is.na(pop_lookup$lon)))

print(pop_lookup %>% select(site, site_code, site_code_meta, study_site_no, site_name, n_ind, bimr_clusters, raw_prefix))

# -----------------------------
# 4) build hierfstat input
# each locus is coded as a 4-digit diploid genotype
# 0000 = missing
# -----------------------------
hf_df <- geno_df %>%
  left_join(pop_lookup %>% select(site_code, site), by = "site_code") %>%
  select(site, starts_with("L"))

for (j in 2:ncol(hf_df)) {
  hf_df[[j]] <- as.integer(hf_df[[j]])
  hf_df[[j]][hf_df[[j]] == 0] <- NA_integer_
}

# -----------------------------
# 5) pairwise FST
# -----------------------------
GIRO_2_fst <- hierfstat::pairwise.WCfst(hf_df)

GIRO_2_fst <- as.matrix(GIRO_2_fst)
diag(GIRO_2_fst) <- 0
GIRO_2_fst[is.na(GIRO_2_fst)] <- 0
GIRO_2_fst[GIRO_2_fst < 0] <- 0

rownames(GIRO_2_fst) <- as.character(pop_lookup$site)
colnames(GIRO_2_fst) <- as.character(pop_lookup$site)

stopifnot(nrow(GIRO_2_fst) == 4)
stopifnot(isTRUE(all.equal(GIRO_2_fst, t(GIRO_2_fst))))

# -----------------------------
# 6) coords df in the SAME ORDER as the FST matrix
# Objective 2 format
# -----------------------------
GIRO_2_coords <- pop_lookup %>%
  transmute(
    site = site,
    lat = lat,
    lon = lon
  )

stopifnot(identical(as.character(GIRO_2_coords$site), rownames(GIRO_2_fst)))

GIRO_2_site_lookup <- pop_lookup %>%
  select(
    site,
    site_code,
    study_site_no,
    site_name,
    stream,
    lat,
    lon,
    confidence,
    n_ind,
    bimr_clusters,
    raw_prefix,
    first_individual
  )

# -----------------------------
# 7) map plot
# labels = Objective-2 site number + site code
# map_data is from ggplot2, not maps
# include USA and Canada, then zoom to point extent
# -----------------------------
world_df <- ggplot2::map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

plot_df <- GIRO_2_site_lookup %>%
  mutate(label = paste0(site, " ", site_code))

x_pad <- max(0.3, diff(range(plot_df$lon)) * 0.10)
y_pad <- max(0.3, diff(range(plot_df$lat)) * 0.10)

p_map <- ggplot() +
  geom_polygon(
    data = world_df,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = plot_df,
    aes(x = lon, y = lat),
    size = 2
  ) +
  geom_text(
    data = plot_df,
    aes(x = lon, y = lat, label = label),
    nudge_y = 0.03,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(min(plot_df$lon) - x_pad, max(plot_df$lon) + x_pad),
    ylim = c(min(plot_df$lat) - y_pad, max(plot_df$lat) + y_pad)
  ) +
  labs(
    title = "GIRO-2 sampling sites",
    subtitle = "Labels show Objective-2 site number and site abbreviation",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw()

print(p_map)

# -----------------------------
# 8) IBD plot
# straight-line distance for QC only
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(GIRO_2_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  fst = GIRO_2_fst[upper.tri(GIRO_2_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "GIRO-2 IBD plot",
    x = "Geographic distance (km)",
    y = "Pairwise FST"
  ) +
  theme_bw()

print(p_ibd)

# -----------------------------
# 9) save outputs
# -----------------------------
save(
  GIRO_2_fst,
  GIRO_2_coords,
  file = file.path(out_dir, "GIRO-2.RData")
)

# optional debug export
# write.csv(GIRO_2_site_lookup, file.path(out_dir, "GIRO-2_site_lookup_debug.csv"), row.names = FALSE)
