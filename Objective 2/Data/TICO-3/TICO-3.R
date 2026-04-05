# ============================================================
# TICO-3 | Loach minnow
# Tiaroga cobitis (formerly Rhinichthys cobitis)
# Pilger et al. 2022, The American Naturalist
# Rebuilt from Tc09-10_dat.gen BIMr input file
# ============================================================

suppressPackageStartupMessages({
  library(hierfstat)
  library(ggplot2)
  library(dplyr)
  library(geosphere)
})

study_code <- "TICO-3"
study_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/TICO-3"
out_dir <- file.path(study_dir, "data")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

# -----------------------------
# 1) read full 16-site metadata
# -----------------------------
metadata_file <- "16site_metadata.csv"

if (!file.exists(metadata_file)) {
  stop("Metadata file not found in TICO-3 directory.")
}

site_metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

required_cols <- c("study_site_no", "site_code", "site_name", "stream", "lat", "lon", "confidence")
stopifnot(all(required_cols %in% names(site_metadata)))

site_metadata$study_site_no <- as.integer(site_metadata$study_site_no)
site_metadata <- site_metadata[, required_cols]

# -----------------------------
# 1b) standardize site metadata to the Figure 1 numbering
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
# 2) read genotype file (BIMr format)
# Tc09-10_dat.gen combines 2009 + 2010 genotypes
# for loach minnow, matching the supplement note
# -----------------------------
geno_file <- "Tc09-10_dat.gen"
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
# 3) recover site-level populations from individual ID prefixes
#
# Supplement Table S.4 indicates loach minnow samples at:
# MB (3), BA (4), RS (5), GV (9)
#
# BIMr grouped these into two migration units: 3-5 and 9.
# For Objective 2 we rebuild the four true site-level pops.
#
# Expected raw prefixes are usually:
# MBTc, BATc, RSTc, GVTc
# The script prints unique_prefixes and stops with the exact
# unmapped prefixes if any differ from expectation.
# -----------------------------
unique_prefixes <- sort(unique(geno_df$raw_prefix))
print(unique_prefixes)

prefix_to_site_code <- c(
  "MBTc" = "MB",
  "BATc" = "BA",
  "RSTc" = "RS",
  "GVTc" = "GV"
)

geno_df$site_code <- unname(prefix_to_site_code[geno_df$raw_prefix])

if (any(is.na(geno_df$site_code))) {
  bad_prefixes <- sort(unique(geno_df$raw_prefix[is.na(geno_df$site_code)]))
  stop(paste("Unmapped prefixes detected:", paste(bad_prefixes, collapse = ", ")))
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

expected_codes <- c("MB", "BA", "RS", "GV")
stopifnot(setequal(pop_lookup$site_code, expected_codes))

site_code_to_study_no <- c(
  "MB" = 3,
  "BA" = 4,
  "RS" = 5,
  "GV" = 9
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

print(pop_lookup %>%
  select(site, site_code, site_code_meta, study_site_no, site_name,
         n_ind, bimr_clusters, raw_prefix))

# -----------------------------
# 4) build hierfstat input
# each locus is a 4-digit diploid genotype
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
TICO_3_fst <- hierfstat::pairwise.WCfst(hf_df)

TICO_3_fst <- as.matrix(TICO_3_fst)
diag(TICO_3_fst) <- 0
TICO_3_fst[is.na(TICO_3_fst)] <- 0
TICO_3_fst[TICO_3_fst < 0] <- 0

rownames(TICO_3_fst) <- as.character(pop_lookup$site)
colnames(TICO_3_fst) <- as.character(pop_lookup$site)

stopifnot(nrow(TICO_3_fst) == 4)
stopifnot(isTRUE(all.equal(TICO_3_fst, t(TICO_3_fst))))

# -----------------------------
# 6) coords df in the SAME ORDER as the FST matrix
# Objective 2 format
# -----------------------------
TICO_3_coords <- pop_lookup %>%
  transmute(
    site = site,
    lat = lat,
    lon = lon
  )

stopifnot(identical(as.character(TICO_3_coords$site), rownames(TICO_3_fst)))

TICO_3_site_lookup <- pop_lookup %>%
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
# -----------------------------
world_df <- ggplot2::map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

plot_df <- TICO_3_site_lookup %>%
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
    title = "TICO-3 sampling sites",
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
  x = as.matrix(TICO_3_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  fst = TICO_3_fst[upper.tri(TICO_3_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "TICO-3 IBD plot",
    x = "Geographic distance (km)",
    y = "Pairwise FST"
  ) +
  theme_bw()

print(p_ibd)

# -----------------------------
# 9) save outputs
# -----------------------------
save(
  TICO_3_fst,
  TICO_3_coords,
  file = file.path(out_dir, "TICO-3.RData")
)

# optional debug export
# write.csv(TICO_3_site_lookup, file.path(out_dir, "TICO-3_site_lookup_debug.csv"), row.names = FALSE)
