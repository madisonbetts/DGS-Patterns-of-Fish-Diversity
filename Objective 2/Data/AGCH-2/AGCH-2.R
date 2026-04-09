# ============================================================
# AGCH-2 | Longfin dace
# Agosia chrysogaster
# Pilger et al. 2022, The American Naturalist
# Updated to match the actual Ac10_dat.gen file structure
# ============================================================

suppressPackageStartupMessages({
  library(hierfstat)
  library(ggplot2)
  library(dplyr)
  library(geosphere)
})

study_code <- "AGCH-2"
study_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/AGCH-2"
out_dir <- file.path(study_dir, "data")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

# -----------------------------
# 1) read full 16-site metadata
# use the user-made metadata csv if present
# -----------------------------
metadata_candidates <- c(
  "16site_metadata.csv",
  file.path("data", "16site_metadata.csv")
)

metadata_file <- metadata_candidates[file.exists(metadata_candidates)][1]

if (length(metadata_file) == 0 || is.na(metadata_file)) {
  stop("Could not find 16site_metadata.csv in the study directory.")
}

site_metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

required_cols <- c("study_site_no", "site_code", "site_name", "stream", "lat", "lon", "confidence")
stopifnot(all(required_cols %in% names(site_metadata)))

site_metadata$study_site_no <- as.integer(site_metadata$study_site_no)

# keep only the columns we actually need
site_metadata <- site_metadata[, required_cols]

# -----------------------------
# 2) read the actual genotype file
# IMPORTANT:
# this is not standard GenePop format
# it is a custom tab-delimited file with a short header:
# [individuals]=
# [populations]=
# [loci]=
# [alleles]=
# [genotypes]=
# then rows:
# individual_id_comma <tab> pop_id <tab> locus1 ... locus10
# -----------------------------
geno_file <- file.path("Agochr", "Ac10_dat.gen")
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

# 1 column for raw individual id, 1 for pop_id, then n_loci_expected loci
stopifnot(all(lengths(split_tab) == (2 + n_loci_expected)))

geno_df <- data.frame(
  individual_raw = vapply(split_tab, `[`, character(1), 1),
  pop_id = as.integer(vapply(split_tab, `[`, character(1), 2)),
  stringsAsFactors = FALSE
)

locus_mat <- do.call(rbind, lapply(split_tab, function(x) x[3:(2 + n_loci_expected)]))
colnames(locus_mat) <- paste0("L", seq_len(n_loci_expected))
geno_df <- cbind(geno_df, as.data.frame(locus_mat, stringsAsFactors = FALSE))

geno_df$individual <- trimws(sub(",$", "", geno_df$individual_raw))
geno_df$individual_prefix <- sub("Ac.*$", "", geno_df$individual)

stopifnot(nrow(geno_df) == n_ind_expected)
stopifnot(length(unique(geno_df$pop_id)) == n_pop_expected)

# -----------------------------
# 3) recover the actual site abbreviations in FILE ORDER
#
# The individual IDs in Ac10_dat.gen do not always use the same
# abbreviation strings as Table S.4 / your site metadata:
#   BC    -> BCA
#   SC    -> SAP
#   TC    -> TURC
#   BlueC -> BLUC
#
# So we build the file-order pop lookup from pop_id and then
# explicitly translate those prefixes to the final site codes.
# -----------------------------
pop_lookup <- geno_df %>%
  group_by(pop_id) %>%
  summarise(
    first_individual = first(individual),
    raw_prefix = first(individual_prefix),
    n_ind = n(),
    .groups = "drop"
  ) %>%
  arrange(pop_id) %>%
  mutate(site = row_number())

prefix_to_site_code <- c(
  "WF" = "WF",
  "MF" = "MF",
  "HB" = "HB",
  "BC" = "BCA",
  "GV" = "GV",
  "SC" = "SAP",
  "TC" = "TURC",
  "GF" = "GF",
  "RS" = "RS",
  "BA" = "BA",
  "MB" = "MB",
  "BlueC" = "BLUC",
  "SSC" = "SSC"
)

pop_lookup$site_code <- unname(prefix_to_site_code[pop_lookup$raw_prefix])

if (any(is.na(pop_lookup$site_code))) {
  stop("At least one raw individual prefix could not be translated to a site_code.")
}

# expected AGCH order from the ACTUAL file
expected_codes <- c("WF","MF","HB","BCA","GV","SAP","TURC","GF","RS","BA","MB","BLUC","SSC")
stopifnot(identical(pop_lookup$site_code, expected_codes))

# join to metadata
pop_lookup <- pop_lookup %>%
  left_join(site_metadata, by = "site_code")

stopifnot(!any(is.na(pop_lookup$study_site_no)))
stopifnot(!any(is.na(pop_lookup$lat)))
stopifnot(!any(is.na(pop_lookup$lon)))

# print for debugging
print(pop_lookup)

# -----------------------------
# 4) build hierfstat input
# each locus is coded as a 4-digit diploid genotype
# 0000 = missing
# -----------------------------
hf_df <- geno_df %>%
  left_join(pop_lookup %>% select(pop_id, site), by = "pop_id") %>%
  select(site, starts_with("L"))

for (j in 2:ncol(hf_df)) {
  hf_df[[j]] <- as.integer(hf_df[[j]])
  hf_df[[j]][hf_df[[j]] == 0] <- NA_integer_
}

# -----------------------------
# 5) pairwise FST
# -----------------------------
AGCH_2_fst <- hierfstat::pairwise.WCfst(hf_df)

AGCH_2_fst <- as.matrix(AGCH_2_fst)
diag(AGCH_2_fst) <- 0
AGCH_2_fst[is.na(AGCH_2_fst)] <- 0
AGCH_2_fst[AGCH_2_fst < 0] <- 0

rownames(AGCH_2_fst) <- as.character(pop_lookup$site)
colnames(AGCH_2_fst) <- as.character(pop_lookup$site)

stopifnot(nrow(AGCH_2_fst) == 13)
stopifnot(isTRUE(all.equal(AGCH_2_fst, t(AGCH_2_fst))))

# -----------------------------
# 6) coords df in the SAME ORDER as the FST matrix
# Objective 2 format
# -----------------------------
AGCH_2_coords <- pop_lookup %>%
  transmute(
    site = site,
    lat = lat,
    lon = lon
  )

stopifnot(identical(as.character(AGCH_2_coords$site), rownames(AGCH_2_fst)))

# optional richer lookup for interpretation
AGCH_2_site_lookup <- pop_lookup %>%
  select(
    site,
    pop_id,
    raw_prefix,
    site_code,
    study_site_no,
    site_name,
    stream,
    lat,
    lon,
    confidence,
    n_ind,
    first_individual
  )

# -----------------------------
# 7) map plot
# labels = file-order site number + site code
# -----------------------------
world_df <- ggplot2::map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

plot_df <- AGCH_2_site_lookup %>%
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
    title = "AGCH-2 sampling sites",
    subtitle = "Labels show file-order site number and site abbreviation",
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
  x = as.matrix(AGCH_2_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  fst = AGCH_2_fst[upper.tri(AGCH_2_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "AGCH-2 IBD plot",
    x = "Geographic distance (km)",
    y = "Pairwise FST"
  ) +
  theme_bw()

print(p_ibd)

# -----------------------------
# 9) save outputs
# save matrix + coords to RData
# -----------------------------
save(
  AGCH_2_fst,
  AGCH_2_coords,
  file = file.path(out_dir, "AGCH-2.RData")
)

# optional debug file
#write.csv(AGCH_2_site_lookup, file.path(out_dir, "AGCH-2_site_lookup_debug.csv"), row.names = FALSE)
