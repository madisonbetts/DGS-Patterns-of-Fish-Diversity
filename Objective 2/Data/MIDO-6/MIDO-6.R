# ============================================================
# MIDO-6 | Smallmouth bass
# Micropterus dolomieu
# Pilger et al. 2022, The American Naturalist
# Rebuilt from Md10_k4_dat.gen BIMr input file
#
# NOTE:
# This study folder did not include a site metadata csv, so this script
# uses an internal site metadata table tied to Figure 1 numbering.
# Coordinates below are best-available georeferenced approximations for
# the Figure 1 study sites used in this dataset.
# ============================================================

suppressPackageStartupMessages({
  library(hierfstat)
  library(ggplot2)
  library(dplyr)
  library(geosphere)
})

study_code <- "MIDO-6"
study_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/MIDO-6"
out_dir <- file.path(study_dir, "data")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

# -----------------------------
# 1) site metadata
# linked to Figure 1 numbering in Pilger et al. 2022
# coordinates are best-available approximations
# -----------------------------
site_metadata <- data.frame(
  study_site_no = 1:16,
  site_code = c("SSC","BLUC","MB","BA","RS","GF","TURC","SAP",
                "GV","LEF","LC","HB","WF","MF","BCA","UEF"),
  site_name = c("Sunset Canal","Blue Creek","Middle Box","Bird Area",
                "Riverside","Gila Farm","Turkey Creek","Sapillo Creek",
                "Grapevine","Lower East Fork","Little Creek","Removal Reach",
                "West Fork","Middle Fork","Black Canyon","Upper East Fork"),
  stream = c("Gila River","Blue Creek","Gila River","Gila River",
             "Gila River","Gila River","Turkey Creek","Sapillo Creek",
             "East Fork Gila River","East Fork Gila River","Little Creek",
             "Gila River","West Fork Gila River","Middle Fork Gila River",
             "Black Canyon","East Fork Gila River"),
  lat = c(32.6491908, 32.8763110, 32.7151752, 32.8369826,
          32.9503780, 33.0468810, 33.0629280, 33.0377740,
          33.1777370, 33.1784200, 33.1920990, 33.2078030,
          33.2289000, 33.2266000, 33.1786840, 33.3020230),
  lon = c(-108.9297125, -108.8375535, -108.7016474, -108.6042210,
          -108.6068090, -108.5286470, -108.5046790, -108.2323510,
          -108.2089570, -108.1973350, -108.2302120, -108.2268550,
          -108.2629000, -108.2415000, -108.0559580, -108.1245170),
  confidence = c("moderate","high","moderate","moderate",
                 "moderate","moderate","moderate","moderate",
                 "moderate","moderate","moderate","high",
                 "high","high","high","moderate"),
  stringsAsFactors = FALSE
)

print(site_metadata)

# -----------------------------
# 2) read genotype file (BIMr format)
# -----------------------------
geno_file <- "Md10_k4_dat.gen"
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
# Actual prefixes in Md10_k4_dat.gen:
# RSMd   -> Riverside (5)
# GFMd   -> Gila Farm (6)
# GVMd   -> Grapevine (9)
# LEFMd  -> Lower East Fork (10)
# HBMd   -> Removal Reach (12)
# MFMd   -> Middle Fork (14)
#
# BIMr k4 combines some of these sites for migration modeling,
# but for Objective 2 we rebuild true site-level populations.
# Figure S.9 / Table S.4 indicate Md sites include 5, 6, 9, 10, 12, 14
# among the locations with sufficient sample sizes used here.
# -----------------------------
unique_prefixes <- sort(unique(geno_df$raw_prefix))
print(unique_prefixes)

prefix_to_site_code <- c(
  "RSMd"  = "RS",
  "GFMd"  = "GF",
  "GVMd"  = "GV",
  "LEFMd" = "LEF",
  "HBMd"  = "HB",
  "MFMd"  = "MF"
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

expected_codes <- c("RS", "GF", "GV", "LEF", "HB", "MF")
stopifnot(setequal(pop_lookup$site_code, expected_codes))

site_code_to_study_no <- c(
  "RS"  = 5,
  "GF"  = 6,
  "GV"  = 9,
  "LEF" = 10,
  "HB"  = 12,
  "MF"  = 14
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
MIDO_6_fst <- hierfstat::pairwise.WCfst(hf_df)

MIDO_6_fst <- as.matrix(MIDO_6_fst)
diag(MIDO_6_fst) <- 0
MIDO_6_fst[is.na(MIDO_6_fst)] <- 0
MIDO_6_fst[MIDO_6_fst < 0] <- 0

rownames(MIDO_6_fst) <- as.character(pop_lookup$site)
colnames(MIDO_6_fst) <- as.character(pop_lookup$site)

stopifnot(nrow(MIDO_6_fst) == 6)
stopifnot(isTRUE(all.equal(MIDO_6_fst, t(MIDO_6_fst))))

# -----------------------------
# 6) coords df in the SAME ORDER as the FST matrix
# -----------------------------
MIDO_6_coords <- pop_lookup %>%
  transmute(
    site = site,
    lat = lat,
    lon = lon
  )

stopifnot(identical(as.character(MIDO_6_coords$site), rownames(MIDO_6_fst)))

MIDO_6_site_lookup <- pop_lookup %>%
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
# Figure 1 site linkage used for XYs
# include USA and Canada, zoom to points
# -----------------------------
world_df <- ggplot2::map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

plot_df <- MIDO_6_site_lookup %>%
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
    title = "MIDO-6 sampling sites",
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
  x = as.matrix(MIDO_6_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  fst = MIDO_6_fst[upper.tri(MIDO_6_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "MIDO-6 IBD plot",
    x = "Geographic distance (km)",
    y = "Pairwise FST"
  ) +
  theme_bw()

print(p_ibd)

# -----------------------------
# 9) save outputs
# -----------------------------
save(
  MIDO_6_fst,
  MIDO_6_coords,
  file = file.path(out_dir, "MIDO-6.RData")
)

# optional debug export
# write.csv(MIDO_6_site_lookup, file.path(out_dir, "MIDO-6_site_lookup_debug.csv"), row.names = FALSE)
