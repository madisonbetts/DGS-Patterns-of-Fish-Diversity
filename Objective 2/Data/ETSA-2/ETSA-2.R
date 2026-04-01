# ========================================
# ETSA-2
# Etheostoma sagitta
# Kentucky Arrow Darter
#
# Primary extraction target:
# Culley et al. 2025 dataset in:
#   ETSA-2/Culley_etal_2025_KAD_LE/
#
# Supporting raw genotype files:
#   LM_geneticdata.gen
#   LM_geneticdata.txt
#
# Supporting older / parallel river_data files:
#   KAD_Watson_GenepopFile_12popRearranged.txt
#   SiteCoords.csv
#   riverdist_correctsitenumbers.csv
#
# This script does the usual Objective 2 pipeline:
#   - read exact site coordinates from the Excel "Site Information" sheet
#   - read the published pairwise FST matrix from the Excel "Pairwise FST" sheet
#   - standardize the objects to numeric site labels 1:n
#   - save ETSA_2_coords and ETSA_2_fst to data/ETSA-2.RData
#   - make a map and Euclidean IBD plot
#
# It also includes lightweight parsing of the .gen / Genepop files so you can
# inspect loci, pop counts, and merge potential across studies.
#
# IMPORTANT MERGE NOTE
# The Culley and Watson datasets use the same 11 microsatellite loci by name:
#   EosC2, EosC3, EosC6, EosC117, EosD10, EosD11, EosD107,
#   EosD116, EosD131, Esc26b, Cv12
#
# However, the allele coding does not look plug-and-play identical across files.
# Example:
#   - Culley LM_geneticdata.gen has EosC3 alleles up to 03 and EosD11 up to 46
#   - Watson Genepop file has EosC3 only up to 02 and EosD11 up to 23
#
# So you should NOT blindly row-bind the two genotype files and recompute FST.
# A combined raw-genotype dataset is only defensible if you can verify that:
#   1) allele binning was done on the same size standard / scoring system
#   2) allele states are directly comparable across studies
#   3) overlapping streams / populations score concordantly
#
# Pragmatically:
#   - combining the published FST matrices is NOT valid because not all among-study
#     pairwise comparisons were estimated
#   - combining the raw Genepop files MAY be possible, but only after explicit
#     harmonization / validation of allele bins
# ========================================

# -----------------------------
# 0) setup
# -----------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(geosphere)
library(maps)

study_code <- "ETSA-2"

base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETSA-2"
culley_dir <- file.path(base_dir, "Culley_etal_2025_KAD_LE")
river_dir  <- file.path(base_dir, "river_data")
data_dir   <- file.path(base_dir, "data")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

excel_file      <- file.path(culley_dir, "Culley_etal_2025_KAD_LE_Sites_Distance_Fst.xlsx")
lm_gen_file     <- file.path(culley_dir, "LM_geneticdata.gen")
lm_txt_file     <- file.path(culley_dir, "LM_geneticdata.txt")
watson_gen_file <- file.path(river_dir, "KAD_Watson_GenepopFile_12popRearranged.txt")
sitecoords_file <- file.path(river_dir, "SiteCoords.csv")
riverdist_file  <- file.path(river_dir, "riverdist_correctsitenumbers.csv")

stopifnot(file.exists(excel_file))

# -----------------------------
# 1) helper: parse Genepop / .gen file
# useful for checking locus names, n pops, and merge potential
# -----------------------------
read_genepop_lines <- function(path) {
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x[x != ""]
}

parse_genepop <- function(path) {
  x <- read_genepop_lines(path)

  title <- x[1]
  i <- 2
  loci <- character(0)

  while (i <= length(x) && toupper(x[i]) != "POP") {
    loci <- c(loci, gsub(",$", "", x[i]))
    i <- i + 1
  }

  if (i > length(x)) stop("No POP marker found in Genepop file: ", path)

  pops <- list()
  current <- character(0)

  for (j in (i + 1):length(x)) {
    if (toupper(x[j]) == "POP") {
      pops[[length(pops) + 1]] <- current
      current <- character(0)
    } else {
      current <- c(current, x[j])
    }
  }
  pops[[length(pops) + 1]] <- current

  list(
    title = title,
    loci = loci,
    pops = pops
  )
}

genepop_alleles <- function(path) {
  gp <- parse_genepop(path)
  out <- setNames(vector("list", length(gp$loci)), gp$loci)

  for (line in unlist(gp$pops)) {
    bits <- strsplit(line, ",")[[1]]
    geno <- trimws(bits[2])
    gvec <- strsplit(geno, "[[:space:]]+")[[1]]

    for (k in seq_along(gp$loci)) {
      g <- gvec[k]
      if (is.na(g) || g == "0000") next
      half <- nchar(g) / 2
      a1 <- substr(g, 1, half)
      a2 <- substr(g, half + 1, nchar(g))
      out[[k]] <- sort(unique(c(out[[k]], a1, a2)))
    }
  }

  out
}

# inspect available genotype files if present
if (file.exists(lm_gen_file)) {
  lm_gp <- parse_genepop(lm_gen_file)
  lm_alleles <- genepop_alleles(lm_gen_file)
} else {
  lm_gp <- NULL
  lm_alleles <- NULL
}

if (file.exists(watson_gen_file)) {
  watson_gp <- parse_genepop(watson_gen_file)
  watson_alleles <- genepop_alleles(watson_gen_file)
} else {
  watson_gp <- NULL
  watson_alleles <- NULL
}

# -----------------------------
# 2) read Excel sheets
# -----------------------------
sheet_names <- excel_sheets(excel_file)
stopifnot(all(c("Site Information", "Distance Matrix (meters)", "Pairwise FST") %in% sheet_names))

site_info <- read_excel(excel_file, sheet = "Site Information")
dist_raw  <- read_excel(excel_file, sheet = "Distance Matrix (meters)")
fst_raw   <- read_excel(excel_file, sheet = "Pairwise FST")

# -----------------------------
# 3) clean site information sheet
# expected columns:
#   Site Name, id, lat, long
# -----------------------------
names(site_info) <- trimws(names(site_info))
site_name_col <- names(site_info)[tolower(names(site_info)) %in% c("site name", "site_name")]
id_col        <- names(site_info)[tolower(names(site_info)) %in% c("id", "site", "site id", "site_id")]
lat_col       <- names(site_info)[grepl("^lat$|latitude", tolower(names(site_info)))]
lon_col       <- names(site_info)[grepl("^long$|^lon$|longitude", tolower(names(site_info)))]

stopifnot(length(site_name_col) >= 1, length(id_col) >= 1, length(lat_col) >= 1, length(lon_col) >= 1)

site_name_col <- site_name_col[1]
id_col        <- id_col[1]
lat_col       <- lat_col[1]
lon_col       <- lon_col[1]

site_info <- site_info %>%
  transmute(
    site_name = as.character(.data[[site_name_col]]),
    site_orig = as.integer(.data[[id_col]]),
    lat = as.numeric(.data[[lat_col]]),
    lon = as.numeric(.data[[lon_col]])
  ) %>%
  filter(!is.na(site_orig)) %>%
  arrange(site_orig)

# -----------------------------
# 4) clean FST sheet
# first column = row IDs
# remaining columns = site IDs as character headers
# -----------------------------
fst_df <- as.data.frame(fst_raw, stringsAsFactors = FALSE)
names(fst_df)[1] <- "row_id"

fst_row_ids <- as.integer(fst_df$row_id)
fst_col_ids <- suppressWarnings(as.integer(names(fst_df)[-1]))

# keep only real matrix columns
keep_cols <- !is.na(fst_col_ids)
fst_col_ids <- fst_col_ids[keep_cols]
fst_vals <- fst_df[, c(TRUE, keep_cols), drop = FALSE]

fst_mat <- as.matrix(fst_vals[, -1, drop = FALSE])
mode(fst_mat) <- "numeric"

rownames(fst_mat) <- fst_row_ids
colnames(fst_mat) <- fst_col_ids

stopifnot(nrow(fst_mat) == ncol(fst_mat))
stopifnot(setequal(as.integer(rownames(fst_mat)), site_info$site_orig))
stopifnot(setequal(as.integer(colnames(fst_mat)), site_info$site_orig))

# reorder to site_info order
fst_mat <- fst_mat[as.character(site_info$site_orig), as.character(site_info$site_orig)]

# enforce symmetry + zero diagonal + no negative values
fst_mat[lower.tri(fst_mat)] <- t(fst_mat)[lower.tri(fst_mat)]
fst_mat[fst_mat < 0] <- 0
diag(fst_mat) <- 0

# -----------------------------
# 5) build final Objective 2 objects
# numeric site labels 1:n
# coords must be: site / lat / lon
# -----------------------------
ETSA_2_coords <- data.frame(
  site = as.character(seq_len(nrow(site_info))),
  lat = site_info$lat,
  lon = site_info$lon,
  stringsAsFactors = FALSE
)

site_lookup <- data.frame(
  site = ETSA_2_coords$site,
  site_orig = site_info$site_orig,
  site_name = site_info$site_name,
  stringsAsFactors = FALSE
)

ETSA_2_fst <- fst_mat
rownames(ETSA_2_fst) <- ETSA_2_coords$site
colnames(ETSA_2_fst) <- ETSA_2_coords$site

# -----------------------------
# 6) optionally read river distance and Watson coords for reference only
# not used in final merged analysis because among-study river distances
# are incomplete
# -----------------------------
if (file.exists(sitecoords_file)) {
  watson_coords <- read.csv(sitecoords_file, stringsAsFactors = FALSE)
} else {
  watson_coords <- NULL
}

if (file.exists(riverdist_file)) {
  watson_riverdist <- read.csv(riverdist_file, stringsAsFactors = FALSE)
} else {
  watson_riverdist <- NULL
}

# -----------------------------
# 7) Euclidean geographic distance matrix
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(ETSA_2_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- ETSA_2_coords$site
colnames(geo_dist_km) <- ETSA_2_coords$site

# -----------------------------
# 8) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1 = rownames(ETSA_2_fst)[row(ETSA_2_fst)[upper.tri(ETSA_2_fst)]],
  site2 = colnames(ETSA_2_fst)[col(ETSA_2_fst)[upper.tri(ETSA_2_fst)]],
  fst = ETSA_2_fst[upper.tri(ETSA_2_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  stringsAsFactors = FALSE
)

ibd_df$site1_name <- site_lookup$site_name[match(ibd_df$site1, site_lookup$site)]
ibd_df$site2_name <- site_lookup$site_name[match(ibd_df$site2, site_lookup$site)]

# -----------------------------
# 9) map of site locations
# tight extent around eastern KY localities
# -----------------------------
state_map <- map_data("state")

lon_rng <- range(ETSA_2_coords$lon, na.rm = TRUE)
lat_rng <- range(ETSA_2_coords$lat, na.rm = TRUE)

x_pad <- max(0.20, diff(lon_rng) * 0.12)
y_pad <- max(0.15, diff(lat_rng) * 0.12)

xlim_use <- c(lon_rng[1] - x_pad, lon_rng[2] + x_pad)
ylim_use <- c(lat_rng[1] - y_pad, lat_rng[2] + y_pad)

plot_sites <- merge(ETSA_2_coords, site_lookup, by = "site", sort = FALSE)

map_plot <- ggplot() +
  geom_polygon(
    data = state_map,
    aes(x = long, y = lat, group = group),
    fill = "grey96",
    color = "grey55",
    linewidth = 0.25
  ) +
  geom_point(
    data = plot_sites,
    aes(x = lon, y = lat),
    size = 3
  ) +
  geom_text(
    data = plot_sites,
    aes(x = lon, y = lat, label = paste0(site, ". ", site_name)),
    nudge_y = 0.012,
    size = 3.2
  ) +
  coord_fixed(
    xlim = xlim_use,
    ylim = ylim_use
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "ETSA-2 sampling locations"
  )

print(map_plot)

ggsave(
  filename = file.path(base_dir, "ETSA-2_map.png"),
  plot = map_plot,
  width = 7.5,
  height = 5.75,
  dpi = 300
)

# -----------------------------
# 10) IBD plot
# -----------------------------
ibd_plot <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.7, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "ETSA-2 isolation by distance"
  )

print(ibd_plot)

ggsave(
  filename = file.path(base_dir, "ETSA-2_IBD.png"),
  plot = ibd_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# -----------------------------
# 11) checks
# -----------------------------
stopifnot(identical(rownames(ETSA_2_fst), ETSA_2_coords$site))
stopifnot(identical(colnames(ETSA_2_fst), ETSA_2_coords$site))
stopifnot(isTRUE(all.equal(ETSA_2_fst, t(ETSA_2_fst))))
stopifnot(isTRUE(all.equal(geo_dist_km, t(geo_dist_km))))

# -----------------------------
# 12) save
# -----------------------------
save(
  ETSA_2_fst,
  ETSA_2_coords,
  file = file.path(data_dir, "ETSA-2.RData")
)

# -----------------------------
# 13) optional console summaries for merge diagnostics
# -----------------------------
cat("\n--- Culley Excel sheets ---\n")
print(sheet_names)

if (!is.null(lm_gp)) {
  cat("\n--- LM_geneticdata.gen summary ---\n")
  cat("n loci:", length(lm_gp$loci), "\n")
  cat("loci:\n")
  print(lm_gp$loci)
  cat("n pops:", length(lm_gp$pops), "\n")
  cat("pop sizes:\n")
  print(sapply(lm_gp$pops, length))
}

if (!is.null(watson_gp)) {
  cat("\n--- Watson Genepop summary ---\n")
  cat("n loci:", length(watson_gp$loci), "\n")
  cat("loci:\n")
  print(watson_gp$loci)
  cat("n pops:", length(watson_gp$pops), "\n")
  cat("pop sizes:\n")
  print(sapply(watson_gp$pops, length))
}

if (!is.null(lm_alleles) && !is.null(watson_alleles)) {
  cat("\n--- Merge diagnostic: allele ranges by locus ---\n")
  merge_diag <- data.frame(
    locus = names(lm_alleles),
    lm_n_alleles = sapply(lm_alleles, length),
    watson_n_alleles = sapply(watson_alleles, length),
    lm_min = sapply(lm_alleles, function(x) min(as.integer(x))),
    lm_max = sapply(lm_alleles, function(x) max(as.integer(x))),
    watson_min = sapply(watson_alleles, function(x) min(as.integer(x))),
    watson_max = sapply(watson_alleles, function(x) max(as.integer(x)))
  )
  print(merge_diag)
}
