# ============================================================
# AMRU-2
# Niagara Falls paper pilot workflow
# Rock bass / Ambloplites rupestris
#
# Updated workflow requested by user:
# - Read VCF
# - Immediately convert VCF -> genlight
# - Trim genlight individual names to the 4-digit tag IDs
# - Read AMRU_xys.xlsx
# - Match spreadsheet Tag # to genlight IDs
# - Store individual Latitude/Longitude as gl@other$latlon
# - Assign populations from spreadsheet "Locality" column ("upper"/"lower")
# - Calculate pairwise FST between upper and lower using gl.fst.pop
# - Build population coordinates as mean XY of each population
# - Plot site map and IBD check
# - Save AMRU_2_fst and AMRU_2_coords to data/AMRU-2.RData
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(vcfR)
  library(adegenet)
  library(dartR)
  library(geosphere)
  library(ggplot2)
  library(maps)
})

# -----------------------------
# 0) paths
# -----------------------------
study_code <- "AMRU-2"
species_code <- "AMRU"

base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/AMRU-2"
base_dir <- path.expand(base_dir)

data_dir <- file.path(base_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

xlsx_file <- file.path(base_dir, paste0(species_code, "_xys.xlsx"))
vcf_file  <- file.path(base_dir, "Ambloplites.vcf")

if (!file.exists(xlsx_file)) {
  stop("Could not find expected XY workbook: ", xlsx_file)
}
if (!file.exists(vcf_file)) {
  stop("Could not find expected VCF: ", vcf_file)
}

cat("Workbook: ", xlsx_file, "\n", sep = "")
cat("VCF: ", vcf_file, "\n", sep = "")

# -----------------------------
# 1) read VCF and immediately convert to genlight
# -----------------------------
vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)
gl <- vcfR::vcfR2genlight(vcf)

# trim genlight individual names to just the 4-digit tag
extract_tag4 <- function(x) {
  m <- regexpr("[0-9]{4}", x)
  out <- ifelse(m > 0, regmatches(x, m), NA_character_)
  as.character(out)
}

trimmed_ids <- extract_tag4(indNames(gl))

if (any(is.na(trimmed_ids))) {
  stop(
    "Could not extract a 4-digit tag ID from some VCF sample names.\n",
    "Problem names: ",
    paste(indNames(gl)[is.na(trimmed_ids)], collapse = ", ")
  )
}

indNames(gl) <- trimmed_ids

if (anyDuplicated(indNames(gl))) {
  stop(
    "Duplicate 4-digit IDs were produced after trimming VCF sample names.\n",
    "Duplicates: ",
    paste(unique(indNames(gl)[duplicated(indNames(gl))]), collapse = ", ")
  )
}

cat("Individuals in genlight: ", nInd(gl), "\n", sep = "")
cat("Loci in genlight: ", nLoc(gl), "\n", sep = "")
cat("Trimmed individual IDs:\n")
print(indNames(gl))

# -----------------------------
# 2) read XY spreadsheet
# -----------------------------
xy_raw <- readxl::read_excel(xlsx_file) |>
  as.data.frame()

names(xy_raw) <- trimws(names(xy_raw))
names(xy_raw) <- gsub("\\s+", "_", names(xy_raw))
names(xy_raw) <- gsub("#", "", names(xy_raw), fixed = TRUE)

tag_col <- names(xy_raw)[tolower(names(xy_raw)) %in% c("tag", "tag_", "tag_number", "tagno", "tagid")]
if (length(tag_col) == 0) {
  stop("Could not identify a Tag column in ", xlsx_file,
       ". Found columns: ", paste(names(xy_raw), collapse = ", "))
}
tag_col <- tag_col[1]

lat_col <- names(xy_raw)[tolower(names(xy_raw)) == "latitude"]
lon_col <- names(xy_raw)[tolower(names(xy_raw)) == "longitude"]
loc_col <- names(xy_raw)[tolower(names(xy_raw)) == "locality"]

if (length(lat_col) == 0 || length(lon_col) == 0) {
  stop("Could not identify Latitude and Longitude columns in ", xlsx_file)
}
if (length(loc_col) == 0) {
  stop("Could not identify a Locality column in ", xlsx_file)
}

xy_meta <- xy_raw |>
  transmute(
    tag = extract_tag4(as.character(.data[[tag_col]])),
    locality = tolower(trimws(as.character(.data[[loc_col[1]]]))),
    lat = as.numeric(.data[[lat_col[1]]]),
    lon = as.numeric(.data[[lon_col[1]]])
  ) |>
  filter(!is.na(tag), !is.na(locality), !is.na(lat), !is.na(lon))

if (!all(xy_meta$locality %in% c("upper", "lower"))) {
  stop(
    "Locality column must contain only 'upper' or 'lower'.\n",
    "Found: ",
    paste(sort(unique(xy_meta$locality)), collapse = ", ")
  )
}

if (anyDuplicated(xy_meta$tag)) {
  dup <- unique(xy_meta$tag[duplicated(xy_meta$tag)])
  stop("Duplicate Tag values found in spreadsheet: ", paste(dup, collapse = ", "))
}

# -----------------------------
# 3) match spreadsheet metadata to genlight IDs
# -----------------------------
match_idx <- match(indNames(gl), xy_meta$tag)

if (any(is.na(match_idx))) {
  stop(
    "Some genlight IDs were not found in spreadsheet Tag values.\n",
    "Missing IDs: ",
    paste(indNames(gl)[is.na(match_idx)], collapse = ", ")
  )
}

sample_meta <- xy_meta[match_idx, , drop = FALSE]
stopifnot(all(sample_meta$tag == indNames(gl)))

# -----------------------------
# 4) add lat/lon matrix to gl@other$latlon
# -----------------------------
gl@other$latlon <- as.matrix(sample_meta[, c("lat", "lon")])
rownames(gl@other$latlon) <- sample_meta$tag
colnames(gl@other$latlon) <- c("lat", "lon")

# also keep full matched metadata for debugging
gl@other$sample_meta <- sample_meta

# -----------------------------
# 5) assign populations from Locality
# -----------------------------
pop(gl) <- as.factor(sample_meta$locality)

cat("Population counts:\n")
print(table(pop(gl)))

# -----------------------------
# 6) calculate pairwise FST between upper/lower
# -----------------------------
fst_out <- dartR::gl.fst.pop(gl)

AMRU_2_fst <- as.matrix(fst_out)
storage.mode(AMRU_2_fst) <- "numeric"

# enforce symmetry explicitly and clean negatives
AMRU_2_fst[is.na(AMRU_2_fst)] <- 0
AMRU_2_fst <- (AMRU_2_fst + t(AMRU_2_fst)) / 2
AMRU_2_fst[AMRU_2_fst < 0] <- 0
diag(AMRU_2_fst) <- 0

# order populations lower then upper if both are present
pop_order <- c("lower", "upper")
pop_order <- pop_order[pop_order %in% rownames(AMRU_2_fst)]
AMRU_2_fst <- AMRU_2_fst[pop_order, pop_order, drop = FALSE]

# rename rows/cols to numeric site IDs for Objective 2 workflow
rownames(AMRU_2_fst) <- colnames(AMRU_2_fst) <- as.character(seq_len(nrow(AMRU_2_fst)))

# -----------------------------
# 7) build coords from mean XY of each pop
# -----------------------------
coords_key <- sample_meta |>
  group_by(locality) |>
  summarise(
    lat = mean(lat, na.rm = TRUE),
    lon = mean(lon, na.rm = TRUE),
    n_ind = n(),
    .groups = "drop"
  ) |>
  mutate(locality = factor(locality, levels = c("lower", "upper"))) |>
  arrange(locality) |>
  mutate(site_id = seq_len(n())) |>
  select(site_id, locality, n_ind, lat, lon)

AMRU_2_coords <- coords_key |>
  select(site_id, lat, lon)

cat("Population centroid coordinates:\n")
print(coords_key)

# -----------------------------
# 8) map plot
# -----------------------------
world_df <- map_data("world") |>
  filter(region %in% c("USA", "Canada"))

coords_plot <- coords_key

xpad <- max(0.1, diff(range(coords_plot$lon)) * 0.40)
ypad <- max(0.1, diff(range(coords_plot$lat)) * 0.40)

p_map <- ggplot() +
  geom_polygon(
    data = world_df,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey65",
    linewidth = 0.2
  ) +
  geom_point(
    data = coords_plot,
    aes(x = lon, y = lat),
    size = 2.8
  ) +
  geom_text(
    data = coords_plot,
    aes(x = lon, y = lat, label = paste0(site_id, ". ", locality)),
    nudge_y = ypad * 0.04,
    size = 3
  ) +
  coord_quickmap(
    xlim = range(coords_plot$lon) + c(-xpad, xpad),
    ylim = range(coords_plot$lat) + c(-ypad, ypad),
    expand = FALSE
  ) +
  theme_classic() +
  labs(
    title = "AMRU-2 sampling populations",
    x = "Longitude",
    y = "Latitude"
  )

print(p_map)

# -----------------------------
# 9) quick IBD plot
# -----------------------------
coord_mat <- as.matrix(AMRU_2_coords[, c("lon", "lat")])
geo_dist_km <- geosphere::distm(coord_mat, fun = geosphere::distHaversine) / 1000

upper_idx <- upper.tri(AMRU_2_fst)

ibd_df <- data.frame(
  fst = AMRU_2_fst[upper_idx],
  dist_km = geo_dist_km[upper_idx]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    title = "AMRU-2 IBD quick-check",
    x = "Great-circle distance (km)",
    y = expression(F[ST])
  )

print(p_ibd)

# -----------------------------
# 10) save outputs
# -----------------------------
save(
  AMRU_2_fst,
  AMRU_2_coords,
  file = file.path(data_dir, "AMRU-2.RData")
)

assign("gl", gl, envir = .GlobalEnv)
assign("AMRU_2_fst", AMRU_2_fst, envir = .GlobalEnv)
assign("AMRU_2_coords", AMRU_2_coords, envir = .GlobalEnv)

cat("\n==================== SUMMARY ====================\n")
cat("Study code: ", study_code, "\n", sep = "")
cat("Workbook: ", xlsx_file, "\n", sep = "")
cat("VCF: ", vcf_file, "\n", sep = "")
cat("Individuals matched: ", nrow(sample_meta), "\n", sep = "")
cat("Populations: ", paste(unique(sample_meta$locality), collapse = ", "), "\n", sep = "")
cat("FST matrix dimensions: ", nrow(AMRU_2_fst), " x ", ncol(AMRU_2_fst), "\n", sep = "")
cat("RData saved to: ", file.path(data_dir, "AMRU-2.RData"), "\n", sep = "")
cat("=================================================\n\n")
