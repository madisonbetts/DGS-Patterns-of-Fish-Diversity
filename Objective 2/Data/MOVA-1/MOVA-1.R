# ============================================================
# MOVA-1
# updated to match MOVA_xys.xlsx columns:
#   tag_, genus, locality, total_reads, latitude, longitude
#
# Key behavior:
# - Read VCF -> genlight
# - Extract numeric tag IDs from VCF sample names
# - Match to spreadsheet tag_ column
# - Assign populations from locality
# - Calculate pairwise FST across retained populations
# - Build coords from site_xys.xlsx if usable, else mean individual XY by locality
# - Save MOVA_1_fst and MOVA_1_coords to data/MOVA-1.RData
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(vcfR)
  library(adegenet)
  library(dartR)
  library(geosphere)
  library(ggplot2)
})

# -----------------------------
# 0) paths
# -----------------------------
study_code <- "MOVA-1"

base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/MOVA-1"
base_dir <- path.expand(base_dir)

data_dir <- file.path(base_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

vcf_file <- file.path(base_dir, "Mox_valenc.vcf")
ind_xlsx <- file.path(base_dir, "MOVA_xys.xlsx")
site_xlsx <- file.path(base_dir, "site_xys.xlsx")

if (!file.exists(vcf_file)) stop("Could not find VCF: ", vcf_file)
if (!file.exists(ind_xlsx)) stop("Could not find metadata workbook: ", ind_xlsx)

# -----------------------------
# helpers
# -----------------------------
clean_names <- function(x) {
  x <- trimws(x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("[#()]", "", x)
  x <- gsub("\\.+", "_", x)
  tolower(x)
}

find_col <- function(df, candidates) {
  nms <- names(df)
  hit <- nms[nms %in% candidates]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

# extract first digit run from IDs / tags
extract_tag <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  m <- regexpr("[0-9]+", x)
  out <- ifelse(m > 0, regmatches(x, m), NA_character_)
  out
}

# -----------------------------
# 1) read VCF -> genlight
# -----------------------------
vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)
gl <- vcfR::vcfR2genlight(vcf)

orig_ids <- indNames(gl)
std_ids <- extract_tag(indNames(gl))

if (all(is.na(std_ids))) {
  stop("Could not extract numeric tag IDs from VCF sample names.")
}

indNames(gl) <- std_ids
gl@other$original_vcf_names <- orig_ids

# drop missing extracted IDs
keep <- !is.na(indNames(gl))
if (!all(keep)) {
  cat("Dropping VCF individuals with no numeric tag extracted:\n")
  print(orig_ids[!keep])
  gl <- gl[keep]
  orig_ids <- orig_ids[keep]
}

# dedupe after trimming
dup_idx <- duplicated(indNames(gl))
if (any(dup_idx)) {
  cat("Duplicate extracted VCF IDs found; keeping first occurrence only:\n")
  print(unique(indNames(gl)[dup_idx]))
  gl <- gl[!dup_idx]
  orig_ids <- orig_ids[!dup_idx]
}

cat("Individuals in genlight after tag extraction: ", nInd(gl), "\n", sep = "")
cat("Loci in genlight: ", nLoc(gl), "\n", sep = "")

# -----------------------------
# 2) read individual metadata workbook
# -----------------------------
ind_raw <- readxl::read_excel(ind_xlsx) |> as.data.frame()
names(ind_raw) <- clean_names(names(ind_raw))

cat("Columns in MOVA_xys.xlsx:\n")
print(names(ind_raw))

id_col <- find_col(ind_raw, c(
  "tag_", "tag", "individual_id", "ind", "id", "sample_id", "sample", "fish_id", "vcf_id", "name"
))
site_col <- find_col(ind_raw, c(
  "locality", "site", "site_name", "population", "pop", "location"
))
lat_col <- find_col(ind_raw, c(
  "latitude", "lat", "y", "site_lat", "site_latitude"
))
lon_col <- find_col(ind_raw, c(
  "longitude", "lon", "long", "x", "site_lon", "site_longitude"
))

if (is.na(id_col)) {
  stop("Could not identify ID/tag column in MOVA_xys.xlsx.\nFound columns: ",
       paste(names(ind_raw), collapse = ", "))
}
if (is.na(site_col)) {
  stop("Could not identify site/locality column in MOVA_xys.xlsx.\nFound columns: ",
       paste(names(ind_raw), collapse = ", "))
}
if (is.na(lat_col) || is.na(lon_col)) {
  stop("Could not identify latitude/longitude columns in MOVA_xys.xlsx.\nFound columns: ",
       paste(names(ind_raw), collapse = ", "))
}

ind_meta <- ind_raw |>
  transmute(
    ind_id = extract_tag(.data[[id_col]]),
    site_name = as.character(.data[[site_col]]),
    lat = as.numeric(.data[[lat_col]]),
    lon = as.numeric(.data[[lon_col]])
  ) |>
  filter(!is.na(ind_id), ind_id != "", !is.na(site_name), !is.na(lat), !is.na(lon))

if (nrow(ind_meta) == 0) {
  stop("No usable rows found in MOVA_xys.xlsx after filtering.")
}

if (anyDuplicated(ind_meta$ind_id)) {
  cat("Duplicate tag IDs in MOVA_xys.xlsx; keeping first occurrence only for:\n")
  print(unique(ind_meta$ind_id[duplicated(ind_meta$ind_id)]))
  ind_meta <- ind_meta |>
    distinct(ind_id, .keep_all = TRUE)
}

# -----------------------------
# 3) match VCF individuals to metadata
# -----------------------------
match_idx <- match(indNames(gl), ind_meta$ind_id)
matched <- !is.na(match_idx)

cat("Matched individuals: ", sum(matched), "\n", sep = "")
cat("Dropped unmatched individuals: ", sum(!matched), "\n", sep = "")

if (sum(matched) == 0) {
  stop("No VCF individuals matched MOVA_xys.xlsx tag values.")
}

if (any(!matched)) {
  dropped_df <- data.frame(
    original_name = gl@other$original_vcf_names[!matched],
    extracted_id = indNames(gl)[!matched],
    stringsAsFactors = FALSE
  )
  cat("Dropped unmatched individuals:\n")
  print(dropped_df)
}

gl <- gl[matched]
meta_matched <- ind_meta[match_idx[matched], , drop = FALSE]
stopifnot(all(indNames(gl) == meta_matched$ind_id))

gl@other$latlon <- as.matrix(meta_matched[, c("lat", "lon")])
rownames(gl@other$latlon) <- meta_matched$ind_id
colnames(gl@other$latlon) <- c("lat", "lon")
gl@other$sample_meta <- meta_matched

# -----------------------------
# 4) assign populations from locality
# -----------------------------
pop(gl) <- as.factor(meta_matched$site_name)

cat("Population counts before singleton removal:\n")
print(table(pop(gl)))

pop_counts <- table(pop(gl))
keep_pops <- names(pop_counts[pop_counts >= 2])

if (length(keep_pops) < 2) {
  stop("Need at least 2 populations with >= 2 individuals each for pairwise FST.")
}

drop_singletons <- !(pop(gl) %in% keep_pops)
if (any(drop_singletons)) {
  cat("Dropping singleton populations:\n")
  print(as.data.frame(pop_counts[pop_counts < 2]))
  gl <- gl[!drop_singletons]
  meta_matched <- meta_matched[!drop_singletons, , drop = FALSE]
  pop(gl) <- droplevels(pop(gl))
  gl@other$latlon <- as.matrix(meta_matched[, c("lat", "lon")])
  rownames(gl@other$latlon) <- meta_matched$ind_id
  colnames(gl@other$latlon) <- c("lat", "lon")
  gl@other$sample_meta <- meta_matched
}

cat("Population counts retained for FST:\n")
print(table(pop(gl)))

# -----------------------------
# 5) calculate pairwise FST
# -----------------------------
fst_out <- dartR::gl.fst.pop(gl)

if (is.list(fst_out) && "Fsts" %in% names(fst_out)) {
  fst_mat <- as.matrix(fst_out$Fsts)
} else {
  fst_mat <- as.matrix(fst_out)
}

storage.mode(fst_mat) <- "numeric"
fst_mat[is.na(fst_mat)] <- 0
fst_mat <- (fst_mat + t(fst_mat)) / 2
fst_mat[fst_mat < 0] <- 0
diag(fst_mat) <- 0

pop_levels <- rownames(fst_mat)
if (is.null(pop_levels)) {
  pop_levels <- levels(pop(gl))
  rownames(fst_mat) <- colnames(fst_mat) <- pop_levels
}

MOVA_1_fst <- fst_mat
rownames(MOVA_1_fst) <- colnames(MOVA_1_fst) <- as.character(seq_len(nrow(MOVA_1_fst)))

# -----------------------------
# 6) build population coordinates
# use site_xys.xlsx if matchable, else mean individual XY by locality
# -----------------------------
coords_key <- NULL

if (file.exists(site_xlsx)) {
  site_raw <- readxl::read_excel(site_xlsx) |> as.data.frame()
  names(site_raw) <- clean_names(names(site_raw))

  site_name_col <- find_col(site_raw, c("site", "site_name", "locality", "population", "pop", "location"))
  site_lat_col  <- find_col(site_raw, c("latitude", "lat", "y", "site_lat", "site_latitude"))
  site_lon_col  <- find_col(site_raw, c("longitude", "lon", "long", "x", "site_lon", "site_longitude"))

  if (!is.na(site_name_col) && !is.na(site_lat_col) && !is.na(site_lon_col)) {
    site_meta <- site_raw |>
      transmute(
        site_name = as.character(.data[[site_name_col]]),
        lat = as.numeric(.data[[site_lat_col]]),
        lon = as.numeric(.data[[site_lon_col]])
      ) |>
      filter(!is.na(site_name), !is.na(lat), !is.na(lon)) |>
      distinct(site_name, .keep_all = TRUE)

    matched_sites <- match(pop_levels, site_meta$site_name)
    if (all(!is.na(matched_sites))) {
      coords_key <- site_meta[matched_sites, , drop = FALSE] |>
        mutate(
          n_ind = as.integer(table(pop(gl))[site_name]),
          site_id = seq_len(n())
        ) |>
        select(site_id, site_name, n_ind, lat, lon)
      cat("Using site_xys.xlsx for population coordinates.\n")
    }
  }
}

if (is.null(coords_key)) {
  coords_key <- meta_matched |>
    group_by(site_name) |>
    summarise(
      lat = mean(lat, na.rm = TRUE),
      lon = mean(lon, na.rm = TRUE),
      n_ind = n(),
      .groups = "drop"
    ) |>
    slice(match(pop_levels, site_name)) |>
    mutate(site_id = seq_len(n())) |>
    select(site_id, site_name, n_ind, lat, lon)

  cat("Using mean individual XY by locality for coordinates.\n")
}

MOVA_1_coords <- coords_key |>
  select(site_id, lat, lon)

cat("Population centroid coordinates:\n")
print(coords_key)

# -----------------------------
# 7) map plot
# -----------------------------
map_df <- map_data("state")

xpad <- max(0.1, diff(range(coords_key$lon)) * 0.30)
ypad <- max(0.1, diff(range(coords_key$lat)) * 0.30)

ggplot() +
  geom_polygon(
    data = map_df,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey65",
    linewidth = 0.2
  ) +
  geom_point(
    data = coords_key,
    aes(x = lon, y = lat),
    size = 2.7
  ) +
  geom_text(
    data = coords_key,
    aes(x = lon, y = lat, label = paste0(site_id, '. ', site_name)),
    nudge_y = ypad * 0.04,
    size = 3
  ) +
  coord_quickmap(
    xlim = range(coords_key$lon) + c(-xpad, xpad),
    ylim = range(coords_key$lat) + c(-ypad, ypad),
    expand = FALSE
  ) +
  theme_classic() +
  labs(
    title = "MOVA-1 sampling populations",
    x = "Longitude",
    y = "Latitude"
  )

# -----------------------------
# 8) quick IBD plot
# -----------------------------
coord_mat <- as.matrix(MOVA_1_coords[, c("lon", "lat")])
geo_dist_km <- geosphere::distm(coord_mat, fun = geosphere::distHaversine) / 1000
upper_idx <- upper.tri(MOVA_1_fst)

ibd_df <- data.frame(
  fst = MOVA_1_fst[upper_idx],
  dist_km = geo_dist_km[upper_idx]
)

ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    title = "MOVA-1 IBD quick-check",
    x = "Great-circle distance (km)",
    y = expression(F[ST])
  )

# -----------------------------
# 9) save outputs
# -----------------------------
save(
  MOVA_1_fst,
  MOVA_1_coords,
  file = file.path(data_dir, "MOVA-1.RData")
)

assign("gl", gl, envir = .GlobalEnv)
assign("MOVA_1_fst", MOVA_1_fst, envir = .GlobalEnv)
assign("MOVA_1_coords", MOVA_1_coords, envir = .GlobalEnv)

cat("\n==================== SUMMARY ====================\n")
cat("Study code: ", study_code, "\n", sep = "")
cat("Workbook: ", ind_xlsx, "\n", sep = "")
cat("VCF: ", vcf_file, "\n", sep = "")
cat("Individuals retained after matching: ", nrow(meta_matched), "\n", sep = "")
cat("Populations retained for FST: ", nrow(coords_key), "\n", sep = "")
cat("FST matrix dimensions: ", nrow(MOVA_1_fst), " x ", ncol(MOVA_1_fst), "\n", sep = "")
cat("RData saved to: ", file.path(data_dir, "MOVA-1.RData"), "\n", sep = "")
cat("=================================================\n\n")
