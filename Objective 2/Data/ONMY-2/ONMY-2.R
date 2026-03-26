############################
## ONMY-all steelhead
## Hand et al. 2016 (Molecular Ecology)
## Neutral SNPs only
##
## This script:
##   1) reads all five neutral genotype sheets from SH_Raw_Genetic.xlsx
##   2) converts the two-column-per-locus GenAlEx-style layout into
##      hierfstat-ready diploid genotype codes
##   3) combines all individuals across all 79 populations
##   4) calculates pairwise Weir & Cockerham FST among all populations
##   5) forces symmetry and sets negative FST to 0
##   6) assigns numeric site IDs 1:79
##   7) builds ONMY_2_coords as site_id / lat / lon only
##   8) plots sites and straight-line IBD
############################

suppressPackageStartupMessages({
  library(readxl)
  library(hierfstat)
  library(geosphere)
  library(ggplot2)
  library(maps)
})

# -----------------------------
# paths
# -----------------------------
xlsx_file <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ONMY-2/doi_10_5061_dryad_pr065__v20151211/DataforDryad/SH_Raw_Genetic.xlsx"

out_dir  <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ONMY-2"
data_dir <- file.path(out_dir, "data")

# -----------------------------
# sheets to combine
# -----------------------------
neutral_sheets <- c(
  "Lo_Col_genos_neutral",
  "Salmon_MPG_genos_neutral",
  "Clearwater_MPG_genos_neutral",
  "Cascades_MPG_genos_neutral",
  "JohnDay_genos_neutral"
)

# -----------------------------
# helper functions
# -----------------------------
fmt_allele <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "<NA>", "NULL")] <- NA_character_
  x
}

make_hf_genotype <- function(a1, a2) {
  a1 <- fmt_allele(a1)
  a2 <- fmt_allele(a2)
  
  bad1 <- !is.na(a1) & !grepl("^[0-9]+$", a1)
  bad2 <- !is.na(a2) & !grepl("^[0-9]+$", a2)
  
  if (any(bad1 | bad2)) {
    stop("Non-numeric allele codes detected in the raw genotype workbook.")
  }
  
  ifelse(
    is.na(a1) | is.na(a2),
    NA_integer_,
    as.integer(sprintf("%03d%03d", as.integer(a1), as.integer(a2)))
  )
}

sanitize_locus_names <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  make.unique(x, sep = "_")
}

read_neutral_sheet <- function(xlsx_file, sheet_name) {
  raw_xlsx <- readxl::read_excel(
    path = xlsx_file,
    sheet = sheet_name,
    col_names = FALSE
  )
  raw_xlsx <- as.data.frame(raw_xlsx, stringsAsFactors = FALSE)
  
  header_row <- raw_xlsx[3, ]
  ind_dat <- raw_xlsx[-c(1, 2, 3), , drop = FALSE]
  rownames(ind_dat) <- NULL
  names(ind_dat)[1:2] <- c("sample_id", "pop_code")
  
  locus_names <- as.character(unlist(header_row[3:ncol(raw_xlsx)]))
  locus_names <- trimws(locus_names)
  
  if (length(locus_names) %% 2 != 0) {
    stop(paste0("Expected an even number of genotype columns in sheet ", sheet_name, "."))
  }
  
  locus_pairs <- split(seq_along(locus_names), ceiling(seq_along(locus_names) / 2))
  
  pair_ok <- vapply(
    locus_pairs,
    function(idx) length(unique(locus_names[idx])) == 1,
    logical(1)
  )
  
  if (!all(pair_ok)) {
    stop(paste0("Not all genotype columns occur as clean allele pairs in sheet ", sheet_name, "."))
  }
  
  paired_names <- sanitize_locus_names(
    vapply(locus_pairs, function(idx) locus_names[idx[1]], character(1))
  )
  
  pop_codes <- unique(ind_dat$pop_code)
  pop_n <- as.integer(table(factor(ind_dat$pop_code, levels = pop_codes)))
  
  list(
    sheet_name  = sheet_name,
    ind_dat     = ind_dat,
    pop_codes   = pop_codes,
    pop_n       = pop_n,
    locus_pairs = locus_pairs,
    locus_names = paired_names
  )
}

sheet_to_hf <- function(sheet_obj, global_site_lookup) {
  ind_dat <- sheet_obj$ind_dat
  locus_pairs <- sheet_obj$locus_pairs
  locus_names <- sheet_obj$locus_names
  
  hf <- data.frame(
    pop = match(ind_dat$pop_code, global_site_lookup$pop_code),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(locus_pairs)) {
    cols <- locus_pairs[[i]] + 2L
    hf[[locus_names[i]]] <- make_hf_genotype(ind_dat[[cols[1]]], ind_dat[[cols[2]]])
  }
  
  hf
}

# -----------------------------
# site metadata in within-sheet order
# -----------------------------
site_meta <- list(
  Lo_Col_genos_neutral = data.frame(
    site_name = c(
      "Cowlitz R.", "Kalama R.", "E. F. Lewis R.", "N. F. Lewis R.",
      "Clackamas R.", "N. F. Eagle Cr.", "Still Cr.", "E. Fork Hood R.",
      "W. Fork Hood R.", "Coweeman R."
    ),
    metapop = "Lower Columbia",
    stringsAsFactors = FALSE
  ),
  Salmon_MPG_genos_neutral = data.frame(
    site_name = c(
      "Hazard Cr.", "Sawtooth Weir", "W. F. Yankee Fork R.", "Morgan Cr.",
      "Pahsimeroi Weir", "Hayden Cr.", "N. F. Salmon R.", "Marsh Cr.",
      "Sulphur Cr.", "Rapid R.", "Pistol Cr.", "Camas Cr.", "Loon Cr.",
      "Upper Big Cr.", "Chamberlain Cr.", "Bargamin Cr.",
      "East Fork S. F. Salmon R.", "Secesh R.", "Lick Cr.",
      "Stolle Meadows", "Boulder Cr.", "Slate Cr.", "Whitebird Cr."
    ),
    metapop = "Salmon",
    stringsAsFactors = FALSE
  ),
  Clearwater_MPG_genos_neutral = data.frame(
    site_name = c(
      "Mission Cr.", "Colt/ Lochsa R.", "Storm/ Lochsa R.", "Crooked Fork/Lochsa R.",
      "Lake/ Lochsa R.", "Fish/ Lochsa R.", "Canyon/ Lochsa R.",
      "Little Clearwater/ Selway R.", "Whitecap/ Selway R.", "Bear/ Selway R.",
      "N. F. Moose/ Selway R.", "Three Links/ Selway R.", "Gedney/ Selway R.",
      "OHara/ Selway R.", "Clear Cr.", "Crooked R.", "Tenmile Cr.",
      "Johns Cr.", "W. F. Potlatch R.", "E. F. Potlatch R.",
      "Big Bear Cr.", "Little Bear Cr."
    ),
    metapop = "Clearwater",
    stringsAsFactors = FALSE
  ),
  Cascades_MPG_genos_neutral = data.frame(
    site_name = c(
      "Shitike Cr.", "Buckhollow Cr.", "Trout Cr.", "Fifteen Cr.",
      "Lower Little Klickitat R.", "Lower Summit Cr.", "Upper Trout Cr.",
      "Deadcanyon Cr.", "Lower White Cr.", "Snyder Cr.", "Surveyor Cr.",
      "Swale Cr.", "Rock Cr.", "Squaw Cr."
    ),
    metapop = "Eastern Cascades",
    stringsAsFactors = FALSE
  ),
  JohnDay_genos_neutral = data.frame(
    site_name = c(
      "Beech Cr.", "Upper Mainstem J.D.", "Baldy Cr.", "Lower Mainstem J.D.",
      "Upper M. F. J.D.", "Granite Cr.", "Middle N. F. J.D.",
      "Big Wall Cr.", "Deer Cr.", "Murderers Cr."
    ),
    metapop = "John Day",
    stringsAsFactors = FALSE
  )
)

# -----------------------------
# read all sheets
# -----------------------------
sheet_objs <- lapply(neutral_sheets, function(sh) read_neutral_sheet(xlsx_file, sh))
names(sheet_objs) <- neutral_sheets

# -----------------------------
# verify same neutral locus panel across sheets
# -----------------------------
ref_loci <- sheet_objs[[1]]$locus_names
same_loci <- vapply(
  sheet_objs,
  function(x) identical(x$locus_names, ref_loci),
  logical(1)
)

if (!all(same_loci)) {
  stop("Neutral locus names are not identical across all sheets.")
}

# -----------------------------
# build one site lookup across all 79 populations
# -----------------------------
site_lookup_list <- vector("list", length(neutral_sheets))

for (i in seq_along(neutral_sheets)) {
  sh <- neutral_sheets[i]
  sh_obj <- sheet_objs[[sh]]
  sh_meta <- site_meta[[sh]]
  
  if (nrow(sh_meta) != length(sh_obj$pop_codes)) {
    stop(paste0("Site metadata length does not match pop code count for sheet ", sh, "."))
  }
  
  site_lookup_list[[i]] <- data.frame(
    sheet_name = sh,
    metapop    = sh_meta$metapop,
    pop_code   = sh_obj$pop_codes,
    site_name  = sh_meta$site_name,
    n_ind      = sh_obj$pop_n,
    stringsAsFactors = FALSE
  )
}

site_lookup <- do.call(rbind, site_lookup_list)
site_lookup$site_id <- seq_len(nrow(site_lookup))
site_lookup <- site_lookup[, c("site_id", "sheet_name", "metapop", "pop_code", "site_name", "n_ind")]

# -----------------------------
# convert each sheet to hierfstat format and combine individuals
# -----------------------------
hf_list <- lapply(sheet_objs, function(x) sheet_to_hf(x, site_lookup))
ONMY_2_hf <- do.call(rbind, hf_list)

keep_loci <- colSums(!is.na(ONMY_2_hf[, -1, drop = FALSE])) > 0
ONMY_2_hf <- ONMY_2_hf[, c(TRUE, keep_loci), drop = FALSE]

# -----------------------------
# pairwise FST across all 79 populations
# -----------------------------
ONMY_2_fst <- hierfstat::pairwise.WCfst(ONMY_2_hf)
ONMY_2_fst <- as.matrix(ONMY_2_fst)

ONMY_2_fst[lower.tri(ONMY_2_fst)] <- t(ONMY_2_fst)[lower.tri(ONMY_2_fst)]
diag(ONMY_2_fst) <- 0
ONMY_2_fst[ONMY_2_fst < 0] <- 0

rownames(ONMY_2_fst) <- site_lookup$site_id
colnames(ONMY_2_fst) <- site_lookup$site_id

# -----------------------------
# hard-coded coordinates in site_lookup order
# final object wanted: ONMY_2_coords = site_id / lat / lon
# -----------------------------
coord_lookup <- data.frame(
  site_id = site_lookup$site_id,
  lat = c(
    46.5026, 46.0449, 45.8655, 45.9516, 45.2417,
    45.3514, 45.3309, 45.5745, 45.6047, 46.1408,
    45.1836, 44.1506, 44.3514, 44.6135, 44.6844,
    44.8616, 45.4094, 44.4493, 44.5526, 44.6790,
    44.7217, 44.8918, 44.5976, 45.1523, 45.4523,
    45.5716, 45.0127, 45.0268, 45.0692, 44.6070,
    45.2019, 45.6380, 45.7523,
    46.3672, 46.4311, 46.4607, 46.5251, 46.4632,
    46.3336, 46.2161, 45.7441, 45.8689, 46.0191,
    46.1634, 46.0981, 46.0583, 46.0809, 46.0486,
    45.8211, 45.8057, 45.8224, 46.8054, 46.7984,
    46.6306, 46.6372,
    44.7615, 45.2507, 44.8217, 45.6251, 45.8434,
    45.9876, 46.0774, 45.9420, 46.0133, 45.8281,
    46.1957, 45.8091, 45.8066, 45.8191,
    44.4733, 44.3945, 44.3629, 44.7358, 44.6195,
    44.8383, 45.0213, 44.8833, 44.1956, 44.2646
  ),
  lon = c(
    -122.5881, -122.8039, -122.7184, -122.5654, -122.2817,
    -122.3840, -121.9158, -121.6271, -121.6335, -122.8536,
    -116.2995, -114.8851, -114.7297, -114.1641, -114.0403,
    -113.6319, -113.9918, -115.2301, -115.2974, -115.1490,
    -115.1488, -114.7222, -114.8123, -115.2975, -114.9310,
    -115.1919, -115.7129, -115.7082, -115.8140, -115.6810,
    -116.3114, -116.2828, -116.3198,
    -116.7360, -114.5395, -114.5467, -114.6786, -114.9965,
    -115.3471, -115.5559, -114.7895, -114.7205, -114.8378,
    -114.9006, -115.0728, -115.3141, -115.5179, -115.7814,
    -115.5272, -115.6833, -115.8887, -116.4182, -116.4194,
    -116.6562, -116.6780,
    -121.2288, -121.0230, -121.0858, -121.0656, -121.0605,
    -121.1255, -121.2122, -121.1439, -121.1500, -121.1724,
    -121.2557, -121.0652, -120.5104, -120.4890,
    -119.0332, -118.5764, -119.7700, -120.3071, -118.5679,
    -118.4770, -118.9905, -119.4121, -119.4716, -119.2857
  )
)

if (nrow(coord_lookup) != nrow(site_lookup)) {
  stop("Coordinate table length does not match site_lookup.")
}

ONMY_2_coords <- merge(site_lookup, coord_lookup, by = "site_id", sort = FALSE)
ONMY_2_coords <- ONMY_2_coords[match(site_lookup$site_id, ONMY_2_coords$site_id), ]
rownames(ONMY_2_coords) <- NULL

ONMY_2_coords <- ONMY_2_coords[, c("site_id", "lat", "lon")]

# -----------------------------
# straight-line distances + IBD
# -----------------------------
coord_mat <- as.matrix(ONMY_2_coords[, c("lon", "lat")])

geo_dist_km <- geosphere::distm(coord_mat, fun = geosphere::distHaversine) / 1000
rownames(geo_dist_km) <- ONMY_2_coords$site_id
colnames(geo_dist_km) <- ONMY_2_coords$site_id

ONMY_2_ibd <- data.frame(
  site1   = rownames(ONMY_2_fst)[row(ONMY_2_fst)[upper.tri(ONMY_2_fst)]],
  site2   = colnames(ONMY_2_fst)[col(ONMY_2_fst)[upper.tri(ONMY_2_fst)]],
  fst     = ONMY_2_fst[upper.tri(ONMY_2_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  stringsAsFactors = FALSE
)

# -----------------------------
# plots
# -----------------------------
usa <- map_data("state")
canada <- map_data("world", region = "Canada")

ggplot() +
  geom_polygon(
    data = usa,
    aes(x = long, y = lat, group = group),
    fill = "gray95",
    color = "gray70",
    linewidth = 0.2
  ) +
  geom_polygon(
    data = canada,
    aes(x = long, y = lat, group = group),
    fill = "gray95",
    color = "gray70",
    linewidth = 0.2
  ) +
  geom_point(
    data = ONMY_2_coords,
    aes(x = lon, y = lat, shape = metapop),
    size = 2.8
  ) +
  geom_text(
    data = ONMY_2_coords,
    aes(x = lon, y = lat, label = site_id),
    size = 2.7,
    nudge_y = 0.18
  ) +
  coord_quickmap(
    xlim = range(ONMY_2_coords$lon, na.rm = TRUE) + c(-2, 2),
    ylim = range(ONMY_2_coords$lat, na.rm = TRUE) + c(-2, 2)
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    shape = "Metapopulation",
    title = "ONMY-2 steelhead sampling sites"
  )


ggplot(ONMY_2_ibd, aes(x = dist_km, y = fst)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Straight-line distance (km)",
    y = expression(F[ST]),
    title = "ONMY-2 isolation by distance"
  )



# -----------------------------
# save objects
# -----------------------------
save(
  ONMY_2_fst,
  ONMY_2_coords,
  #ONMY_2_coords,
  #ONMY_2_ibd,
  file = file.path(data_dir, "ONMY-2.RData")
)

# -----------------------------
# quick output
# -----------------------------
cat("\nFinished combined neutral extraction across all five sheets.\n")
cat("Total populations:", nrow(ONMY_2_coords), "\n")
cat("Total individuals:", sum(ONMY_2_coords$n_ind), "\n")
cat("Total loci:", ncol(ONMY_2_hf) - 1, "\n\n")

print(head(ONMY_2_coords, 10))
print(dim(ONMY_2_fst))
print(round(ONMY_2_fst[1:min(10, nrow(ONMY_2_fst)), 1:min(10, ncol(ONMY_2_fst))], 5))
