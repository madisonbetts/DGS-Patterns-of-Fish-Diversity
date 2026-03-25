############################
## ONMY-3 Eastern Cascades steelhead
## Hand et al. 2016 (Molecular Ecology)
## Neutral SNPs only
##
## This script:
##   1) reads the raw genotype workbook
##   2) extracts the sheet: Cascades_MPG_genos_neutral
##   3) converts the two-column-per-locus genalex-style layout into
##      hierfstat-ready diploid genotype codes
##   4) calculates pairwise Weir & Cockerham FST among populations
##   5) forces symmetry and sets negative FST to 0
##   6) renames rows/cols with numeric site IDs
##   7) builds a matching site dataframe
##   8) plots sites and an IBD relationship using straight-line distances
############################

suppressPackageStartupMessages({
  library(readxl)
  library(hierfstat)
  library(geosphere)
  library(ggplot2)
  library(maps)
})

# -----------------------------
# file path + sheet
# -----------------------------
xlsx_file <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ONMY-2/doi_10_5061_dryad_pr065__v20151211/DataforDryad/SH_Raw_Genetic.xlsx"
sheet_name <- "Cascades_MPG_genos_neutral"

# -----------------------------
# read raw sheet (no header)
# -----------------------------
raw_xlsx <- readxl::read_excel(
  path = xlsx_file,
  sheet = sheet_name,
  col_names = FALSE
)
raw_xlsx <- as.data.frame(raw_xlsx, stringsAsFactors = FALSE)

# row 3 contains the usable header
header_row <- raw_xlsx[3, ]
ind_dat <- raw_xlsx[-c(1, 2, 3), , drop = FALSE]
rownames(ind_dat) <- NULL
names(ind_dat)[1:2] <- c("sample_id", "pop_code")

# -----------------------------
# parse loci (2 columns per locus)
# -----------------------------
locus_names <- as.character(unlist(header_row[3:ncol(raw_xlsx)]))
locus_names <- trimws(locus_names)

if (length(locus_names) %% 2 != 0) {
  stop("Expected an even number of genotype columns after sample_id and pop_code.")
}

locus_pairs <- split(seq_along(locus_names), ceiling(seq_along(locus_names) / 2))
paired_names <- vapply(locus_pairs, function(idx) locus_names[idx[1]], character(1))

if (!all(vapply(locus_pairs, function(idx) length(unique(locus_names[idx])) == 1, logical(1)))) {
  stop("Not all genotype columns occur as clean allele pairs with the same locus name.")
}

# sanitize locus names for safe dataframe handling
paired_names <- gsub("[[:space:]]+", "_", paired_names)
paired_names <- gsub("[^A-Za-z0-9_]", "_", paired_names)
paired_names <- make.unique(paired_names, sep = "_")

fmt_allele <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "<NA>", "NULL")] <- NA_character_
  x
}

# Convert two 3-digit alleles into a single 6-digit diploid genotype code
# expected by hierfstat, e.g. 101 and 103 -> 101103
make_hf_genotype <- function(a1, a2) {
  a1 <- fmt_allele(a1)
  a2 <- fmt_allele(a2)

  bad1 <- !is.na(a1) & !grepl("^[0-9]+$", a1)
  bad2 <- !is.na(a2) & !grepl("^[0-9]+$", a2)
  if (any(bad1 | bad2)) {
    stop("Non-numeric allele codes detected in the raw genotype sheet.")
  }

  out <- ifelse(
    is.na(a1) | is.na(a2),
    NA_integer_,
    as.integer(sprintf("%03d%03d", as.integer(a1), as.integer(a2)))
  )
  out
}

# -----------------------------
# population order present in the xlsx sheet
# -----------------------------
pop_codes <- unique(ind_dat$pop_code)

# -----------------------------
# site lookup
# -----------------------------
site_lookup <- data.frame(
  pop_code = pop_codes,
  site_id  = seq_along(pop_codes),
  site_name = c(
    "Shitike Cr.",
    "Buckhollow Cr.",
    "Trout Cr.",
    "Fifteen Cr.",
    "Lower Little Klickitat R.",
    "Lower Summit Cr.",
    "Upper Trout Cr.",
    "Deadcanyon Cr.",
    "Lower White Cr.",
    "Snyder Cr.",
    "Surveyor Cr.",
    "Swale Cr.",
    "Rock Cr.",
    "Squaw Cr."
  ),
  tributary_region = c(
    "Deschutes",
    "Deschutes",
    "Deschutes",
    "Fifteen",
    "Little Klickitat",
    "Klickitat",
    "Klickitat",
    "Klickitat",
    "Klickitat",
    "Klickitat",
    "Klickitat",
    "Klickitat",
    "Rock",
    "Rock"
  ),
  lat = c(
    44.76178,
    45.12350,
    45.26900,
    45.60300,
    45.84279,
    45.99300,
    45.99100,
    45.94234,
    45.97600,
    45.95400,
    45.96200,
    45.91900,
    45.83700,
    45.86500
  ),
  lon = c(
    -121.25697,
    -120.85800,
    -120.83300,
    -121.07100,
    -121.05917,
    -121.09100,
    -121.28600,
    -121.14146,
    -121.23200,
    -121.18800,
    -121.39500,
    -121.01200,
    -120.88400,
    -120.82200
  ),
  stringsAsFactors = FALSE
)

if (nrow(site_lookup) != length(pop_codes)) {
  stop("site_lookup row count does not match the number of populations in the xlsx sheet.")
}

# -----------------------------
# build hierfstat input directly from the raw alleles
# -----------------------------
ONMY_2_hf <- data.frame(
  pop = match(ind_dat$pop_code, site_lookup$pop_code),
  stringsAsFactors = FALSE
)

for (i in seq_along(locus_pairs)) {
  cols <- locus_pairs[[i]] + 2L
  ONMY_2_hf[[paired_names[i]]] <- make_hf_genotype(ind_dat[[cols[1]]], ind_dat[[cols[2]]])
}

# drop loci that are entirely missing
keep_loci <- colSums(!is.na(ONMY_2_hf[, -1, drop = FALSE])) > 0
ONMY_2_hf <- ONMY_2_hf[, c(TRUE, keep_loci), drop = FALSE]

# -----------------------------
# pairwise FST
# -----------------------------
ONMY_2_fst <- hierfstat::pairwise.WCfst(ONMY_2_hf)
ONMY_2_fst <- as.matrix(ONMY_2_fst)

# force symmetry + cleanup
ONMY_2_fst[lower.tri(ONMY_2_fst)] <- t(ONMY_2_fst)[lower.tri(ONMY_2_fst)]
diag(ONMY_2_fst) <- 0
ONMY_2_fst[ONMY_2_fst < 0] <- 0

# reorder to site lookup / numeric ids
ONMY_2_fst <- ONMY_2_fst[as.character(site_lookup$site_id), as.character(site_lookup$site_id)]
rownames(ONMY_2_fst) <- site_lookup$site_id
colnames(ONMY_2_fst) <- site_lookup$site_id

# -----------------------------
# coords dataframe
# -----------------------------
ONMY_2_coords <- site_lookup[, c("site_id", "site_name", "tributary_region", "lat", "lon")]

# -----------------------------
# straight-line distance matrix (km)
# -----------------------------
xy <- as.matrix(ONMY_2_coords[, c("lon", "lat")])
ONMY_2_geo_km <- geosphere::distm(xy, fun = geosphere::distHaversine) / 1000
rownames(ONMY_2_geo_km) <- ONMY_2_coords$site_id
colnames(ONMY_2_geo_km) <- ONMY_2_coords$site_id

# -----------------------------
# pairwise dataframe for IBD
# -----------------------------
ONMY_2_ibd <- data.frame(
  site1   = rownames(ONMY_2_fst)[row(ONMY_2_fst)[upper.tri(ONMY_2_fst)]],
  site2   = colnames(ONMY_2_fst)[col(ONMY_2_fst)[upper.tri(ONMY_2_fst)]],
  fst     = ONMY_2_fst[upper.tri(ONMY_2_fst)],
  dist_km = ONMY_2_geo_km[upper.tri(ONMY_2_geo_km)],
  stringsAsFactors = FALSE
)

# -----------------------------
# save data objects
# -----------------------------
save(
  ONMY_2_fst,
  ONMY_2_coords,
  file =  "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ONMY-2/data/ONMY-2.RData"
)

# -----------------------------
# plot sites
# -----------------------------
world_map <- map_data("world")

ggplot() +
  geom_polygon(
    data = subset(world_map, region %in% c("USA", "Canada")),
    aes(x = long, y = lat, group = group),
    fill = "grey95", color = "grey65", linewidth = 0.25
  ) +
  geom_point(
    data = ONMY_2_coords,
    aes(x = lon, y = lat),
    size = 2.4, alpha = 0.95
  ) +
  geom_text(
    data = ONMY_2_coords,
    aes(x = lon, y = lat, label = site_id),
    nudge_y = 0.08, size = 3
  ) +
  coord_quickmap(
    xlim = c(min(ONMY_2_coords$lon) - 0.6, max(ONMY_2_coords$lon) + 0.6),
    ylim = c(min(ONMY_2_coords$lat) - 0.07, max(ONMY_2_coords$lat) + 0.07)
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "ONMY-3 Eastern Cascades steelhead sampling sites"
  )

# -----------------------------
# IBD plot (straight-line distance)
# -----------------------------
ggplot(ONMY_2_ibd, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.4, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  theme_classic() +
  labs(
    x = "Straight-line distance (km)",
    y = expression(F[ST]),
    title = "ONMY-3 isolation by distance"
  )
