# -----------------------------
# COAS-1 prickly sculpin (Cottus asper)
# microsatellite FST matrix + site coordinates + map + IBD plot
# -----------------------------

library(ggplot2)
library(dplyr)
library(geosphere)
library(adegenet)
library(hierfstat)
library(maps)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/COAS-1"
geno_file <- file.path(save_dir, "Final_Genotypes_Prickly.csv")

if (!file.exists(geno_file)) {
  stop("Could not find Final_Genotypes_Prickly.csv in COAS-1 folder.")
}

# -----------------------------
# 1) read + parse the genotype file
# file is structured as Genepop-like blocks:
#   row 1 title
#   row 2 npops
#   row 3 nloci
#   row 4 locus names in alternating columns
#   then repeated blocks beginning with 'pop = <site>'
# -----------------------------
raw <- read.csv(
  geno_file,
  header = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "", "?")
)

# locus names are on row 4 (index 4 in R), odd-numbered genotype columns
loci <- raw[4, seq(2, ncol(raw), by = 2)] |> unlist(use.names = FALSE)
loci <- loci[!is.na(loci)]
loci <- trimws(loci)

pop_start <- grep("^pop\\s*=", raw[[1]], ignore.case = TRUE)
if (length(pop_start) == 0) {
  stop("Could not find any 'pop =' blocks in Final_Genotypes_Prickly.csv.")
}

parse_pop_block <- function(start_row, end_row, site_name, loci, raw_df) {
  block <- raw_df[(start_row + 1):(end_row - 1), , drop = FALSE]
  block <- block[!is.na(block[[1]]) & block[[1]] != "", , drop = FALSE]

  if (nrow(block) == 0) return(NULL)

  out <- data.frame(
    sample_id = trimws(block[[1]]),
    site_name = site_name,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(loci)) {
    c1 <- 2 * i
    c2 <- c1 + 1
    out[[paste0(loci[i], "_1")]] <- suppressWarnings(as.integer(block[[c1]]))
    out[[paste0(loci[i], "_2")]] <- suppressWarnings(as.integer(block[[c2]]))
  }

  out
}

pop_names <- trimws(sub("^pop\\s*=\\s*", "", raw[pop_start, 1], ignore.case = TRUE))
pop_end <- c(pop_start[-1], nrow(raw) + 1)

coas1_raw <- do.call(
  rbind,
  lapply(seq_along(pop_start), function(i) {
    parse_pop_block(
      start_row = pop_start[i],
      end_row   = pop_end[i],
      site_name = pop_names[i],
      loci      = loci,
      raw_df    = raw
    )
  })
)

if (is.null(coas1_raw) || nrow(coas1_raw) == 0) {
  stop("No genotype rows were parsed from Final_Genotypes_Prickly.csv.")
}

# duplicated sample IDs exist in the source file (e.g. 16_019, 24_057)
# so create unique internal individual IDs while preserving original IDs
coas1_raw$sample_id_orig <- coas1_raw$sample_id
coas1_raw$sample_id <- make.unique(coas1_raw$sample_id, sep = "_dup")


# preserve file order as site order
site_levels <- unique(coas1_raw$site_name)
coas1_raw$site_name <- factor(coas1_raw$site_name, levels = site_levels)
site_key <- data.frame(
  site = as.character(seq_along(site_levels)),
  site_name = site_levels,
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) site coordinates
# best-available approximate site coordinates assembled from the paper map,
# named-feature localities, and river/estuary reference points.
# these are treated as approximate XYs for Objective 2 IBD use.
# -----------------------------
COAS_1_coords <- data.frame(
  site = as.character(1:28),
  site_name = c(
    "Winchuck_Lagoon",
    "Smith",
    "Wilson",
    "Redwood",
    "Big _Lagoon",
    "Mad",
    "Freshwater",
    "Van_Duzen",
    "Cottoneva",
    "Navaro",
    "Garcia",
    "Russian_AC",
    "Estero_Americano",
    "Santa_Clara",
    "Napa",
    "SuisunBay",
    "Sacramento_Red_Bluff",
    "Sacramento_RM205.5",
    "ClearL",
    "Bear",
    "Feather",
    "American",
    "Mokelumne",
    "Calaveras",
    "Stanislaus",
    "Merced",
    "San_Joaquin",
    "Kings"
  ),
  lat = c(
    42.0192,
    41.9280,
    41.7730,
    41.3010,
    41.1614,
    40.9150,
    40.7868,
    40.5370,
    39.7360,
    39.1921,
    38.8780,
    38.5150,
    38.2955,
    34.2789,
    38.3700,
    38.0700,
    40.1544,
    39.7898,
    39.0340,
    39.0090,
    38.7810,
    38.5974,
    38.0960,
    37.9669,
    37.7900,
    37.3020,
    36.8220,
    36.0620
  ),
  lon = c(
    -124.1056,
    -124.1470,
    -124.2130,
    -124.0420,
    -124.1158,
    -124.1090,
    -124.0970,
    -124.1480,
    -123.8292,
    -123.7611,
    -123.6310,
    -123.0580,
    -123.0025,
    -119.1421,
    -122.3000,
    -122.0700,
    -122.2020,
    -122.0500,
    -122.8070,
    -122.7430,
    -121.6180,
    -121.5080,
    -121.5700,
    -121.3677,
    -120.7600,
    -120.4750,
    -120.9650,
    -119.6030
  ),
  stringsAsFactors = FALSE
)

# reorder to exactly match genotype file order
COAS_1_coords <- COAS_1_coords[match(site_key$site_name, COAS_1_coords$site_name), ]
COAS_1_coords$site <- as.character(seq_len(nrow(COAS_1_coords)))

if (any(is.na(COAS_1_coords$site_name))) {
  stop("Some COAS-1 site coordinates did not match parsed genotype site names.")
}

# -----------------------------
# 3) map of sampling locations
# plotted with USA + Canada for context
# -----------------------------
map_df <- map_data("world")
map_df <- subset(map_df, region %in% c("USA", "Canada"))

p_sites <- ggplot() +
  geom_polygon(
    data = map_df,
    aes(x = long, y = lat, group = group),
    fill = "grey92",
    color = "grey65",
    linewidth = 0.2
  ) +
  geom_point(
    data = COAS_1_coords,
    aes(x = lon, y = lat),
    size = 2.6
  ) +
  geom_text(
    data = COAS_1_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.22,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(-126, -118),
    ylim = c(33, 43)
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "COAS-1 sampling sites"
  )

print(p_sites)

# -----------------------------
# 4) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(COAS_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- COAS_1_coords$site
colnames(geo_dist_km) <- COAS_1_coords$site

# -----------------------------
# 5) build genind object correctly from paired allele columns
# adegenet wants each locus as e.g. 272/281 when ncode = 3 and sep = "/"
# -----------------------------
collapse_for_genind <- function(a1, a2) {
  ifelse(
    is.na(a1) | is.na(a2),
    NA_character_,
    paste0(sprintf("%03d", as.integer(a1)), "/", sprintf("%03d", as.integer(a2)))
  )
}

geno_df <- as.data.frame(
  setNames(
    lapply(loci, function(loc) {
      collapse_for_genind(
        coas1_raw[[paste0(loc, "_1")]],
        coas1_raw[[paste0(loc, "_2")]]
      )
    }),
    loci
  ),
  stringsAsFactors = FALSE
)

COAS_1_gen <- adegenet::df2genind(
  X = geno_df,
  ploidy = 2,
  sep = "/",
  ncode = 3,
  ind.names = coas1_raw$sample_id,
  pop = factor(match(coas1_raw$site_name, site_key$site_name), levels = seq_len(nrow(site_key))),
  NA.char = NA,
  type = "codom"
)

# -----------------------------
# 6) build hierfstat dataframe directly from paired allele columns
# avoids downstream rowname issues in genind2hierfstat for this file format
# each diploid genotype is stored as 6-digit allele codes, e.g. 272281
# -----------------------------
collapse_for_hierfstat <- function(a1, a2) {
  out <- ifelse(
    is.na(a1) | is.na(a2),
    NA_character_,
    paste0(sprintf("%03d", as.integer(a1)), sprintf("%03d", as.integer(a2)))
  )
  as.integer(out)
}

COAS_1_hier <- data.frame(
  pop = match(coas1_raw$site_name, site_key$site_name),
  setNames(
    lapply(loci, function(loc) {
      collapse_for_hierfstat(
        coas1_raw[[paste0(loc, "_1")]],
        coas1_raw[[paste0(loc, "_2")]]
      )
    }),
    loci
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rownames(COAS_1_hier) <- coas1_raw$sample_id

# -----------------------------
# 7) pairwise FST matrix
# Weir & Cockerham FST from microsatellite genotypes
# -----------------------------
COAS_1_fst <- as.matrix(hierfstat::pairwise.WCfst(COAS_1_hier))

# enforce site order 1:28 and symmetry
site_ids <- as.character(seq_len(nrow(COAS_1_coords)))
COAS_1_fst <- COAS_1_fst[site_ids, site_ids]
diag(COAS_1_fst) <- 0
COAS_1_fst[COAS_1_fst < 0] <- 0
COAS_1_fst <- (COAS_1_fst + t(COAS_1_fst)) / 2

rownames(COAS_1_fst) <- colnames(COAS_1_fst) <- site_ids

# -----------------------------
# 8) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(COAS_1_fst)[row(COAS_1_fst)[upper.tri(COAS_1_fst)]],
  site2   = colnames(COAS_1_fst)[col(COAS_1_fst)[upper.tri(COAS_1_fst)]],
  fst     = COAS_1_fst[upper.tri(COAS_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 9) IBD plot
# -----------------------------
p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.6, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "COAS-1 isolation by distance"
  )

print(p_ibd)

# -----------------------------
# 10) quick checks
# -----------------------------
stopifnot(identical(rownames(COAS_1_fst), COAS_1_coords$site))
stopifnot(identical(colnames(COAS_1_fst), COAS_1_coords$site))
stopifnot(isTRUE(all.equal(COAS_1_fst, t(COAS_1_fst))))

# -----------------------------
# 11) save outputs
# -----------------------------
out_dir <- file.path(save_dir, "data")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

write.csv(
  COAS_1_coords,
  file.path(out_dir, "COAS-1_coords.csv"),
  row.names = FALSE
)

write.csv(
  COAS_1_fst,
  file.path(out_dir, "COAS-1_fst_matrix.csv"),
  row.names = TRUE
)

save(
  COAS_1_fst,
  COAS_1_coords,
  file = file.path(out_dir, "COAS-1.RData")
)
