# -----------------------------
# SACO-4 Bull trout (Salvelinus confluentus)
# Hargrove et al. 2025
# Snake River / Hells Canyon baseline
#
# Workflow:
# - read the FishGen GenePop export with assignPOP::read.Genepop()
# - link original POP blocks to the 32 pooled populations in Table 1 / Fig. 1
# - reconstruct diploid multilocus genotypes from the assignPOP allele-dosage matrix
# - compute pairwise Weir & Cockerham FST from the pooled SNP genotypes
# - force symmetry and set negative FST values to 0
# - create matching site coordinate and IBD objects
# -----------------------------

library(assignPOP)
library(adegenet)
library(hierfstat)
library(geosphere)
library(ggplot2)
library(maps)

# -----------------------------
# paths
# -----------------------------
genepop_file <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SACO-4/GenePop.gen.txt"
out_rdata    <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SACO-4/data/SACO-4.RData"

# -----------------------------
# 1) pooled populations in Table 1 / Fig. 1 order
# -----------------------------
SACO_4_coords <- data.frame(
  site = 1:32,
  #reporting_unit = c(
  #  rep("Grande Ronde R", 15),
  #  "Imnaha R",
  #  rep("Pine Cr", 9),
  #  rep("Salmon R", 3),
  #  "Indian Cr",
  #  rep("Wildhorse R", 3)
  #),
  #stream = c(
  #  "Butte Cr",
  #  "Wenaha R",
  #  "Lookingglass Cr",
  #  "Minam R",
  #  "Deer Cr",
  #  "Indian Cr",
  #  "Bear Cr",
  #  "Lostine R",
  #  "Hurricane Cr",
  #  "Little Minam R",
  #  "North Fork Catherine Cr",
  #  "Elk Cr",
  #  "Limber Jim Cr",
  #  "Indiana Cr",
  #  "Clear Cr",
  #  "Imnaha R",
  #  "Elk Cr",
  #  "Aspen Cr",
  #  "Clear Cr",
  #  "Meadow Cr",
  #  "Trail Cr",
  #  "East Pine Cr",
   # "Pine Cr",
    #"Middle Fork Pine Cr",
#    "East Fork Pine Cr",
#    "Rapid R",
#    "4th of July Cr",
#    "East Fork Salmon R",
#    "Indian Cr",
#    "Bear Cr",
#    "Wesley Cr",
#    "Crooked R"
 # ),
  lat = c(
    45.7900, 45.9510, 45.5430, 45.6150, 45.4780,
    45.3900, 45.4480, 45.2930, 45.2480, 45.3340,
    45.2550, 45.2140, 44.9950, 44.9880, 44.9720,
    45.1780,
    44.9920, 44.9700, 44.8640, 44.9190, 44.9810,
    44.9030, 44.8330, 44.9700, 44.8480,
    44.8860, 45.0710, 44.3530,
    45.0150, 44.9340, 45.0080, 44.8520
  ),
  lon = c(
    -117.6960, -117.7941, -117.9830, -117.7220, -117.5310,
    -117.8650, -117.6600, -117.3990, -117.3080, -117.5840,
    -117.6840, -117.3560, -118.1700, -118.2360, -118.1970,
    -116.9590,
    -116.7700, -116.7480, -116.8090, -116.8750, -116.9370,
    -116.8170, -116.8910, -116.9320, -116.7730,
    -116.3420, -114.8490, -114.7440,
    -116.5480, -116.5150, -116.4560, -116.4040
  ),
  stringsAsFactors = FALSE
)

stopifnot(nrow(SACO_4_coords) == 32)
stopifnot(identical(SACO_4_coords$site, 1:32))

# -----------------------------
# 2) map raw GenePop POP blocks to pooled populations
# -----------------------------
collection_map <- data.frame(
  raw_collection = c(
    "ScoBCWE95C",
    "ScoNFWE03C", "ScoNFWE04C", "ScoSFWE03C", "ScoSFWE95C", "ScoWENA20C",
    "ScoLKGC03C",
    "ScoMINA18C", "ScoMINA19C", "ScoMINA20C",
    "ScoDECW03C",
    "ScoINDG95C",
    "ScoBCWL95C",
    "ScoLSTN03C", "ScoLSTN18C", "ScoLSTN19C", "ScoLSTN20C", "ScoLSTN95C",
    "ScoLSTW18C", "ScoLSTW19C", "ScoLSTW20C",
    "ScoHURR03C", "ScoHURR95C",
    "ScoLMIN03C", "ScoLMIN95C",
    "ScoNCAT95C",
    "ScoELKM95C",
    "ScoLJIM03C", "ScoLJIM95C",
    "ScoICGR03C",
    "ScoCCGR03C", "ScoCCGR95C",
    "ScoIMNA11C", "ScoIMNA12C", "ScoIMNA12C_1",
    "ScoELKF15C",
    "ScoASPN15C",
    "MixCLCR13C",
    "ScoMCCC16C",
    "ScoTRCC16C",
    "ScoEPNC16C",
    "ScoPINC14C",
    "ScoMFPC14C",
    "ScoEFPC14C",
    "ScoRRHW04C",
    "Sco4JUL04C",
    "ScoEFSW08C", "ScoEFSW09C", "ScoEFSW10C", "ScoEFSW11C", "ScoEFSW12C", "ScoEFSW13C", "ScoEFSW14C",
    "MixINDN08C", "ScoINDN11C", "ScoINDN19C",
    "HybBEWR17C",
    "HybWSLY17C",
    "HybCRWH18C"
  ),
  site = c(
    1,
    2,2,2,2,2,
    3,
    4,4,4,
    5,
    6,
    7,
    8,8,8,8,8,
    8,8,8,
    9,9,
    10,10,
    11,
    12,
    13,13,
    14,
    15,15,
    16,16,16,
    17,
    18,
    19,
    20,
    21,
    22,
    23,
    24,
    25,
    26,
    27,
    28,28,28,28,28,28,28,
    29,29,29,
    30,
    31,
    32
  ),
  stringsAsFactors = FALSE
)

collection_map$population <- sprintf("%02d", collection_map$site)
collection_map$stream <- SACO_4_coords$stream[match(collection_map$site, SACO_4_coords$site)]
collection_map$reporting_unit <- SACO_4_coords$reporting_unit[match(collection_map$site, SACO_4_coords$site)]

stopifnot(length(unique(collection_map$raw_collection)) == nrow(collection_map))
stopifnot(all(sort(unique(collection_map$site)) == 1:32))

# -----------------------------
# 3) read GenePop and relabel individuals to pooled populations
# -----------------------------
stopifnot(file.exists(genepop_file))
SACO_4 <- assignPOP::read.Genepop(genepop_file)

raw_collection <- SACO_4$SampleID
raw_collection <- gsub('^"|"$', "", raw_collection)
raw_collection <- trimws(raw_collection)
raw_collection <- sub("-.*", "", raw_collection)

unknown_labels <- sort(setdiff(unique(raw_collection), collection_map$raw_collection))
if (length(unknown_labels) > 0) {
  stop(
    "These raw collections are missing from collection_map: ",
    paste(unknown_labels, collapse = ", ")
  )
}

raw_to_site <- setNames(collection_map$site, collection_map$raw_collection)
SACO_4_site_id <- unname(raw_to_site[raw_collection])
SACO_4_pop_lab <- sprintf("%02d", SACO_4_site_id)

# -----------------------------
# 4) reconstruct diploid genotypes from assignPOP allele-dosage matrix
# assignPOP returns one column per allele with values 0 / 0.5 / 1
# -----------------------------
dm <- SACO_4$DataMatrix
if (!is.data.frame(dm)) dm <- as.data.frame(dm)

# drop any first column that is a population label rather than a locus/allele column
allele_cols <- names(dm)
allele_cols_trim <- trimws(allele_cols)
is_allele_col <- grepl("_[0-9]+$", allele_cols_trim)

if (!all(is_allele_col)) {
  dm <- dm[, is_allele_col, drop = FALSE]
  allele_cols <- names(dm)
  allele_cols_trim <- trimws(allele_cols)
}

stopifnot(all(grepl("_[0-9]+$", allele_cols_trim)))

# clean individual names
ind_names <- gsub('^"|"$', "", SACO_4$SampleID)
ind_names <- make.unique(trimws(ind_names))
rownames(dm) <- ind_names

# locus names and allele codes
locus_id <- sub("_[0-9]+$", "", allele_cols_trim)
allele_code <- sub("^.*_([0-9]+)$", "\\1", allele_cols_trim)

locus_index <- split(seq_along(locus_id), locus_id)

collapse_locus <- function(mat, codes) {
  apply(mat, 1, function(v) {
    v <- suppressWarnings(as.numeric(v))
    v[is.na(v)] <- 0
    counts <- round(v * 2)
    counts[counts < 0] <- 0

    if (sum(counts) == 0) {
      return(NA_character_)
    }

    # try to coerce unexpected dosage totals back to a diploid call
    if (sum(counts) != 2) {
      if (sum(v) > 0) {
        scaled <- (v / sum(v)) * 2
        counts <- floor(scaled + 1e-8)
        remainder <- 2 - sum(counts)
        if (remainder > 0) {
          add_idx <- order(scaled - counts, decreasing = TRUE)[seq_len(remainder)]
          counts[add_idx] <- counts[add_idx] + 1
        }
      }
    }

    if (sum(counts) == 1) {
      counts[which.max(counts)] <- 2
    }

    if (sum(counts) != 2) {
      return(NA_character_)
    }

    geno <- rep(codes, counts)
    if (length(geno) != 2) {
      return(NA_character_)
    }

    geno <- sort(geno)
    paste0(geno, collapse = "")
  })
}

geno_df <- as.data.frame(
  lapply(locus_index, function(idx) collapse_locus(as.matrix(dm[, idx, drop = FALSE]), allele_code[idx])),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# preserve original locus order from assignPOP
locus_order <- unique(locus_id)
geno_df <- geno_df[, locus_order, drop = FALSE]
rownames(geno_df) <- ind_names

# keep only loci with at least some called genotypes
keep_loci <- vapply(geno_df, function(x) any(!is.na(x)), logical(1))
geno_df <- geno_df[, keep_loci, drop = FALSE]

# -----------------------------
# 5) convert to genind and compute pairwise WC FST
# use unique numeric site IDs as population labels to avoid duplicate names
# -----------------------------
SACO_4_gen <- adegenet::df2genind(
  X = geno_df,
  ploidy = 2,
  ncode = 3,
  ind.names = rownames(geno_df),
  pop = factor(SACO_4_pop_lab, levels = sprintf("%02d", SACO_4_coords$site)),
  NA.char = NA,
  type = "codom"
)

SACO_4_hf <- hierfstat::genind2hierfstat(SACO_4_gen)
SACO_4_fst <- hierfstat::pairwise.WCfst(SACO_4_hf)

ord <- sprintf("%02d", SACO_4_coords$site)
SACO_4_fst <- SACO_4_fst[ord, ord]
rownames(SACO_4_fst) <- colnames(SACO_4_fst) <- as.character(SACO_4_coords$site)

SACO_4_fst[lower.tri(SACO_4_fst)] <- t(SACO_4_fst)[lower.tri(SACO_4_fst)]
SACO_4_fst[SACO_4_fst < 0] <- 0
diag(SACO_4_fst) <- 0

stopifnot(identical(dim(SACO_4_fst), c(32L, 32L)))
stopifnot(isTRUE(all.equal(SACO_4_fst, t(SACO_4_fst), check.attributes = FALSE)))

# -----------------------------
# 6) matching coordinate dataframe
# -----------------------------
SACO_4_coords <- SACO_4_coords[, c("site", "lat", "lon")]
stopifnot(identical(as.character(SACO_4_coords$site), rownames(SACO_4_fst)))

# -----------------------------
# 7) geographic distances (km)
# -----------------------------
coords <- SACO_4_coords[, c("lon", "lat")]
geo_dist_km <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
rownames(geo_dist_km) <- colnames(geo_dist_km) <- SACO_4_coords$site

# -----------------------------
# 8) IBD dataframe
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(SACO_4_fst)[row(SACO_4_fst)[upper.tri(SACO_4_fst)]],
  site2   = colnames(SACO_4_fst)[col(SACO_4_fst)[upper.tri(SACO_4_fst)]],
  fst     = SACO_4_fst[upper.tri(SACO_4_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 9) map of sampling sites
# -----------------------------
world_map <- map_data("world")

map_df <- SACO_4_coords
map_df$label <- as.character(map_df$site)

p_map <- ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "grey92",
    color = "grey55",
    linewidth = 0.2
  ) +
  geom_point(
    data = map_df,
    aes(x = lon, y = lat),
    color = "red3",
    size = 3
  ) +
  geom_text(
    data = map_df,
    aes(x = lon, y = lat, label = label),
    nudge_y = 0.15,
    size = 3
  ) +
  coord_fixed(
    ratio = 1.3,
    xlim = range(map_df$lon) + c(-1.0, 1.0),
    ylim = range(map_df$lat) + c(-0.7, 0.7)
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "SACO-4 sampling sites"
  )

# -----------------------------
# 10) IBD plot
# -----------------------------
p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Geographic distance (km)",
    y = expression(F[ST]),
    title = "SACO-4 isolation by distance"
  )

# -----------------------------
# 11) save RData
# -----------------------------
save(
  SACO_4_fst,
  SACO_4_coords,
  file = out_rdata
)

# -----------------------------
# 12) print outputs
# -----------------------------
print(SACO_4_coords)
print(round(SACO_4_fst, 4))
print(head(ibd_df))
print(p_map)
print(p_ibd)
