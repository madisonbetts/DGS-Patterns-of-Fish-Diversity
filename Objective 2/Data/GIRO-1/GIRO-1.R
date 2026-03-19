# --------------
# roundtail chub
# GIRO-1
# --------------

library(adegenet)
library(hierfstat)
library(dplyr)
library(sf)
library(ggplot2)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/GIRO-1/"

geno_path <- file.path(save_dir, "Gila_robusta_microsat_genotypes.csv")

# ------------------------------
# read genotypes
# ------------------------------
geno <- read.csv(geno_path, stringsAsFactors = FALSE, check.names = FALSE)

# first col = sample IDs
sample_ids <- geno[[1]]

# everything else = allele columns
geno_loci <- geno[, -1]

# ------------------------------
# convert 0 -> NA
# ------------------------------
geno_loci[geno_loci == 0] <- NA
geno_loci[geno_loci == "0"] <- NA

# make sure all allele columns are character
geno_loci[] <- lapply(geno_loci, as.character)

# ------------------------------
# collapse paired allele columns
# e.g. 222_A1 + 222_A2 -> "148/150"
# ------------------------------
locus_names <- unique(sub("_(A1|A2)$", "", names(geno_loci)))

geno_joined <- data.frame(row.names = sample_ids)

for (loc in locus_names) {
  a1_col <- paste0(loc, "_A1")
  a2_col <- paste0(loc, "_A2")
  
  a1 <- geno_loci[[a1_col]]
  a2 <- geno_loci[[a2_col]]
  
  both_na <- is.na(a1) & is.na(a2)
  joined <- ifelse(both_na, NA, paste(a1, a2, sep = "/"))
  
  geno_joined[[loc]] <- joined
}

# inspect
head(geno_joined)

# ------------------------------
# convert to genind
# ------------------------------
gi <- df2genind(
  geno_joined,
  ploidy = 2,
  sep = "/",
  ind.names = sample_ids,
  NA.char = "NA"
)

# ------------------------------
# assign populations from sample IDs
# extracts leading letters before the first number
# works for 3-char and longer acronyms
# ------------------------------
pop(gi) <- factor(toupper(sub("^([A-Za-z]+).*", "\\1", sample_ids)))

table(pop(gi))

# ------------------------------
# hierfstat + pairwise FST
# ------------------------------
hf <- genind2hierfstat(gi)
GIRO_1_fst <- as.matrix(pairwise.WCfst(hf))
diag(GIRO_1_fst) <- 0

rownames(GIRO_1_fst) <- toupper(rownames(GIRO_1_fst))
colnames(GIRO_1_fst) <- toupper(colnames(GIRO_1_fst))

round(GIRO_1_fst, 4)

# ------------------------------
# population metadata
# coordinates are UTM from Table 1
# keep only site + coordinates
# ------------------------------
all_sites <- data.frame(
  site = c(
    "SIL", "BLU", "BON", "DIX", "EFE", "HCN", "LEG", "NMFKS", "TURNM", "UEG",
    "BLK", "CHR", "MAR", "ROC", "SPRSA", "TON",
    "CC", "SAB", "SHY",
    "ARA", "BAS", "ODN", "RDF", "TURAZ",
    "EVR", "FOS", "LSALT", "SPRVE", "VDP", "WAK", "WCL", "WVW",
    "BOL", "TRT"
  ),
  # group = c(
  #   "Agua Fria River basin",
  #   rep("Gila River sub-basin", 9),
  #   rep("Salt River sub-basin", 6),
  #   rep("Santa Cruz River sub-basin", 3),
  #   rep("San Pedro River sub-basin", 5),
  #   rep("Verde River sub-basin", 8),
  #   rep("Bill Williams River basin", 2)
  # ),
  # site_name = c(
  #   "Silver Creek, Yavapai Co., AZ",
  #   "Blue River, Gila Co., AZ",
  #   "Bonita Creek, Graham Co., AZ",
  #   "Dix Creek, Greenlee Co., AZ",
  #   "East Fork Eagle Creek, Greenlee Co., AZ",
  #   "Harden-Cienega Creek, Greenlee Co., AZ",
  #   "Eagle Creek, lower, Greenlee Co., AZ",
  #   "East, Middle and West Forks Gila River, Catron Co., NM",
  #   "Turkey Creek, Grant Co., NM",
  #   "Eagle Creek, upper, Greenlee Co., AZ",
  #   "Black River, Greenlee Co., AZ",
  #   "Cherry Creek, Gila Co., AZ",
  #   "Marsh Creek, Gila Co., AZ",
  #   "Rock Creek, Gila Co., AZ",
  #   "Spring Creek, Gila Co., AZ",
  #   "Tonto Creek, Gila Co., AZ",
  #   "Cienega Creek, Pima Co., AZ",
  #   "Sabino Canyon, Pima Co., AZ",
  #   "Sheehy Spring, Santa Cruz, AZ",
  #   "Aravaipa Creek, Pinal Co., AZ",
  #   "Bass Canyon, Cochise Co., AZ",
  #   "O'Donnell Canyon, Santa Cruz Co., AZ",
  #   "Redfield Canyon, Pima Co., AZ",
  #   "Turkey Creek, Santa Cruz Co., AZ",
  #   "East Verde River, Gila Co., AZ",
  #   "Fossil Creek, Yavapai Co., AZ",
  #   "Canal downstream from confluence Salt and Verde rivers, Maricopa Co., AZ",
  #   "Spring Creek, Yavapai Co., AZ",
  #   "Verde River, Perkinsville, Yavapai Co., AZ",
  #   "Walker Creek, Yavapai Co., AZ",
  #   "West Clear Creek, Yavapai Co., AZ",
  #   "Williamson Valley Wash, Yavapai Co., AZ",
  #   "Boulder Creek, Yavapai Co., AZ",
  #   "Trout Creek, Mohave Co., AZ"
  # ),
  utm_zone = "12S",
  easting = c(
    414762, 566895, 637253, 671832, 643185, 673853, 648752, 758619, 734293, 641212,
    639590, 515036, 497689, 489588, 495954, 491190,
    540260, 519663, 540028,
    551774, 571383, 544969, 562895, 546154,
    465504, 447363, 438441, 415889, 391039, 435967, 436195, 364924,
    302702, 275797
  ),
  northing = c(
    3793564, 3708078, 3647293, 3675085, 3707050, 3674850, 3651019, 3677745, 3662829, 3704936,
    3724720, 3740368, 3803361, 3759443, 3765606, 3786183,
    3524841, 3577319, 3470448,
    3640104, 3579679, 3491596, 3588722, 3498280,
    3794077, 3861552, 3712009, 3853475, 3862667, 3833635, 3822262, 3822260,
    3834064, 3875303
  ),
  # species = c(
  #   "intermedia", "intermedia", "intermedia", "intermedia", "robusta", "intermedia", "robusta", "nigra", "nigra", "intermedia",
  #   "robusta", "robusta", "nigra", "nigra", "intermedia", "nigra",
  #   "intermedia", "intermedia", "intermedia",
  #   "robusta", "intermedia", "intermedia", "intermedia", "intermedia",
  #   "nigra", "nigra", "robusta", "intermedia", "robusta", "robusta", "robusta", "intermedia",
  #   "robusta", "robusta"
  # ),
  stringsAsFactors = FALSE
)

# ------------------------------
# convert UTM (zone 12N) -> lon/lat
# ------------------------------
utm_sf <- st_as_sf(
  all_sites,
  coords = c("easting", "northing"),
  crs = 32612
)

lonlat_sf <- st_transform(utm_sf, crs = 4326)
coords_ll <- st_coordinates(lonlat_sf)

all_sites$lon <- coords_ll[, 1]
all_sites$lat <- coords_ll[, 2]

# keep only site + decimal degree coordinates
all_sites <- all_sites[, c("site", "lat", "lon")]

# ------------------------------
# match metadata to FST matrix order
# use only site for formatting
# ------------------------------
GIRO_1_coords <- all_sites[match(rownames(GIRO_1_fst), all_sites$site), ]

if (any(is.na(GIRO_1_coords$site))) {
  missing_sites <- rownames(GIRO_1_fst)[is.na(GIRO_1_coords$site)]
  stop("Missing metadata for these populations: ", paste(missing_sites, collapse = ", "))
}

# ------------------------------
# quick checks
# ------------------------------
stopifnot(identical(rownames(GIRO_1_fst), GIRO_1_coords$site))
stopifnot(identical(colnames(GIRO_1_fst), GIRO_1_coords$site))
stopifnot(isTRUE(all.equal(GIRO_1_fst, t(GIRO_1_fst))))

# ------------------------------
# map of sampling locations
# ------------------------------
ggplot(GIRO_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.12, size = 3.5) +
  coord_equal() +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "GIRO-1 sampling locations"
  )

# ------------------------------
# IBD plot
# straight-line distance vs FST
# ------------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(GIRO_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- GIRO_1_coords$site
colnames(geo_dist_km) <- GIRO_1_coords$site

ibd_df <- data.frame(
  site1 = rownames(GIRO_1_fst)[row(GIRO_1_fst)[upper.tri(GIRO_1_fst)]],
  site2 = colnames(GIRO_1_fst)[col(GIRO_1_fst)[upper.tri(GIRO_1_fst)]],
  fst = GIRO_1_fst[upper.tri(GIRO_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

ggplot(ibd_df, aes(dist_km, fst)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(
    x = "Straight-line distance (km)",
    y = expression(F[ST]),
    title = "GIRO-1 isolation by distance"
  )

# ------------------------------
# save RData
# ------------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  GIRO_1_fst,
  GIRO_1_coords,
  file = file.path(out_dir, "GIRO-1.RData")
)