# -----------------------------
# SACO-4 Bull trout (Salvelinus confluentus)
# Hargrove et al. 2025
# Snake River / Hells Canyon baseline
#
# Workflow:
# - read the FishGen GenePop export
# - link original POP blocks to the 32 pooled populations in Table 1 / Fig. 1
# - use feature-anchor coordinates for each pooled population
# - compute pairwise Weir & Cockerham FST from the raw SNP genotypes
# - force symmetry and set negative FST values to 0
# - create a matching site coordinate dataframe and IBD objects
#
# NOTE:
# Coordinates below are feature-anchor coordinates tied to the named stream /
# river location shown in Fig. 1 and Table 1. They are intended for Objective 2
# mapping + IBD workflows and can be refined later if a higher-precision site
# source becomes available.
# -----------------------------





SACO_4 = read.Genepop("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SACO-4/GenePop.gen.txt")



library(adegenet)
library(HWxtest)
library(hierfstat)
library(geosphere)
library(ggplot2)
library(maps)

# -----------------------------
# paths
# -----------------------------
genepop_file <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SACO-4/GenePop.gen"
out_rdata    <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SACO-4/data/SACO-4.RData"

# -----------------------------
# 1) pooled populations in Table 1 / Fig. 1 order
# -----------------------------
SACO_4_sites <- data.frame(
  site = 1:32,
  reporting_unit = c(
    rep("Grande Ronde R", 15),
    "Imnaha R",
    rep("Pine Cr", 9),
    rep("Salmon R", 3),
    "Indian Cr",
    rep("Wildhorse R", 3)
  ),
  stream = c(
    "Butte Cr",
    "Wenaha R",
    "Lookingglass Cr",
    "Minam R",
    "Deer Cr",
    "Indian Cr",
    "Bear Cr",
    "Lostine R",
    "Hurricane Cr",
    "Little Minam R",
    "North Fork Catherine Cr",
    "Elk Cr",
    "Limber Jim Cr",
    "Indiana Cr",
    "Clear Cr",
    "Imnaha R",
    "Elk Cr",
    "Aspen Cr",
    "Clear Cr",
    "Meadow Cr",
    "Trail Cr",
    "East Pine Cr",
    "Pine Cr",
    "Middle Fork Pine Cr",
    "East Fork Pine Cr",
    "Rapid R",
    "4th of July Cr",
    "East Fork Salmon R",
    "Indian Cr",
    "Bear Cr",
    "Wesley Cr",
    "Crooked R"
  ),
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

stopifnot(nrow(SACO_4_sites) == 32)
stopifnot(identical(SACO_4_sites$site, 1:32))

# -----------------------------
# 2) map raw GenePop POP blocks to pooled populations
# pooled exactly as reported in Table 1
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

collection_map$population <- SACO_4_sites$stream[match(collection_map$site, SACO_4_sites$site)]
collection_map$reporting_unit <- SACO_4_sites$reporting_unit[match(collection_map$site, SACO_4_sites$site)]

stopifnot(length(unique(collection_map$raw_collection)) == nrow(collection_map))
stopifnot(all(sort(unique(collection_map$site)) == 1:32))

# -----------------------------
# 3) read GenePop and relabel individuals to pooled populations
# HWxtest::genepop.to.genind works around the .gen-only expectation
# in adegenet::read.genepop for GenePop files saved as .txt
# -----------------------------
stopifnot(file.exists(genepop_file))
SACO_4_gen <- HWxtest::genepop.to.genind(genepop_file, ncode = 3)

raw_collection <- sub("-.*", "", adegenet::indNames(SACO_4_gen))

stopifnot(all(raw_collection %in% collection_map$raw_collection))

raw_to_site <- setNames(collection_map$site, collection_map$raw_collection)
raw_to_pop  <- setNames(collection_map$population, collection_map$raw_collection)

SACO_4_site_id <- raw_to_site[raw_collection]
SACO_4_pop_lab <- raw_to_pop[raw_collection]

adegenet::pop(SACO_4_gen) <- factor(SACO_4_pop_lab, levels = SACO_4_sites$stream)

# -----------------------------
# 4) pairwise Weir & Cockerham FST from pooled populations
# -----------------------------
SACO_4_hf  <- hierfstat::genind2hierfstat(SACO_4_gen)
SACO_4_fst_named <- hierfstat::pairwise.WCfst(SACO_4_hf)

# reorder to Table 1 / Fig. 1 order and rename to numeric site IDs
ord <- SACO_4_sites$stream
SACO_4_fst <- SACO_4_fst_named[ord, ord]
rownames(SACO_4_fst) <- colnames(SACO_4_fst) <- as.character(SACO_4_sites$site)

# force symmetry + clean
SACO_4_fst[lower.tri(SACO_4_fst)] <- t(SACO_4_fst)[lower.tri(SACO_4_fst)]
SACO_4_fst[SACO_4_fst < 0] <- 0
diag(SACO_4_fst) <- 0

stopifnot(identical(dim(SACO_4_fst), c(32L, 32L)))
stopifnot(isTRUE(all.equal(SACO_4_fst, t(SACO_4_fst), check.attributes = FALSE)))

# -----------------------------
# 5) matching coordinate dataframe
# keep only numeric site ID + lat/lon for saved workflow object
# -----------------------------
SACO_4_coords <- SACO_4_sites[, c("site", "lat", "lon")]

stopifnot(identical(as.character(SACO_4_coords$site), rownames(SACO_4_fst)))

# -----------------------------
# 6) geographic distances (km)
# -----------------------------
coords <- SACO_4_coords[, c("lon", "lat")]
geo_dist_km <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
rownames(geo_dist_km) <- colnames(geo_dist_km) <- SACO_4_coords$site

# -----------------------------
# 7) IBD dataframe
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(SACO_4_fst)[row(SACO_4_fst)[upper.tri(SACO_4_fst)]],
  site2   = colnames(SACO_4_fst)[col(SACO_4_fst)[upper.tri(SACO_4_fst)]],
  fst     = SACO_4_fst[upper.tri(SACO_4_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 8) map of sampling sites
# -----------------------------
world_map <- map_data("world")

map_df <- SACO_4_sites
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
# 9) IBD plot
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
# 10) save RData
# -----------------------------
save(
  SACO_4_fst,
  SACO_4_coords,
  file = out_rdata
)

# -----------------------------
# 11) print outputs
# -----------------------------
print(SACO_4_sites)
print(collection_map)
print(round(SACO_4_fst, 4))
print(head(ibd_df))
print(p_map)
print(p_ibd)
