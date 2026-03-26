# -----------------------------
# NOHU-1 Olympic mudminnow
# Novumbra hubbsi
# Stack et al. 2014 / associated microsatellite workbook
# pairwise Weir & Cockerham FST + site coordinates + map + IBD plot
# -----------------------------

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(geosphere)
  library(adegenet)
  library(hierfstat)
  library(maps)
})

# -----------------------------
# paths
# -----------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/NOHU-1"
geno_file <- file.path(base_dir, "Mudminnow Genetic Data 8-26-13.xlsx")
out_dir   <- file.path(base_dir, "data")

if (!file.exists(geno_file)) {
  stop("Could not find 'Mudminnow Genetic Data 8-26-13.xlsx' in the NOHU-1 folder.")
}
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# -----------------------------
# 1) read genotype workbook
# -----------------------------
mud_raw <- readxl::read_excel(
  geno_file,
  sheet = "Mudminnow Genotypes",
  skip = 2
) %>%
  dplyr::filter(!is.na(`Analysis ID`)) %>%
  dplyr::mutate(
    `Analysis ID`    = as.integer(`Analysis ID`),
    `Collection #`   = as.character(`Collection #`),
    `Individual name` = as.character(`Individual name`),
    `Site Name`      = as.character(`Site Name`)
  )

site_key <- readxl::read_excel(
  geno_file,
  sheet = "Site Summary"
) %>%
  dplyr::rename(
    site_name = `Site Name`,
    site_id   = `Analysis ID`
  ) %>%
  dplyr::mutate(site_id = as.integer(site_id)) %>%
  dplyr::arrange(site_id)

# -----------------------------
# 2) site coordinates
# best-available approximate coordinates from site names / mapped localities
# paper + workbook do not provide a tabulated lat/lon table
# site numbering follows Analysis ID 1:23 exactly
# -----------------------------
NOHU_1_coords <- data.frame(
  site = as.character(1:23),
  site_name = c(
    "Green Cove 1",
    "Green Cove 2",
    "Green Cove 3",
    "Woodward Creek",
    "Spurgeon Creek",
    "Hopkins Ditch 1 (Salmon Creek 1)",
    "Hopkins Ditch 2 (Salmon Creek 2)",
    "Hopkins Ditch 3 (Salmon Creek 3)",
    "South Hanaford Cr",
    "Adna (Stearns) Wetland",
    "Chehalis Oxbow Lake",
    "Satsop Slough",
    "People's Creek",
    "Cherry Creek",
    "EF Issaquah Creek",
    "Gillis Slough Pond",
    "Connor Creek",
    "Ditch along Hwy 109",
    "Upper Cook Creek",
    "Whale Creek",
    "Steamboat Creek Bog",
    "James Pond",
    "Lake Ozette"
  ),
  lat = c(
    47.0620961, 47.0620853, 47.0614197, 47.1164850, 46.9512087,
    46.9383150, 46.9421400, 46.9459622, 46.6951033, 46.6381590,
    46.6565000, 46.9735000, 47.7939883, 47.7684000, 47.5362121,
    47.0424640, 47.0911969, 47.0675000, 47.3672500, 47.4898032,
    47.7200000, 47.9194440, 48.1144750
  ),
  lon = c(
    -122.9375614, -122.9378308, -122.9393707, -122.8631941, -122.8487461,
    -122.9396633, -122.9501570, -122.9606514, -122.8367914, -123.0192967,
    -123.0950000, -123.4870000, -121.9842903, -121.9603000, -122.0398415,
    -124.0465670, -124.1751764, -124.1110000, -124.0264170, -124.3424115,
    -124.3100000, -124.6000000, -124.6007890
  ),
  stringsAsFactors = FALSE
)

# enforce exact site-name order from workbook
NOHU_1_coords <- NOHU_1_coords %>%
  dplyr::arrange(as.integer(site))
NOHU_1_coords$site_name <- site_key$site_name[match(NOHU_1_coords$site, site_key$site_id)]

# -----------------------------
# 3) build genotype tables
# workbook stores diploid alleles in separate *_A and *_B columns
# for adegenet::df2genind, each locus must look like "201/205"
# for hierfstat::pairwise.WCfst, each locus must look like 201205
# -----------------------------
geno_cols <- grep("^N[hH]ub[0-9]+_[AB]$", names(mud_raw), value = TRUE)
loci <- unique(sub("_[AB]$", "", geno_cols))

sample_names <- ifelse(
  !is.na(mud_raw$`Individual name`) & mud_raw$`Individual name` != "",
  paste0(mud_raw$`Analysis ID`, "_", mud_raw$`Individual name`),
  paste0(mud_raw$`Analysis ID`, "_", mud_raw$`Collection #`, "_", seq_len(nrow(mud_raw)))
)

if (any(duplicated(sample_names))) {
  sample_names <- make.unique(sample_names)
}

collapse_for_genind <- function(a1, a2) {
  ifelse(
    is.na(a1) | is.na(a2),
    NA_character_,
    paste0(
      sprintf("%03d", as.integer(a1)),
      "/",
      sprintf("%03d", as.integer(a2))
    )
  )
}

collapse_for_hierfstat <- function(a1, a2) {
  out <- ifelse(
    is.na(a1) | is.na(a2),
    NA_character_,
    paste0(
      sprintf("%03d", as.integer(a1)),
      sprintf("%03d", as.integer(a2))
    )
  )
  as.integer(out)
}

geno_df <- as.data.frame(
  setNames(
    lapply(loci, function(loc) {
      collapse_for_genind(
        mud_raw[[paste0(loc, "_A")]],
        mud_raw[[paste0(loc, "_B")]]
      )
    }),
    loci
  ),
  stringsAsFactors = FALSE
)
rownames(geno_df) <- sample_names

NOHU_1_gen <- adegenet::df2genind(
  X = geno_df,
  ploidy = 2,
  sep = "/",
  ncode = 3,
  ind.names = sample_names,
  pop = factor(mud_raw$`Analysis ID`, levels = 1:23),
  NA.char = NA,
  type = "codom"
)

NOHU_1_hier <- data.frame(
  pop = mud_raw$`Analysis ID`,
  setNames(
    lapply(loci, function(loc) {
      collapse_for_hierfstat(
        mud_raw[[paste0(loc, "_A")]],
        mud_raw[[paste0(loc, "_B")]]
      )
    }),
    loci
  ),
  check.names = FALSE
)
rownames(NOHU_1_hier) <- sample_names

# -----------------------------
# 4) pairwise FST matrix
# -----------------------------
NOHU_1_fst <- hierfstat::pairwise.WCfst(NOHU_1_hier)
NOHU_1_fst <- as.matrix(NOHU_1_fst)

# reorder to exact site order 1:23
NOHU_1_fst <- NOHU_1_fst[as.character(1:23), as.character(1:23), drop = FALSE]

# clean matrix to project standard
diag(NOHU_1_fst) <- 0
NOHU_1_fst[NOHU_1_fst < 0] <- 0
NOHU_1_fst[lower.tri(NOHU_1_fst)] <- t(NOHU_1_fst)[lower.tri(NOHU_1_fst)]
rownames(NOHU_1_fst) <- as.character(1:23)
colnames(NOHU_1_fst) <- as.character(1:23)

# -----------------------------
# 5) geographic distance matrix
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(NOHU_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- NOHU_1_coords$site
colnames(geo_dist_km) <- NOHU_1_coords$site

# -----------------------------
# 6) pairwise dataframe for IBD
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(NOHU_1_fst)[row(NOHU_1_fst)[upper.tri(NOHU_1_fst)]],
  site2   = colnames(NOHU_1_fst)[col(NOHU_1_fst)[upper.tri(NOHU_1_fst)]],
  fst     = NOHU_1_fst[upper.tri(NOHU_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 7) sampling-site map
# -----------------------------
world_map <- ggplot2::map_data("world")
world_map_sub <- dplyr::filter(world_map, region %in% c("USA", "Canada"))

ggplot() +
  geom_polygon(
    data = world_map_sub,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = NOHU_1_coords,
    aes(x = lon, y = lat),
    size = 2.5
  ) +
  geom_text(
    data = NOHU_1_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.10,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(-125.3, -121.3),
    ylim = c(46.3, 48.8)
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "NOHU-1 sampling sites"
  )


# -----------------------------
# 8) IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "NOHU-1 isolation by distance"
  )



# -----------------------------
# 9) save outputs
# -----------------------------
save(
  NOHU_1_fst,
  NOHU_1_coords,
  #geo_dist_km,
  #ibd_df,
  file = file.path(out_dir, "NOHU-1.RData")
)
# -----------------------------
# 10) checks
# -----------------------------
stopifnot(identical(rownames(NOHU_1_fst), NOHU_1_coords$site))
stopifnot(identical(colnames(NOHU_1_fst), NOHU_1_coords$site))
stopifnot(isTRUE(all.equal(NOHU_1_fst, t(NOHU_1_fst))))
