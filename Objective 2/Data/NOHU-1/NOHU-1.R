# -----------------------------
# NOHU-1 Olympic mudminnow
# Novumbra hubbsi
# pairwise FST + site coordinates + map + IBD plot
# -----------------------------

library(readxl)
library(ggplot2)
library(dplyr)
library(geosphere)
library(adegenet)
library(hierfstat)
library(maps)

# ------------------------------
# paths
# ------------------------------
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
raw_sheet <- readxl::read_excel(
  geno_file,
  sheet = "Mudminnow Genotypes",
  col_names = FALSE
)

hdr <- as.character(unlist(raw_sheet[3, ]))
mud_raw <- raw_sheet[-c(1:3), seq_along(hdr)]
names(mud_raw) <- hdr

mud_raw <- mud_raw |>
  dplyr::filter(!is.na(`Analysis ID`)) |>
  dplyr::mutate(
    `Analysis ID` = as.integer(`Analysis ID`),
    `Collection #` = as.character(`Collection #`),
    `Individual name` = as.character(`Individual name`),
    `Site Name` = as.character(`Site Name`)
  )

site_key <- readxl::read_excel(geno_file, sheet = "Site Summary") |>
  dplyr::rename(site_name = `Site Name`, site_id = `Analysis ID`) |>
  dplyr::mutate(site_id = as.integer(site_id))

# -----------------------------
# 2) best-available site coordinates
#
# explicit collection coordinates were not tabulated in the paper
# or genotype workbook. these XYs are best-available approximations
# assembled from:
#   - named-site coordinates when available from external sources
#   - mapped named features/drainages
#   - paper Figure 1 placement for sites lacking explicit XYs
#
# site order matches Table 1 in DeHaan et al. 2014
# -----------------------------
NOHU_1_coords <- data.frame(
  site = as.character(1:23),
  site_name = c(
    "Green Cove 1",
    "Green Cove 2",
    "Green Cove 3",
    "Woodard Creek",
    "Spurgeon Creek",
    "Hopkins Ditch 1",
    "Hopkins Ditch 2",
    "Hopkins Ditch 3",
    "South Hanaford Creek",
    "Adna Wetland",
    "Chehalis Oxbow Lake",
    "Satsop Slough",
    "Peoples Creek",
    "Cherry Creek",
    "E.F. Issaquah Creek",
    "Gillis Slough Pond",
    "Connor Creek",
    "Ditch along Hwy 109",
    "Upper Cook Creek",
    "N.F. Whale Creek",
    "Steamboat Creek Bog",
    "James Pond",
    "Lake Ozette Pond"
  ),
  lat = c(
    47.0620961,
    47.0620853,
    47.0614197,
    47.1164850,
    46.9512087,
    46.9383150,
    46.9421400,
    46.9459622,
    46.6951033,
    46.6381590,
    46.6565000,
    46.9735000,
    47.7939883,
    47.7684000,
    47.5362121,
    47.0424640,
    47.0911969,
    47.0675000,
    47.3672500,
    47.4898032,
    47.7200000,
    47.9194440,
    48.1144750
  ),
  lon = c(
    -122.9375614,
    -122.9378308,
    -122.9393707,
    -122.8631941,
    -122.8487461,
    -122.9396633,
    -122.9501570,
    -122.9606514,
    -122.8367914,
    -123.0192967,
    -123.0950000,
    -123.4870000,
    -121.9842903,
    -121.9603000,
    -122.0398415,
    -124.0465670,
    -124.1751764,
    -124.1110000,
    -124.0264170,
    -124.3424115,
    -124.3100000,
    -124.6000000,
    -124.6007890
  ),
  coord_basis = c(
    rep("named sampling-site coordinates from Green Cove WDFW report", 3),
    "named stream coordinate; pond placement approximate within Woodard Creek drainage",
    "named stream coordinate; site placement approximate within Spurgeon Creek drainage",
    "Hopkins Ditch upstream extent from WA DNR naming proposal",
    "midpoint interpolated along Hopkins Ditch between upstream extent and mouth",
    "Hopkins Ditch mouth at Salmon Creek / Jones Road from WA DNR naming proposal",
    "USGS/WQP monitoring coordinate for South Hanaford Creek near Kopiah",
    "named stream coordinate; wetland placement approximate near Stearns Creek / Adna",
    "paper Figure 1 placement in lower Chehalis floodplain; oxbow centroid approximate",
    "paper Figure 1 placement near lower Satsop / Chehalis floodplain; slough centroid approximate",
    "named stream coordinate for Peoples Creek",
    "lower Cherry Creek drainage approximation using mapped mouth coordinate",
    "named stream coordinate for East Fork Issaquah Creek",
    "named slough coordinate",
    "named stream coordinate",
    "paper Figure 1 placement for roadside ditch along SR 109 between Connor and Gillis",
    "named stream coordinate for Cook Creek / upper reach approximation",
    "named stream coordinate",
    "paper Figure 1 placement; bog centroid approximate in Steamboat Creek drainage",
    "James Pond near Mora Ranger Station",
    "Swan Bay / Lake Ozette beaver-pond wetland approximation"
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 3) build a clean genotype table
# -----------------------------
geno_cols <- grep("^N[hH]ub[0-9]+_[AB]$", names(mud_raw), value = TRUE)

# preserve locus order as they appear in workbook
loci <- unique(sub("_[AB]$", "", geno_cols))

collapse_alleles <- function(a1, a2) {
  ifelse(
    is.na(a1) | is.na(a2),
    NA_character_,
    sprintf("%03d%03d", as.integer(a1), as.integer(a2))
  )
}

geno_df <- as.data.frame(
  setNames(
    lapply(loci, function(loc) {
      collapse_alleles(
        mud_raw[[paste0(loc, "_A")]],
        mud_raw[[paste0(loc, "_B")]]
      )
    }),
    loci
  ),
  stringsAsFactors = FALSE
)

sample_names <- ifelse(
  !is.na(mud_raw$`Individual name`) & mud_raw$`Individual name` != "",
  mud_raw$`Individual name`,
  paste0(mud_raw$`Analysis ID`, "_", mud_raw$`Collection #`)
)

NOHU_1_gen <- adegenet::df2genind(
  X = geno_df,
  ploidy = 2,
  sep = "",
  ncode = 3,
  ind.names = sample_names,
  pop = factor(mud_raw$`Analysis ID`, levels = 1:23),
  NA.char = NA,
  type = "codom"
)

# -----------------------------
# 4) pairwise FST matrix
# -----------------------------
NOHU_1_hier <- hierfstat::genind2hierfstat(NOHU_1_gen)
NOHU_1_fst <- as.matrix(hierfstat::pairwise.WCfst(NOHU_1_hier))

NOHU_1_fst <- NOHU_1_fst[as.character(1:23), as.character(1:23)]
diag(NOHU_1_fst) <- 0
NOHU_1_fst[NOHU_1_fst < 0] <- 0
NOHU_1_fst[lower.tri(NOHU_1_fst)] <- t(NOHU_1_fst)[lower.tri(NOHU_1_fst)]

rownames(NOHU_1_fst) <- colnames(NOHU_1_fst) <- as.character(1:23)

# -----------------------------
# 5) geographic distance matrix
# straight-line distance in km
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
# 7) map of sites
# includes USA + Canada context
# -----------------------------
world_map <- ggplot2::map_data("world")
world_map_sub <- dplyr::filter(world_map, region %in% c("USA", "Canada"))

site_map <- ggplot() +
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
    nudge_y = 0.12,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(-125.2, -121.3),
    ylim = c(46.3, 48.7)
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "NOHU-1 sampling sites"
  )

print(site_map)

# -----------------------------
# 8) IBD plot
# -----------------------------
ibd_plot <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "NOHU-1 isolation by distance"
  )

print(ibd_plot)

# -----------------------------
# 9) save outputs
# -----------------------------
write.csv(
  NOHU_1_coords,
  file.path(out_dir, "NOHU-1_coords.csv"),
  row.names = FALSE
)

write.csv(
  NOHU_1_fst,
  file.path(out_dir, "NOHU-1_fst_matrix.csv"),
  row.names = TRUE
)

save(
  NOHU_1_fst,
  NOHU_1_coords,
  geo_dist_km,
  ibd_df,
  file = file.path(out_dir, "NOHU-1.RData")
)

# -----------------------------
# 10) checks
# -----------------------------
stopifnot(identical(rownames(NOHU_1_fst), NOHU_1_coords$site))
stopifnot(identical(colnames(NOHU_1_fst), NOHU_1_coords$site))
stopifnot(isTRUE(all.equal(NOHU_1_fst, t(NOHU_1_fst))))
