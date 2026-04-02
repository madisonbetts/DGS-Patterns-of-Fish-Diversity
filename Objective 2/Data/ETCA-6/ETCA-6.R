# -----------------------------
# ETCA-6 rainbow darter
# Etheostoma caeruleum
# published pairwise FST + site coordinates
# with VCF-linked population assignments for QC / downstream work
# -----------------------------

library(vcfR)
library(dartR)
library(adegenet)
library(ggplot2)
library(ggrepel)
library(geosphere)
library(maps)
library(readxl)

# -----------------------------
# 0) paths
# -----------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETCA-6/doi_10_5061_dryad_4b8gtht9v__v20200928"
out_dir  <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETCA-6/data"

vcf_path <- file.path(base_dir, "VCF_RDarters.vcf")
ibd_path <- file.path(base_dir, "IBD_RDarters.xlsx")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# -----------------------------
# 1) read VCF and link individuals to sites
# -----------------------------
ETCA_6_gl <- vcfR::read.vcfR(vcf_path, verbose = FALSE) |>
  vcfR::vcfR2genlight()

# site codes are the alphabetic prefixes in sample names
# e.g. KFB1.1.sort -> KFB
ETCA_6_ind_sites <- sub("^([A-Z]+).*", "\\1", ETCA_6_gl@ind.names)

site_codes <- c("KFB", "UAC", "LGC", "UGC", "GLD", "UWC", "LWC", "MKR")

stopifnot(all(ETCA_6_ind_sites %in% site_codes))

ETCA_6_gl@pop <- factor(ETCA_6_ind_sites, levels = site_codes)
ETCA_6_gl@other$ind_sites <- ETCA_6_ind_sites

# quick check of individuals per site
ETCA_6_site_n <- table(ETCA_6_gl@pop)
print(ETCA_6_site_n)

# -----------------------------
# 2) site coordinates
# Table 1 from Oliveira et al. 2020
# numeric site order used for final matrix / coords df:
# 1 KFB
# 2 UAC
# 3 LGC
# 4 UGC
# 5 GLD
# 6 UWC
# 7 LWC
# 8 MKR
# -----------------------------
ETCA_6_sites <- data.frame(
  site = as.character(1:8),
  site_code = site_codes,
  site_name = c(
    "Kellogg Forest",
    "Upper Augusta Creek",
    "Lower Gull Creek",
    "Upper Gull Creek",
    "Gull Lake",
    "Upper Wabascon Creek",
    "Lower Wabascon Creek",
    "Kalamazoo River"
  ),
  lat = c(
    42.3633,
    42.4176,
    42.3686,
    42.4268,
    42.4063,
    42.3623,
    42.3541,
    42.3241
  ),
  lon = c(
    -85.3533,
    -85.3524,
    -85.4037,
    -85.4281,
    -85.4041,
    -85.2359,
    -85.2492,
    -85.3584
  ),
  stringsAsFactors = FALSE
)

# final coords df saved to RData
ETCA_6_coords <- ETCA_6_sites[, c("site", "lat", "lon")]

# lookup vector for converting site codes to numeric site ids
site_id_lookup <- setNames(ETCA_6_sites$site, ETCA_6_sites$site_code)

# -----------------------------
# 3) map
# include US + Canada, zoom to point extent
# plot only in RStudio; do not save
# -----------------------------
usa <- map_data("state")
can <- map_data("world", region = "Canada")

xpad <- 1.0
ypad <- 0.8

xlim_map <- c(
  min(ETCA_6_sites$lon) - xpad,
  max(ETCA_6_sites$lon) + xpad
)

ylim_map <- c(
  min(ETCA_6_sites$lat) - ypad,
  max(ETCA_6_sites$lat) + ypad
)

ggplot() +
  geom_polygon(
    data = usa,
    aes(x = long, y = lat, group = group),
    fill = "grey92", color = "black", linewidth = 0.25
  ) +
  geom_polygon(
    data = can,
    aes(x = long, y = lat, group = group),
    fill = "grey86", color = "black", linewidth = 0.25
  ) +
  geom_point(
    data = ETCA_6_sites,
    aes(x = lon, y = lat),
    size = 2.8
  ) +
  ggrepel::geom_text_repel(
    data = ETCA_6_sites,
    aes(x = lon, y = lat, label = paste0(site, ": ", site_code)),
    size = 3.3,
    seed = 1,
    min.segment.length = 0
  ) +
  coord_fixed(
    xlim = xlim_map,
    ylim = ylim_map
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "ETCA-6 sampling sites"
  )

# -----------------------------
# 4) published pairwise FST matrix
# source: IBD_RDarters.xlsx from Dryad
# use published pairwise FST for final saved matrix
# -----------------------------
ibd_df <- readxl::read_excel(ibd_path)

# split comparisons into site codes
pair_mat <- do.call(rbind, strsplit(ibd_df$Comparison, "-", fixed = TRUE))
ibd_df$site1_code <- pair_mat[, 1]
ibd_df$site2_code <- pair_mat[, 2]

stopifnot(all(ibd_df$site1_code %in% site_codes))
stopifnot(all(ibd_df$site2_code %in% site_codes))

# initialize symmetric matrix
ETCA_6_fst <- matrix(
  0,
  nrow = length(site_codes),
  ncol = length(site_codes),
  dimnames = list(ETCA_6_sites$site, ETCA_6_sites$site)
)

for (i in seq_len(nrow(ibd_df))) {
  s1 <- site_id_lookup[[ibd_df$site1_code[i]]]
  s2 <- site_id_lookup[[ibd_df$site2_code[i]]]
  fst_val <- as.numeric(ibd_df$Fst[i])

  ETCA_6_fst[s1, s2] <- fst_val
  ETCA_6_fst[s2, s1] <- fst_val
}

diag(ETCA_6_fst) <- 0
ETCA_6_fst[ETCA_6_fst < 0] <- 0

# -----------------------------
# 5) geographic distance matrix from site coordinates
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(ETCA_6_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- ETCA_6_coords$site
colnames(geo_dist_km) <- ETCA_6_coords$site

# -----------------------------
# 6) IBD plot
# use published FST + coordinate-based Euclidean distance
# -----------------------------
ETCA_6_ibd_df <- data.frame(
  site1   = rownames(ETCA_6_fst)[row(ETCA_6_fst)[upper.tri(ETCA_6_fst)]],
  site2   = colnames(ETCA_6_fst)[col(ETCA_6_fst)[upper.tri(ETCA_6_fst)]],
  fst     = ETCA_6_fst[upper.tri(ETCA_6_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  stringsAsFactors = FALSE
)

ggplot(ETCA_6_ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.8, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "ETCA-6 isolation by distance"
  )

# -----------------------------
# 7) optional QC:
# compare published FST to dartR-derived site-level FST
# this block does not affect saved ETCA_6_fst
# uncomment if you want to recalculate pairwise FST from the VCF
# -----------------------------
# ETCA_6_fst_calc <- dartR::gl.fst.pop(ETCA_6_gl, nboots = 0, plot.out = FALSE)
# if (is.list(ETCA_6_fst_calc) && "Fsts" %in% names(ETCA_6_fst_calc)) {
#   ETCA_6_fst_calc <- ETCA_6_fst_calc$Fsts
# }
# ETCA_6_fst_calc <- as.matrix(ETCA_6_fst_calc)
# rownames(ETCA_6_fst_calc) <- colnames(ETCA_6_fst_calc) <- site_codes
# ETCA_6_fst_calc[ETCA_6_fst_calc < 0] <- 0
# print(round(ETCA_6_fst_calc, 6))

# -----------------------------
# 8) quick checks
# -----------------------------
stopifnot(nrow(ETCA_6_coords) == 8)
stopifnot(identical(rownames(ETCA_6_fst), ETCA_6_coords$site))
stopifnot(identical(colnames(ETCA_6_fst), ETCA_6_coords$site))
stopifnot(isTRUE(all.equal(ETCA_6_fst, t(ETCA_6_fst))))
stopifnot(all(diag(ETCA_6_fst) == 0))

# -----------------------------
# 9) save RData
# -----------------------------
save(
  ETCA_6_fst,
  ETCA_6_coords,
  file = file.path(out_dir, "ETCA-6.RData")
)
