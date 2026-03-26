############################
## SAFO-5 Brook Trout
## Stack et al. 2025/2026
##
## This script:
##   1) reads the uploaded pairwise FST matrix (.rds)
##   2) checks and reorders sites to the published code order
##   3) assigns numeric site IDs 1:22
##   4) builds SAFO_5_coords as site_id / lat / lon
##   5) plots site locations and straight-line IBD
##   6) saves SAFO_5_fst, SAFO_5_coords, SAFO_5_sites, SAFO_5_ibd
##
## IMPORTANT:
## The coordinates below are approximate site proxies inferred from
## the published site map and watershed landmarks in the upper Cache
## la Poudre / Long Draw region. They are suitable for a first-pass
## Euclidean IBD workflow, but they are not exact field GPS points.
############################

suppressPackageStartupMessages({
  library(geosphere)
  library(ggplot2)
  library(maps)
})

# -----------------------------
# paths
# -----------------------------
fst_path <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SAFO-5/doi_10_5061_dryad_nk98sf84b__v20250923/Stack_et_al_data/stack_et_al_Fst_matrix.rds"

out_dir  <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SAFO-5"
data_dir <- file.path(out_dir, "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# expected site order from the matrix
# -----------------------------
site_codes <- c(
  "UCO", "MCO", "LCO",
  "UW", "MW", "LW",
  "MUM", "HAZ1", "HAG3", "HAG2", "HAG1",
  "UCH", "MCH", "LCH",
  "UPOU", "PBC", "MPOU", "POU1",
  "ULPP", "LLPP", "CLPA", "CLPB"
)

# -----------------------------
# read fst matrix
# -----------------------------
fst_mat <- readRDS(fst_path)
fst_mat <- as.matrix(fst_mat)

if (!all(site_codes %in% rownames(fst_mat))) {
  stop("Not all expected site codes are present in rownames(fst_mat).")
}
if (!all(site_codes %in% colnames(fst_mat))) {
  stop("Not all expected site codes are present in colnames(fst_mat).")
}

fst_mat <- fst_mat[site_codes, site_codes]

diag(fst_mat) <- 0
fst_mat[fst_mat < 0] <- 0
fst_mat[lower.tri(fst_mat)] <- t(fst_mat)[lower.tri(fst_mat)]

# -----------------------------
# site table
# -----------------------------
SAFO_5_sites <- data.frame(
  site_id = seq_along(site_codes),
  site_code = site_codes,
  stringsAsFactors = FALSE
)

# -----------------------------
# approximate coordinates
# -----------------------------
coord_lookup <- data.frame(
  site_code = site_codes,
  lat = c(
    40.5231, 40.5203, 40.5166,
    40.4818, 40.4910, 40.5001,
    40.5260, 40.5127, 40.5212, 40.5186, 40.5131,
    40.4680, 40.4762, 40.4864,
    40.4696, 40.4938, 40.4836, 40.5041,
    40.5056, 40.5111, 40.5166, 40.5196
  ),
  lon = c(
    -105.7489, -105.7366, -105.7227,
    -105.7542, -105.7516, -105.7487,
    -105.6660, -105.6882, -105.6752, -105.6861, -105.6942,
    -105.7000, -105.6978, -105.7032,
    -105.7276, -105.7097, -105.7162, -105.7087,
    -105.7687, -105.7521, -105.7441, -105.7394
  ),
  stringsAsFactors = FALSE
)

SAFO_5_sites <- merge(
  SAFO_5_sites,
  coord_lookup,
  by = "site_code",
  all.x = TRUE,
  sort = FALSE
)

SAFO_5_sites <- SAFO_5_sites[match(site_codes, SAFO_5_sites$site_code), ]
rownames(SAFO_5_sites) <- NULL

SAFO_5_coords <- SAFO_5_sites[, c("site_id", "lat", "lon")]

# -----------------------------
# numeric fst object
# -----------------------------
SAFO_5_fst <- fst_mat
rownames(SAFO_5_fst) <- SAFO_5_sites$site_id
colnames(SAFO_5_fst) <- SAFO_5_sites$site_id

# -----------------------------
# straight-line distances
# -----------------------------
coord_mat <- as.matrix(SAFO_5_coords[, c("lon", "lat")])

geo_dist_km <- geosphere::distm(coord_mat, fun = geosphere::distHaversine) / 1000
rownames(geo_dist_km) <- SAFO_5_coords$site_id
colnames(geo_dist_km) <- SAFO_5_coords$site_id

# -----------------------------
# IBD dataframe
# -----------------------------
SAFO_5_ibd <- data.frame(
  site1   = rownames(SAFO_5_fst)[row(SAFO_5_fst)[upper.tri(SAFO_5_fst)]],
  site2   = colnames(SAFO_5_fst)[col(SAFO_5_fst)[upper.tri(SAFO_5_fst)]],
  fst     = SAFO_5_fst[upper.tri(SAFO_5_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  stringsAsFactors = FALSE
)

# -----------------------------
# map
# -----------------------------
usa <- map_data("state")

ggplot() +
  geom_polygon(
    data = usa,
    aes(x = long, y = lat, group = group),
    fill = "gray95",
    color = "gray70",
    linewidth = 0.2
  ) +
  geom_point(
    data = SAFO_5_sites,
    aes(x = lon, y = lat),
    size = 2.8
  ) +
  geom_text(
    data = SAFO_5_sites,
    aes(x = lon, y = lat, label = site_id),
    size = 2.8,
    nudge_y = 0.01
  ) +
  coord_quickmap(
    xlim = range(SAFO_5_sites$lon) + c(-0.03, 0.03),
    ylim = range(SAFO_5_sites$lat) + c(-0.03, 0.03)
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "SAFO-5 Brook Trout sampling sites"
  )


# -----------------------------
# IBD plot
# -----------------------------
ggplot(SAFO_5_ibd, aes(x = dist_km, y = fst)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Straight-line distance (km)",
    y = expression(F[ST]),
    title = "SAFO-5 isolation by distance"
  )


# -----------------------------
# save objects
# -----------------------------
save(
  SAFO_5_fst,
  SAFO_5_coords,
  #SAFO_5_sites,
  #SAFO_5_ibd,
  file = file.path(data_dir, "SAFO-5.RData")
)

# -----------------------------
# quick output
# -----------------------------
cat("\nFinished SAFO-5 processing.\n")
cat("Sites:", nrow(SAFO_5_coords), "\n")
cat("Dimensions of FST matrix:", dim(SAFO_5_fst)[1], "x", dim(SAFO_5_fst)[2], "\n\n")

print(SAFO_5_coords)
print(round(SAFO_5_fst[1:10, 1:10], 5))
