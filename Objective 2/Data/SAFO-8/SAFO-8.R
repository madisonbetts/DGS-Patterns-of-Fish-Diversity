# ============================================================
# SAFO-8
# White, Hanks, and Wagner (2020)
# Brook trout (Salvelinus fontinalis)
# Loyalsock Creek watershed, Pennsylvania, USA
#
# Notes
# - The README in the supplied BGR materials states that D.rds is
#   the FST matrix and covariates.rds contains the covariate matrices.
# - The provided BGR model code reads D.rds directly and fits the
#   model with df = 12 microsatellite loci and 33 observed nodes.
# - The paper states the analysis included 33 sites and used pairwise
#   FST as the genetic distance response.
# - This workflow therefore uses D.rds directly as the published
#   33 x 33 pairwise FST matrix for Objective 2 extraction.
# - Site labels and approximate coordinates below are reconstructed
#   from the labeled sampling map in Figure 4A and watershed
#   geography. They are best-available inferred coordinates, not
#   original field GPS coordinates.
# - If you later recover exact coordinates, update only the coords
#   object while keeping site order identical to the FST matrix.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(geosphere)
  library(maps)
  library(grid)
})

# -----------------------------
# 0) paths
# -----------------------------
study_code <- "SAFO-8"

base_dir <- "~/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SAFO-8"
base_dir <- path.expand(base_dir)

data_dir <- file.path(base_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# expected supplemental directory provided by user
supp_dir <- file.path(base_dir, "slw361-BGR_Model-783dcbe")

# locate D.rds robustly
d_candidates <- c(
  file.path(supp_dir, "D.rds"),
  file.path(base_dir, "D.rds")
)
d_file <- d_candidates[file.exists(d_candidates)][1]

if (is.na(d_file)) {
  stop(
    "Could not find D.rds. Checked:\n",
    paste(" -", d_candidates, collapse = "\n")
  )
}

# -----------------------------
# 1) read the published FST matrix
# -----------------------------
fst_raw <- readRDS(d_file)

if (!is.matrix(fst_raw) && !is.data.frame(fst_raw)) {
  stop("D.rds did not read in as a matrix / data.frame.")
}

fst_raw <- as.matrix(fst_raw)

if (nrow(fst_raw) != 33 || ncol(fst_raw) != 33) {
  stop(
    "Expected a 33 x 33 FST matrix in D.rds; got ",
    nrow(fst_raw), " x ", ncol(fst_raw)
  )
}

orig_node_ids <- rownames(fst_raw)
if (is.null(orig_node_ids)) {
  orig_node_ids <- paste0("node_", seq_len(nrow(fst_raw)))
}

# -----------------------------
# 2) site order
# -----------------------------
# Order below follows the supplied SAFO-8 templates / BGR materials.
# Original D.rds rownames correspond to observed node IDs in the SSEN:
# X103, X72:X74, X76:X102, X104, X105

site_key <- tibble::tribble(
  ~site_id, ~node_id, ~site_name,
   1, "X103", "FLAG",
   2, "X72",  "CONK",
   3, "X73",  "MILA",
   4, "X74",  "BEAR",
   5, "X76",  "USCO",
   6, "X77",  "DSCO",
   7, "X78",  "DPOL",
   8, "X79",  "UPOL",
   9, "X80",  "DSHA",
  10, "X81",  "USHA",
  11, "X82",  "DOUB",
  12, "X83",  "DSEA",
  13, "X84",  "USEA",
  14, "X85",  "YELL",
  15, "X86",  "ROCK",
  16, "X87",  "LEV",
  17, "X88",  "STRB",
  18, "X89",  "SCAR",
  19, "X90",  "UNT",
  20, "X91",  "LICK",
  21, "X92",  "USWE",
  22, "X93",  "DSWE",
  23, "X94",  "DRHO",
  24, "X95",  "SWAM",
  25, "X96",  "MIHI",
  26, "X97",  "HUCK",
  27, "X98",  "BRUN",
  28, "X99",  "GRAN",
  29, "X100", "SNAK",
  30, "X101", "SSR",
  31, "X102", "RED",
  32, "X104", "DSLB",
  33, "X105", "JACO"
)

if (!all(site_key$node_id == orig_node_ids)) {
  warning(
    "The inferred node_id -> site_name mapping does not perfectly match ",
    "the row order in D.rds. Check `site_key` if anything looks off."
  )
}

# -----------------------------
# 3) approximate site coordinates
# -----------------------------
# Best-available approximate centroids inferred from the labeled
# sampling map in Figure 4A plus watershed geography.
# These are intentionally transparent approximations.

SAFO_8_coords_full <- tibble::tribble(
  ~site_id, ~node_id, ~site_name, ~lat,    ~lon,
   1, "X103", "FLAG", 41.559, -76.739,
   2, "X72",  "CONK", 41.485, -76.915,
   3, "X73",  "MILA", 41.510, -76.941,
   4, "X74",  "BEAR", 41.523, -76.978,
   5, "X76",  "USCO", 41.504, -77.013,
   6, "X77",  "DSCO", 41.497, -77.020,
   7, "X78",  "DPOL", 41.492, -77.028,
   8, "X79",  "UPOL", 41.484, -77.034,
   9, "X80",  "DSHA", 41.486, -77.061,
  10, "X81",  "USHA", 41.475, -77.074,
  11, "X82",  "DOUB", 41.495, -77.126,
  12, "X83",  "DSEA", 41.487, -77.111,
  13, "X84",  "USEA", 41.478, -77.105,
  14, "X85",  "YELL", 41.540, -77.081,
  15, "X86",  "ROCK", 41.533, -77.016,
  16, "X87",  "LEV",  41.578, -77.084,
  17, "X88",  "STRB", 41.567, -77.077,
  18, "X89",  "SCAR", 41.514, -77.154,
  19, "X90",  "UNT",  41.593, -77.192,
  20, "X91",  "LICK", 41.544, -77.225,
  21, "X92",  "USWE", 41.563, -77.437,
  22, "X93",  "DSWE", 41.533, -77.470,
  23, "X94",  "DRHO", 41.516, -77.326,
  24, "X95",  "SWAM", 41.492, -77.377,
  25, "X96",  "MIHI", 41.471, -77.414,
  26, "X97",  "HUCK", 41.417, -77.240,
  27, "X98",  "BRUN", 41.400, -77.136,
  28, "X99",  "GRAN", 41.357, -77.288,
  29, "X100", "SNAK", 41.362, -77.345,
  30, "X101", "SSR",  41.343, -77.385,
  31, "X102", "RED",  41.330, -77.414,
  32, "X104", "DSLB", 41.315, -77.451,
  33, "X105", "JACO", 41.354, -77.628
)

stopifnot(nrow(SAFO_8_coords_full) == 33)

# -----------------------------
# 4) final FST object
# -----------------------------
SAFO_8_fst <- fst_raw
storage.mode(SAFO_8_fst) <- "numeric"
SAFO_8_fst[is.na(SAFO_8_fst)] <- 0
SAFO_8_fst[SAFO_8_fst < 0] <- 0
diag(SAFO_8_fst) <- 0

# numeric site IDs as row/col names
rownames(SAFO_8_fst) <- colnames(SAFO_8_fst) <- as.character(site_key$site_id)

# final saved coords object: matching site_id order with lat/lon
SAFO_8_coords <- SAFO_8_coords_full |>
  select(site_id, lat, lon)

# labeled version retained only for plotting / QC
SAFO_8_coords_plot <- SAFO_8_coords_full

# -----------------------------
# 5) map of sampling locations
# -----------------------------
# user requested US + Canada context and zoom to points
world_df <- map_data("world") |>
  filter(region %in% c("USA", "Canada"))

xpad <- max(0.35, diff(range(SAFO_8_coords_plot$lon)) * 0.12)
ypad <- max(0.20, diff(range(SAFO_8_coords_plot$lat)) * 0.12)

xlim_use <- range(SAFO_8_coords_plot$lon) + c(-xpad, xpad)
ylim_use <- range(SAFO_8_coords_plot$lat) + c(-ypad, ypad)

p_map <- ggplot() +
  geom_polygon(
    data = world_df,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey65",
    linewidth = 0.2
  ) +
  geom_point(
    data = SAFO_8_coords_plot,
    aes(x = lon, y = lat),
    size = 2.1
  ) +
  geom_text(
    data = SAFO_8_coords_plot,
    aes(x = lon, y = lat, label = site_name),
    nudge_y = 0.010,
    size = 2.6
  ) +
  coord_quickmap(
    xlim = xlim_use,
    ylim = ylim_use,
    expand = FALSE
  ) +
  theme_classic() +
  labs(
    title = "SAFO-8 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )

print(p_map)

# -----------------------------
# 6) IBD plot
# -----------------------------
# quick-check plot using great-circle distance because river distance
# is not provided here as a simple site x site matrix for the
# Objective 2 extraction objects

coord_mat <- as.matrix(SAFO_8_coords[, c("lon", "lat")])
geo_dist_km <- geosphere::distm(coord_mat, fun = geosphere::distHaversine) / 1000

upper_idx <- upper.tri(SAFO_8_fst)

ibd_df <- tibble(
  fst = SAFO_8_fst[upper_idx],
  dist_km = geo_dist_km[upper_idx]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    title = "SAFO-8 IBD quick-check",
    x = "Great-circle distance (km)",
    y = expression(F[ST])
  )

print(p_ibd)

# -----------------------------
# 7) save RData
# -----------------------------
save(
  SAFO_8_fst,
  SAFO_8_coords,
  file = file.path(data_dir, "SAFO-8.RData")
)

# -----------------------------
# 8) optional assignment to global env
# -----------------------------
assign("SAFO_8_fst", SAFO_8_fst, envir = .GlobalEnv)
assign("SAFO_8_coords", SAFO_8_coords, envir = .GlobalEnv)

cat("\n==================== SUMMARY ====================\n")
cat("Study code: ", study_code, "\n", sep = "")
cat("D.rds source: ", d_file, "\n", sep = "")
cat("Sites: ", nrow(SAFO_8_coords), "\n", sep = "")
cat("FST matrix dimensions: ", nrow(SAFO_8_fst), " x ", ncol(SAFO_8_fst), "\n", sep = "")
cat("RData saved to: ", file.path(data_dir, "SAFO-8.RData"), "\n", sep = "")
cat("=================================================\n\n")
