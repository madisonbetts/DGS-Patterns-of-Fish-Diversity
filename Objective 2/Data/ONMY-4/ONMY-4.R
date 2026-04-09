#-----------------
# ONMY-4
# Oncorhynchus mykiss
# South Fork Eel River tributaries
#-----------------

# This workflow uses the called genotype matrix from Kelson et al. (2020)
# and retains the named, georeferenceable site groups present directly in
# sample IDs:
# Barnwell, EelSF, ElderA, ElderB, Fox, Misery, Paralyze
#
# NOTE:
# The genotype readme states that exact sample_ID locations can be merged
# with "platemaps_omy5genotypes.csv". That file was not available locally
# here, so this workflow only uses named sample groups that are already
# explicit in the uploaded genotype matrix.
#
# Coordinates are best-available estimates from:
# - named stream mouth / confluence locations from South Fork Eel watershed
#   sources and study map context
# - reach centroids for ElderA / ElderB inferred from Figure 1 of Kelson et al. (2020)

library(dplyr)
library(ggplot2)
library(geosphere)
library(hierfstat)
library(readr)
library(maps)

# -----------------------------
# 1) read called genotype matrix
# -----------------------------
geno_raw <- read_csv("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ONMY-4/calledgenos_minInd50_95p.csv", show_col_types = FALSE)

# -----------------------------
# 2) identify sample IDs
# pattern in this file:
# col 1 = Chr
# col 2 = Pos
# then per individual:
# numeric genotype, sample_ID, posterior probability
# -----------------------------
all_numeric_cols <- seq(3, ncol(geno_raw), by = 3)
all_id_cols      <- seq(4, ncol(geno_raw), by = 3)

sample_ids <- colnames(geno_raw)[all_id_cols]
sample_prefix <- sub("_.*", "", sample_ids)

# retain only named, georeferenceable groups
keep_groups <- c("Barnwell", "EelSF", "ElderA", "ElderB", "Fox", "Misery", "Paralyze")
keep_idx <- which(sample_prefix %in% keep_groups)

keep_numeric_cols <- all_numeric_cols[keep_idx]
keep_id_cols      <- all_id_cols[keep_idx]
keep_ids          <- sample_ids[keep_idx]
keep_pops         <- sample_prefix[keep_idx]

# -----------------------------
# 3) build individual x locus genotype matrix
# numeric coding in file:
# -1 = no data
#  0 = hom major
#  1 = het
#  2 = hom minor
# convert to hierfstat-style diploid genotype codes:
# 11, 12, 22
# -----------------------------
geno_num <- as.data.frame(geno_raw[, keep_numeric_cols, drop = FALSE])
geno_num_t <- as.data.frame(t(as.matrix(geno_num)))
colnames(geno_num_t) <- paste0(geno_raw$Chr, "_", geno_raw$Pos)
rownames(geno_num_t) <- keep_ids

geno_num_t[geno_num_t == -1] <- NA

recode_snp <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  out[x == 0] <- 11
  out[x == 1] <- 12
  out[x == 2] <- 22
  out
}

geno_hf <- as.data.frame(lapply(geno_num_t, recode_snp))
geno_hf <- cbind(pop = as.numeric(factor(keep_pops, levels = keep_groups)), geno_hf)

# -----------------------------
# 4) pairwise FST
# -----------------------------
ONMY_4_fst <- pairwise.WCfst(geno_hf)
ONMY_4_fst <- as.matrix(ONMY_4_fst)

# force order to match keep_groups
ONMY_4_fst <- ONMY_4_fst[keep_groups, keep_groups]
ONMY_4_fst[is.na(ONMY_4_fst)] <- 0
ONMY_4_fst[ONMY_4_fst < 0] <- 0
diag(ONMY_4_fst) <- 0

# replace names with numeric site IDs
rownames(ONMY_4_fst) <- colnames(ONMY_4_fst) <- as.character(1:nrow(ONMY_4_fst))

# -----------------------------
# 5) coordinates
# site IDs
# 1 Barnwell Creek
# 2 South Fork Eel River mainstem
# 3 Elder Creek above waterfall
# 4 Elder Creek below waterfall
# 5 Fox Creek
# 6 Misery Creek
# 7 Paralyze Creek
# -----------------------------
ONMY_4_coords <- data.frame(
  site_id = 1:7,
  lat = c(
    39.7444018,  # 1 Barnwell Creek mouth
    39.7331716,  # 2 South Fork Eel mainstem near Fox/Elder confluences
    39.7218000,  # 3 ElderA reach centroid (above waterfall)
    39.7248000,  # 4 ElderB reach centroid (below waterfall)
    39.7389965,  # 5 Fox Creek mouth
    39.7181895,  # 6 Misery Creek confluence
    39.7170000   # 7 Paralyze Creek confluence / lower reach proxy
  ),
  lon = c(
    -123.6420360, # 1 Barnwell Creek
    -123.6369130, # 2 South Fork Eel mainstem
    -123.6225000, # 3 ElderA
    -123.6360000, # 4 ElderB
    -123.6287163, # 5 Fox Creek
    -123.6106415, # 6 Misery Creek
    -123.6055000  # 7 Paralyze Creek
  )
)

# -----------------------------
# 6) plotting labels only
# -----------------------------
plot_labs <- c(
  "1. Barnwell Creek",
  "2. South Fork Eel",
  "3. Elder above waterfall",
  "4. Elder below waterfall",
  "5. Fox Creek",
  "6. Misery Creek",
  "7. Paralyze Creek"
)

plot_df <- ONMY_4_coords
plot_df$label <- plot_labs

# -----------------------------
# 7) site map
# -----------------------------
usa <- ggplot2::map_data("world", region = "USA")
can <- ggplot2::map_data("world", region = "Canada")
world_map <- bind_rows(usa, can)

ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "grey90", color = "grey50", linewidth = 0.2
  ) +
  geom_point(
    data = plot_df,
    aes(x = lon, y = lat),
    size = 3
  ) +
  geom_text(
    data = plot_df,
    aes(x = lon, y = lat, label = label),
    hjust = 0, nudge_x = 0.01, size = 3
  ) +
  coord_fixed(1.25, xlim = c(-123.70, -123.56), ylim = c(39.70, 39.76)) +
  theme_bw() +
  labs(
    title = "ONMY-4 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )

# -----------------------------
# 8) pairwise geographic distance
# -----------------------------
coords_mat <- as.matrix(ONMY_4_coords[, c("lon", "lat")])
geo_dist_km <- geosphere::distm(coords_mat, fun = geosphere::distGeo) / 1000
rownames(geo_dist_km) <- colnames(geo_dist_km) <- as.character(ONMY_4_coords$site_id)

# -----------------------------
# 9) IBD dataframe
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(ONMY_4_fst)[row(ONMY_4_fst)[upper.tri(ONMY_4_fst)]],
  site2   = colnames(ONMY_4_fst)[col(ONMY_4_fst)[upper.tri(ONMY_4_fst)]],
  fst     = ONMY_4_fst[upper.tri(ONMY_4_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 10) IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(
    title = "ONMY-4 isolation by distance",
    x = "Geographic distance (km)",
    y = expression(F[ST])
  )

# -----------------------------
# 11) save RData
# -----------------------------
save(
  ONMY_4_fst,
  ONMY_4_coords,
  file = file.path("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ONMY-4/data/ONMY-4.RData")
)
