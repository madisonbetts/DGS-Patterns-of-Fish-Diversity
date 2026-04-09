# ============================================================
# SAFO-13 | Salvelinus fontinalis
# Michaelides et al. 2025
#
# Objective 2 workflow
# - read BT_ALL.vcf as genlight
# - link individuals to populations from BT_ALL_popmap.txt
# - scalp site coordinates directly from Table 1 in the paper
#   using the exact decimals shown there (no rounding)
# - attach individual lat/lon to BT_snps@other$latlon
# - subset to inland wild Rhode Island sites for FST / IBD
#   to match the paper's riverscape workflow
# - calculate pairwise FST with gl.fst.pop() using no extra arguments
# - clamp negative FST to 0
# - plot sites and QC IBD in RStudio only
# - save SAFO_13_fst and SAFO_13_coords to data/SAFO-13.RData
# ============================================================

suppressPackageStartupMessages({
  library(vcfR)
  library(dartR)
  library(dplyr)
  library(ggplot2)
  library(geosphere)
})

study_code <- "SAFO-13"
study_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SAFO-13"
out_dir <- file.path(study_dir, "data")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

# -----------------------------
# 1) read VCF as genlight
# -----------------------------
BT_snps <- vcfR2genlight(
  read.vcfR("28628804/BT_ALL.vcf")
)

# -----------------------------
# 2) read popmap and match to individuals
# -----------------------------
popmap <- read.table("BT_ALL_popmap.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(popmap) <- c("ind_name", "pop")

popmap <- popmap[match(indNames(BT_snps), popmap$ind_name), ]

if (any(is.na(popmap$ind_name))) {
  stop("Some genlight individuals were not matched in BT_ALL_popmap.txt")
}
if (!all(popmap$ind_name == indNames(BT_snps))) {
  stop("Popmap ordering does not match BT_snps individual names")
}

BT_snps@pop <- as.factor(popmap$pop)
BT_snps@other$ind_sites <- popmap$pop

# -----------------------------
# 3) exact coordinates from Table 1 in Michaelides et al. 2025
# hatchery strains do not have coordinates in Table 1
# so individual latlon are attached for the 16 georeferenced wild/coastal sites
# -----------------------------
coords_lookup <- data.frame(
  pop = c(
    "AFP", "BAK", "BHU", "BRH", "BRO", "BRP", "BRU", "FRC",
    "MBC", "MBP", "SBD", "SHB", "UTA", "UTB", "UTB2", "WR"
  ),
  site_name = c(
    "Acid Factory Brook",
    "Baker Brook",
    "Breakheart Brook above Breakheart Pond",
    "Beaver River at Hillsdale Road",
    "Beaver River at Old Mountain Road",
    "Beaver River at Punchbowl Trail",
    "Brushy Brook",
    "Falls River",
    "Meadow Brook Carolina Management Area",
    "Meadow Brook at Pine Hill Road",
    "Smelt/Brown Brook below Rum Pond",
    "Shelter Harbour Brook",
    "Unnamed tributary to Wood River below Alton Dam",
    "Unnamed tributary to Breakheart Pond (above dam)",
    "Unnamed tributary to Breakheart Pond (below dam)",
    "Wood River upstream of Frying Pan Pond"
  ),
  lon = c(
    -71.721953,
    -71.693075,
    -71.697832,
    -71.640094,
    -71.640852,
    -71.641560,
    -71.736417,
    -71.735044,
    -71.689147,
    -71.688469,
    -71.512712,
    -71.738029,
    -71.723342,
    -71.699931,
    -71.700343,
    -71.714549
  ),
  lat = c(
    41.633769,
    41.542195,
    41.604672,
    41.522293,
    41.538449,
    41.512152,
    41.528783,
    41.588722,
    41.457545,
    41.471734,
    41.414437,
    41.341279,
    41.437903,
    41.595471,
    41.595870,
    41.565066
  ),
  stringsAsFactors = FALSE
)

ind_latlon <- popmap %>%
  left_join(coords_lookup[, c("pop", "lat", "lon")], by = "pop")

BT_snps@other$latlon <- ind_latlon[, c("lat", "lon")]
rownames(BT_snps@other$latlon) <- ind_latlon$ind_name

# -----------------------------
# 4) subset to inland RI wild sites for FST / IBD
# matches BT_Inland logic in the paper:
# exclude hatcheries LFA, LFR and coastal sites SHB, SBD
# -----------------------------
exclude_pops <- c("LFA", "LFR", "SHB", "SBD")
keep_ind <- !(BT_snps@pop %in% exclude_pops)

BT_inland <- BT_snps[keep_ind, ]
BT_inland@pop <- droplevels(BT_inland@pop)
BT_inland@other$ind_sites <- BT_inland@pop
BT_inland@other$latlon <- BT_snps@other$latlon[keep_ind, , drop = FALSE]

if (any(is.na(BT_inland@other$latlon[, 1])) || any(is.na(BT_inland@other$latlon[, 2]))) {
  stop("Some inland individuals are missing coordinates after pop join")
}

# -----------------------------
# 5) pairwise FST from scratch using gl.fst.pop()
# -----------------------------
fst_raw <- gl.fst.pop(BT_inland)

SAFO_13_fst <- as.matrix(fst_raw)

# handle NA + negatives first
SAFO_13_fst[is.na(SAFO_13_fst)] <- 0
SAFO_13_fst[SAFO_13_fst < 0] <- 0

# enforce symmetry (CRITICAL)
SAFO_13_fst <- pmax(SAFO_13_fst, t(SAFO_13_fst))

# ensure diagonal = 0
diag(SAFO_13_fst) <- 0

# -----------------------------
# 6) coords df aligned to FST order
# -----------------------------
fst_pops <- rownames(SAFO_13_fst)

coords_ordered <- coords_lookup %>%
  filter(pop %in% fst_pops) %>%
  slice(match(fst_pops, pop))

if (!all(coords_ordered$pop == fst_pops)) {
  stop("Coordinate order does not match FST population order")
}

coords_ordered$site <- seq_len(nrow(coords_ordered))

rownames(SAFO_13_fst) <- coords_ordered$site
colnames(SAFO_13_fst) <- coords_ordered$site

SAFO_13_coords <- coords_ordered %>%
  select(site, lat, lon)

SAFO_13_coords_check <- coords_ordered %>%
  select(site, pop, site_name, lat, lon)

# -----------------------------
# 7) map plot
# include US + Canada and zoom to point extent
# -----------------------------
world_df <- ggplot2::map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

plot_df <- SAFO_13_coords_check %>%
  mutate(label = paste0(site, " ", pop))

x_pad <- max(0.15, diff(range(plot_df$lon)) * 0.12)
y_pad <- max(0.10, diff(range(plot_df$lat)) * 0.12)

p_map <- ggplot() +
  geom_polygon(
    data = world_df,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = plot_df,
    aes(x = lon, y = lat),
    size = 2
  ) +
  geom_text(
    data = plot_df,
    aes(x = lon, y = lat, label = label),
    nudge_y = 0.01,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(min(plot_df$lon) - x_pad, max(plot_df$lon) + x_pad),
    ylim = c(min(plot_df$lat) - y_pad, max(plot_df$lat) + y_pad)
  ) +
  labs(
    title = "SAFO-13 inland sampling sites",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw()

print(p_map)

# -----------------------------
# 8) QC IBD plot
# straight-line distance only for QC
# use linearized FST to match paper convention
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(SAFO_13_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

#fst_linearized <- SAFO_13_fst / (1 - SAFO_13_fst)
#diag(fst_linearized) <- 0

ibd_df <- data.frame(
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  fst_lin = SAFO_13_fst[upper.tri(SAFO_13_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst_lin)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "SAFO-13 inland IBD plot",
    x = "Geographic distance (km)",
    y = "FST"
  ) +
  theme_bw()

print(p_ibd)

# -----------------------------
# 9) save outputs
# -----------------------------
save(
  SAFO_13_fst,
  SAFO_13_coords,
  file = file.path(out_dir, "SAFO-13.RData")
)

write.csv(
  SAFO_13_coords_check,
  file = file.path(out_dir, "SAFO-13_coords_check.csv"),
  row.names = FALSE
)