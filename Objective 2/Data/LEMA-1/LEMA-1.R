# -----------------------------
# LEMA-1 bluegill
# genepop -> genlight -> downsample SNPs -> dartR FST
# attach lat/lon to @other$latlon
# minimal map + gl.ibd plot
# save LEMA_1_fst and LEMA_1_coords
# -----------------------------

library(adegenet)
library(dplyr)
library(dartR)
library(ggplot2)
library(ggrepel)
library(maps)

# -----------------------------
# paths
# -----------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/LEMA-1"

gen_file <- file.path(base_dir, "genotypes", "blgl.neutral_thin.gen")
popmap_file <- file.path(base_dir, "popmap_filtered_blgl.txt")

# target number of SNPs
target_snps <- 7000

# reproducible downsampling
set.seed(123)

# -----------------------------
# 1) read genepop + popmap
# -----------------------------
LEMA_1_gen <- read.genepop(gen_file, ncode = 3)

LEMA_1_pops <- read.table(
  popmap_file,
  header = TRUE,
  stringsAsFactors = FALSE
)

LEMA_1_pops <- LEMA_1_pops[, 1:2]
names(LEMA_1_pops) <- c("sample_id", "pop")

LEMA_1_pops$sample_id <- trimws(LEMA_1_pops$sample_id)
LEMA_1_pops$pop <- trimws(LEMA_1_pops$pop)

# -----------------------------
# 2) match genepop individuals to popmap
# -----------------------------
LEMA_1_id_df <- data.frame(
  sample_id = trimws(indNames(LEMA_1_gen)),
  stringsAsFactors = FALSE
) |>
  left_join(LEMA_1_pops, by = "sample_id")

if (any(is.na(LEMA_1_id_df$pop))) {
  stop(
    "Some genepop individual IDs did not match the popmap:\n",
    paste(LEMA_1_id_df$sample_id[is.na(LEMA_1_id_df$pop)], collapse = ", ")
  )
}

pop(LEMA_1_gen) <- as.factor(LEMA_1_id_df$pop)

cat("\nPopulation counts:\n")
print(table(pop(LEMA_1_gen)))

# expected pops:
# LG_IL OR_MO Pool04_MN Pool08_WI Pool13_IA Pool26_IL

# -----------------------------
# 3) define population-level coordinates
# previously approximated reach midpoints
# -----------------------------
reach_key <- data.frame(
  pop = c("Pool04_MN", "Pool08_WI", "Pool13_IA", "Pool26_IL", "OR_MO", "LG_IL"),
  site = 1:6,
  lat = c(
    44.4494,  # Pool 4 / Lake City, MN
    43.8014,  # Pool 8 / La Crosse, WI
    42.2586,  # Pool 13 / Bellevue, IA
    38.8906,  # Pool 26 / Alton, IL
    37.3059,  # Open River / Cape Girardeau, MO
    40.3003   # La Grange / Havana, IL
  ),
  lon = c(
    -92.2660,
    -91.2396,
    -90.4229,
    -90.1843,
    -89.5181,
    -90.0609
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 4) convert genind -> genlight
# -----------------------------
LEMA_1_gl <- gi2gl(LEMA_1_gen)
pop(LEMA_1_gl) <- pop(LEMA_1_gen)

# -----------------------------
# 5) downsample to target SNP count
# -----------------------------
n_loci_total <- nLoc(LEMA_1_gl)
cat("\nTotal loci before downsampling:", n_loci_total, "\n")

if (n_loci_total > target_snps) {
  keep_idx <- sort(sample(seq_len(n_loci_total), target_snps))
  LEMA_1_gl <- LEMA_1_gl[, keep_idx]
  cat("Downsampled to", nLoc(LEMA_1_gl), "SNPs\n")
} else {
  cat("Dataset has", n_loci_total, "SNPs; keeping all loci\n")
}

# -----------------------------
# 6) attach individual-level coordinates to genlight @other$latlon
# each individual gets the lat/lon of its assigned population
# -----------------------------
ind_pops <- as.character(pop(LEMA_1_gl))

latlon_df <- reach_key[match(ind_pops, reach_key$pop), c("lon", "lat")]
rownames(latlon_df) <- indNames(LEMA_1_gl)

if (any(is.na(latlon_df$lon)) || any(is.na(latlon_df$lat))) {
  stop("Could not match all individuals to reach_key coordinates.")
}

LEMA_1_gl@other$latlon <- latlon_df

# -----------------------------
# 7) compute pairwise FST with dartR
# -----------------------------
LEMA_1_fst_dist <- gl.fst.pop(
  x = LEMA_1_gl,
  nboots = 1,
  nclusters = 1,
  verbose = 2
)

LEMA_1_fst_raw <- as.matrix(LEMA_1_fst_dist)

# enforce desired population order
pop_order <- reach_key$pop

missing_pops <- setdiff(pop_order, rownames(LEMA_1_fst_raw))
if (length(missing_pops) > 0) {
  stop(
    "These populations were missing from the FST output: ",
    paste(missing_pops, collapse = ", ")
  )
}

LEMA_1_fst_raw <- LEMA_1_fst_raw[pop_order, pop_order, drop = FALSE]

# -----------------------------
# 8) rebuild full symmetric matrix with zeros on diagonal
# -----------------------------
LEMA_1_fst <- LEMA_1_fst_raw

# fill upper triangle from lower triangle
LEMA_1_fst[upper.tri(LEMA_1_fst)] <-
  t(LEMA_1_fst_raw)[upper.tri(LEMA_1_fst_raw)]

# clean up
LEMA_1_fst[LEMA_1_fst < 0] <- 0
diag(LEMA_1_fst) <- 0

# relabel rows/cols
rownames(LEMA_1_fst) <- reach_key$site
colnames(LEMA_1_fst) <- reach_key$site

# -----------------------------
# 9) coords dataframe for saved object
# -----------------------------
LEMA_1_coords <- reach_key[, c("site", "lat", "lon")]

# helper object for plotting
LEMA_1_coords_plot <- reach_key[, c("site", "pop", "lat", "lon")]

# -----------------------------
# 10) minimal map
# states only, no Canada
# -----------------------------
us_map <- map_data("state")

x_pad <- 1.5
y_pad <- 1.2

xlim_use <- range(LEMA_1_coords_plot$lon) + c(-x_pad, x_pad)
ylim_use <- range(LEMA_1_coords_plot$lat) + c(-y_pad, y_pad)

p_map <- ggplot() +
  geom_polygon(
    data = us_map,
    aes(x = long, y = lat, group = group),
    fill = "gray97",
    color = "gray75",
    linewidth = 0.2
  ) +
  geom_point(
    data = LEMA_1_coords_plot,
    aes(x = lon, y = lat),
    size = 2.5
  ) +
  geom_text_repel(
    data = LEMA_1_coords_plot,
    aes(x = lon, y = lat, label = paste0(site, ": ", pop)),
    size = 3,
    max.overlaps = 100
  ) +
  coord_quickmap(
    xlim = xlim_use,
    ylim = ylim_use,
    expand = FALSE
  ) +
  theme_classic() +
  labs(
    title = "LEMA-1 sampling reaches",
    x = "Longitude",
    y = "Latitude"
  )

print(p_map)

# -----------------------------
# 11) IBD with dartR gl.ibd
# uses coordinates stored in @other$latlon
# -----------------------------
LEMA_1_ibd <- gl.ibd(
  x = LEMA_1_gl,
  Dgen = as.dist(LEMA_1_fst_raw),
  coordinates = "latlon",
  plot.out = TRUE,
  verbose = 2
)

# -----------------------------
# 12) inspect outputs
# -----------------------------
print(LEMA_1_coords)
print(round(LEMA_1_fst, 5))

cat("\nSite key:\n")
cat("1 = Pool04_MN\n")
cat("2 = Pool08_WI\n")
cat("3 = Pool13_IA\n")
cat("4 = Pool26_IL\n")
cat("5 = OR_MO\n")
cat("6 = LG_IL\n")

# -----------------------------
# 13) save
# -----------------------------
save(
  LEMA_1_fst,
  LEMA_1_coords,
  file = file.path(base_dir, "data", "LEMA-1.RData")
)