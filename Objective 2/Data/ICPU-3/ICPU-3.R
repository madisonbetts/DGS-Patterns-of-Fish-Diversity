# -----------------------------
# ICPU-3 channel catfish
# genepop -> genlight -> downsample SNPs -> dartR FST
# attach lat/lon to @other$latlon
# minimal map + gl.ibd plot
# save ICPU_3_fst and ICPU_3_coords
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
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ICPU-3"

gen_file <- file.path(base_dir, "genotypes", "cncf.neutral_thin.gen")
popmap_file <- file.path(base_dir, "popmap_filtered_cncf.txt")

# target number of SNPs
target_snps <- 7000

# reproducible downsampling
set.seed(123)

# -----------------------------
# 1) read genepop + popmap
# -----------------------------
ICPU_3_gen <- read.genepop(gen_file, ncode = 3)

ICPU_3_pops <- read.table(
  popmap_file,
  header = TRUE,
  stringsAsFactors = FALSE
)

ICPU_3_pops <- ICPU_3_pops[, 1:2]
names(ICPU_3_pops) <- c("sample_id", "pop")

ICPU_3_pops$sample_id <- trimws(ICPU_3_pops$sample_id)
ICPU_3_pops$pop <- trimws(ICPU_3_pops$pop)

# -----------------------------
# 2) match genepop individuals to popmap
# -----------------------------
ICPU_3_id_df <- data.frame(
  sample_id = trimws(indNames(ICPU_3_gen)),
  stringsAsFactors = FALSE
) |>
  left_join(ICPU_3_pops, by = "sample_id")

if (any(is.na(ICPU_3_id_df$pop))) {
  stop(
    "Some genepop individual IDs did not match the popmap:\n",
    paste(ICPU_3_id_df$sample_id[is.na(ICPU_3_id_df$pop)], collapse = ", ")
  )
}

pop(ICPU_3_gen) <- as.factor(ICPU_3_id_df$pop)

cat("\nPopulation counts:\n")
print(table(pop(ICPU_3_gen)))

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
ICPU_3_gl <- gi2gl(ICPU_3_gen)
pop(ICPU_3_gl) <- pop(ICPU_3_gen)

# -----------------------------
# 5) downsample to target SNP count
# -----------------------------
n_loci_total <- nLoc(ICPU_3_gl)
cat("\nTotal loci before downsampling:", n_loci_total, "\n")

if (n_loci_total > target_snps) {
  keep_idx <- sort(sample(seq_len(n_loci_total), target_snps))
  ICPU_3_gl <- ICPU_3_gl[, keep_idx]
  cat("Downsampled to", nLoc(ICPU_3_gl), "SNPs\n")
} else {
  cat("Dataset has", n_loci_total, "SNPs; keeping all loci\n")
}

# -----------------------------
# 6) attach individual-level coordinates to genlight @other$latlon
# each individual gets the lat/lon of its assigned population
# -----------------------------
ind_pops <- as.character(pop(ICPU_3_gl))

latlon_df <- reach_key[match(ind_pops, reach_key$pop), c("lon", "lat")]
rownames(latlon_df) <- indNames(ICPU_3_gl)

if (any(is.na(latlon_df$lon)) || any(is.na(latlon_df$lat))) {
  stop("Could not match all individuals to reach_key coordinates.")
}

ICPU_3_gl@other$latlon <- latlon_df

# -----------------------------
# 7) compute pairwise FST with dartR
# -----------------------------
ICPU_3_fst_dist <- gl.fst.pop(
  x = ICPU_3_gl,
  nboots = 1,
  nclusters = 1,
  verbose = 2
)

ICPU_3_fst_raw <- as.matrix(ICPU_3_fst_dist)

# enforce desired population order
pop_order <- reach_key$pop

missing_pops <- setdiff(pop_order, rownames(ICPU_3_fst_raw))
if (length(missing_pops) > 0) {
  stop(
    "These populations were missing from the FST output: ",
    paste(missing_pops, collapse = ", ")
  )
}

ICPU_3_fst_raw <- ICPU_3_fst_raw[pop_order, pop_order, drop = FALSE]

# -----------------------------
# 8) rebuild full symmetric matrix with zeros on diagonal
# -----------------------------
ICPU_3_fst <- ICPU_3_fst_raw

# fill upper triangle from lower triangle
ICPU_3_fst[upper.tri(ICPU_3_fst)] <-
  t(ICPU_3_fst_raw)[upper.tri(ICPU_3_fst_raw)]

# clean up
ICPU_3_fst[ICPU_3_fst < 0] <- 0
diag(ICPU_3_fst) <- 0

# relabel rows/cols
rownames(ICPU_3_fst) <- reach_key$site
colnames(ICPU_3_fst) <- reach_key$site

# -----------------------------
# 9) coords dataframe for saved object
# -----------------------------
ICPU_3_coords <- reach_key[, c("site", "lat", "lon")]

# helper object for plotting
ICPU_3_coords_plot <- reach_key[, c("site", "pop", "lat", "lon")]

# -----------------------------
# 10) minimal map
# -----------------------------
us_map <- map_data("state")

x_pad <- 1.5
y_pad <- 1.2

xlim_use <- range(ICPU_3_coords_plot$lon) + c(-x_pad, x_pad)
ylim_use <- range(ICPU_3_coords_plot$lat) + c(-y_pad, y_pad)

p_map <- ggplot() +
  geom_polygon(
    data = us_map,
    aes(x = long, y = lat, group = group),
    fill = "gray97",
    color = "gray75",
    linewidth = 0.2
  ) +
  geom_point(
    data = ICPU_3_coords_plot,
    aes(x = lon, y = lat),
    size = 2.5
  ) +
  geom_text_repel(
    data = ICPU_3_coords_plot,
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
    title = "ICPU-3 sampling reaches",
    x = "Longitude",
    y = "Latitude"
  )

print(p_map)

# -----------------------------
# 11) IBD with dartR gl.ibd
# uses coordinates stored in @other$latlon
# -----------------------------
ICPU_3_ibd <- gl.ibd(
  x = ICPU_3_gl,
  Dgen = as.dist(ICPU_3_fst_raw),
  coordinates = "latlon",
  plot.out = TRUE,
  verbose = 2
)

# -----------------------------
# 12) inspect outputs
# -----------------------------
print(ICPU_3_coords)
print(round(ICPU_3_fst, 5))

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
  ICPU_3_fst,
  ICPU_3_coords,
  file = file.path(base_dir, "data", "ICPU-3.RData")
)