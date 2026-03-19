# -------------------------
# rainbow darter
# ETCA-1 Luiken et al.
# site metadata + FST workflow
# -------------------------

library(adegenet)
library(magrittr)
library(dartR)
library(ggplot2)

# -------------------------
# 0) coordinates
# -------------------------
ETCA_1_coords <- data.frame(
  site = c("V01", "V02", "V03", "V04", "V05",
           "M01", "M02", "M03", "M04", "M05"),
  lat = c(42.78801, 42.81858, 42.84436, 42.83861, 42.86248,
          37.63617, 37.65924, 37.70125, 37.74572, 37.95130),
  lon = c(-91.88205, -91.87847, -91.79318, -91.77200, -91.76365,
          -91.41185, -91.41496, -91.44697, -91.43385, -91.50868),
  stringsAsFactors = FALSE
)

ETCA_1_coords

ggplot(ETCA_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.08, size = 4) +
  coord_fixed() +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude")

# -------------------------
# 1) genotype file
# -------------------------
file_path <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETCA-1/doi_10_5061_dryad_6m905qg17__v20221129/Combined/populations.snps.gen"

# -------------------------
# 2) individual -> site crosswalk
# -------------------------
pop_key <- c(
  Ec601 = "V05", Ec602 = "V05", Ec603 = "V05", Ec605 = "V05", Ec606 = "V05",
  Ec608 = "V05", Ec609 = "V05", Ec610 = "V05", Ec611 = "V05", Ec612 = "V05",
  Ec613 = "V03", Ec614 = "V03", Ec615 = "V03", Ec616 = "V03", Ec617 = "V03",
  Ec619 = "V03", Ec620 = "V03", Ec621 = "V03", Ec622 = "V03", Ec624 = "V03",
  Ec626 = "V02", Ec627 = "V02", Ec628 = "V02", Ec629 = "V02", Ec630 = "V02",
  Ec631 = "V02", Ec633 = "V02", Ec634 = "V02", Ec635 = "V02", Ec636 = "V02",
  Ec637 = "V01", Ec638 = "V01", Ec639 = "V01", Ec641 = "V01", Ec642 = "V01",
  Ec644 = "V01", Ec645 = "V01", Ec646 = "V01", Ec647 = "V01", Ec648 = "V01",
  Ec649 = "V04", Ec650 = "V04", Ec651 = "V04", Ec652 = "V04", Ec654 = "V04",
  Ec655 = "V04", Ec656 = "V04", Ec657 = "V04", Ec658 = "V04", Ec659 = "V04",
  Ec685 = "M05", Ec686 = "M05", Ec687 = "M05", Ec689 = "M05", Ec690 = "M05",
  Ec691 = "M05", Ec692 = "M05", Ec693 = "M05", Ec694 = "M05", Ec696 = "M05",
  Ec697 = "M01", Ec698 = "M01", Ec699 = "M01", Ec700 = "M01", Ec702 = "M01",
  Ec703 = "M01", Ec704 = "M01", Ec705 = "M01", Ec706 = "M01", Ec707 = "M01",
  Ec709 = "M02", Ec710 = "M02", Ec712 = "M02", Ec713 = "M02", Ec714 = "M02",
  Ec715 = "M02", Ec716 = "M02", Ec718 = "M02", Ec719 = "M02", Ec720 = "M02",
  Ec721 = "M03", Ec723 = "M03", Ec725 = "M03", Ec726 = "M03", Ec727 = "M03",
  Ec728 = "M03", Ec729 = "M03", Ec730 = "M03", Ec731 = "M03", Ec732 = "M03",
  Ec733 = "M04", Ec734 = "M04", Ec735 = "M04", Ec736 = "M04", Ec738 = "M04",
  Ec739 = "M04", Ec740 = "M04", Ec741 = "M04", Ec743 = "M04", Ec744 = "M04"
)

# -------------------------
# 3) read, filter, assign pops
# -------------------------
ETCAE_gl <-
  read.genepop(file_path, ncode = 2) %>%
  gi2gl() %>%
  gl.filter.monomorphs() %>%
  gl.filter.callrate(method = "loc", threshold = 0.9) %>%
  gl.filter.callrate(method = "ind", threshold = 0.5)

pop(ETCAE_gl) <- factor(pop_key[indNames(ETCAE_gl)], levels = ETCA_1_coords$site)

if (any(is.na(pop(ETCAE_gl)))) {
  stop("Some individuals were not matched in pop_key.")
}

table(pop(ETCAE_gl))

# -------------------------
# 4) attach coordinates to individuals
# -------------------------
ETCAE_gl@other$latlon <- as.matrix(
  ETCA_1_coords[match(pop(ETCAE_gl), ETCA_1_coords$site), c("lat", "lon")]
)
rownames(ETCAE_gl@other$latlon) <- indNames(ETCAE_gl)
colnames(ETCAE_gl@other$latlon) <- c("lat", "lon")

head(ETCAE_gl@other$latlon)

# -------------------------
# 5) pairwise FST from dartR
# -------------------------
ETCA_1_fst <- gl.fst.pop(ETCAE_gl)

if (inherits(ETCA_1_fst, "dist")) {
  ETCA_1_fst <- as.matrix(ETCA_1_fst)
}

# fill missing triangle if dartR returned only one half
ETCA_1_fst[is.na(ETCA_1_fst)] <- t(ETCA_1_fst)[is.na(ETCA_1_fst)]

# reorder to match coordinate table
ETCA_1_fst <- ETCA_1_fst[
  match(ETCA_1_coords$site, rownames(ETCA_1_fst)),
  match(ETCA_1_coords$site, colnames(ETCA_1_fst))
]

diag(ETCA_1_fst) <- 0

stopifnot(identical(rownames(ETCA_1_fst), ETCA_1_coords$site))
stopifnot(identical(colnames(ETCA_1_fst), ETCA_1_coords$site))
stopifnot(isTRUE(all.equal(ETCA_1_fst, t(ETCA_1_fst))))

ETCA_1_fst

# -------------------------
# 6) optional IBD plot
# -------------------------
site_dist <- as.matrix(dist(ETCA_1_coords[, c("lon", "lat")]))

ibd_df <- data.frame(
  dist = site_dist[upper.tri(site_dist)],
  fst  = ETCA_1_fst[upper.tri(ETCA_1_fst)]
)

plot(
  ibd_df$dist,
  ibd_df$fst,
  xlab = "Euclidean distance (decimal-degree units)",
  ylab = "Pairwise FST",
  pch = 21
)

# -------------------------
# 7) save RData
# -------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETCA-1/data"

if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

save(
  ETCA_1_fst,
  ETCA_1_coords,
  file = file.path(save_dir, "ETCA-1.RData")
)