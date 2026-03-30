# -----------------------------
# SAFO-6
# Salvelinus fontinalis (brook trout)
# Whiteley et al. 2016 (West Brook system)
# -----------------------------

library(dplyr)
library(hierfstat)
library(geosphere)
library(ggplot2)

# -----------------------------
# paths
# -----------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SAFO-6"
setwd(file.path(base_dir, "doi_10_5061_dryad_s4td3__v20161215"))

geno_file <- "WB+2001-2009+dryad.txt"

if (!file.exists(geno_file)) {
  stop("Genotype file not found in working directory: ", getwd())
}

# -----------------------------
# 1) read + clean genotype data
# -----------------------------
dat <- read.table(
  geno_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dat <- dat %>%
  filter(Species == "BKT")

# -----------------------------
# 2) define section-level pops
# -----------------------------
dat <- dat %>%
  mutate(
    pop = paste(River, Section, sep = "_")
  )

# -----------------------------
# 3) inspect cohorts and choose one
# section-level pops must have >= min_n individuals
# -----------------------------
min_n <- 20

cohort_pop_summary <- dat %>%
  group_by(Cohort, pop) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Cohort) %>%
  summarise(
    n_pops_ge_min = sum(n >= min_n),
    total_n = sum(n),
    .groups = "drop"
  ) %>%
  arrange(desc(n_pops_ge_min), desc(total_n), Cohort)

print(cohort_pop_summary)

# if you want to force a cohort, set it here; otherwise leave as NA
target_cohort <- NA

if (is.na(target_cohort)) {
  target_cohort <- cohort_pop_summary$Cohort[1]
}

message("Using cohort: ", target_cohort)

dat_sub <- dat %>%
  filter(Cohort == target_cohort)

pop_sizes <- table(dat_sub$pop)
keep_pops <- names(pop_sizes[pop_sizes >= min_n])

dat_sub <- dat_sub %>%
  filter(pop %in% keep_pops)

if (nrow(dat_sub) == 0) {
  stop("No rows left after cohort and population-size filtering. Check cohort_pop_summary above.")
}

if (length(unique(dat_sub$pop)) < 2) {
  stop("Need at least 2 populations after filtering; currently have ", length(unique(dat_sub$pop)))
}

message("Rows in dat_sub: ", nrow(dat_sub))
message("Retained pops: ", length(unique(dat_sub$pop)))
print(summary(table(dat_sub$pop)))

# -----------------------------
# 4) identify microsatellite loci
# -----------------------------
loci_cols <- grep("\\.(1|2)$", names(dat_sub), value = TRUE)

if (length(loci_cols) != 24) {
  stop("Expected 24 allele columns (12 loci x 2), found: ", length(loci_cols))
}

locus_bases <- unique(sub("\\.(1|2)$", "", loci_cols))

# -----------------------------
# 5) build hierfstat input directly
# hierfstat wants one column per locus coded like 128132
# -----------------------------
make_hf_locus <- function(df, locus) {
  a1 <- df[[paste0(locus, ".1")]]
  a2 <- df[[paste0(locus, ".2")]]
  
  a1 <- as.character(a1)
  a2 <- as.character(a2)
  
  miss <- is.na(a1) | is.na(a2) | a1 == "" | a2 == "" | a1 == "NA" | a2 == "NA"
  
  vals <- c(a1[!miss], a2[!miss])
  width <- max(nchar(vals), na.rm = TRUE)
  
  a1_pad <- ifelse(miss, NA, sprintf(paste0("%0", width, "d"), as.integer(a1)))
  a2_pad <- ifelse(miss, NA, sprintf(paste0("%0", width, "d"), as.integer(a2)))
  
  out <- ifelse(miss, NA, as.numeric(paste0(a1_pad, a2_pad)))
  out
}

pop_levels <- sort(unique(dat_sub$pop))

hf_dat <- data.frame(
  pop = as.numeric(factor(dat_sub$pop, levels = pop_levels))
)

for (loc in locus_bases) {
  hf_dat[[loc]] <- make_hf_locus(dat_sub, loc)
}

# optional: drop loci that are entirely NA
all_na_loci <- names(hf_dat)[colSums(!is.na(hf_dat)) == 0]
if (length(all_na_loci) > 0) {
  hf_dat <- hf_dat[, !names(hf_dat) %in% all_na_loci, drop = FALSE]
}

# -----------------------------
# 6) FST matrix
# -----------------------------
fst_mat <- pairwise.WCfst(hf_dat)
fst_mat <- as.matrix(fst_mat)

rownames(fst_mat) <- colnames(fst_mat) <- pop_levels

fst_mat[is.na(fst_mat)] <- 0
fst_mat[fst_mat < 0] <- 0
diag(fst_mat) <- 0

# -----------------------------
# 7) coordinates
# reconstruct linear stream positions
# from Section x 20 m
# -----------------------------
coords_df <- dat_sub %>%
  group_by(pop, River, Section) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(River, Section) %>%
  mutate(
    site_id = 1:n(),
    dist_m = Section * 20
  )

stream_centers <- data.frame(
  River = c("WEST.BROOK", "WB.JIMMY", "WB.MITCHELL", "WB.OBEAR"),
  lat0  = c(42.4400, 42.4450, 42.4430, 42.4420),
  lon0  = c(-72.6700, -72.6650, -72.6680, -72.6660)
)

coords_df <- coords_df %>%
  left_join(stream_centers, by = "River") %>%
  mutate(
    lat = lat0 + (dist_m / 1000) * 0.005,
    lon = lon0
  ) %>%
  select(site_id, pop, lat, lon)

# -----------------------------
# 8) match FST order
# -----------------------------
common <- intersect(rownames(fst_mat), coords_df$pop)

fst_mat <- fst_mat[common, common, drop = FALSE]

coords_df <- coords_df %>%
  filter(pop %in% common) %>%
  arrange(match(pop, common))

rownames(fst_mat) <- colnames(fst_mat) <- as.character(coords_df$site_id)

# -----------------------------
# 9) geographic distance
# -----------------------------
coords_mat <- as.matrix(coords_df[, c("lon", "lat")])

geo_dist <- geosphere::distm(coords_mat, fun = geosphere::distGeo) / 1000
rownames(geo_dist) <- colnames(geo_dist) <- as.character(coords_df$site_id)

# -----------------------------
# 10) IBD
# -----------------------------
ibd_df <- data.frame(
  fst  = fst_mat[upper.tri(fst_mat)],
  dist = geo_dist[upper.tri(geo_dist)]
)

ibd_plot <- ggplot(ibd_df, aes(dist, fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(
    title = paste0("SAFO-6 IBD (West Brook), cohort ", target_cohort),
    x = "Distance (km)",
    y = "FST"
  )

print(ibd_plot)

# -----------------------------
# 11) save
# -----------------------------
SAFO_6_fst <- fst_mat
SAFO_6_coords <- coords_df[, c("site_id", "lat", "lon")]

save(
  SAFO_6_fst,
  SAFO_6_coords,
  file = file.path(base_dir, "data", "SAFO-6.RData")
)