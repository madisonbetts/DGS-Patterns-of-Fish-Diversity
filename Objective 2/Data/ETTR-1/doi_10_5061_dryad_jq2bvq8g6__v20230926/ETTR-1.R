# -----------------------------
# ETTR-1
# Etheostoma trisella
# trispot darter
# -----------------------------

# -----------------------------
# 0) setup
# -----------------------------
library(dplyr)
library(ggplot2)
library(geosphere)
library(hierfstat)

base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETTR-1"

data_dir <- file.path(base_dir, "data")
dir.create(data_dir, showWarnings = FALSE)

setwd(file.path(
  base_dir,
  "doi_10_5061_dryad_jq2bvq8g6__v20230926"
))

geno_file <- "Etheostoma_trisella_SNP_report_2023.csv"

if (!file.exists(geno_file)) {
  stop("Genotype file not found: ", geno_file)
}

# -----------------------------
# 1) read genotype data
# real header begins after 6 lines
# -----------------------------
raw <- read.csv(
  geno_file,
  skip = 6,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# debug
stopifnot(nrow(raw) > 0)

# -----------------------------
# 2) identify metadata + sample genotype columns
# DArT marker metadata occupy cols 1:18
# sample genotype columns begin at col 19
# -----------------------------
meta_cols   <- colnames(raw)[1:18]
sample_cols <- colnames(raw)[19:ncol(raw)]

# debug
print(meta_cols)
print(head(sample_cols))

# -----------------------------
# 3) parse sample metadata from sample column names
# sample names look like:
# 1694_Little Canoe Creek
# 1711_Coosawattee River
# -----------------------------
sample_info <- data.frame(
  sample_col = sample_cols,
  stringsAsFactors = FALSE
) %>%
  mutate(
    sample_id   = sub("_.*$", "", sample_col),
    population  = sub("^[^_]+_", "", sample_col)
  )

print(table(sample_info$population))

# -----------------------------
# 4) transpose to sample x locus matrix
# loci are rows, samples are columns in raw
# -----------------------------
geno_mat <- t(as.matrix(raw[, sample_cols, drop = FALSE]))
geno_mat <- as.data.frame(geno_mat, stringsAsFactors = FALSE)

# use AlleleID if present, otherwise fall back to first metadata column
if ("AlleleID" %in% colnames(raw)) {
  locus_names <- raw$AlleleID
} else {
  locus_names <- raw[[1]]
}

colnames(geno_mat) <- make.unique(as.character(locus_names))
rownames(geno_mat) <- sample_cols

# -----------------------------
# 5) convert genotype matrix to numeric
# only sample genotype block is converted
# -----------------------------
geno_mat[] <- lapply(geno_mat, function(x) as.numeric(as.character(x)))

# bind metadata
dat <- cbind(sample_info, geno_mat)

# -----------------------------
# 6) define populations
# -----------------------------
dat$pop <- as.factor(dat$population)

print(table(dat$pop))

# -----------------------------
# 7) build hierfstat dataframe
# first column = numeric population code
# rest = SNP loci
# -----------------------------
snp_cols <- colnames(geno_mat)

pop_numeric <- as.numeric(dat$pop)

hf_dat <- data.frame(
  pop = pop_numeric,
  dat[, snp_cols, drop = FALSE]
)

# remove loci that are all NA
keep_loci <- colSums(!is.na(hf_dat[, -1, drop = FALSE])) > 0
hf_dat <- data.frame(
  pop = hf_dat$pop,
  hf_dat[, -1, drop = FALSE][, keep_loci, drop = FALSE]
)

# debug
stopifnot(nrow(hf_dat) > 0)
stopifnot(ncol(hf_dat) > 10)

# -----------------------------
# 8) compute pairwise FST
# -----------------------------
fst_res <- pairwise.WCfst(hf_dat)
fst_res <- as.matrix(fst_res)

pop_levels <- levels(dat$pop)
rownames(fst_res) <- colnames(fst_res) <- pop_levels

fst_res[is.na(fst_res)] <- 0
fst_res[fst_res < 0] <- 0
diag(fst_res) <- 0

print(fst_res)

# -----------------------------
# 9) build coordinate dataframe
# best-available creek-system centroids from paper map
# creek systems used in SNP analyses:
# Little Canoe Creek, Ballplay Creek, Coahulla Creek,
# Mill Creek, Coosawattee River
# -----------------------------
coords_lookup <- data.frame(
  pop = c(
    "Little Canoe Creek",
    "Ballplay Creek",
    "Coahulla Creek",
    "Mill Creek",
    "Coosawattee River"
  ),
  lat = c(
    33.9395,  # Little Canoe Creek
    34.0905,  # Ballplay Creek
    34.8410,  # Coahulla Creek
    34.8920,  # Mill Creek
    34.6750   # Coosawattee River
  ),
  lon = c(
    -86.0220, # Little Canoe Creek
    -85.7890, # Ballplay Creek
    -84.8960, # Coahulla Creek
    -84.8270, # Mill Creek
    -84.7380  # Coosawattee River
  ),
  stringsAsFactors = FALSE
)

coords_df <- data.frame(
  pop = pop_levels,
  stringsAsFactors = FALSE
) %>%
  left_join(coords_lookup, by = "pop") %>%
  mutate(site = seq_len(n())) %>%
  select(site, pop, lat, lon)

print(coords_df)

# -----------------------------
# 10) geographic distance matrix
# -----------------------------
coords_mat <- as.matrix(coords_df[, c("lon", "lat")])
dist_mat <- geosphere::distm(coords_mat, fun = geosphere::distHaversine) / 1000

rownames(dist_mat) <- colnames(dist_mat) <- as.character(coords_df$site)

# -----------------------------
# 11) IBD dataframe
# -----------------------------
fst <- fst_res[coords_df$pop, coords_df$pop, drop = FALSE]
rownames(fst) <- colnames(fst) <- as.character(coords_df$site)

upper_idx <- upper.tri(fst)

ibd_df <- data.frame(
  fst  = fst[upper_idx],
  dist = dist_mat[upper_idx]
)

# -----------------------------
# 12) IBD plot
# -----------------------------
p_ibd <- ggplot(ibd_df, aes(x = dist, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  theme_bw() +
  labs(
    title = "ETTR-1 isolation by distance",
    x = "Geographic distance (km)",
    y = "FST"
  )

print(p_ibd)

# -----------------------------
# 13) map of sites
# -----------------------------
p_map <- ggplot(coords_df, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.03) +
  theme_bw() +
  labs(
    title = "ETTR-1 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )

print(p_map)

# -----------------------------
# 14) format outputs
# -----------------------------
ETTR_1_fst <- fst
ETTR_1_coords <- coords_df[, c("site", "lat", "lon")]

# -----------------------------
# 15) save
# -----------------------------
save(
  ETTR_1_fst,
  ETTR_1_coords,
  file = file.path(data_dir, "ETTR-1.RData")
)