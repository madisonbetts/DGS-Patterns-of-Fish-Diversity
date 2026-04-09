# -----------------------------
# SAFO-9 | Tennessee Brook Trout
# Salvelinus fontinalis
# Hargrove et al. 2022, Conservation Genetics
#
# Workflow:
# - read multilocus microsatellite genotypes from Supplemental File 1
# - read site coordinates from Supplemental File 2
# - link individuals to populations using Collection ID parsed from sample IDs
# - calculate pairwise Weir & Cockerham FST with hierfstat
# - set negative FST values to 0
# - build matching coords dataframe with numeric site IDs
# - plot sites and an IBD scatter in RStudio
# - save SAFO_9_fst and SAFO_9_coords to data/SAFO-9.RData
# -----------------------------

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(adegenet)
library(hierfstat)
library(geosphere)
library(ggplot2)
library(maps)

# -----------------------------
# 0) paths
# -----------------------------
geno_file <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SAFO-9/10592_2021_1404_MOESM1_ESM.xlsx"
site_file <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SAFO-9/10592_2021_1404_MOESM2_ESM.xlsx"

out_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SAFO-9/data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) read coordinates / site metadata
# Supplemental File 2
# keep only site, lat, lon in final coords df
# -----------------------------
site_df <- read_excel(site_file, sheet = 1, skip = 1) %>%
  transmute(
    site_name = `Site name`,
    collection_id = `Collection ID`,
    lat = as.numeric(Latitude),
    lon = as.numeric(Longitude)
  ) %>%
  mutate(
    collection_id = paste0("TN", as.integer(str_extract(collection_id, "\\d+")))
  ) %>%
  arrange(as.integer(str_extract(collection_id, "\\d+"))) %>%
  mutate(site = seq_len(n()))

# optional filter if you ever want to exclude Rock Creek because the paper
# notes it represented a single full-sibling family:
# site_df <- site_df %>% filter(collection_id != "TN33") %>% mutate(site = seq_len(n()))

# coords object for saving
SAFO_9_coords <- site_df %>%
  select(site, lat, lon)

# -----------------------------
# 2) read multilocus genotypes
# Supplemental File 1 is in GenAlEx-style layout
# -----------------------------
raw_geno <- read_excel(geno_file, sheet = 1, skip = 2)

# keep sample id, population label, and allele columns only
geno <- raw_geno %>%
  select(1:26) %>%
  rename(sample = 1, pop_name = 2) %>%
  mutate(
    collection_id = paste0("TN", as.integer(str_extract(sample, "(?<=TN)0*\\d+(?=-)")))
  ) %>%
  filter(!is.na(collection_id)) %>%
  inner_join(site_df %>% select(site, collection_id, site_name), by = "collection_id")

# check that all populations in coords are represented in genotype data
missing_pops <- setdiff(site_df$collection_id, unique(geno$collection_id))
if (length(missing_pops) > 0) {
  stop("Some populations in the site metadata were not found in the genotype file: ",
       paste(missing_pops, collapse = ", "))
}

# -----------------------------
# 3) build a genind object
# convert each allele pair to 'a/b' character format
# -----------------------------
allele_cols <- names(geno)[3:26]
loci <- unique(str_replace(allele_cols, " [12]$", ""))

geno_loci <- vector("list", length(loci))
names(geno_loci) <- loci

for (loc in loci) {
  a1 <- geno[[paste0(loc, " 1")]]
  a2 <- geno[[paste0(loc, " 2")]]
  geno_loci[[loc]] <- ifelse(is.na(a1) | is.na(a2), NA, paste0(a1, "/", a2))
}

geno_loci <- as.data.frame(geno_loci, stringsAsFactors = FALSE)

genind_obj <- df2genind(
  X = geno_loci,
  sep = "/",
  ncode = NULL,
  ploidy = 2,
  pop = as.factor(geno$site),
  ind.names = geno$sample
)

# -----------------------------
# 4) pairwise FST
# Weir & Cockerham pairwise FST via hierfstat
# -----------------------------
hf_dat <- genind2hierfstat(genind_obj)

fst_mat <- pairwise.WCfst(hf_dat)

# enforce full matrix / correct ordering
site_order <- as.character(site_df$site)
fst_mat <- as.matrix(fst_mat)[site_order, site_order, drop = FALSE]
diag(fst_mat) <- 0
fst_mat[is.na(fst_mat)] <- 0
fst_mat[fst_mat < 0] <- 0

rownames(fst_mat) <- site_order
colnames(fst_mat) <- site_order

SAFO_9_fst <- fst_mat

# -----------------------------
# 5) quick site lookup for plotting labels
# -----------------------------
site_lookup <- site_df %>%
  select(site, site_name, collection_id, lat, lon)

# -----------------------------
# 6) map plot
# Plot US and Canada and zoom to the extent of the points
# -----------------------------
world_df <- map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

x_pad <- max(0.5, diff(range(site_lookup$lon)) * 0.15)
y_pad <- max(0.5, diff(range(site_lookup$lat)) * 0.15)

ggplot() +
  geom_polygon(
    data = world_df,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = site_lookup,
    aes(x = lon, y = lat),
    size = 2
  ) +
  geom_text(
    data = site_lookup,
    aes(x = lon, y = lat, label = site),
    nudge_y = y_pad * 0.08,
    size = 2.5
  ) +
  coord_quickmap(
    xlim = c(min(site_lookup$lon) - x_pad, max(site_lookup$lon) + x_pad),
    ylim = c(min(site_lookup$lat) - y_pad, max(site_lookup$lat) + y_pad)
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()

# -----------------------------
# 7) IBD plot
# straight-line geographic distance as pilot distance layer
# -----------------------------
geo_dist <- geosphere::distm(
  x = as.matrix(site_lookup[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist[upper.tri(geo_dist)],
  fst = SAFO_9_fst[upper.tri(SAFO_9_fst)]
)

ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 1.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Geographic distance (km)", y = "Pairwise FST") +
  theme_bw()

# -----------------------------
# 8) save objects
# -----------------------------
save(
  SAFO_9_fst,
  SAFO_9_coords,
  file = file.path(out_dir, "SAFO-9.RData")
)
