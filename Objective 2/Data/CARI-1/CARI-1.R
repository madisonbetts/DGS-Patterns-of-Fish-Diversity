# -------------------------
# Klamath smallscale sucker
# site metadata + FST workflow
# genotype file contains 14 populations,
# corresponding to the 14 sites sampled in 2013
# -------------------------
library(ggplot2)
library(adegenet)
library(hierfstat)
library(dplyr)
library(stringr)

# -------------------------
# 0) 2013 site metadata only
#    these are the sites represented
#    in the genotype file
# -------------------------
kss_sites <- data.frame(
  site = c("N1", "MN2", "MN3", "MN4",
           "L5", "L6", "L7", "L8",
           "L9", "L10", "L11",
           "S4", "S3", "S2"),
  year = 2013,
  n = c(77, 37, 96, 66,
        55, 54, 79, 96,
        16, 53, 20,
        23, 32, 42),
  lat = c(41.97773, 41.84289, 41.84417, 41.84199,
          41.80046, 41.80988, 41.81364, 41.81574,
          41.85460, 41.87092, 41.87801,
          41.73729, 41.73560, 41.70272),
  lon = c(-123.96172, -123.99501, -124.00387, -124.01427,
          -124.05446, -124.08224, -124.08655, -124.09945,
          -124.12184, -124.12592, -124.13308,
          -123.98290, -123.98335, -123.93401),
  stringsAsFactors = FALSE
)

kss_sites

# sanity plot
ggplot(kss_sites, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.01, size = 4) +
  coord_fixed() +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude")

# -------------------------
# 1) genotype file path
# -------------------------
file_path <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/CARI-1/10592_2024_1651_MOESM2_ESM.txt"

# -------------------------
# 2) read raw lines
# -------------------------
x <- readLines(file_path)

# -------------------------
# 3) find population blocks
# -------------------------
sample_start <- grep('^\\s*SampleName = ', x)
sample_names <- str_match(x[sample_start], 'SampleName = "(Pop[0-9]+)"')[, 2]
sample_end <- c(sample_start[-1] - 1, length(x))

# -------------------------
# 4) hardcoded crosswalk:
#    genotype populations -> 2013 sites
#    matched by sample order and N
# -------------------------
pop_crosswalk <- data.frame(
  pop = c("Pop001","Pop002","Pop003","Pop004","Pop005","Pop006","Pop007",
          "Pop008","Pop009","Pop010","Pop011","Pop012","Pop013","Pop014"),
  site = c("N1","MN2","MN3","MN4","L5","L6","L7",
           "L8","L9","L10","L11","S4","S3","S2"),
  expected_n = c(77, 37, 96, 66, 55, 54, 79,
                 96, 16, 53, 20, 23, 32, 42),
  stringsAsFactors = FALSE
)

pop_crosswalk

# -------------------------
# 5) parse each individual
#    each genotype spans 2 lines
# -------------------------
parsed_list <- list()

for(i in seq_along(sample_start)) {
  
  block <- x[sample_start[i]:sample_end[i]]
  pop_i <- sample_names[i]
  
  line1_idx <- grep('^\\s*[0-9]{2}_SMI-', block)
  
  for(j in line1_idx) {
    
    l1 <- str_squish(block[j])
    l2 <- str_squish(block[j + 1])
    
    tok1 <- str_split(l1, "\\s+")[[1]]
    tok2 <- str_split(l2, "\\s+")[[1]]
    
    ind_id <- tok1[1]
    site_code <- str_extract(ind_id, "(?<=SMI-)[A-Z0-9]+")
    
    a1 <- tok1[-c(1, 2)]
    a2 <- tok2
    
    if(length(a1) != length(a2)) {
      stop(paste("Allele count mismatch for", ind_id))
    }
    
    locus_names <- paste0("L", seq_along(a1))
    
    geno_strings <- ifelse(
      a1 == "?" | a2 == "?",
      NA,
      paste(a1, a2, sep = "/")
    )
    
    df_i <- data.frame(
      pop = pop_i,
      ind = ind_id,
      site_code = site_code,
      stringsAsFactors = FALSE
    )
    
    df_i[locus_names] <- as.list(geno_strings)
    parsed_list[[length(parsed_list) + 1]] <- df_i
  }
}

geno_df <- bind_rows(parsed_list)

# -------------------------
# 6) check sample sizes
# -------------------------
pop_counts <- geno_df %>%
  count(pop, name = "observed_n") %>%
  left_join(pop_crosswalk, by = "pop")

pop_counts

stopifnot(all(pop_counts$observed_n == pop_counts$expected_n))

# -------------------------
# 7) assign actual site names
# -------------------------
geno_df <- geno_df %>%
  left_join(pop_crosswalk[, c("pop", "site")], by = "pop")

# quick check
geno_df %>%
  count(pop, site, site_code)

# -------------------------
# 8) convert to genind
#    use actual site labels
# -------------------------
loc_cols <- grep("^L[0-9]+$", names(geno_df), value = TRUE)

genind_obj <- df2genind(
  X = geno_df[, loc_cols],
  ploidy = 2,
  sep = "/",
  pop = as.factor(geno_df$site),
  ind.names = geno_df$ind,
  NA.char = NA
)

genind_obj

# -------------------------
# 9) pairwise FST
#    Weir & Cockerham
# -------------------------
hf_obj <- genind2hierfstat(genind_obj)

pairwise_fst <- pairwise.WCfst(hf_obj)
fst_mat <- as.matrix(pairwise_fst)
diag(fst_mat) <- 0

fst_mat

# -------------------------
# 10) reorder site metadata
#     to match fst matrix
# -------------------------
kss_sites_fst <- kss_sites[match(rownames(fst_mat), kss_sites$site), ]

stopifnot(identical(rownames(fst_mat), kss_sites_fst$site))

kss_sites_fst

# -------------------------
# 11) optional:
#     Euclidean geographic
#     distance matrix in km
# -------------------------
site_dist_km <- as.matrix(
  dist(kss_sites_fst[, c("lon", "lat")])
)

rownames(site_dist_km) <- kss_sites_fst$site
colnames(site_dist_km) <- kss_sites_fst$site

site_dist_km

# -------------------------
# 12) optional IBD plot
# -------------------------
ibd_df <- data.frame(
  dist = site_dist_km[upper.tri(site_dist_km)],
  fst  = fst_mat[upper.tri(fst_mat)]
)

plot(
  ibd_df$dist,
  ibd_df$fst,
  xlab = "Euclidean distance (decimal-degree units)",
  ylab = "Pairwise FST",
  pch = 21
)


# -------------------------
# 13) prepare outputs
# -------------------------

# rename objects to requested names
CARI_1_fst <- fst_mat

CARI_1_coords <- kss_sites_fst[, c("site", "lat", "lon")]

# -------------------------
# 14) save RData
# -------------------------

save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/CARI-1/data"

if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

save(
  CARI_1_fst,
  CARI_1_coords,
  file = file.path(save_dir, "CARI-1.RData")
)