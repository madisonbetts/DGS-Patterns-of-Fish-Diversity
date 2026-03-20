# ------------------
# NOJE-1
# Notropis jemezanus
# Rio Grande shiner
# ddRAD Pecos sites only
# ------------------

library(vcfR)
library(adegenet)
library(dartR)
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(magrittr)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/NOJE-1"

vcf_path <- file.path(
  save_dir,
  "doi_10_5061_dryad_4b8gthtb7__v20201014",
  "Final.Nj.mac.hap.vcf"
)

meta_path <- file.path(save_dir, "doi_10_5061_dryad_4b8gthtb7__v20201014//Dryad_Not_jem_Mac_aest_usats_mtDNA.xlsx")
xy_path   <- file.path(save_dir, "N. jemezanus_xys.csv")

# -----------------------------
# site key
# -----------------------------
# 1 = Willow Creek
# 2 = Gasline
# 3 = US Hwy 70
# 4 = Scout Camp
# 5 = US Hwy 380
# 6 = BLM ACEC

# -----------------------------
# 0) site coordinates
# note: csv column is `name`, not `site`
# -----------------------------
NOJE_1_coords <- read.csv(xy_path, stringsAsFactors = FALSE) %>%
  mutate(
    lon = as.numeric(str_match(WKT, "POINT \\(([-0-9.]+) ([-0-9.]+)\\)")[, 2]),
    lat = as.numeric(str_match(WKT, "POINT \\(([-0-9.]+) ([-0-9.]+)\\)")[, 3]),
    site = dplyr::case_when(
      name == "Willow Creek" ~ "1",
      name == "Gasline"      ~ "2",
      name == "US Hwy 70"    ~ "3",
      name == "Scout Camp"   ~ "4",
      name == "US Hwy 380"   ~ "5",
      name == "BLM ACEC"     ~ "6",
      TRUE                   ~ NA_character_
    )
  ) %>%
  select(site, lat, lon) %>%
  mutate(site = factor(site, levels = as.character(1:6))) %>%
  arrange(site) %>%
  mutate(site = as.character(site))

stopifnot(!any(is.na(NOJE_1_coords$site)))

# -----------------------------
# 1) ddRAD metadata
# keep Pecos individuals only
# -----------------------------
nj_meta <- read_excel(meta_path, sheet = "Not_jem_ddRAD_samples") %>%
  rename(ind = INDV_Index) %>%
  filter(River == "Pecos") %>%
  mutate(
    site_num = unname(c(
      "Willow"  = "1",
      "Gasline" = "2",
      "Hwy 70"  = "3",
      "Scout"   = "4",
      "Hwy 82"  = "5",   # typo in source; corresponds to US Hwy 380
      "BLMACEC" = "6"
    )[Site])
  )

stopifnot(!any(is.na(nj_meta$site_num)))

# -----------------------------
# 2) read vcf -> genlight
# -----------------------------
n_jemezanus <- vcfR2genlight(read.vcfR(vcf_path))

# -----------------------------
# 3) keep only individuals in metadata
# and align metadata to genlight order
# -----------------------------
n_jemezanus <- n_jemezanus[n_jemezanus@ind.names %in% nj_meta$ind]

nj_meta <- nj_meta %>%
  slice(match(n_jemezanus@ind.names, ind))

stopifnot(identical(n_jemezanus@ind.names, nj_meta$ind))

# assign populations
pop(n_jemezanus) <- factor(nj_meta$site_num, levels = as.character(1:6))

# -----------------------------
# 4) assign individual coordinates to genlight
# gl.ibd wants coords for EACH INDIVIDUAL in @other$latlon
# -----------------------------
n_jemezanus@other$latlon <- nj_meta %>%
  transmute(site = site_num) %>%
  left_join(NOJE_1_coords, by = "site") %>%
  select(lon, lat) %>%
  as.data.frame()

rownames(n_jemezanus@other$latlon) <- indNames(n_jemezanus)

stopifnot(nrow(n_jemezanus@other$latlon) == nInd(n_jemezanus))
stopifnot(!any(is.na(n_jemezanus@other$latlon$lon)))
stopifnot(!any(is.na(n_jemezanus@other$latlon$lat)))

# -----------------------------
# 5) map of sampling locations
# -----------------------------
ggplot(NOJE_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.03, size = 4) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "NOJE-1 sampling locations"
  )

# -----------------------------
# 6) isolation by distance
# uses individual coordinates stored in @other$latlon
# -----------------------------
gl.ibd(
  x = n_jemezanus,
  plot.out = TRUE
)

# -----------------------------
# 7) pairwise FST matrix
# set negative FST values to 0
# -----------------------------
NOJE_1_fst <- gl.fst.pop(n_jemezanus)
NOJE_1_fst[NOJE_1_fst < 0] <- 0

# -----------------------------
# fill upper triangle + zero diagonal
# -----------------------------
NOJE_1_fst <- NOJE_1_fst %>%
  { 
    m <- as.matrix(.)
    
    # mirror lower -> upper
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    
    # set diagonal to 0
    diag(m) <- 0
    
    m
  }

# -----------------------------
# 8) save RData
# final objects go in the data subfolder
# -----------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  NOJE_1_coords,
  NOJE_1_fst,
  file = file.path(out_dir, "NOJE-1.RData")
)