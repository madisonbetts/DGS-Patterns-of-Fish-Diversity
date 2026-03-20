# ------------------
# MAAE-1
# Macrhybopsis aestivalis
# Speckled chub
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
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/MAAE-1"

vcf_path <- file.path(
  save_dir,
  "doi_10_5061_dryad_4b8gthtb7__v20201014",
  "Final.Ma.hap.vcf"
)

meta_path <- file.path(
  save_dir,
  "doi_10_5061_dryad_4b8gthtb7__v20201014",
  "Dryad_Not_jem_Mac_aest_usats_mtDNA.xlsx"
)

xy_path <- file.path(save_dir, "M_aestivalis_xys.csv")

# -----------------------------
# site key
# -----------------------------
# 1 = Willow Creek
# 2 = Six Mile Draw
# 3 = Crockett Draw
# 4 = US Hwy 70
# 5 = US Hwy 380

# -----------------------------
# 0) site coordinates
# note: csv column is `site`
# -----------------------------
MAAE_1_coords <- read.csv(xy_path, stringsAsFactors = FALSE) %>%
  mutate(
    site = dplyr::case_when(
      site == "Willow Creek"  ~ "1",
      site == "Six Mile Draw" ~ "2",
      site == "Crockett Draw" ~ "3",
      site == "US Hwy 70"     ~ "4",
      site == "US Hwy 380"    ~ "5",
      TRUE                    ~ NA_character_
    )
  ) %>%
  filter(!is.na(site)) %>%
  select(site, lat, lon) %>%
  mutate(site = factor(site, levels = as.character(1:5))) %>%
  arrange(site) %>%
  mutate(site = as.character(site))

stopifnot(!any(is.na(MAAE_1_coords$site)))

# -----------------------------
# 1) ddRAD metadata
# keep Pecos individuals only
# -----------------------------
ma_meta <- read_excel(meta_path, sheet = "Mac_aest_ddRAD_samples") %>%
  rename(ind = INDV) %>%
  mutate(
    `Collection Site` = trimws(`Collection Site`),
    site_num = unname(c(
      "Willow"   = "1",
      "Six Mile" = "2",
      "Crockett" = "3",
      "Hwy 82"   = "4",  # assumed typo in source; corresponds to US Hwy 70
      "Hwy 380"  = "5"
    )[`Collection Site`])
  ) %>%
  filter(River_Year == "Pecos 2017")

stopifnot(!any(is.na(ma_meta$site_num)))

# -----------------------------
# 2) read vcf -> genlight
# -----------------------------
m_aestivalis <- vcfR2genlight(read.vcfR(vcf_path))

# -----------------------------
# 3) keep only individuals in metadata
# and align metadata to genlight order
# -----------------------------
m_aestivalis <- m_aestivalis[m_aestivalis@ind.names %in% ma_meta$ind]

ma_meta <- ma_meta %>%
  slice(match(m_aestivalis@ind.names, ind))

stopifnot(identical(m_aestivalis@ind.names, ma_meta$ind))

# assign populations
pop(m_aestivalis) <- factor(ma_meta$site_num, levels = as.character(1:5))

# -----------------------------
# 4) assign individual coordinates to genlight
# gl.ibd wants coords for EACH INDIVIDUAL in @other$latlon
# -----------------------------
m_aestivalis@other$latlon <- ma_meta %>%
  transmute(site = site_num) %>%
  left_join(MAAE_1_coords, by = "site") %>%
  select(lon, lat) %>%
  as.data.frame()

rownames(m_aestivalis@other$latlon) <- indNames(m_aestivalis)

stopifnot(nrow(m_aestivalis@other$latlon) == nInd(m_aestivalis))
stopifnot(!any(is.na(m_aestivalis@other$latlon$lon)))
stopifnot(!any(is.na(m_aestivalis@other$latlon$lat)))

# -----------------------------
# 5) map of sampling locations
# -----------------------------
ggplot(MAAE_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.03, size = 4) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "MAAE-1 sampling locations"
  )

# -----------------------------
# 6) isolation by distance
# uses individual coordinates stored in @other$latlon
# -----------------------------
gl.ibd(
  x = m_aestivalis,
  plot.out = TRUE
)

# -----------------------------
# 7) pairwise FST matrix
# set negative FST values to 0
# -----------------------------
MAAE_1_fst <- gl.fst.pop(m_aestivalis)
MAAE_1_fst[MAAE_1_fst < 0] <- 0

# -----------------------------
# fill upper triangle + zero diagonal
# -----------------------------
MAAE_1_fst <- MAAE_1_fst %>%
  {
    m <- as.matrix(.)
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
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
  MAAE_1_coords,
  MAAE_1_fst,
  file = file.path(out_dir, "MAAE-1.RData")
)