# -----------------------------
# CASA-1 Santa Ana sucker
# Richmond et al. 2018 / USGS data release
# Microsatellite-based pairwise FST + site coordinates + map + IBD plot
# -----------------------------

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(hierfstat)
library(geosphere)
library(ggplot2)
library(maps)

# -----------------------------
# paths
# -----------------------------
setwd("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/CASA-1/Data Package")

full_file <- "Full dataset_final.csv"
pop_file  <- "Population identifiers.csv"
xy_file   <- "CASA sampling points.csv"

# -----------------------------
# helper to parse CONVERT-style microsat file
# -----------------------------
parse_convert_microsat <- function(file) {
  raw <- read_csv(file, col_names = FALSE, show_col_types = FALSE, na = c("", "NA"))
  raw <- as.data.frame(raw, stringsAsFactors = FALSE)

  # locus names are on row 3 in every other column starting at V2
  locus_row <- raw[3, ]
  locus_names <- as.character(unlist(locus_row[-1]))
  locus_names <- locus_names[seq(1, length(locus_names), by = 2)]
  locus_names <- locus_names[!is.na(locus_names)]

  out_list <- list()
  current_pop <- NA_character_

  for (i in seq_len(nrow(raw))) {
    first_val <- raw[i, 1][[1]]

    if (is.na(first_val)) next

    first_val <- trimws(as.character(first_val))

    if (first_val == "" || grepl("^npops=", first_val) || grepl("^nloci=", first_val)) next

    if (grepl("^pop\\s*=", first_val)) {
      current_pop <- trimws(sub("^pop\\s*=", "", first_val))
      next
    }

    # skip locus-name row and any explanatory rows
    if (first_val %in% locus_names || grepl("^row ", first_val)) next

    geno_vals <- as.character(unlist(raw[i, -1]))
    geno_vals <- geno_vals[seq_len(length(locus_names) * 2)]
    geno_vals[geno_vals == "?"] <- NA

    out_list[[length(out_list) + 1]] <- data.frame(
      indiv = first_val,
      pop_code = current_pop,
      matrix(geno_vals, nrow = 1, dimnames = list(NULL, as.vector(rbind(paste0(locus_names, ".1"), paste0(locus_names, ".2"))))),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  bind_rows(out_list)
}

# -----------------------------
# read and parse genotype data
# -----------------------------
microsat_df <- parse_convert_microsat(full_file)

# convert allele columns to numeric where possible
allele_cols <- setdiff(names(microsat_df), c("indiv", "pop_code"))
microsat_df[allele_cols] <- lapply(microsat_df[allele_cols], function(x) as.numeric(x))

# -----------------------------
# population names
# -----------------------------
pop_key <- read_csv(pop_file, show_col_types = FALSE) %>%
  rename(pop_code = 1, location = 2) %>%
  filter(!is.na(pop_code))

# -----------------------------
# generalized sampling coordinates
# exact GPS were withheld in the XML metadata, so use the generalized
# sampling points supplied in CASA sampling points.csv and collapse them
# to one representative XY per sampled population
# -----------------------------
xy_raw <- read_csv(xy_file, show_col_types = FALSE) %>%
  rename(site_no = 1, lon = 2, lat = 3, pop_code = 4, location = 5)

xy_df <- xy_raw %>%
  group_by(pop_code) %>%
  summarise(
    lon = mean(lon, na.rm = TRUE),
    lat = mean(lat, na.rm = TRUE),
    .groups = "drop"
  )

# add San Dimas Canyon approximate coordinates from Fig. 1 because the
# generalized points file has no SDC coordinate even though SDC is a sampled population.
# this estimate places SDC east of the main San Gabriel River localities,
# consistent with the study map and label position.
if (!"SDC" %in% xy_df$pop_code) {
  xy_df <- bind_rows(
    xy_df,
    data.frame(
      pop_code = "SDC",
      lon = -117.9070,
      lat = 34.1400
    )
  )
}

# -----------------------------
# final site order
# keep same biological order as paper / dataset
# -----------------------------
site_key <- tibble(
  site = 1:8,
  pop_code = c("LSCR", "SFC", "VAL", "BTC", "HAC", "SGR", "SDC", "SAR")
) %>%
  left_join(pop_key, by = "pop_code") %>%
  left_join(xy_df, by = "pop_code")

CASA_1_sites <- site_key %>%
  select(site, lat, lon)

# -----------------------------
# pairwise FST from microsatellite genotypes
# hierfstat expects first column = population code numeric
# -----------------------------
microsat_hf <- microsat_df %>%
  left_join(site_key %>% select(site, pop_code), by = "pop_code") %>%
  select(site, all_of(allele_cols))

# remove any individuals missing site assignment
microsat_hf <- microsat_hf %>% filter(!is.na(site))

CASA_1_fst <- pairwise.WCfst(microsat_hf)
CASA_1_fst[CASA_1_fst < 0] <- 0
CASA_1_fst <- as.matrix(CASA_1_fst)
rownames(CASA_1_fst) <- colnames(CASA_1_fst) <- site_key$site

# -----------------------------
# pairwise geographic distances for IBD
# -----------------------------
coords <- CASA_1_sites %>% select(lon, lat)
geo_dist_km <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
rownames(geo_dist_km) <- colnames(geo_dist_km) <- CASA_1_sites$site

ibd_df <- data.frame(
  site1   = rownames(CASA_1_fst)[row(CASA_1_fst)[upper.tri(CASA_1_fst)]],
  site2   = colnames(CASA_1_fst)[col(CASA_1_fst)[upper.tri(CASA_1_fst)]],
  fst     = CASA_1_fst[upper.tri(CASA_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# map of sampling sites
# -----------------------------
world_map <- map_data("world")
usa_map   <- subset(world_map, region %in% c("USA", "Canada", "Mexico"))

map_df <- site_key %>%
  mutate(label = as.character(site))

p_map <- ggplot() +
  geom_polygon(
    data = usa_map,
    aes(x = long, y = lat, group = group),
    fill = "grey92",
    color = "grey55",
    linewidth = 0.2
  ) +
  geom_point(
    data = map_df,
    aes(x = lon, y = lat),
    color = "red3",
    size = 3
  ) +
  geom_text(
    data = map_df,
    aes(x = lon, y = lat, label = label),
    nudge_y = 0.07,
    size = 3.5
  ) +
  coord_fixed(
    ratio = 1.3,
    xlim = range(map_df$lon) + c(-0.5, 0.5),
    ylim = range(map_df$lat) + c(-0.5, 0.5)
  )
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "CASA-1 sampling sites"
  )

# -----------------------------
# IBD plot
# -----------------------------
p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Geographic distance (km)",
    y = expression(F[ST]),
    title = "CASA-1 isolation by distance"
  )


# -----------------------------
# save RData
# -----------------------------
diag(CASA_1_fst) <- 0 # make diags zero

save(
  CASA_1_fst,
  CASA_1_sites,
  file = "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/CASA-1/data/CASA-1.RData"
)

# -----------------------------
# print outputs
# -----------------------------
print(CASA_1_sites)
print(round(CASA_1_fst, 4))
print(head(ibd_df))
