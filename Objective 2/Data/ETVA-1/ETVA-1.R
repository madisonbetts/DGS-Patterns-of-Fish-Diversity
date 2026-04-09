# -----------------------------
# ETVA-1 | Variegate darter
# Etheostoma variatum
# Argentina et al. 2018, Freshwater Biology
#
# Workflow:
# - use Appendix S1 coordinates from the paper as the georeferencing source
# - read individual microsatellite genotypes from AllVDmicrosats_NoHybrids_M13Corrr_20180118.xlsx
# - link spreadsheet locality codes to paper site numbers
# - allow inferred many-to-one joins where spreadsheet localities were split across batches
# - calculate pairwise Weir & Cockerham FST with hierfstat
# - set negative FST values to 0
# - build matching numeric-site coordinate dataframe
# - plot sites and IBD in RStudio
# - save ETVA_1_fst and ETVA_1_coords to data/ETVA-1.RData
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

study_code <- "ETVA-1"

# -----------------------------
# 0) paths
# -----------------------------
geno_file <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETVA-1/AllVDmicrosats_NoHybrids_M13Corrr_20180118.xlsx"

out_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETVA-1/data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) paper Appendix S1 site table
# use the site coordinates reported in the paper / appendix
# site 22-24 were pooled in the paper analyses and are treated as one site group here
# -----------------------------
paper_sites <- tibble::tribble(
  ~site, ~paper_label, ~site_name, ~lat, ~lon, ~paper_n,
   1, "1",     "Bell Run at State Highway 44", 41.97970, -78.23660, 12,
   2, "2",     "Allegheny River downstream of Union Street bridge, Olean, NY", 42.07052, -78.42996, 20,
   3, "3",     "Stillwater Creek downstream of Bacon Road near Jamestown, PA", 42.04372, -79.23126, 20,
   4, "4",     "Brokenstraw Creek, Buckaloons Recreation Area", 41.83528, -79.25975, 30,
   5, "5",     "French Creek at NY and PA state line", 42.01589, -79.50438, 19,
   6, "6",     "South Branch French Creek, VFW Park, Union City, PA", 41.90269, -79.85466, 33,
   7, "7",     "South Branch French Creek at Bridge Street", 41.89770, -79.83500, 8,
   8, "8",     "Mahoning Creek at State Route 1025", 40.94880, -79.31980, 8,
   9, "9",     "Cheat River at Holly Meadows Road", 39.12250, -79.67510, 4,
  10, "10",    "Middle Fork Beaver Creek at State Route 30", 40.77480, -80.77980, 4,
  11, "11",    "Back Fork Elk River near Webster Springs, WV", 38.47905, -80.38910, 23,
  12, "12",    "Birch River near Glendon, WV", 38.59257, -80.87021, 21,
  13, "13",    "Dry Fork downstream of Crane Creek confluence", 37.41460, -81.78312, 31,
  14, "14",    "Tug Fork at Horse Creek confluence", 37.46841, -81.86341, 17,
  15, "15",    "Dismal River", 37.24127, -81.87148, 2,
  16, "16",    "Dismal River off County Road 638", 37.24624, -82.02314, 2,
  17, "17",    "Slate Creek in Grundy, VA near confluence with Levisa Fork", 37.27739, -82.09370, 15,
  18, "18",    "Levisa Fork 500 m downstream of Highway 83", 37.23333, -82.10251, 19,
  19, "19",    "Levisa Fork upstream of Vansant, VA", 37.23682, -82.06055, 10,
  20, "20",    "Levisa Fork at County Road 604 near Harman Junction, VA", 37.30709, -82.15599, 33,
  21, "21",    "Levisa Fork upstream of Fishtrap Lake Wildlife Management Area", 37.39688, -82.25520, 40,
  22, "22_24", "Levisa Fork downstream of Fishtrap Lake near Shelbiana, KY", 37.41162, -82.42752, 57,
  25, "25",    "Levisa Fork near Pikeville, KY", 37.46396, -82.52574, 10,
  26, "26",    "Russell Fork in Regina, KY", 37.36608, -82.41172, 43,
  27, "27",    "Tygarts Creek at State Routes 7 and 784", 38.46640, -83.05080, 15,
  28, "28",    "Big Darby Creek, Pleasant Valley Metropark, Orient, OH", 39.80716, -83.15501, 34,
  29, "29",    "Salt Creek at Narrows Road", 39.34880, -82.67760, 9,
  30, "30",    "Little Miami River, State Route 3, Kings Mills", 39.37889, -84.25159, 30,
  31, "31",    "Salt Creek at Bull Fork Road", 39.40330, -85.20640, 5,
  32, "32",    "Whitewater River at St. Peters Road", 39.31440, -84.90500, 5,
  33, "33",    "Middle Fork Red River at State Route 77", 37.81510, -83.71840, 7,
  34, "34",    "Redbird River at Eriline Road", 37.18310, -83.60030, 10
)

# -----------------------------
# 2) spreadsheet locality -> paper site mapping
# these are the genotype locality prefixes after stripping trailing sample numbers
# several were split across spreadsheet batches and must be merged to a single paper site
# -----------------------------
site_map <- tibble::tribble(
  ~legacy_loc, ~site,
  "BellRun",   1,
  "Allegh",    2,
  "Stillwtr",  3,
  "brkStrw",   4,
  "French",    5,
  "PAFrench",  6,
  "41French",  7,
  "Mahoning",  8,
  "Cheat",     9,
  "MidlFkB",  10,
  "BFElk",    11,
  "Birch",    12,
  "DryF",     13,
  "TugHorse", 14,
  "Dismal",   15,
  "Ltwelve",  16,
  "Slate",    17,
  "SlateA",   17,
  "LFouteen", 18,
  "Lfive",    19,
  "LEight",   20,
  "LevWMA",   21,
  "LKenOne",  22,
  "LFPike",   25,
  "RussOne",  26,
  "Tygarts",  27,
  "BigDarby", 28,
  "SaltCrk",  29,
  "LtMiami",  30,
  "34Salt",   31,
  "Whitewatr",32,
  "RedR",     33,
  "Redbird",  34
) %>%
  left_join(paper_sites %>% select(site, paper_label, site_name, lat, lon, paper_n), by = "site")

# -----------------------------
# 3) read genotype file
# -----------------------------
raw_geno <- read_excel(geno_file, sheet = 1)

geno <- raw_geno %>%
  rename(sample = sampleID, site_sample_no = SiteSampleNo) %>%
  mutate(
    legacy_loc = str_remove(site_sample_no, "\\d+$")
  ) %>%
  left_join(site_map, by = "legacy_loc")

unmatched_locs <- geno %>%
  filter(is.na(site)) %>%
  distinct(legacy_loc) %>%
  pull(legacy_loc)

if (length(unmatched_locs) > 0) {
  stop("Some spreadsheet localities were not mapped to paper sites: ",
       paste(unmatched_locs, collapse = ", "))
}

# -----------------------------
# 4) count checks
# compare genotype totals against expected paper counts after aggregating
# split spreadsheet batches that map to the same paper site are summed first
# -----------------------------
count_check <- geno %>%
  count(site, name = "geno_n") %>%
  left_join(paper_sites %>% select(site, paper_n), by = "site") %>%
  arrange(site)

if (any(is.na(count_check$paper_n))) {
  stop("At least one mapped paper site is missing an expected paper sample size.")
}

bad_counts <- count_check %>%
  filter(geno_n != paper_n)

if (nrow(bad_counts) > 0) {
  stop(
    paste0(
      "Mapped genotype totals do not match expected paper sample sizes for: ",
      paste0("site ", bad_counts$site, " (geno=", bad_counts$geno_n,
             ", paper=", bad_counts$paper_n, ")", collapse = "; ")
    )
  )
}

message("All mapped genotype totals match the paper sample sizes after aggregation.")

# optional per-prefix check as a warning only
prefix_check <- geno %>%
  count(legacy_loc, site, name = "geno_n") %>%
  left_join(site_map %>% distinct(legacy_loc, site, paper_n), by = c("legacy_loc", "site")) %>%
  mutate(prefix_matches_paper = geno_n == paper_n)

if (any(!prefix_check$prefix_matches_paper)) {
  warning(
    "Some spreadsheet locality prefixes are split across the same paper site. ",
    "This is expected for combined sites such as Slate + SlateA and should not be treated as an error."
  )
}

# -----------------------------
# 5) georeferenced individual-level genotype table
# keep this object so individuals are tied to paper site coords
# -----------------------------
ETVA_1_geno_georef <- geno %>%
  select(sample, site_sample_no, legacy_loc, site, paper_label, site_name, lat, lon, everything())

# -----------------------------
# 6) build genind object
# sanitize locus names so df2genind does not choke on punctuation
# -----------------------------
meta_cols <- c("sample", "site_sample_no", "legacy_loc", "site", "paper_label", "site_name", "lat", "lon")
allele_cols <- setdiff(names(ETVA_1_geno_georef), meta_cols)

# keep only columns that look like allele columns
allele_cols <- allele_cols[str_detect(allele_cols, " - 1$| - 2$|- 1$|- 2$| 1$| 2$")]

raw_loci <- unique(str_replace(allele_cols, "( - 1| - 2|- 1|- 2| 1| 2)$", ""))
clean_loci <- raw_loci %>%
  str_replace_all("[^[:alnum:]_]", "_") %>%
  str_replace_all("_+", "_")

if (any(duplicated(clean_loci))) {
  dupes <- unique(clean_loci[duplicated(clean_loci)])
  stop("After cleaning locus names, duplicates were created: ", paste(dupes, collapse = ", "))
}

locus_key <- tibble(raw = raw_loci, clean = clean_loci)

geno_loci <- vector("list", length(clean_loci))
names(geno_loci) <- clean_loci

for (i in seq_len(nrow(locus_key))) {
  raw_loc <- locus_key$raw[i]
  clean_loc <- locus_key$clean[i]

  a1_col <- allele_cols[str_detect(allele_cols, paste0("^", stringr::fixed(raw_loc), "( - 1|- 1| 1)$"))]
  a2_col <- allele_cols[str_detect(allele_cols, paste0("^", stringr::fixed(raw_loc), "( - 2|- 2| 2)$"))]

  if (length(a1_col) != 1 || length(a2_col) != 1) {
    stop("Could not uniquely identify paired allele columns for locus: ", raw_loc)
  }

  a1 <- ETVA_1_geno_georef[[a1_col]]
  a2 <- ETVA_1_geno_georef[[a2_col]]

  geno_loci[[clean_loc]] <- ifelse(
    is.na(a1) | is.na(a2),
    NA,
    paste0(a1, "/", a2)
  )
}

geno_loci <- as.data.frame(geno_loci, stringsAsFactors = FALSE, check.names = FALSE)

genind_obj <- df2genind(
  X = geno_loci,
  sep = "/",
  ploidy = 2,
  pop = as.factor(ETVA_1_geno_georef$site),
  ind.names = ETVA_1_geno_georef$sample
)

# -----------------------------
# 7) pairwise FST
# Weir & Cockerham pairwise FST via hierfstat
# -----------------------------
hf_dat <- genind2hierfstat(genind_obj)
fst_mat <- pairwise.WCfst(hf_dat)

site_order <- sort(unique(ETVA_1_geno_georef$site))
site_order <- as.character(site_order)

fst_mat <- as.matrix(fst_mat)[site_order, site_order, drop = FALSE]
diag(fst_mat) <- 0
fst_mat[is.na(fst_mat)] <- 0
fst_mat[fst_mat < 0] <- 0
rownames(fst_mat) <- site_order
colnames(fst_mat) <- site_order

ETVA_1_fst <- fst_mat

# -----------------------------
# 8) coords dataframe in matching order
# -----------------------------
ETVA_1_coords <- paper_sites %>%
  filter(site %in% as.integer(site_order)) %>%
  arrange(match(site, as.integer(site_order))) %>%
  transmute(site = site, lat = lat, lon = lon)

stopifnot(identical(as.character(ETVA_1_coords$site), rownames(ETVA_1_fst)))

# -----------------------------
# 9) map plot
# Plot US and Canada and zoom to the extent of the points
# -----------------------------
world_df <- map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

site_lookup <- paper_sites %>%
  filter(site %in% ETVA_1_coords$site) %>%
  arrange(match(site, ETVA_1_coords$site))

x_pad <- max(0.5, diff(range(site_lookup$lon)) * 0.12)
y_pad <- max(0.5, diff(range(site_lookup$lat)) * 0.12)

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
    nudge_y = y_pad * 0.04,
    size = 2.7
  ) +
  coord_quickmap(
    xlim = c(min(site_lookup$lon) - x_pad, max(site_lookup$lon) + x_pad),
    ylim = c(min(site_lookup$lat) - y_pad, max(site_lookup$lat) + y_pad)
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()

# -----------------------------
# 10) IBD plot
# straight-line geographic distance as pilot distance layer
# -----------------------------
geo_dist <- geosphere::distm(
  x = as.matrix(ETVA_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist[upper.tri(geo_dist)],
  fst = ETVA_1_fst[upper.tri(ETVA_1_fst)]
)

ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 1.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Geographic distance (km)", y = "Pairwise FST") +
  theme_bw()

# -----------------------------
# 11) save objects
# -----------------------------
save(
  ETVA_1_fst,
  ETVA_1_coords,
  file = file.path(out_dir, "ETVA-1.RData")
)
