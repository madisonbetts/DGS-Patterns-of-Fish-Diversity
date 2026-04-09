
# ============================================================
# NOTA-1 | Caddo Madtom
# Noturus taylori
# McCall & Fluker 2020, Conservation Genetics
#
# Objective 2 workflow
# - transcribes microsatellite pairwise FST from Table 3
# - builds a 5 x 5 symmetric FST matrix using the microsatellite
#   FST values (lower triangle in Table 3), not mtDNA PhiST
# - reads user-updated grouped coordinates from a csv when present
# - links grouped coordinates to the published FST populations using
#   the site_name text rather than assuming the csv row order is the
#   same as the FST table order
# - plots shown in RStudio only; not saved
# - saves NOTA_1_fst and NOTA_1_coords to data/NOTA-1.RData
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(geosphere)
  library(readr)
  library(stringr)
})

study_code <- "NOTA-1"
study_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/NOTA-1"
out_dir <- file.path(study_dir, "data")

if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

# -----------------------------
# 1) published FST population order
# Table 1 / Table 3:
# 1 = Caddo River (E, G, H)
# 2 = Lick Creek (F)
# 3 = South Fork Caddo River (I, J)
# 4 = Ouachita River (A, B, C)
# 5 = South Fork Ouachita River (D)
# -----------------------------
fst_groups <- data.frame(
  site = 1:5,
  site_name = c(
    "Caddo River",
    "Lick Creek",
    "South Fork Caddo River",
    "Ouachita River",
    "South Fork Ouachita River"
  ),
  site_abbr = c("CAD", "LIC", "SFC", "OUA", "SFO"),
  map_group = c("(E/G/H)", "(F)", "(I/J)", "(A/B/C)", "(D)"),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) grouped coordinates
# Prefer the user-updated csv when present. Match rows to the FST
# groups using the text in site_name so the workflow is robust even
# when the csv row order or numeric ids differ from the FST order.
# -----------------------------
coord_candidates <- c(
  "NOTA-1_coords.csv",
  "NOTA-1- NOTA-1_coords.csv.csv"
)

coord_file <- coord_candidates[file.exists(coord_candidates)][1]

if (!is.na(coord_file) && nzchar(coord_file)) {
  raw_coords <- suppressMessages(readr::read_csv(coord_file, show_col_types = FALSE))
} else {
  # fallback coordinates from the previous draft workflow
  raw_coords <- data.frame(
    site_name = c(
      "Ouachita River - Oden Access (A/B/C)",
      "South Fork Ouachita River (D)",
      "Caddo River (E/G/H)",
      "Lick Creek (F)",
      "South Fork Caddo River (I/J)"
    ),
    lat = c(34.6065899, 34.5084005, 34.4608793, 34.4732828, 34.3630233),
    lon = c(-93.8211198, -93.8188255, -93.6946497, -93.8099403, -93.7659144),
    stringsAsFactors = FALSE
  )
}

required_cols <- c("site_name", "lat", "lon")
stopifnot(all(required_cols %in% names(raw_coords)))

coord_lookup <- raw_coords %>%
  transmute(
    source_name = site_name,
    lat = as.numeric(lat),
    lon = as.numeric(lon),
    fst_name = case_when(
      str_detect(str_to_lower(site_name), "south fork ouachita") ~ "South Fork Ouachita River",
      str_detect(str_to_lower(site_name), "ouachita river") ~ "Ouachita River",
      str_detect(str_to_lower(site_name), "south fork caddo") ~ "South Fork Caddo River",
      str_detect(str_to_lower(site_name), "lick creek") ~ "Lick Creek",
      str_detect(str_to_lower(site_name), "caddo") ~ "Caddo River",
      TRUE ~ NA_character_
    )
  )

if (any(is.na(coord_lookup$fst_name))) {
  stop("At least one coordinate row could not be matched to an FST group.")
}

if (anyDuplicated(coord_lookup$fst_name)) {
  dupes <- coord_lookup$fst_name[duplicated(coord_lookup$fst_name)]
  stop(paste("Duplicate coordinate matches found for:", paste(unique(dupes), collapse = ", ")))
}

NOTA_1_coords <- fst_groups %>%
  left_join(coord_lookup, by = c("site_name" = "fst_name")) %>%
  select(site, site_name, site_abbr, map_group, source_name, lat, lon)

stopifnot(!any(is.na(NOTA_1_coords$lat)))
stopifnot(!any(is.na(NOTA_1_coords$lon)))
stopifnot(identical(NOTA_1_coords$site, 1:5))

print(NOTA_1_coords)

# -----------------------------
# 3) pairwise FST matrix
# transcribed from Table 3, microsatellite FST values
# (lower triangle in the published table)
# order:
# 1 Caddo River
# 2 Lick Creek
# 3 South Fork Caddo River
# 4 Ouachita River
# 5 South Fork Ouachita River
# -----------------------------
site_ids <- as.character(NOTA_1_coords$site)

NOTA_1_fst <- matrix(0, nrow = 5, ncol = 5)
rownames(NOTA_1_fst) <- site_ids
colnames(NOTA_1_fst) <- site_ids

lower_vals <- list(
  c(2, 1, 0.066),
  c(3, 1, 0.068),
  c(3, 2, 0.141),
  c(4, 1, 0.227),
  c(4, 2, 0.275),
  c(4, 3, 0.231),
  c(5, 1, 0.305),
  c(5, 2, 0.357),
  c(5, 3, 0.312),
  c(5, 4, 0.069)
)

for (x in lower_vals) {
  i <- x[1]
  j <- x[2]
  v <- x[3]
  NOTA_1_fst[i, j] <- v
  NOTA_1_fst[j, i] <- v
}

diag(NOTA_1_fst) <- 0
NOTA_1_fst[NOTA_1_fst < 0] <- 0

# -----------------------------
# 4) map plot
# map_data is from ggplot2, not maps
# include USA and Canada, then zoom to point extent
# -----------------------------
world_df <- ggplot2::map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

plot_df <- NOTA_1_coords %>%
  mutate(label = paste0(site, " ", site_abbr))

x_pad <- max(0.25, diff(range(plot_df$lon)) * 0.15)
y_pad <- max(0.25, diff(range(plot_df$lat)) * 0.15)

p_map <- ggplot() +
  geom_polygon(
    data = world_df,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = plot_df,
    aes(x = lon, y = lat),
    size = 2
  ) +
  geom_text(
    data = plot_df,
    aes(x = lon, y = lat, label = label),
    nudge_y = 0.02,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(min(plot_df$lon) - x_pad, max(plot_df$lon) + x_pad),
    ylim = c(min(plot_df$lat) - y_pad, max(plot_df$lat) + y_pad)
  ) +
  labs(
    title = "NOTA-1 sampling sites",
    subtitle = "Grouped localities matched to the published FST populations",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw()

print(p_map)

# -----------------------------
# 5) IBD plot
# straight-line distance for QC only
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(NOTA_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  fst = NOTA_1_fst[upper.tri(NOTA_1_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "NOTA-1 IBD plot",
    x = "Geographic distance (km)",
    y = "Pairwise FST"
  ) +
  theme_bw()

print(p_ibd)

# -----------------------------
# 6) save outputs
# save the Objective 2 objects only
# -----------------------------
NOTA_1_coords <- NOTA_1_coords %>%
  select(site, lat, lon)

save(
  NOTA_1_fst,
  NOTA_1_coords,
  file = file.path(out_dir, "NOTA-1.RData")
)
