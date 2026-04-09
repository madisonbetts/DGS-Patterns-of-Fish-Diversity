# ============================================================
# ESMA-2 | Muskellunge
# Esox masquinongy
# White et al. 2018, Transactions of the American Fisheries Society
#
# Objective 2 workflow
# - pairwise FST values transcribed from Appendix Table A.2
# - site coordinates hardcoded from user-corrected Google plotting
# - plots are shown in RStudio only; not saved
# - saves ESMA_2_fst, ESMA_2_coords, and ESMA_2_coords_check
#   to data/ESMA-2.RData
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(geosphere)
})

study_code <- "ESMA-2"
study_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ESMA-2"
out_dir <- file.path(study_dir, "data")

if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

# -----------------------------
# 1) site metadata
# order follows Table 1 / Appendix Table A.2:
# BER, LMR, SUN, BR, ELK, MIC, NKR, NYH, KH, LH, MCH
#
# coordinates hardcoded from the user-corrected ESMA-2 csv
# -----------------------------
ESMA_2_coords_check <- data.frame(
  site_id = 1:11,
  site = c("BER", "LMR", "SUN", "BR", "ELK", "MIC", "NKR", "NYH", "KH", "LH", "MCH"),
  site_name = c(
    "Berlin Lake, Portage County, Ohio",
    "Little Muskingum River, Washington County, Ohio",
    "Sunfish Creek, Pike County, Ohio",
    "Buckhannon River, Upshur County, West Virginia",
    "Elk River, Kanawha County, West Virginia",
    "Middle Island Creek, Pleasants County, West Virginia",
    "New and Kanawha rivers, Kanawha County, West Virginia",
    "Chautauqua State Fish Hatchery, New York",
    "Kincaid State Fish Hatchery, Ohio",
    "London State Fish Hatchery, Ohio",
    "Minor Clark State Fish Hatchery, Kentucky"
  ),
  lat = c(
    41.0120000,
    39.6100925,
    39.0155930,
    38.9720000,
    38.4997557,
    39.4397877,
    38.1595557,
    42.1590000,
    39.1007291,
    39.8967994,
    38.1136375
  ),
  lon = c(
    -81.0160000,
    -81.1183578,
    -83.0945000,
    -80.2300000,
    -81.2088848,
    -81.0119047,
    -81.2016310,
    -79.4300000,
    -83.2624998,
    -83.5105147,
    -83.5267855
  ),
  stringsAsFactors = FALSE
)

print(ESMA_2_coords_check)

# final Objective 2 coords df should only have site, lat, lon
ESMA_2_coords <- ESMA_2_coords_check %>%
  transmute(site = site_id, lat = lat, lon = lon)

# -----------------------------
# 2) pairwise FST matrix
# values transcribed from Appendix Table A.2
# order:
# BER, LMR, SUN, BR, ELK, MIC, NKR, NYH, KH, LH, MCH
# -----------------------------
site_ids <- as.character(ESMA_2_coords$site)

ESMA_2_fst <- matrix(0, nrow = 11, ncol = 11)
rownames(ESMA_2_fst) <- site_ids
colnames(ESMA_2_fst) <- site_ids

# fill lower triangle from published table, then mirror
lower_vals <- list(
  c(2,1,0.055),
  c(3,1,0.039), c(3,2,0.065),
  c(4,1,0.224), c(4,2,0.208), c(4,3,0.132),
  c(5,1,0.119), c(5,2,0.114), c(5,3,0.057), c(5,4,0.100),
  c(6,1,0.146), c(6,2,0.124), c(6,3,0.113), c(6,4,0.113), c(6,5,0.078),
  c(7,1,0.167), c(7,2,0.153), c(7,3,0.082), c(7,4,0.064), c(7,5,0.059), c(7,6,0.094),
  c(8,1,0.113), c(8,2,0.091), c(8,3,0.098), c(8,4,0.198), c(8,5,0.102), c(8,6,0.107), c(8,7,0.136),
  c(9,1,0.074), c(9,2,0.093), c(9,3,0.167), c(9,4,0.248), c(9,5,0.157), c(9,6,0.174), c(9,7,0.179), c(9,8,0.138),
  c(10,1,0.073), c(10,2,0.083), c(10,3,0.050), c(10,4,0.230), c(10,5,0.146), c(10,6,0.149), c(10,7,0.171), c(10,8,0.115), c(10,9,0.028),
  c(11,1,0.126), c(11,2,0.123), c(11,3,0.082), c(11,4,0.244), c(11,5,0.139), c(11,6,0.150), c(11,7,0.174), c(11,8,0.111), c(11,9,0.099), c(11,10,0.084)
)

for (x in lower_vals) {
  i <- x[1]; j <- x[2]; v <- x[3]
  ESMA_2_fst[i, j] <- v
  ESMA_2_fst[j, i] <- v
}

diag(ESMA_2_fst) <- 0
ESMA_2_fst[is.na(ESMA_2_fst)] <- 0
ESMA_2_fst[ESMA_2_fst < 0] <- 0

stopifnot(identical(as.character(ESMA_2_coords$site), rownames(ESMA_2_fst)))

# -----------------------------
# 3) map plot
# include USA and Canada, then zoom to point extent
# labels use site_id + site code for QC
# -----------------------------
world_df <- ggplot2::map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

plot_df <- ESMA_2_coords_check %>%
  mutate(label = paste0(site_id, " ", site))

x_pad <- max(1.2, diff(range(plot_df$lon)) * 0.08)
y_pad <- max(0.8, diff(range(plot_df$lat)) * 0.10)

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
    nudge_y = 0.15,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(min(plot_df$lon) - x_pad, max(plot_df$lon) + x_pad),
    ylim = c(min(plot_df$lat) - y_pad, max(plot_df$lat) + y_pad)
  ) +
  labs(
    title = "ESMA-2 sampling sites",
    subtitle = "Hardcoded user-corrected coordinates",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw()

print(p_map)

# -----------------------------
# 4) IBD plot
# straight-line distance for QC only
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(ESMA_2_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  fst = ESMA_2_fst[upper.tri(ESMA_2_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "ESMA-2 IBD plot",
    x = "Geographic distance (km)",
    y = "Pairwise FST"
  ) +
  theme_bw()

print(p_ibd)

# -----------------------------
# 5) save outputs
# -----------------------------
save(
  ESMA_2_fst,
  ESMA_2_coords,
  ESMA_2_coords_check,
  file = file.path(out_dir, "ESMA-2.RData")
)

write.csv(
  ESMA_2_coords_check,
  file = file.path(out_dir, "ESMA-2_coords_check.csv"),
  row.names = FALSE
)
