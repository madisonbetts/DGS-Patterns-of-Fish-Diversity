# ========================================
# SAFO-7
# Salvelinus fontinalis
# Brook Trout
#
# Beer et al. 2019
# Pairwise FST + GPS extraction workflow
#
# Reads:
#   NY_BKT_FST.xlsx
#   NY_BKT_GPS.xlsx
#
# Writes:
#   SAFO_7_fst
#   SAFO_7_coords
#
# to:
#   /Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/
#   Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/
#   DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SAFO-7/data/SAFO-7.RData
# ========================================

# -----------------------------
# 0) setup
# -----------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(geosphere)
library(maps)

study_code <- "SAFO-7"

base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SAFO-7"
data_dir <- file.path(base_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

fst_file <- file.path(base_dir, "NY_BKT_FST.xlsx")
gps_file <- file.path(base_dir, "NY_BKT_GPS.xlsx")

stopifnot(file.exists(fst_file))
stopifnot(file.exists(gps_file))

# -----------------------------
# 1) read and clean FST matrix
# file structure:
# row 1 = title
# row 2 = column IDs
# rows 3:80 = FST matrix (78 x 78)
# row 82 onward = p-values block that we do NOT want
# -----------------------------
fst_raw <- read_excel(
  fst_file,
  sheet = 1,
  col_names = FALSE
)

# row labels begin in row 3 until the first blank row
fst_row_ids_all <- fst_raw[[1]][3:nrow(fst_raw)]
first_blank <- which(is.na(fst_row_ids_all) | trimws(as.character(fst_row_ids_all)) == "")[1]

if (is.na(first_blank)) {
  stop("Could not detect the end of the FST block.")
}

n_sites_total <- first_blank - 1

fst_row_ids <- trimws(as.character(fst_raw[[1]][3:(2 + n_sites_total)]))
fst_col_ids <- trimws(as.character(unlist(fst_raw[2, 2:(1 + n_sites_total)])))

fst_vals <- as.data.frame(fst_raw[3:(2 + n_sites_total), 2:(1 + n_sites_total)])
fst_vals[] <- lapply(fst_vals, function(x) as.numeric(as.character(x)))

fst_mat <- as.matrix(fst_vals)
rownames(fst_mat) <- fst_row_ids
colnames(fst_mat) <- fst_col_ids

stopifnot(nrow(fst_mat) == ncol(fst_mat))
stopifnot(identical(rownames(fst_mat), colnames(fst_mat)))

# lower triangle is filled in the file; mirror it to upper
fst_mat[upper.tri(fst_mat)] <- t(fst_mat)[upper.tri(fst_mat)]

fst_mat[fst_mat < 0] <- 0
diag(fst_mat) <- 0

# -----------------------------
# 2) read and clean GPS table
# expected columns in uploaded file:
# Latitude, Longitude, Lab ID, Field ID Number, N, Basin,
# NYS_HUC_12, Location, Survey_Date
# -----------------------------
gps_raw <- read_excel(
  gps_file,
  sheet = 1
)

names(gps_raw) <- trimws(names(gps_raw))

lat_col <- names(gps_raw)[grepl("^latitude$|lat", tolower(names(gps_raw)))]
lon_col <- names(gps_raw)[grepl("^longitude$|^long$|lon", tolower(names(gps_raw)))]
lab_col <- names(gps_raw)[tolower(names(gps_raw)) %in% c("lab id", "lab_id", "labid")]

stopifnot(length(lat_col) >= 1)
stopifnot(length(lon_col) >= 1)
stopifnot(length(lab_col) >= 1)

lat_col <- lat_col[1]
lon_col <- lon_col[1]
lab_col <- lab_col[1]

gps <- gps_raw %>%
  transmute(
    site_code = trimws(as.character(.data[[lab_col]])),
    lat = as.numeric(.data[[lat_col]]),
    lon = as.numeric(.data[[lon_col]])
  ) %>%
  filter(!is.na(site_code), site_code != "") %>%
  distinct(site_code, .keep_all = TRUE)

# -----------------------------
# 3) keep only wild sites that have GPS
# FST file includes 78 IDs:
# 75 wild sites + ROM, OSW, TYL hatchery strains
# GPS file includes the 75 wild sites only
# -----------------------------
wild_ids <- rownames(fst_mat)[rownames(fst_mat) %in% gps$site_code]

if (length(wild_ids) == 0) {
  stop("No overlap between FST row names and GPS Lab ID values.")
}

fst_mat <- fst_mat[wild_ids, wild_ids, drop = FALSE]

gps <- gps %>%
  filter(site_code %in% wild_ids)

gps <- gps[match(wild_ids, gps$site_code), , drop = FALSE]

stopifnot(length(wild_ids) == nrow(gps))
stopifnot(all(wild_ids == gps$site_code))

# -----------------------------
# 4) build final objects
# fst gets numeric row/col names 1:n
# coords gets site / lat / lon
# -----------------------------
SAFO_7_coords <- data.frame(
  site = seq_len(nrow(gps)),
  lat = gps$lat,
  lon = gps$lon,
  stringsAsFactors = FALSE
)

SAFO_7_fst <- fst_mat
rownames(SAFO_7_fst) <- SAFO_7_coords$site
colnames(SAFO_7_fst) <- SAFO_7_coords$site

stopifnot(nrow(SAFO_7_fst) == ncol(SAFO_7_fst))
stopifnot(nrow(SAFO_7_fst) == nrow(SAFO_7_coords))

# -----------------------------
# 5) geographic distance matrix
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(SAFO_7_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- SAFO_7_coords$site
colnames(geo_dist_km) <- SAFO_7_coords$site

# -----------------------------
# 6) pairwise dataframe for IBD
# -----------------------------
ibd_df <- data.frame(
  site1 = rownames(SAFO_7_fst)[row(SAFO_7_fst)[upper.tri(SAFO_7_fst)]],
  site2 = colnames(SAFO_7_fst)[col(SAFO_7_fst)[upper.tri(SAFO_7_fst)]],
  fst = SAFO_7_fst[upper.tri(SAFO_7_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  stringsAsFactors = FALSE
)

# -----------------------------
# 7) sampling map
# tight extent around western NY points
# -----------------------------
state_map <- map_data("state")

lon_rng <- range(SAFO_7_coords$lon, na.rm = TRUE)
lat_rng <- range(SAFO_7_coords$lat, na.rm = TRUE)

x_pad <- max(0.35, diff(lon_rng) * 0.08)
y_pad <- max(0.25, diff(lat_rng) * 0.08)

xlim_use <- c(lon_rng[1] - x_pad, lon_rng[2] + x_pad)
ylim_use <- c(lat_rng[1] - y_pad, lat_rng[2] + y_pad)

# plot sites
ggplot() +
  geom_polygon(
    data = state_map,
    aes(x = long, y = lat, group = group),
    fill = "grey96",
    color = "grey55",
    linewidth = 0.25
  ) +
  geom_point(
    data = SAFO_7_coords,
    aes(x = lon, y = lat),
    size = 2.5
  ) +
  coord_fixed(
    xlim = xlim_use,
    ylim = ylim_use
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "SAFO-7 sampling locations"
  )


# -----------------------------
# 8) IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.2, alpha = 0.75) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "SAFO-7 isolation by distance"
  )


# -----------------------------
# 9) checks
# -----------------------------
stopifnot(identical(rownames(SAFO_7_fst), as.character(SAFO_7_coords$site)))
stopifnot(identical(colnames(SAFO_7_fst), as.character(SAFO_7_coords$site)))
stopifnot(isTRUE(all.equal(SAFO_7_fst, t(SAFO_7_fst))))
stopifnot(isTRUE(all.equal(geo_dist_km, t(geo_dist_km))))

# -----------------------------
# 10) save RData
# -----------------------------
save(
  SAFO_7_fst,
  SAFO_7_coords,
  file = file.path(data_dir, "SAFO-7.RData")
)
