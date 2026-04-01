# ========================================
# PYOL-1
# Pylodictis olivaris
# Flathead catfish
# Four collections from the Schuylkill and Susquehanna rivers
# Build collection centroids from provided genotype dataset,
# calculate pairwise FST from microsatellite genotypes,
# map sites, make IBD plot, and save PYOL-1.RData
# ========================================

# -----------------------------
# 0) setup
# -----------------------------
library(dplyr)
library(ggplot2)
library(geosphere)
library(adegenet)
library(hierfstat)
library(maps)

study_code <- "PYOL-1"

base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/PYOL-1"
data_dir <- file.path(base_dir, "data")
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) locate and read genotype data
# -----------------------------
possible_geno_files <- c(
  file.path(base_dir, "DataRelease.csv"),
  file.path(getwd(), "DataRelease.csv")
)

geno_file <- possible_geno_files[file.exists(possible_geno_files)][1]
stopifnot(!is.na(geno_file))

dat <- read.csv(
  geno_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(nrow(dat) > 0)

# CSV has duplicate column names for diploid microsatellite alleles
# make them unique before any dplyr operations
colnames(dat) <- make.unique(colnames(dat), sep = "_")

# -----------------------------
# 2) site order and labels
# article + metadata indicate four collections:
# Schuylkill, Lower_Susquehanna, Middle_Susquehanna, Upper_Susquehanna
# -----------------------------
site_order <- c(
  "Schuylkill",
  "Lower_Susquehanna",
  "Middle_Susquehanna",
  "Upper_Susquehanna"
)

site_labels <- c(
  "Schuylkill",
  "Lower Susquehanna",
  "Middle Susquehanna",
  "Upper Susquehanna"
)

stopifnot(all(site_order %in% unique(dat$`Sample Location`)))

# -----------------------------
# 3) clean and order data
# -----------------------------
dat <- dat %>%
  mutate(
    sample_id = `Sample Info`,
    sample_location = factor(`Sample Location`, levels = site_order),
    Lat = as.numeric(Lat),
    Long = as.numeric(Long)
  ) %>%
  arrange(sample_location, sample_id)

stopifnot(!any(is.na(dat$sample_location)))

# -----------------------------
# 4) derive one centroid per collection
# using mean coordinates of all fish in each collection
# -----------------------------
site_centroids <- dat %>%
  group_by(sample_location) %>%
  summarise(
    lat = mean(Lat, na.rm = TRUE),
    lon = mean(Long, na.rm = TRUE),
    n_ind = n(),
    .groups = "drop"
  ) %>%
  arrange(match(sample_location, site_order))

PYOL_1_coords <- data.frame(
  site_id = as.character(seq_len(nrow(site_centroids))),
  lat = site_centroids$lat,
  lon = site_centroids$lon,
  stringsAsFactors = FALSE
)

site_lookup <- data.frame(
  site_id = PYOL_1_coords$site_id,
  site_name = site_labels,
  sample_location = site_order,
  n_ind = site_centroids$n_ind,
  stringsAsFactors = FALSE
)

stopifnot(nrow(PYOL_1_coords) == 4)

# -----------------------------
# 5) build diploid microsatellite table for df2genind
# each locus is stored in two columns after make.unique(),
# e.g. GY113J02 and GY113J02_1
# 000 denotes missing and should be NA
# -----------------------------
locus_names <- c(
  "GY113J02",
  "GY047K03",
  "IpCG0071",
  "71-75",
  "Ip077",
  "IpCG00189",
  "Ip271",
  "Ip357",
  "Ip365",
  "Ip372",
  "Ip554",
  "Ip591",
  "IpCG0001"
)

expected_cols <- c(rbind(locus_names, paste0(locus_names, "_1")))
stopifnot(all(expected_cols %in% colnames(dat)))

geno_list <- lapply(locus_names, function(loc) {
  a1 <- dat[[loc]]
  a2 <- dat[[paste0(loc, "_1")]]

  a1 <- ifelse(a1 %in% c("000", 0, "0", ""), NA, as.character(a1))
  a2 <- ifelse(a2 %in% c("000", 0, "0", ""), NA, as.character(a2))

  ifelse(is.na(a1) | is.na(a2), NA, paste(a1, a2, sep = "/"))
})

geno_df <- as.data.frame(geno_list, stringsAsFactors = FALSE)
colnames(geno_df) <- paste0("L", seq_along(locus_names))
rownames(geno_df) <- dat$sample_id

# -----------------------------
# 6) convert to genind and calculate pairwise FST
# Weir & Cockerham pairwise FST via hierfstat
# -----------------------------
pyol_gen <- adegenet::df2genind(
  X = geno_df,
  sep = "/",
  ploidy = 2,
  ind.names = dat$sample_id,
  pop = dat$sample_location,
  NA.char = NA,
  type = "codom"
)

pyol_hf <- hierfstat::genind2hierfstat(pyol_gen)
pyol_fst_named <- hierfstat::pairwise.WCfst(pyol_hf)

# reorder to match desired site order
pyol_fst_named <- pyol_fst_named[site_order, site_order]

# force any negative FST to zero
pyol_fst_named[is.na(pyol_fst_named)] <- 0
pyol_fst_named[pyol_fst_named < 0] <- 0
diag(pyol_fst_named) <- 0

PYOL_1_fst <- as.matrix(pyol_fst_named)
rownames(PYOL_1_fst) <- PYOL_1_coords$site_id
colnames(PYOL_1_fst) <- PYOL_1_coords$site_id

# -----------------------------
# 7) geographic distance matrix
# straight-line distance among collection centroids
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(PYOL_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- PYOL_1_coords$site_id
colnames(geo_dist_km) <- PYOL_1_coords$site_id

# -----------------------------
# 8) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1 = rownames(PYOL_1_fst)[row(PYOL_1_fst)[upper.tri(PYOL_1_fst)]],
  site2 = colnames(PYOL_1_fst)[col(PYOL_1_fst)[upper.tri(PYOL_1_fst)]],
  fst = PYOL_1_fst[upper.tri(PYOL_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  stringsAsFactors = FALSE
)

ibd_df$site1_name <- site_lookup$site_name[match(ibd_df$site1, site_lookup$site_id)]
ibd_df$site2_name <- site_lookup$site_name[match(ibd_df$site2, site_lookup$site_id)]

# -----------------------------
# 9) map of sampling locations
# plot to the extent of the points and include states
# -----------------------------
state_map <- map_data("state")
canada_map <- subset(map_data("world"), region == "Canada")

plot_sites <- cbind(PYOL_1_coords, site_lookup[c("site_name", "n_ind")])

lon_rng <- range(plot_sites$lon, na.rm = TRUE)
lat_rng <- range(plot_sites$lat, na.rm = TRUE)

x_pad <- max(0.75, diff(lon_rng) * 0.8)
y_pad <- max(0.60, diff(lat_rng) * 1.2)

xlim_use <- c(lon_rng[1] - x_pad, lon_rng[2] + x_pad)
ylim_use <- c(lat_rng[1] - y_pad, lat_rng[2] + y_pad)

ggplot() +
  geom_polygon(
    data = canada_map,
    aes(x = long, y = lat, group = group),
    fill = "grey97",
    color = "grey65",
    linewidth = 0.25
  ) +
  geom_polygon(
    data = state_map,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey55",
    linewidth = 0.25
  ) +
  geom_point(
    data = plot_sites,
    aes(x = lon, y = lat),
    size = 3
  ) +
  geom_text(
    data = plot_sites,
    aes(x = lon, y = lat, label = paste0(site_id, ". ", site_name)),
    nudge_y = 0.08,
    size = 3.5
  ) +
  coord_fixed(
    xlim = xlim_use,
    ylim = ylim_use
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "PYOL-1 sampling locations"
  )


# -----------------------------
# 10) IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "PYOL-1 isolation by distance"
  )

# -----------------------------
# 11) quick checks
# -----------------------------
stopifnot(identical(rownames(PYOL_1_fst), PYOL_1_coords$site_id))
stopifnot(identical(colnames(PYOL_1_fst), PYOL_1_coords$site_id))
stopifnot(isTRUE(all.equal(PYOL_1_fst, t(PYOL_1_fst))))
stopifnot(isTRUE(all.equal(geo_dist_km, t(geo_dist_km))))

# -----------------------------
# 12) save RData
# -----------------------------
save(
  PYOL_1_fst,
  PYOL_1_coords,
  file = file.path(data_dir, "PYOL-1.RData")
)
