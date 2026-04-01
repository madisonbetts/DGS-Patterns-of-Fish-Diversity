# ========================================
# LEME-1
# Lepomis megalotis
# longear sunfish
#
# Smith 2016 thesis:
# "Low-Head Dams on Habitat, Fish Functional Guilds and Genetic Structuring in a Midwestern River System"
#
# Exact site coordinates from Appendix Table I
# Pairwise LOS FST values from Table 2.5
# ========================================

library(ggplot2)
library(geosphere)
library(maps)

study_code <- "LEME-1"

base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/LEME-1"
data_dir <- file.path(base_dir, "data")

dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) exact site coordinates
# Appendix Table I
# Site order matches Table 2.5 exactly
# -----------------------------
LEME_1_coords <- data.frame(
  site_id = as.character(1:12),
  lat = c(
    40.117908, 40.121214, 40.122522, 40.117217, 40.120531, 40.122811,
    40.120542, 40.122756, 40.124717, 40.128669, 40.132097, 40.129844
  ),
  lon = c(
    -87.629033, -87.631983, -87.633308, -87.650892, -87.653156, -87.660347,
    -87.641531, -87.639333, -87.639153, -87.642339, -87.647511, -87.650286
  ),
  stringsAsFactors = FALSE
)

site_lookup <- data.frame(
  site_id = as.character(1:12),
  site_code = c(
    "V_BD1", "V_BD2", "V_P1", "V_P2", "V_R1", "V_R2",
    "NF_BD1", "NF_BD2", "NF_P1", "NF_P2", "NF_R1", "NF_R2"
  ),
  river = c(
    rep("Vermilion", 6),
    rep("North Fork Vermilion", 6)
  ),
  location = c(
    "below dam", "below dam", "pool", "pool", "river", "river",
    "below dam", "below dam", "pool", "pool", "river", "river"
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) helper
# fill symmetric matrix from lower triangle values
# -----------------------------
fill_sym_from_lower <- function(pops, vals, diag_val = 0) {
  n <- length(pops)
  stopifnot(length(vals) == n * (n - 1) / 2)

  mat <- matrix(
    0,
    nrow = n,
    ncol = n,
    dimnames = list(pops, pops)
  )

  k <- 1
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      mat[i, j] <- vals[k]
      mat[j, i] <- vals[k]
      k <- k + 1
    }
  }

  diag(mat) <- diag_val
  mat
}

# -----------------------------
# 3) pairwise FST matrix
# LOS values from Table 2.5
# negatives set to 0 per workflow
# order:
# V_BD1, V_BD2, V_P1, V_P2, V_R1, V_R2,
# NF_BD1, NF_BD2, NF_P1, NF_P2, NF_R1, NF_R2
# -----------------------------
fst_vals <- c(
  -0.0032,
  -0.0015,  0.0057,
  -0.0040, -0.0028,  0.0083,
  -0.0064, -0.0020,  0.0046, -0.0057,
  -0.0021, -0.0023,  0.0043, -0.0008, -0.0049,
  -0.0015,  0.0101,  0.0058,  0.0012,  0.0013,  0.0025,
  -0.0025,  0.0001, -0.0011,  0.0023, -0.0020, -0.0034,  0.0009,
   0.0007,  0.0076,  0.0154,  0.0028,  0.0043,  0.0028,  0.0090,  0.0029,
  -0.0026, -0.0010,  0.0076, -0.0023, -0.0006, -0.0001,  0.0079,  0.0022, -0.0003,
  -0.0028,  0.0078, -0.0029,  0.0041, -0.0007, -0.0017, -0.0039,  0.0012,  0.0084,  0.0053,
  -0.0052,  0.0066,  0.0106,  0.0025,  0.0015,  0.0038,  0.0018,  0.0036, -0.0011, -0.0031,  0.0052
)

LEME_1_fst <- fill_sym_from_lower(
  pops = LEME_1_coords$site_id,
  vals = fst_vals,
  diag_val = 0
)

LEME_1_fst[LEME_1_fst < 0] <- 0
diag(LEME_1_fst) <- 0

# -----------------------------
# 4) geographic distance matrix
# straight-line distance among exact site coordinates
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(LEME_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- LEME_1_coords$site_id
colnames(geo_dist_km) <- LEME_1_coords$site_id

# -----------------------------
# 5) IBD dataframe
# -----------------------------
ibd_df <- data.frame(
  site1 = rownames(LEME_1_fst)[row(LEME_1_fst)[upper.tri(LEME_1_fst)]],
  site2 = colnames(LEME_1_fst)[col(LEME_1_fst)[upper.tri(LEME_1_fst)]],
  fst = LEME_1_fst[upper.tri(LEME_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  stringsAsFactors = FALSE
)

ibd_df$site1_code <- site_lookup$site_code[match(ibd_df$site1, site_lookup$site_id)]
ibd_df$site2_code <- site_lookup$site_code[match(ibd_df$site2, site_lookup$site_id)]

# -----------------------------
# 6) map
# tight local extent around exact sites
# -----------------------------
state_map <- map_data("state")
plot_sites <- merge(LEME_1_coords, site_lookup, by = "site_id", sort = FALSE)

lon_rng <- range(plot_sites$lon, na.rm = TRUE)
lat_rng <- range(plot_sites$lat, na.rm = TRUE)

x_pad <- max(0.02, diff(lon_rng) * 0.20)
y_pad <- max(0.01, diff(lat_rng) * 0.25)

xlim_use <- c(lon_rng[1] - x_pad, lon_rng[2] + x_pad)
ylim_use <- c(lat_rng[1] - y_pad, lat_rng[2] + y_pad)

ggplot() +
  geom_polygon(
    data = state_map,
    aes(x = long, y = lat, group = group),
    fill = "grey96",
    color = "grey55",
    linewidth = 0.25
  ) +
  geom_point(
    data = plot_sites,
    aes(x = lon, y = lat),
    size = 2.8
  ) +
  geom_text(
    data = plot_sites,
    aes(x = lon, y = lat, label = paste0(site_id, ". ", site_code)),
    nudge_y = 0.0018,
    size = 3.1
  ) +
  coord_fixed(
    xlim = xlim_use,
    ylim = ylim_use
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "LEME-1 sampling locations"
  )


# -----------------------------
# 7) IBD plot
# -----------------------------
 ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.8, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "LEME-1 isolation by distance"
  )



# -----------------------------
# 8) checks + save
# -----------------------------
stopifnot(identical(rownames(LEME_1_fst), LEME_1_coords$site_id))
stopifnot(identical(colnames(LEME_1_fst), LEME_1_coords$site_id))
stopifnot(isTRUE(all.equal(LEME_1_fst, t(LEME_1_fst))))
stopifnot(isTRUE(all.equal(geo_dist_km, t(geo_dist_km))))

save(
  LEME_1_fst,
  LEME_1_coords,
  file = file.path(data_dir, "LEME-1.RData")
)
