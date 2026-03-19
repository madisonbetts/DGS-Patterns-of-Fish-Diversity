# -----------------------------
# PEFL-1 yellow perch
# site coordinates + site-level FST matrix + IBD plot
# -----------------------------

library(ggplot2)
library(sf)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/PEFL-1/"

# -------------------------
# 0) site coordinates
# -------------------------
PEFL_1_coords <- data.frame(
  site = c("QUE","GTL","EUG","LSPNE","LSPSE","LSPNO","LSPSO",
           "IMOI","CTRN","CTRS","BOU","LSLN","LSLS","LDMT","LSFN","LSFS"),
  lat = c(
    46.813900, 46.386500, 46.282377, 46.258361,
    46.214667, 46.225250, 46.137444, 46.118972,
    45.886889, 45.870139, 45.595167, 45.412811,
    45.365284, 45.453167, 45.102306, 45.060361
  ),
  lon = c(
    -71.208000, -72.370000, -72.701890, -72.827667,
    -72.671278, -72.902528, -72.849250, -72.957639,
    -73.266139, -73.229806, -73.465222, -73.870654,
    -73.794181, -74.086694, -74.505611, -74.455833
  ),
  stringsAsFactors = FALSE
)

# optional map
ggplot(PEFL_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.05, size = 4) +
  coord_fixed() +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude")

# -------------------------
# 1) population labels
# -------------------------
pops <- c(
  "QUE (P03)", "QUE (P04)", "QUE (A03)", "QUE (A04)",
  "GTL (03)", "GTL (04)", "EUG (04)",
  "LSPNE (03)", "LSPNE (04)",
  "LSPSE (03)", "LSPSE (04)",
  "LSPNO (03)", "LSPNO (04)",
  "LSPSO (03)", "LSPSO (04)",
  "IMOI (03)",
  "CTRN (03)", "CTRN (04)",
  "CTRS (03)", "CTRS (04)",
  "BOU (03)",
  "LSLN (03)",
  "LSLS (03)", "LSLS (04)",
  "LDMT (03)", "LDMT (04)",
  "LSFN (03)",
  "LSFS (03)", "LSFS (04)"
)

# -------------------------
# 2) helper: fill symmetric matrix
# from row-wise upper triangle
# -------------------------
fill_sym_from_upper <- function(pops, vals, diag_val = 0) {
  n <- length(pops)
  stopifnot(length(vals) == n * (n - 1) / 2)
  
  mat <- matrix(0, nrow = n, ncol = n, dimnames = list(pops, pops))
  
  k <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      mat[i, j] <- vals[k]
      mat[j, i] <- vals[k]
      k <- k + 1
    }
  }
  
  diag(mat) <- diag_val
  mat
}

# -------------------------
# 3) population-level FST matrix
# values entered row-wise from the upper triangle
# -------------------------
fst_vals <- c(
  # QUE (P03)
  0.00873, 0.00179, 0.01927, 0.00205, 0.00346, 0.00374,
  0.00097, 0.00492, 0.00189, 0.00308,
  0.00134, 0.00425, 0.00571, 0.00255, 0.00553,
  0.02022, 0.01759, 0.02468, 0.01481,
  0.03127, 0.03557, 0.03297, 0.03096,
  0.01938, 0.02868, 0.07377, 0.07291, 0.06218,
  
  # QUE (P04)
  0.00920, 0.01045, 0.01176, 0.01113, 0.01814,
  0.01647, 0.00895, 0.01207, 0.01584,
  0.01290, 0.02078, 0.02925, 0.02395, 0.01884,
  0.04068, 0.04218, 0.04410, 0.03431,
  0.05559, 0.06317, 0.05774, 0.05277,
  0.03488, 0.04936, 0.10653, 0.11067, 0.09247,
  
  # QUE (A03)
  0.01870, 0.00063, 0.00804, 0.00616,
  0.00305, 0.00095, 0.00069, 0.00070,
  0.00260, 0.00498, 0.00597, 0.00344, 0.00531,
  0.02369, 0.01772, 0.02627, 0.01999,
  0.03067, 0.03543, 0.03956, 0.03984,
  0.01816, 0.02765, 0.08243, 0.08074, 0.06868,
  
  # QUE (A04)
  0.01471, 0.01826, 0.00114,
  0.02322, 0.02354, 0.01426, 0.01172,
  0.01952, 0.00978, 0.00735, 0.01196, 0.01595,
  0.03632, 0.00974, 0.03540, 0.02238,
  0.03766, 0.04437, 0.05014, 0.04672,
  0.02598, 0.02362, 0.09156, 0.09517, 0.08170,
  
  # GTL (03)
  0.00129, 0.00131,
  0.00400, 0.00622, 0.00199, -0.00005,
  0.00114, 0.00412, 0.00728, 0.00444, 0.00469,
  0.02219, 0.01212, 0.01687, 0.01612,
  0.02922, 0.03392, 0.03248, 0.03952,
  0.01782, 0.02225, 0.08218, 0.08322, 0.06924,
  
  # GTL (04)
  0.00457,
  0.00895, 0.01047, 0.00697, 0.00206,
  0.00207, 0.00531, 0.01418, 0.00798, 0.00383,
  0.01607, 0.00876, 0.01387, 0.00699,
  0.02136, 0.02834, 0.02553, 0.02519,
  0.01622, 0.01276, 0.07193, 0.07600, 0.06177,
  
  # EUG (04)
  0.00645, 0.01234, 0.00536, 0.00328,
  0.00421, 0.00391, 0.00827, 0.00681, 0.00697,
  0.02642, 0.00778, 0.01844, 0.00950,
  0.02312, 0.02575, 0.02717, 0.03465,
  0.01428, 0.02154, 0.07783, 0.06648, 0.05897,
  
  # LSPNE (03)
  0.00545, 0.00453, 0.00094,
  0.00034, 0.00415, 0.00306, 0.00643, 0.00516,
  0.02636, 0.02206, 0.02591, 0.01732,
  0.03384, 0.03626, 0.04020, 0.04024,
  0.02454, 0.02910, 0.08364, 0.08495, 0.07032,
  
  # LSPNE (04)
  0.00259, -0.00025,
  0.00552, 0.00363, 0.00673, 0.00444, 0.00761,
  0.02299, 0.02522, 0.02963, 0.01933,
  0.03558, 0.04176, 0.04251, 0.04483,
  0.02237, 0.03337, 0.09248, 0.09223, 0.07726,
  
  # LSPSE (03)
  0.00246, 0.00408, 0.00248, 0.00382, 0.00109,
  0.00378, 0.02683, 0.01816, 0.02691, 0.01912,
  0.03632, 0.04646, 0.04162, 0.04496,
  0.02807, 0.03120, 0.09669, 0.09464, 0.08130,
  
  # LSPSE (04)
  0.00055, -0.00305, 0.00133, 0.00063, 0.00445,
  0.01750, 0.01104, 0.01757, 0.01018,
  0.02215, 0.02842, 0.02901, 0.03371,
  0.01668, 0.01939, 0.07739, 0.07767, 0.06487,
  
  # LSPNO (03)
  0.00397, 0.00580, 0.00250, 0.00250,
  0.01657, 0.01329, 0.01917, 0.01162,
  0.02304, 0.02909, 0.02890, 0.03179,
  0.01885, 0.02198, 0.07002, 0.07357, 0.05834,
  
  # LSPNO (04)
  0.00149, -0.00117, 0.00541,
  0.02147, 0.01482, 0.02002, 0.01163,
  0.02542, 0.03382, 0.03065, 0.03551,
  0.01840, 0.02191, 0.08439, 0.08423, 0.07345,
  
  # LSPSO (03)
  0.00365, 0.00613,
  0.02453, 0.01698, 0.02719, 0.01901,
  0.02973, 0.03479, 0.03695, 0.04340,
  0.02505, 0.02769, 0.08141, 0.07979, 0.07048,
  
  # LSPSO (04)
  0.00847,
  0.02014, 0.01569, 0.02204, 0.01175,
  0.02641, 0.03383, 0.03181, 0.03472,
  0.02011, 0.02834, 0.07897, 0.07891, 0.06666,
  
  # IMOI (03)
  0.02060, 0.01534, 0.02416, 0.01487,
  0.02412, 0.02979, 0.03542, 0.03417,
  0.01865, 0.01269, 0.08116, 0.08139, 0.07111,
  
  # CTRN (03)
  0.00887, 0.00425, 0.00717,
  0.00148, 0.01250, 0.00543, 0.00977,
  0.01252, 0.00751, 0.03238, 0.03971, 0.02588,
  
  # CTRN (04)
  0.00571, 0.00210,
  0.00481, 0.01039, 0.01151, 0.01286,
  0.00867, 0.00419, 0.04662, 0.04459, 0.03725,
  
  # CTRS (03)
  0.00563,
  0.00248, 0.01065, 0.00253, 0.01597,
  0.01481, 0.00929, 0.03959, 0.04063, 0.02961,
  
  # CTRS (04)
  0.00665, 0.00844, 0.01152, 0.01214,
  0.00554, 0.01068, 0.05456, 0.04867, 0.03942,
  
  # BOU (03)
  0.00000, -0.00046, 0.00615,
  0.00876, 0.00404, 0.02509, 0.02306, 0.01793,
  
  # LSLN (03)
  0.01109, 0.01330,
  0.00592, 0.01053, 0.03213, 0.01939, 0.02007,
  
  # LSLS (03)
  0.01147,
  0.01974, 0.01804, 0.02592, 0.02379, 0.01970,
  
  # LSLS (04)
  0.01898, 0.01679, 0.03101, 0.03254, 0.02749,
  
  # LDMT (03)
  0.01216, 0.05175, 0.04720, 0.03814,
  
  # LDMT (04)
  0.04837, 0.05280, 0.04353,
  
  # LSFN (03)
  0.01137, 0.00506,
  
  # LSFS (03)
  0.00279
)

fst_mat <- fill_sym_from_upper(pops, fst_vals, diag_val = 0)
fst_mat[fst_mat < 0] <- 0

# -------------------------
# 4) project site coordinates
# once into meters
# -------------------------
coords_sf <- st_as_sf(PEFL_1_coords, coords = c("lon", "lat"), crs = 4326)
coords_sf <- st_transform(coords_sf, 32198)   # NAD83 / Quebec Lambert

coords_xy <- cbind(PEFL_1_coords["site"], st_coordinates(coords_sf))

# -------------------------
# 5) population-level distance matrix (km)
# -------------------------
pop_meta <- data.frame(
  pop  = pops,
  site = sub(" \\([^)]*\\)$", "", pops),
  stringsAsFactors = FALSE
)

pop_coords <- merge(pop_meta, coords_xy, by = "site", all.x = TRUE, sort = FALSE)
pop_coords <- pop_coords[match(pops, pop_coords$pop), ]

dist_mat <- as.matrix(dist(pop_coords[, c("X", "Y")])) / 1000
rownames(dist_mat) <- pop_coords$pop
colnames(dist_mat) <- pop_coords$pop

# -------------------------
# 6) collapse population FST
# to mean site-pair values
# -------------------------
site_ids <- sub(" \\([^)]*\\)$", "", rownames(fst_mat))
site_levels <- unique(site_ids)

PEFL_1_fst <- matrix(
  0,
  nrow = length(site_levels),
  ncol = length(site_levels),
  dimnames = list(site_levels, site_levels)
)

for (i in seq_along(site_levels)) {
  for (j in i:length(site_levels)) {
    
    idx_i <- which(site_ids == site_levels[i])
    idx_j <- which(site_ids == site_levels[j])
    
    if (i == j) {
      PEFL_1_fst[i, j] <- 0
    } else {
      PEFL_1_fst[i, j] <- mean(fst_mat[idx_i, idx_j], na.rm = TRUE)
      PEFL_1_fst[j, i] <- PEFL_1_fst[i, j]
    }
  }
}

PEFL_1_fst[PEFL_1_fst < 0] <- 0
diag(PEFL_1_fst) <- 0

# enforce site order and keep coordinates
site_order <- rownames(PEFL_1_fst)
PEFL_1_coords <- PEFL_1_coords[match(site_order, PEFL_1_coords$site), c("site", "lat", "lon")]

# -------------------------
# 7) site-level distance matrix (km)
# -------------------------
site_xy <- coords_xy[match(site_order, coords_xy$site), ]

site_dist <- as.matrix(dist(site_xy[, c("X", "Y")])) / 1000
rownames(site_dist) <- site_xy$site
colnames(site_dist) <- site_xy$site

# -------------------------
# 8) site-level IBD data
# -------------------------
ibd_df <- data.frame(
  distance_km = site_dist[upper.tri(site_dist)],
  mean_fst    = PEFL_1_fst[upper.tri(PEFL_1_fst)]
)

# -------------------------
# 9) final site-level IBD plot
# -------------------------
ggplot(ibd_df, aes(x = distance_km, y = mean_fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Distance (km)",
    y = "Mean FST",
    title = "PEFL-1 isolation by distance"
  )

# -------------------------
# 10) quick integrity checks
# -------------------------
stopifnot(
  nrow(PEFL_1_fst) == 16,
  ncol(PEFL_1_fst) == 16,
  identical(rownames(PEFL_1_fst), colnames(PEFL_1_fst)),
  identical(rownames(PEFL_1_fst), PEFL_1_coords$site),
  isTRUE(all.equal(PEFL_1_fst, t(PEFL_1_fst)))
)

# -------------------------
# 11) save objects
# -------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  PEFL_1_fst,
  PEFL_1_coords,
  file = file.path(out_dir, "PEFL-1.RData")
)