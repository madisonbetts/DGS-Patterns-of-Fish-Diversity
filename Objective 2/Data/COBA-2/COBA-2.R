# -----------------------------
# COBA-2 mottled sculpin
# Cottus bairdii
# Homola et al. 2016
# 16 sites, 6 microsatellite loci
# pairwise FST from Table III (below diagonal)
# site XYs are best-inference representative stream points
# from named streams + Fig. 1
# -----------------------------

library(ggplot2)
library(ggrepel)
library(dplyr)
library(maps)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/COBA-2"
data_dir <- file.path(save_dir, "data")

# -----------------------------
# 1) coordinates
# final saved object only has site, lat, lon
# site order follows Table II / Table III exactly:
# SCG AC LC PC RC SCK 7MC SLV SB WC BC BHC CeC CC KC SCW
# -----------------------------
COBA_2_coords <- data.frame(
  site = as.character(1:16),
  lat = c(
    42.946136,  # 1  SCG  Sand Creek (Grand basin; representative stream point)
    42.333096,  # 2  AC   Augusta Creek
    42.310000,  # 3  LC   Lee Creek (best-inference from Fig. 1 + watershed position)
    42.295318,  # 4  PC   Portage Creek
    42.292819,  # 5  RC   Rice Creek
    42.592252,  # 6  SCK  Sand Creek (Kalamazoo basin representative point)
    42.360600,  # 7  7MC  Sevenmile Creek
    42.414481,  # 8  SLV  Silver Creek
    42.356428,  # 9  SB   Spring Brook
    42.227250,  # 10 WC   Wilder Creek
    43.420855,  # 11 BC   Bigelow Creek
    43.791687,  # 12 BHC  Buckhorn Creek
    43.772510,  # 13 CeC  Cedar Creek
    41.891100,  # 14 CC   Curtis Creek
    43.530570,  # 15 KC   Knutson Creek
    43.456680   # 16 SCW  Sand Creek (White basin)
  ),
  lon = c(
    -85.852260, # 1
    -85.348333, # 2
    -85.650000, # 3
    -85.573062, # 4
    -84.906925, # 5
    -85.977532, # 6
    -85.299170, # 7
    -85.594734, # 8
    -85.576953, # 9
    -84.907806, # 10
    -85.771994, # 11
    -85.500325, # 12
    -85.983950, # 13
    -85.760800, # 14
    -86.191180, # 15
    -86.262570  # 16
  ),
  stringsAsFactors = FALSE
)

# plotting helper with stream codes
COBA_2_coords_plot <- data.frame(
  site = as.character(1:16),
  code = c("SCG","AC","LC","PC","RC","SCK","7MC","SLV","SB","WC","BC","BHC","CeC","CC","KC","SCW"),
  lat = COBA_2_coords$lat,
  lon = COBA_2_coords$lon,
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) helper: fill symmetric matrix from lower triangle
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
# Table III values entered from LOWER triangle row-wise
# order:
# SCG AC LC PC RC SCK 7MC SLV SB WC BC BHC CeC CC KC SCW
# -----------------------------
fst_vals <- c(
  # AC
  0.222,
  
  # LC
  0.322, 0.334,
  
  # PC
  0.122, 0.239, 0.342,
  
  # RC
  0.248, 0.152, 0.224, 0.314,
  
  # SCK
  0.159, 0.255, 0.342, 0.127, 0.310,
  
  # 7MC
  0.209, 0.104, 0.243, 0.221, 0.094, 0.243,
  
  # SLV
  0.268, 0.265, 0.039, 0.291, 0.112, 0.324, 0.173,
  
  # SB
  0.252, 0.253, 0.129, 0.270, 0.192, 0.311, 0.170, 0.082,
  
  # WC
  0.135, 0.148, 0.161, 0.138, 0.108, 0.136, 0.093, 0.090, 0.115,
  
  # BC
  0.108, 0.100, 0.168, 0.156, 0.138, 0.156, 0.127, 0.134, 0.137, 0.081,
  
  # BHC
  0.314, 0.112, 0.548, 0.383, 0.311, 0.364, 0.207, 0.481, 0.432, 0.291, 0.263,
  
  # CeC
  0.183, 0.180, 0.240, 0.215, 0.253, 0.294, 0.212, 0.201, 0.156, 0.166, 0.110, 0.406,
  
  # CC
  0.278, 0.422, 0.373, 0.303, 0.276, 0.242, 0.276, 0.306, 0.320, 0.172, 0.319, 0.533, 0.402,
  
  # KC
  0.132, 0.167, 0.239, 0.154, 0.249, 0.197, 0.215, 0.235, 0.206, 0.164, 0.089, 0.355, 0.055, 0.358,
  
  # SCW
  0.266, 0.256, 0.278, 0.316, 0.229, 0.323, 0.267, 0.242, 0.260, 0.222, 0.185, 0.459, 0.198, 0.385, 0.144
)

COBA_2_fst <- fill_sym_from_lower(
  pops = COBA_2_coords$site,
  vals = fst_vals,
  diag_val = 0
)

COBA_2_fst[COBA_2_fst < 0] <- 0
diag(COBA_2_fst) <- 0

# -----------------------------
# 4) inspect outputs
# -----------------------------
print(COBA_2_coords)
print(round(COBA_2_fst, 3))

cat("\nSite key:\n")
cat("1  = SCG\n")
cat("2  = AC\n")
cat("3  = LC\n")
cat("4  = PC\n")
cat("5  = RC\n")
cat("6  = SCK\n")
cat("7  = 7MC\n")
cat("8  = SLV\n")
cat("9  = SB\n")
cat("10 = WC\n")
cat("11 = BC\n")
cat("12 = BHC\n")
cat("13 = CeC\n")
cat("14 = CC\n")
cat("15 = KC\n")
cat("16 = SCW\n")

# -----------------------------
# 5) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(COBA_2_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- COBA_2_coords$site
colnames(geo_dist_km) <- COBA_2_coords$site

# -----------------------------
# 6) map of sampling locations
# US + Canada
# -----------------------------
world_map <- map_data("world") |>
  dplyr::filter(region %in% c("USA", "Canada"))

states_map <- map_data("state")

x_pad <- 2.5
y_pad <- 1.8

xlim_use <- range(COBA_2_coords_plot$lon) + c(-x_pad, x_pad)
ylim_use <- range(COBA_2_coords_plot$lat) + c(-y_pad, y_pad)

p_map <- ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "gray98",
    color = "gray80",
    linewidth = 0.2
  ) +
  geom_polygon(
    data = states_map,
    aes(x = long, y = lat, group = group),
    fill = "gray95",
    color = "gray70",
    linewidth = 0.2
  ) +
  geom_point(
    data = COBA_2_coords_plot,
    aes(x = lon, y = lat),
    size = 2.8
  ) +
  geom_text_repel(
    data = COBA_2_coords_plot,
    aes(x = lon, y = lat, label = paste0(site, ": ", code)),
    size = 3.0,
    max.overlaps = 100
  ) +
  coord_quickmap(
    xlim = xlim_use,
    ylim = ylim_use,
    expand = FALSE
  ) +
  theme_classic() +
  labs(
    title = "COBA-2 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )

print(p_map)

# -----------------------------
# 7) IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(COBA_2_fst)[row(COBA_2_fst)[upper.tri(COBA_2_fst)]],
  site2   = colnames(COBA_2_fst)[col(COBA_2_fst)[upper.tri(COBA_2_fst)]],
  fst     = COBA_2_fst[upper.tri(COBA_2_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    title = "COBA-2 isolation by distance",
    x = "Euclidean distance among site centroids (km)",
    y = expression(F[ST])
  )

print(p_ibd)

# -----------------------------
# 8) save
# -----------------------------
save(
  COBA_2_fst,
  COBA_2_coords,
  file = file.path(data_dir, "COBA-2.RData")
)