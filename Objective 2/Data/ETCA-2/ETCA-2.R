# -----------------------------
# Rainbow darter eastern Iowa
# approximate site coordinates + FST workflow
# -----------------------------
library(ggplot2)
library(ggrepel)
library(geosphere)

# -----------------------------
# 0) site metadata
# -----------------------------
etca_2_sites <- data.frame(
  site = c("UI4M","UI3BG","UI1BC","UI2HR",
           "TK1WC","TK3FA","TK2SC","VTVS",
           "MRNL","MRNM","MRBP",
           "OCNA","RCSO","EOC"),
  lat = c(
    43.602, 43.531, 43.555, 43.528,
    43.454, 43.281, 43.404, 43.094,
    42.888, 42.776, 42.891,
    43.529, 43.381, 42.355
  ),
  lon = c(
    -91.983, -91.957, -91.664, -92.204,
    -92.080, -91.785, -91.973, -91.880,
    -91.755, -91.669, -91.657,
    -92.919, -92.812, -91.824
  ),
  stringsAsFactors = FALSE
)

etca_2_sites

# sanity plot
ggplot(etca_2_sites, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = site),
    size = 4,
    box.padding = 0.3,
    point.padding = 0.3
  ) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "Estimated sampling locations for rainbow darter"
  )

# -----------------------------
# 1) pairwise FST matrix
# values transcribed from lower triangle in figure
# -----------------------------
sites <- etca_2_sites$site

fst <- matrix(
  0,
  nrow = length(sites),
  ncol = length(sites),
  dimnames = list(sites, sites)
)

fst["UI3BG", "UI4M"] <- 0.003

fst["UI1BC", c("UI4M","UI3BG")] <- c(0.006, 0.009)

fst["UI2HR", c("UI4M","UI3BG","UI1BC")] <- c(0.023, 0.011, 0.012)

fst["TK1WC", c("UI4M","UI3BG","UI1BC","UI2HR")] <- c(0.022, 0.015, 0.016, 0.031)

fst["TK3FA", c("UI4M","UI3BG","UI1BC","UI2HR","TK1WC")] <- c(0.011, 0.008, 0.014, 0.026, 0.014)

fst["TK2SC", c("UI4M","UI3BG","UI1BC","UI2HR","TK1WC","TK3FA")] <- c(0.005, 0.005, 0.010, 0.025, 0.013, 0.002)

fst["VTVS", c("UI4M","UI3BG","UI1BC","UI2HR","TK1WC","TK3FA","TK2SC")] <- c(0.009, 0.005, 0.008, 0.011, 0.019, 0.012, 0.006)

fst["MRNL", c("UI4M","UI3BG","UI1BC","UI2HR","TK1WC","TK3FA","TK2SC","VTVS")] <- c(0.020, 0.020, 0.017, 0.032, 0.017, 0.014, 0.009, 0.020)

fst["MRNM", c("UI4M","UI3BG","UI1BC","UI2HR","TK1WC","TK3FA","TK2SC","VTVS","MRNL")] <- c(0.015, 0.020, 0.014, 0.024, 0.027, 0.011, 0.013, 0.014, -0.002)

fst["MRBP", c("UI4M","UI3BG","UI1BC","UI2HR","TK1WC","TK3FA","TK2SC","VTVS","MRNL","MRNM")] <- c(0.020, 0.026, 0.019, 0.007, 0.022, 0.024, 0.016, 0.010, 0.005, -0.0007)

fst["OCNA", c("UI4M","UI3BG","UI1BC","UI2HR","TK1WC","TK3FA","TK2SC","VTVS","MRNL","MRNM","MRBP")] <- c(0.010, 0.009, 0.006, 0.024, 0.022, 0.004, 0.012, 0.012, 0.023, 0.0111, 0.024)

fst["RCSO", c("UI4M","UI3BG","UI1BC","UI2HR","TK1WC","TK3FA","TK2SC","VTVS","MRNL","MRNM","MRBP","OCNA")] <- c(0.019, 0.008, 0.006, 0.014, 0.021, 0.006, 0.016, 0.016, 0.022, 0.008, 0.019, -0.004)

fst["EOC", c("UI4M","UI3BG","UI1BC","UI2HR","TK1WC","TK3FA","TK2SC","VTVS","MRNL","MRNM","MRBP","OCNA","RCSO")] <- c(0.010, 0.013, 0.017, 0.017, 0.031, 0.016, 0.015, 0.016, 0.033, 0.0148, 0.018, 0.0057, 0.0006)

fst[upper.tri(fst)] <- t(fst)[upper.tri(fst)]

# set negative FST values to 0
fst[fst < 0] <- 0

fst

# -----------------------------
# 2) reorder site metadata
#    to match fst matrix
# -----------------------------
etca_2_sites_fst <- etca_2_sites[match(rownames(fst), etca_2_sites$site), ]

stopifnot(identical(rownames(fst), etca_2_sites_fst$site))

etca_2_sites_fst

# -----------------------------
# 3) optional:
#    geographic distance matrix in km
# -----------------------------
coords <- etca_2_sites_fst[, c("lon", "lat")]

site_dist_km <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
site_dist_km <- as.matrix(site_dist_km)

rownames(site_dist_km) <- etca_2_sites_fst$site
colnames(site_dist_km) <- etca_2_sites_fst$site

round(site_dist_km, 2)

# -----------------------------
# 4) optional IBD plot
# -----------------------------
ibd_df <- data.frame(
  dist = site_dist_km[upper.tri(site_dist_km)],
  fst  = fst[upper.tri(fst)]
)

plot(
  ibd_df$dist,
  ibd_df$fst,
  xlab = "Geographic distance (km)",
  ylab = "Pairwise FST",
  pch = 21
)

# optional ggplot version
ggplot(ibd_df, aes(x = dist, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Geographic distance (km)",
    y = expression(F[ST]),
    title = "Rainbow darter: geographic distance vs. FST"
  )

# -----------------------------
# 5) prepare outputs
# -----------------------------
ETCA_2_fst <- fst

ETCA_2_coords <- etca_2_sites_fst[, c("site", "lat", "lon")]

# -----------------------------
# 6) save RData
# -----------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETCA-2/data"

if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

save(
  ETCA_2_fst,
  ETCA_2_coords,
  file = file.path(save_dir, "ETCA-2.RData")
)

# -----------------------------
# 7) quick checks
# -----------------------------
stopifnot(identical(rownames(ETCA_2_fst), ETCA_2_coords$site))
stopifnot(identical(colnames(ETCA_2_fst), ETCA_2_coords$site))
stopifnot(isTRUE(all.equal(ETCA_2_fst, t(ETCA_2_fst))))