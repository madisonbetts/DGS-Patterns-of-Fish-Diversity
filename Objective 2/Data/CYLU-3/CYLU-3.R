# -----------------------------
# CYLU-3 | Red shiner
# Cyprinella lutrensis
# Franssen 2012, Freshwater Biology 57:155-165
# Lake Texoma reservoir-altered stream network
# -----------------------------

# Workflow notes:
# - Site coordinates are taken directly from Table 1.
# - Pairwise FST values are taken directly from Table 3.
# - Negative FST values are set to 0, per workflow.
# - The two Brier Creek Cove collections (2008 and 2009) are treated as
#   separate samples because they are reported separately in Table 3.
# - Coordinates dataframe contains only: site, lat, lon.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(geosphere)
  library(maps)
})

# -----------------------------
# 1) site coordinates
# Table 1 from Franssen (2012)
# site numbers here match the FST matrix order below
# -----------------------------

site_lookup <- tibble(
  site = 1:13,
  sample_name = c(
    "Whiskey Creek",
    "Red River",
    "Coffee Pot Creek",
    "Walnut Bayou Creek",
    "Bills Creek",
    "Hickory Creek",
    "Brier Creek Cove 2008",
    "Brier Creek Cove 2009",
    "Brier Creek Station 6",
    "Brier Creek Station 5",
    "Biostation Shore",
    "Glasses Creek Cove",
    "Little Glasses Creek"
  ),
  lat = c(
    34.1529,
    33.9377,
    33.9403,
    33.9166,
    33.8970,
    34.0377,
    33.9248,
    33.9248,
    33.9539,
    33.9990,
    33.8794,
    34.0281,
    34.0256
  ),
  lon = c(
    -98.1565,
    -97.7320,
    -97.4518,
    -97.2824,
    -97.1577,
    -97.1434,
    -96.8654,
    -96.8654,
    -96.8422,
    -96.8286,
    -96.8021,
    -96.6611,
    -96.6983
  )
)

CYLU_3_coords <- site_lookup %>%
  select(site, lat, lon)

# -----------------------------
# 2) pairwise FST matrix
# Table 3 from Franssen (2012)
# lower triangle entered row-by-row in site order above
# -----------------------------

lower_vals <- list(
  c(),
  c(0.006),
  c(0.012, 0.001),
  c(0.011, 0.001, 0.003),
  c(0.009, 0.005, -0.008, 0.000),
  c(0.038, 0.020, 0.006, 0.009, 0.005),
  c(0.025, 0.010, 0.009, -0.001, 0.002, -0.003),
  c(0.015, 0.017, 0.003, 0.005, 0.001, -0.001, 0.001),
  c(0.008, 0.000, -0.016, -0.002, -0.013, -0.006, 0.002, 0.004),
  c(0.055, 0.040, 0.034, 0.034, 0.032, 0.030, 0.039, 0.034, 0.026),
  c(0.006, -0.006, -0.008, -0.001, 0.002, 0.005, 0.009, 0.006, -0.016, 0.023),
  c(0.078, 0.042, 0.016, 0.027, 0.026, 0.008, 0.020, 0.029, 0.019, 0.048, 0.035),
  c(0.067, 0.037, 0.016, 0.031, 0.021, 0.011, 0.027, 0.028, 0.002, 0.033, 0.020, -0.001)
)

n_sites <- nrow(CYLU_3_coords)
fst <- matrix(0, n_sites, n_sites)
rownames(fst) <- CYLU_3_coords$site
colnames(fst) <- CYLU_3_coords$site

for (i in seq_len(n_sites)) {
  if (length(lower_vals[[i]]) > 0) {
    fst[i, seq_len(i - 1)] <- lower_vals[[i]]
  }
}

fst <- fst + t(fst)
diag(fst) <- 0
fst[fst < 0] <- 0

CYLU_3_fst <- fst

# -----------------------------
# 3) quick checks
# -----------------------------

stopifnot(
  nrow(CYLU_3_fst) == nrow(CYLU_3_coords),
  all(rownames(CYLU_3_fst) == CYLU_3_coords$site),
  all(colnames(CYLU_3_fst) == CYLU_3_coords$site),
  all(CYLU_3_fst == t(CYLU_3_fst)),
  all(diag(CYLU_3_fst) == 0),
  !any(is.na(CYLU_3_fst))
)

# -----------------------------
# 4) map plot (US + Canada context; zoom to points)
# not saved, just plotted in RStudio
# -----------------------------

world_map <- map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

x_pad <- diff(range(CYLU_3_coords$lon)) * 0.18
y_pad <- diff(range(CYLU_3_coords$lat)) * 0.18
if (x_pad == 0) x_pad <- 0.5
if (y_pad == 0) y_pad <- 0.5

p_map <- ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "grey95", color = "grey70", linewidth = 0.2
  ) +
  geom_point(
    data = CYLU_3_coords,
    aes(x = lon, y = lat),
    size = 2.3
  ) +
  geom_text(
    data = CYLU_3_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.03,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(min(CYLU_3_coords$lon) - x_pad, max(CYLU_3_coords$lon) + x_pad),
    ylim = c(min(CYLU_3_coords$lat) - y_pad, max(CYLU_3_coords$lat) + y_pad)
  ) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude", title = "CYLU-3 sampling sites")

print(p_map)

# -----------------------------
# 5) IBD plot
# using geographic distance among site coordinates
# not saved, just plotted in RStudio
# -----------------------------

geo_dist_km <- geosphere::distm(
  x = CYLU_3_coords[, c("lon", "lat")],
  fun = geosphere::distHaversine
) / 1000

ibd_df <- tibble(
  site1 = row(CYLU_3_fst)[lower.tri(CYLU_3_fst)],
  site2 = col(CYLU_3_fst)[lower.tri(CYLU_3_fst)],
  fst = CYLU_3_fst[lower.tri(CYLU_3_fst)],
  geo_km = geo_dist_km[lower.tri(geo_dist_km)]
)

p_ibd <- ggplot(ibd_df, aes(x = geo_km, y = fst)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  theme_bw() +
  labs(
    x = "Geographic distance (km)",
    y = expression(italic(F)[ST]),
    title = "CYLU-3 IBD"
  )

print(p_ibd)

# -----------------------------
# 6) save objects
# -----------------------------

out_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/CYLU-3/data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save(
  CYLU_3_fst,
  CYLU_3_coords,
  file = file.path(out_dir, "CYLU-3.RData")
)
