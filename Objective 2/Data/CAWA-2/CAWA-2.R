# ============================================================
# CAWA-2 | Warner sucker
# Catostomus warnerensis
# DeHaan et al. 2017, Ecology of Freshwater Fish
# DOI: 10.1111/eff.12305
# ============================================================
# Workflow notes
# - pairwise FST values are transcribed directly from Table 4
# - coordinates below are best-available approximate mid-creek / mid-canal
#   georeferences inferred from Figure 1
# - these are broad Objective 2 site estimates, not exact capture points
# - negative FST values are set to 0
# - plots are shown in RStudio and are not saved
# ============================================================

library(ggplot2)
library(dplyr)
library(geosphere)
library(maps)

study_code <- "CAWA-2"
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"
study_dir <- file.path(base_dir, study_code)
out_dir <- file.path(study_dir, "data")

if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

# -----------------------------
# 1) site key
# order follows Table 4
# -----------------------------
site_key <- data.frame(
  site = 1:5,
  site_code = c(
    "Honey Creek",
    "Snyder Creek",
    "Deep Creek",
    "Twentymile Creek",
    "Summer Lake Irrigation Canal"
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) site coordinates
# hard-coded from user-mapped CAWA-2 site csv
# -----------------------------
site_lookup <- data.frame(
  site = 1:5,
  site_code = site_key$site_code,
  lat = c(
    42.4277670,
    42.4779718,
    42.2190014,
    42.0502558,
    43.0120000
  ),
  lon = c(
    -120.0140863,
    -120.0733845,
    -120.1535843,
    -119.9676383,
    -120.7900000
  ),
  stringsAsFactors = FALSE
)

CAWA_2_coords <- site_lookup %>%
  select(site, lat, lon)

# -----------------------------
# 3) pairwise FST matrix
# values below diagonal in Table 4
# -----------------------------
CAWA_2_fst <- matrix(0, nrow = 5, ncol = 5)

CAWA_2_fst[2,1] <- 0.092
CAWA_2_fst[3,1] <- 0.082
CAWA_2_fst[3,2] <- 0.114
CAWA_2_fst[4,1] <- 0.180
CAWA_2_fst[4,2] <- 0.175
CAWA_2_fst[4,3] <- 0.095
CAWA_2_fst[5,1] <- 0.212
CAWA_2_fst[5,2] <- 0.211
CAWA_2_fst[5,3] <- 0.104
CAWA_2_fst[5,4] <- 0.139

CAWA_2_fst[upper.tri(CAWA_2_fst)] <- t(CAWA_2_fst)[upper.tri(CAWA_2_fst)]
diag(CAWA_2_fst) <- 0
CAWA_2_fst[CAWA_2_fst < 0] <- 0

rownames(CAWA_2_fst) <- as.character(site_key$site)
colnames(CAWA_2_fst) <- as.character(site_key$site)

stopifnot(isTRUE(all.equal(CAWA_2_fst, t(CAWA_2_fst))))

# -----------------------------
# 4) map plot
# -----------------------------
world_df <- map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

x_pad <- max(0.5, diff(range(site_lookup$lon)) * 0.12)
y_pad <- max(0.5, diff(range(site_lookup$lat)) * 0.12)

p_map <- ggplot() +
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
    nudge_y = y_pad * 0.03,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(min(site_lookup$lon) - x_pad, max(site_lookup$lon) + x_pad),
    ylim = c(min(site_lookup$lat) - y_pad, max(site_lookup$lat) + y_pad)
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "CAWA-2 sampling sites",
    subtitle = "User-mapped site coordinates"
  ) +
  theme_bw()

print(p_map)

# -----------------------------
# 5) IBD plot
# -----------------------------
geo_dist <- geosphere::distm(
  x = as.matrix(site_lookup[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist[upper.tri(geo_dist)],
  fst = CAWA_2_fst[upper.tri(CAWA_2_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Geographic distance (km)",
    y = "Pairwise FST",
    title = "CAWA-2 IBD plot"
  ) +
  theme_bw()

print(p_ibd)

# -----------------------------
# 6) save objects
# -----------------------------
save(
  CAWA_2_fst,
  CAWA_2_coords,
  file = file.path(out_dir, "CAWA-2.RData")
)
