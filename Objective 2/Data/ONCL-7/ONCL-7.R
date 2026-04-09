# ============================================================
# ONCL-7 | Coastal cutthroat trout
# Oncorhynchus clarki clarki
# Wofford et al. 2005. Ecological Applications 15:628-637
#
# Corrected Objective 2 extraction workflow
# - uses the user-updated coordinates from ONCL-7_coords.csv
# - collapses T4-0 and T4-1 to one site ("T4") by averaging their
#   pairwise FST values to all other populations
# - uses the mean XY of T4-0 and T4-1 for the collapsed T4 site
# - keeps MS4 as MS4-1 because the paper removed age-0 fish from MS4
# - plots shown in RStudio only; not saved
# - saves ONCL_7_fst and ONCL_7_coords to data/ONCL-7.RData
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(geosphere)
})

study_code <- "ONCL-7"
study_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ONCL-7"
out_dir <- file.path(study_dir, "data")

if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(study_dir)

# -----------------------------
# 1) original 11 x 11 pairwise FST matrix from Table 1
# order in paper:
# MS0, MS1, T1, T2, T3, MS2, T4-0, T4-1, UT4, MS3, MS4-1
# -----------------------------
orig_names <- c("MS0", "MS1", "T1", "T2", "T3", "MS2", "T4-0", "T4-1", "UT4", "MS3", "MS4-1")

fst11 <- matrix(0, nrow = 11, ncol = 11, dimnames = list(orig_names, orig_names))

upper_vals <- list(
  c(1,2,0.04), c(1,3,0.19), c(1,4,0.06), c(1,5,0.03), c(1,6,0.07),
  c(1,7,0.10), c(1,8,0.10), c(1,9,0.19), c(1,10,0.06), c(1,11,0.09),
  c(2,3,0.22), c(2,4,0.03), c(2,5,0.01), c(2,6,0.04), c(2,7,0.08),
  c(2,8,0.08), c(2,9,0.19), c(2,10,0.05), c(2,11,0.09),
  c(3,4,0.25), c(3,5,0.24), c(3,6,0.26), c(3,7,0.37), c(3,8,0.33),
  c(3,9,0.39), c(3,10,0.23), c(3,11,0.32),
  c(4,5,0.04), c(4,6,0.08), c(4,7,0.11), c(4,8,0.10), c(4,9,0.21),
  c(4,10,0.04), c(4,11,0.10),
  c(5,6,0.05), c(5,7,0.09), c(5,8,0.09), c(5,9,0.21), c(5,10,0.04),
  c(5,11,0.06),
  c(6,7,0.05), c(6,8,0.03), c(6,9,0.15), c(6,10,0.03), c(6,11,0.06),
  c(7,8,0.12), c(7,9,0.20), c(7,10,0.09), c(7,11,0.15),
  c(8,9,0.20), c(8,10,0.07), c(8,11,0.06),
  c(9,10,0.12), c(9,11,0.17),
  c(10,11,0.06)
)

for (x in upper_vals) {
  i <- x[1]; j <- x[2]; v <- x[3]
  fst11[i, j] <- v
  fst11[j, i] <- v
}
diag(fst11) <- 0

# -----------------------------
# 2) collapse T4-0 and T4-1 to one T4 site
#
# Paper handling:
# - T4 was one sampling location that was split into age-0 and age-1+
#   groups because those age classes differed at multiple loci.
#
# Meta-analysis handling here:
# - collapse to a single mapped site T4
# - average FST(T4-0, x) and FST(T4-1, x) for every other site x
# -----------------------------
final_names <- c("MS0", "MS1", "T1", "T2", "T3", "MS2", "T4", "UT4", "MS3", "MS4-1")
ONCL_7_fst <- matrix(0, nrow = 10, ncol = 10, dimnames = list(as.character(1:10), as.character(1:10)))

carry <- fst11[c("MS0", "MS1", "T1", "T2", "T3", "MS2", "UT4", "MS3", "MS4-1"),
               c("MS0", "MS1", "T1", "T2", "T3", "MS2", "UT4", "MS3", "MS4-1")]

ONCL_7_fst[c(1,2,3,4,5,6,8,9,10), c(1,2,3,4,5,6,8,9,10)] <- carry

t4_means <- c(
  MS0   = mean(c(fst11["T4-0","MS0"],   fst11["T4-1","MS0"])),
  MS1   = mean(c(fst11["T4-0","MS1"],   fst11["T4-1","MS1"])),
  T1    = mean(c(fst11["T4-0","T1"],    fst11["T4-1","T1"])),
  T2    = mean(c(fst11["T4-0","T2"],    fst11["T4-1","T2"])),
  T3    = mean(c(fst11["T4-0","T3"],    fst11["T4-1","T3"])),
  MS2   = mean(c(fst11["T4-0","MS2"],   fst11["T4-1","MS2"])),
  UT4   = mean(c(fst11["T4-0","UT4"],   fst11["T4-1","UT4"])),
  MS3   = mean(c(fst11["T4-0","MS3"],   fst11["T4-1","MS3"])),
  MS4_1 = mean(c(fst11["T4-0","MS4-1"], fst11["T4-1","MS4-1"]))
)

put_sym <- function(mat, i, j, value) {
  mat[i, j] <- value
  mat[j, i] <- value
  mat
}

for (pair in list(
  c(7,1,t4_means["MS0"]),
  c(7,2,t4_means["MS1"]),
  c(7,3,t4_means["T1"]),
  c(7,4,t4_means["T2"]),
  c(7,5,t4_means["T3"]),
  c(7,6,t4_means["MS2"]),
  c(7,8,t4_means["UT4"]),
  c(7,9,t4_means["MS3"]),
  c(7,10,t4_means["MS4_1"])
)) {
  ONCL_7_fst <- put_sym(ONCL_7_fst, pair[1], pair[2], pair[3])
}

diag(ONCL_7_fst) <- 0
ONCL_7_fst[ONCL_7_fst < 0] <- 0

# -----------------------------
# 3) final site coordinates
# these are the updated user-provided coordinates
# T4 = mean of T4-0 and T4-1 coordinates
# -----------------------------
ONCL_7_coords <- data.frame(
  site = 1:10,
  site_abbr = final_names,
  lat = c(
    43.5552992, # MS0
    43.5531056, # MS1
    43.5583207, # T1
    43.5475059, # T2
    43.5531422, # T3
    43.5435478, # MS2
    mean(c(43.5402545, 43.5400670)), # T4
    43.5327188, # UT4
    43.5390342, # MS3
    43.5362672  # MS4-1
  ),
  lon = c(
    -123.7137362, # MS0
    -123.7070670, # MS1
    -123.6956159, # T1
    -123.7057134, # T2
    -123.6822118, # T3
    -123.6810677, # MS2
    mean(c(-123.6824978, -123.6824127)), # T4
    -123.6863834, # UT4
    -123.6748643, # MS3
    -123.6656124  # MS4-1
  ),
  stringsAsFactors = FALSE
)

ONCL_7_lookup <- ONCL_7_coords %>%
  mutate(label = paste0(site, " ", site_abbr))

# -----------------------------
# 4) map plot
# -----------------------------
world_df <- ggplot2::map_data("world") %>%
  filter(region %in% c("USA", "Canada"))

x_pad <- max(0.05, diff(range(ONCL_7_lookup$lon)) * 0.12)
y_pad <- max(0.05, diff(range(ONCL_7_lookup$lat)) * 0.12)

p_map <- ggplot() +
  geom_polygon(
    data = world_df,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = ONCL_7_lookup,
    aes(x = lon, y = lat),
    size = 2
  ) +
  geom_text(
    data = ONCL_7_lookup,
    aes(x = lon, y = lat, label = label),
    nudge_y = 0.0010,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(min(ONCL_7_lookup$lon) - x_pad, max(ONCL_7_lookup$lon) + x_pad),
    ylim = c(min(ONCL_7_lookup$lat) - y_pad, max(ONCL_7_lookup$lat) + y_pad)
  ) +
  labs(
    title = "ONCL-7 sampling sites",
    subtitle = "Updated user coordinates; T4 collapsed to one mapped site",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw()

print(p_map)

# -----------------------------
# 5) IBD plot
# -----------------------------
geo_dist_km <- geosphere::distm(
  x = as.matrix(ONCL_7_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

ibd_df <- data.frame(
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  fst = ONCL_7_fst[upper.tri(ONCL_7_fst)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "ONCL-7 IBD plot",
    x = "Geographic distance (km)",
    y = "Pairwise FST"
  ) +
  theme_bw()

print(p_ibd)

# -----------------------------
# 6) save outputs
# -----------------------------
save(
  ONCL_7_fst,
  ONCL_7_coords,
  file = file.path(out_dir, "ONCL-7.RData")
)

write.csv(ONCL_7_coords, file.path(out_dir, "ONCL-7_coords_final.csv"), row.names = FALSE)
