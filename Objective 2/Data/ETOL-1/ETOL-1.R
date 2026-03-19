# -----------------------------
# ETOL-1 tessellated darter
# site coordinates + FST matrix
# -----------------------------

library(tidyverse)
library(geosphere)
library(ggplot2)

# directory where everything is
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETOL-1/"

# -----------------------------
# 0) site metadata
# remove Canadian Missisquoi sites up front
# -----------------------------
ETOL_1_coords <- data.frame(
  site = c(
    "ChIS","ChML1","ChMalE1","ChMalE2","ChMalW","ChML2",
    "IBADAF","IBADBF","IBBDBF",
    "LADAF1","LADAF2","LBDAF1","LBDAF2","LBDAF3","LBFBD1","LBFBD2"
  ),
  #site_name = c(
  #  "Lake Champlain – Inland Sea",
  #  "Lake Champlain – Main Lake 1 (Law Island)",
  #  "Lake Champlain – E. Malletts Bay 1",
  #  "Lake Champlain – E. Malletts Bay 2",
  #  "Lake Champlain – W. Malletts Bay",
  #  "Lake Champlain – Main Lake 2 (Sunset Island)",
  #  "Indian Brook Above Dam Above Fall Line",
  #  "Indian Brook Above Dam Below Fall Line",
  #  "Indian Brook Below Dam Below Fall Line",
  #  "Lewis Creek Above Dam Above Fall line 1",
  #  "Lewis Creek Above Dam Above Fall line 2",
  #  "Lewis Creek Below Dam Above Fall Line 1",
  #  "Lewis Creek Below Dam Above Fall Line 2",
  #  "Lewis Creek Below Dam Above Fall Line 3",
  #  "Lewis Below Dam Below Fall Line 1",
  #  "Lewis Below Dam Below Fall Line 2"
  #),
  lat = c(
    44.63813, 44.56186, 44.62641, 44.63031, 44.55770, 44.56235,
    44.53585, 44.53799, 44.54860,
    44.28339, 44.28750, 44.26557, 44.27860, 44.28031, 44.24549, 44.25964
  ),
  lon = c(
    -73.26140, -73.32148, -73.25384, -73.26620, -73.30088, -73.31015,
    -73.14530, -73.14684, -73.16579,
    -73.17473, -73.17302, -73.20415, -73.17858, -73.17396, -73.24727, -73.21330
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) read and clean pairwise table
# -----------------------------
remove_sites <- c("MissAD", "MissBD")

df <- read.csv(
  file.path(save_dir, "ETOL_fst_gst.csv"),
  stringsAsFactors = FALSE,
  check.names = TRUE
) %>%
  filter(!is.na(Comparison), Comparison != "", grepl(" vs\\. ", Comparison)) %>%
  separate(Comparison, into = c("site1", "site2"), sep = " vs\\. ") %>%
  mutate(across(c(site1, site2), trimws)) %>%
  filter(!(site1 %in% remove_sites | site2 %in% remove_sites))

# -----------------------------
# 2) build symmetric FST matrix
# -----------------------------
ETOL_1_fst <- matrix(
  0,
  nrow = nrow(ETOL_1_coords),
  ncol = nrow(ETOL_1_coords),
  dimnames = list(ETOL_1_coords$site, ETOL_1_coords$site)
)

i1 <- match(df$site1, ETOL_1_coords$site)
i2 <- match(df$site2, ETOL_1_coords$site)

bad <- which(is.na(i1) | is.na(i2))
if (length(bad) > 0) {
  print(df[bad, c("site1", "site2")])
  stop("Some site codes do not match ETOL_1_coords$site.")
}

for (i in seq_len(nrow(df))) {
  ETOL_1_fst[i1[i], i2[i]] <- df$FST[i]
  ETOL_1_fst[i2[i], i1[i]] <- df$FST[i]
}

diag(ETOL_1_fst) <- 0
ETOL_1_fst[ETOL_1_fst < 0] <- 0

# -----------------------------
# IBD plot: straight-line distance vs FST
# -----------------------------
dist_km <- geosphere::distm(
  as.matrix(ETOL_1_coords[, c("lon", "lat")]),
  fun = geosphere::distGeo
) / 1000

ibd_df <- data.frame(
  dist = dist_km[upper.tri(dist_km)],
  fst  = ETOL_1_fst[upper.tri(ETOL_1_fst)]
)

ggplot(ibd_df, aes(x = dist, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Geographic distance (km)",
    y = expression(F[ST]),
    title = "IBD: straight-line distance vs FST"
  )

# -----------------------------
# 3) quick checks
# -----------------------------
stopifnot(identical(rownames(ETOL_1_fst), ETOL_1_coords$site))
stopifnot(identical(colnames(ETOL_1_fst), ETOL_1_coords$site))
stopifnot(isTRUE(all.equal(ETOL_1_fst, t(ETOL_1_fst))))

# -----------------------------
# 4) save RData
# -----------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  ETOL_1_fst,
  ETOL_1_coords,
  file = file.path(out_dir, "ETOL-1.RData")
)