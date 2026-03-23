# -------------------------------
# Global IBD plot across all datasets
# -------------------------------
library(terra)
library(ggplot2)

base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"

# -------------------------------
# find all RData files in /data subfolders
# -------------------------------
rdata_files <- list.files(
  base_dir,
  pattern = "\\.RData$",
  recursive = TRUE,
  full.names = TRUE
)

# keep only those in "data" folders
rdata_files <- rdata_files[grepl("/data/", rdata_files)]

# -------------------------------
# loop + accumulate IBD points
# -------------------------------
all_ibd <- data.frame(dist_km = numeric(), fst = numeric())

for (f in rdata_files) {
  
  e <- new.env()
  load(f, envir = e)
  
  # grab fst + coords objects dynamically
  fst_obj    <- ls(e)[grepl("_fst$", ls(e))]
  coords_obj <- ls(e)[grepl("_coords$", ls(e))]
  
  if (length(fst_obj) == 0 | length(coords_obj) == 0) next
  
  fst    <- e[[fst_obj[1]]]
  coords <- e[[coords_obj[1]]]
  
  # ensure ordering matches
  coords <- coords[match(rownames(fst), coords$site), ]
  
  # skip if mismatch
  if (any(is.na(coords$lat))) next
  
  # distance matrix (km)
  d <- terra::distance(
    vect(coords[, c("lon", "lat")], crs = "EPSG:4326")
  ) / 1000
  
  d <- as.matrix(d)
  
  # extract upper triangle
  all_ibd <- rbind(
    all_ibd,
    data.frame(
      dist_km = d[upper.tri(d)],
      fst     = fst[upper.tri(fst)]
    )
  )
}

# -------------------------------
# clean + bounds
# -------------------------------
all_ibd <- all_ibd[!is.na(all_ibd$fst) & !is.na(all_ibd$dist_km), ]

xmax <- max(all_ibd$dist_km, na.rm = TRUE)
ymax <- max(all_ibd$fst, na.rm = TRUE)

# -------------------------------
# plot
# -------------------------------
ggplot(all_ibd, aes(x = dist_km, y = fst)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax)) +
  labs(
    x = "Geographic distance (km)",
    y = expression(F[ST]),
    title = "Global Isolation by Distance (All Studies)"
  )