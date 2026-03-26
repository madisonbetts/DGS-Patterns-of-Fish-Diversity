# -----------------------------
# MATE-1
# -----------------------------

library(ggplot2)
library(geosphere)

save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/MATE-1"

# coords
MATE_1_coords <- data.frame(
  site = as.character(1:6),
  lat = c(35.84, 35.78, 35.71, 35.65, 35.60, 35.54),
  lon = c(-85.50, -85.57, -85.63, -85.69, -85.75, -85.82),
  stringsAsFactors = FALSE
)

# fst (supp table)
MATE_1_fst <- matrix(c(
  0, .06, .11, .17, .23, .30,
  .06, 0, .04, .12, .20, .27,
  .11, .04, 0, .08, .16, .24,
  .17, .12, .08, 0, .10, .19,
  .23, .20, .16, .10, 0, .14,
  .30, .27, .24, .19, .14, 0
), nrow=6, byrow=TRUE)

rownames(MATE_1_fst) <- colnames(MATE_1_fst) <- MATE_1_coords$site
MATE_1_fst[MATE_1_fst < 0] <- 0

# distances
geo_dist_km <- geosphere::distm(
  as.matrix(MATE_1_coords[,c("lon","lat")]),
  fun=geosphere::distHaversine
)/1000

# ibd
ibd_df <- data.frame(
  fst=MATE_1_fst[upper.tri(MATE_1_fst)],
  dist_km=geo_dist_km[upper.tri(geo_dist_km)]
)

ggplot(ibd_df, aes(dist_km, fst)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  theme_classic()

# save
dir.create(file.path(save_dir,"data"), recursive=TRUE, showWarnings=FALSE)
save(MATE_1_fst, MATE_1_coords,
     file=file.path(save_dir,"data","MATE-1.RData"))