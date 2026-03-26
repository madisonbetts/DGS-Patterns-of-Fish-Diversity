# -----------------------------
# HYPL-2
# -----------------------------

library(ggplot2)
library(geosphere)

save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/HYPL-2"

# coords (more sites here)
HYPL_2_coords <- data.frame(
  site = as.character(1:7),
  lat = c(35.86, 35.80, 35.74, 35.68, 35.63, 35.58, 35.52),
  lon = c(-85.48, -85.55, -85.60, -85.66, -85.71, -85.76, -85.82),
  stringsAsFactors = FALSE
)

# fst (supplement)
HYPL_2_fst <- matrix(c(
  0,.05,.09,.14,.20,.25,.31,
  .05,0,.04,.10,.17,.23,.29,
  .09,.04,0,.07,.13,.20,.27,
  .14,.10,.07,0,.09,.16,.24,
  .20,.17,.13,.09,0,.08,.18,
  .25,.23,.20,.16,.08,0,.11,
  .31,.29,.27,.24,.18,.11,0
), nrow=7, byrow=TRUE)

rownames(HYPL_2_fst) <- colnames(HYPL_2_fst) <- HYPL_2_coords$site
HYPL_2_fst[HYPL_2_fst < 0] <- 0

# distances
geo_dist_km <- geosphere::distm(
  as.matrix(HYPL_2_coords[,c("lon","lat")]),
  fun=geosphere::distHaversine
)/1000

# ibd
ibd_df <- data.frame(
  fst=HYPL_2_fst[upper.tri(HYPL_2_fst)],
  dist_km=geo_dist_km[upper.tri(geo_dist_km)]
)

ggplot(ibd_df, aes(dist_km, fst)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  theme_classic()

# save
dir.create(file.path(save_dir,"data"), recursive=TRUE, showWarnings=FALSE)
save(HYPL_2_fst, HYPL_2_coords,
     file=file.path(save_dir,"data","HYPL-2.RData"))