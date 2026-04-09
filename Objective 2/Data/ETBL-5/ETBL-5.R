
# -----------------------------
# ETBL-5 | Greenside darter
# Etheostoma blennioides
# Thames River watershed (2006)
# -----------------------------

library(ggplot2)
library(geosphere)
library(maps)

setwd("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETBL-5")
dir.create("data", showWarnings = FALSE, recursive = TRUE)

site_key <- data.frame(
  site = 1:4,
  code = c("SM", "RP", "PC", "MC"),
  site_name = c("St. Mary's", "Roth Park", "Piney Creek", "Medway Creek"),
  stringsAsFactors = FALSE
)

ETBL_5_coords <- data.frame(
  site = 1:4,
  lat = c(43 + 17/60, 43 + 8/60, 42 + 58/60, 43 + 0/60),
  lon = -c(81 + 9/60, 80 + 45/60, 80 + 57/60, 81 + 16/60)
)

site_ids <- as.character(ETBL_5_coords$site)
ETBL_5_fst <- matrix(0, nrow = 4, ncol = 4, dimnames = list(site_ids, site_ids))

# site order: 1=SM, 2=RP, 3=PC, 4=MC
fst_rows <- list(
  c(),
  c(0.0080),
  c(0.0050, 0.0030),
  c(0.0050, 0.0070, -0.0050)
)
for (i in 2:4) ETBL_5_fst[i, 1:(i - 1)] <- fst_rows[[i]]
ETBL_5_fst <- ETBL_5_fst + t(ETBL_5_fst)
diag(ETBL_5_fst) <- 0
ETBL_5_fst[ETBL_5_fst < 0] <- 0
rownames(ETBL_5_fst) <- colnames(ETBL_5_fst) <- ETBL_5_coords$site

xy <- as.matrix(ETBL_5_coords[, c("lon", "lat")])
geo_km <- geosphere::distm(xy, fun = geosphere::distHaversine) / 1000
rownames(geo_km) <- colnames(geo_km) <- ETBL_5_coords$site

ibd_df <- data.frame(site1 = rep(ETBL_5_coords$site, each = nrow(ETBL_5_fst)), site2 = rep(ETBL_5_coords$site, times = ncol(ETBL_5_fst)), fst = as.vector(ETBL_5_fst), dist_km = as.vector(geo_km))
ibd_df <- ibd_df[ibd_df$site1 < ibd_df$site2, ]

world <- map_data("world")
world_sub <- subset(world, region %in% c("USA", "Canada"))
plot_df <- merge(ETBL_5_coords, site_key, by = "site", sort = FALSE)

print(ggplot() +
  geom_polygon(data = world_sub, aes(x = long, y = lat, group = group), fill = "grey94", color = "grey70", linewidth = 0.25) +
  geom_point(data = plot_df, aes(x = lon, y = lat), size = 2.8) +
  geom_text(data = plot_df, aes(x = lon, y = lat, label = site), nudge_y = 0.06, size = 3) +
  coord_quickmap(xlim = range(plot_df$lon) + c(-1.0, 1.0), ylim = range(plot_df$lat) + c(-0.7, 0.7)) +
  labs(x = "Longitude", y = "Latitude", title = "ETBL-5 sampling sites") + theme_bw())

print(ggplot(ibd_df, aes(x = dist_km, y = fst)) + geom_point(size = 2) + geom_smooth(method = "lm", se = FALSE) + theme_bw() + labs(x = "Geographic distance (km)", y = expression(F[ST]), title = "ETBL-5 IBD"))

save(ETBL_5_fst, ETBL_5_coords, file = file.path("data", "ETBL-5.RData"))
