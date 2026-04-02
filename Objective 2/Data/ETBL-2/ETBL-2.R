# -----------------------------
# ETBL-2 | Greenside darter
# Etheostoma blennioides
# Grand River watershed (2005)
# Beneteau et al. 2009 Conserv Genet / thesis cross-check
# -----------------------------

library(ggplot2)
library(geosphere)
library(maps)

setwd("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETBL-2")
dir.create("data", showWarnings = FALSE, recursive = TRUE)

site_key <- data.frame(
  site = 1:7,
  code = c("BP", "SJ", "T", "W", "D", "FP", "MW"),
  site_name = c("Bean Park", "St. Jacob's", "Trussler", "Woodlawn", "Doon", "Freeport", "Mannheim Weir"),
  stringsAsFactors = FALSE
)

ETBL_2_coords <- data.frame(
  site = 1:7,
  lat = c(43 + 10/60, 43 + 32/60, 43 + 17/60, 43 + 34/60, 43 + 23/60, 43 + 24/60, 43 + 25/60),
  lon = -c(80 + 22/60, 80 + 33/60, 80 + 28/60, 80 + 16/60, 80 + 23/60, 80 + 26/60, 80 + 25/60)
)

site_ids <- as.character(ETBL_2_coords$site)
ETBL_2_fst <- matrix(0, nrow = 7, ncol = 7, dimnames = list(site_ids, site_ids))

# site order: 1=BP, 2=SJ, 3=T, 4=W, 5=D, 6=FP, 7=MW
fst_rows <- list(
  c(),
  c(0.0037),
  c(0.0019, 0.0091),
  c(0.0115, 0.0238, 0.0137),
  c(0.0002, 0.0041, 0.0045, 0.0069),
  c(0.0037, 0.0035, 0.0067, 0.0164, 0.0038),
  c(0.0084, 0.0029, 0.0113, 0.0224, 0.0126, 0.0006)
)
for (i in 2:7) ETBL_2_fst[i, 1:(i - 1)] <- fst_rows[[i]]
ETBL_2_fst <- ETBL_2_fst + t(ETBL_2_fst)
diag(ETBL_2_fst) <- 0
ETBL_2_fst[ETBL_2_fst < 0] <- 0
rownames(ETBL_2_fst) <- colnames(ETBL_2_fst) <- ETBL_2_coords$site

xy <- as.matrix(ETBL_2_coords[, c("lon", "lat")])
geo_km <- geosphere::distm(xy, fun = geosphere::distHaversine) / 1000
rownames(geo_km) <- colnames(geo_km) <- ETBL_2_coords$site

ibd_df <- data.frame(site1 = rep(ETBL_2_coords$site, each = nrow(ETBL_2_fst)), site2 = rep(ETBL_2_coords$site, times = ncol(ETBL_2_fst)), fst = as.vector(ETBL_2_fst), dist_km = as.vector(geo_km))
ibd_df <- ibd_df[ibd_df$site1 < ibd_df$site2, ]

world <- map_data("world")
world_sub <- subset(world, region %in% c("USA", "Canada"))
plot_df <- merge(ETBL_2_coords, site_key, by = "site", sort = FALSE)

print(ggplot() +
  geom_polygon(data = world_sub, aes(x = long, y = lat, group = group), fill = "grey94", color = "grey70", linewidth = 0.25) +
  geom_point(data = plot_df, aes(x = lon, y = lat), size = 2.8) +
  geom_text(data = plot_df, aes(x = lon, y = lat, label = site), nudge_y = 0.06, size = 3) +
  coord_quickmap(xlim = range(plot_df$lon) + c(-1.0, 1.0), ylim = range(plot_df$lat) + c(-0.8, 0.8)) +
  labs(x = "Longitude", y = "Latitude", title = "ETBL-2 sampling sites") + theme_bw())

print(ggplot(ibd_df, aes(x = dist_km, y = fst)) + geom_point(size = 2) + geom_smooth(method = "lm", se = FALSE) + theme_bw() + labs(x = "Geographic distance (km)", y = expression(F[ST]), title = "ETBL-2 IBD"))

save(ETBL_2_fst, ETBL_2_coords, file = file.path("data", "ETBL-2.RData"))
