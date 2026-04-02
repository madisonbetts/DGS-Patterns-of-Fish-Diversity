
# -----------------------------
# ETBL-4 | Greenside darter
# Etheostoma blennioides
# Ausable River watershed
# temporal replicates + added 2006 site
# -----------------------------

library(ggplot2)
library(geosphere)
library(maps)

setwd("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETBL-4")
dir.create("data", showWarnings = FALSE, recursive = TRUE)

site_key <- data.frame(
  site = 1:7,
  code = c("LP05", "LP06", "B05", "B06", "HH05", "HH06", "S06"),
  site_name = c("Lyon's Park 2005", "Lyon's Park 2006", "Brinsley 2005", "Brinsley 2006", "Hungry Hollow 2005", "Hungry Hollow 2006", "Springbank 2006"),
  stringsAsFactors = FALSE
)

ETBL_4_coords <- data.frame(
  site = 1:7,
  lat = c(43 + 8/60, 43 + 8/60, 43 + 12/60, 43 + 12/60, 43 + 4/60, 43 + 4/60, 43 + 3/60),
  lon = -c(81 + 32/60, 81 + 32/60, 81 + 31/60, 81 + 31/60, 81 + 47/60, 81 + 47/60, 81 + 41/60)
)

site_ids <- as.character(ETBL_4_coords$site)
ETBL_4_fst <- matrix(0, nrow = 7, ncol = 7, dimnames = list(site_ids, site_ids))

# site order: 1=LP05, 2=LP06, 3=B05, 4=B06, 5=HH05, 6=HH06, 7=S06
fst_rows <- list(
  c(),
  c(0.0030),
  c(0.0073, 0.0063),
  c(0.0101, 0.0061, 0.0047),
  c(0.0125, 0.0127, 0.0037, 0.0123),
  c(0.0076, 0.0077, 0.0063, 0.0061, 0.0051),
  c(0.0044, 0.0023, -0.0048, -0.0026, -0.0006, -0.0055)
)
for (i in 2:7) ETBL_4_fst[i, 1:(i - 1)] <- fst_rows[[i]]
ETBL_4_fst <- ETBL_4_fst + t(ETBL_4_fst)
diag(ETBL_4_fst) <- 0
ETBL_4_fst[ETBL_4_fst < 0] <- 0
rownames(ETBL_4_fst) <- colnames(ETBL_4_fst) <- ETBL_4_coords$site

xy <- as.matrix(ETBL_4_coords[, c("lon", "lat")])
geo_km <- geosphere::distm(xy, fun = geosphere::distHaversine) / 1000
rownames(geo_km) <- colnames(geo_km) <- ETBL_4_coords$site

ibd_df <- data.frame(site1 = rep(ETBL_4_coords$site, each = nrow(ETBL_4_fst)), site2 = rep(ETBL_4_coords$site, times = ncol(ETBL_4_fst)), fst = as.vector(ETBL_4_fst), dist_km = as.vector(geo_km))
ibd_df <- ibd_df[ibd_df$site1 < ibd_df$site2, ]

world <- map_data("world")
world_sub <- subset(world, region %in% c("USA", "Canada"))
plot_df <- merge(ETBL_4_coords, site_key, by = "site", sort = FALSE)

print(ggplot() +
  geom_polygon(data = world_sub, aes(x = long, y = lat, group = group), fill = "grey94", color = "grey70", linewidth = 0.25) +
  geom_point(data = plot_df, aes(x = lon, y = lat), size = 2.8) +
  geom_text(data = plot_df, aes(x = lon, y = lat, label = site), nudge_y = 0.06, size = 3) +
  coord_quickmap(xlim = range(plot_df$lon) + c(-0.9, 0.9), ylim = range(plot_df$lat) + c(-0.7, 0.7)) +
  labs(x = "Longitude", y = "Latitude", title = "ETBL-4 sampling sites") + theme_bw())

print(ggplot(ibd_df, aes(x = dist_km, y = fst)) + geom_point(size = 2) + geom_smooth(method = "lm", se = FALSE) + theme_bw() + labs(x = "Geographic distance (km)", y = expression(F[ST]), title = "ETBL-4 IBD"))

save(ETBL_4_fst, ETBL_4_coords, file = file.path("data", "ETBL-4.RData"))
