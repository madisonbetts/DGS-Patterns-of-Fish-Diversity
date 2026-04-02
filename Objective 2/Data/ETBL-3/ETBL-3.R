
# -----------------------------
# ETBL-3 | Greenside darter
# Etheostoma blennioides
# Sydenham River watershed
# temporal replicates + added 2006 site
# -----------------------------

library(ggplot2)
library(geosphere)
library(maps)

setwd("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETBL-3")
dir.create("data", showWarnings = FALSE, recursive = TRUE)

site_key <- data.frame(
  site = 1:8,
  code = c("F05", "F06", "C05", "C06", "PD05", "PD06", "CM05", "CS06"),
  site_name = c("Florence 2005", "Florence 2006", "Carolinian 2005", "Carolinian 2006", "Petrolia Dam 2005", "Petrolia Dam 2006", "Cider Mill 2005", "Coldstream 2006"),
  stringsAsFactors = FALSE
)

ETBL_3_coords <- data.frame(
  site = 1:8,
  lat = c(42 + 38/60, 42 + 38/60, 42 + 46/60, 42 + 46/60, 42 + 52/60, 42 + 52/60, 42 + 36/60, 43 + 0/60),
  lon = -c(82 + 0/60, 82 + 0/60, 81 + 50/60, 81 + 50/60, 82 + 8/60, 82 + 8/60, 82 + 4/60, 81 + 30/60)
)

site_ids <- as.character(ETBL_3_coords$site)
ETBL_3_fst <- matrix(0, nrow = 8, ncol = 8, dimnames = list(site_ids, site_ids))

# site order: 1=F05, 2=F06, 3=C05, 4=C06, 5=PD05, 6=PD06, 7=CM05, 8=CS06
fst_rows <- list(
  c(),
  c(0.0547),
  c(0.0117, 0.0218),
  c(0.0533, 0.0025, 0.0237),
  c(0.0426, 0.0070, 0.0115, 0.0115),
  c(0.0669, 0.0152, 0.0302, 0.0080, 0.0232),
  c(0.0049, 0.0533, 0.0172, 0.0581, 0.0506, 0.0743),
  c(0.0708, 0.0155, 0.0365, 0.0125, 0.0270, 0.0052, 0.0794)
)
for (i in 2:8) ETBL_3_fst[i, 1:(i - 1)] <- fst_rows[[i]]
ETBL_3_fst <- ETBL_3_fst + t(ETBL_3_fst)
diag(ETBL_3_fst) <- 0
ETBL_3_fst[ETBL_3_fst < 0] <- 0
rownames(ETBL_3_fst) <- colnames(ETBL_3_fst) <- ETBL_3_coords$site

xy <- as.matrix(ETBL_3_coords[, c("lon", "lat")])
geo_km <- geosphere::distm(xy, fun = geosphere::distHaversine) / 1000
rownames(geo_km) <- colnames(geo_km) <- ETBL_3_coords$site

ibd_df <- data.frame(site1 = rep(ETBL_3_coords$site, each = nrow(ETBL_3_fst)), site2 = rep(ETBL_3_coords$site, times = ncol(ETBL_3_fst)), fst = as.vector(ETBL_3_fst), dist_km = as.vector(geo_km))
ibd_df <- ibd_df[ibd_df$site1 < ibd_df$site2, ]

world <- map_data("world")
world_sub <- subset(world, region %in% c("USA", "Canada"))
plot_df <- merge(ETBL_3_coords, site_key, by = "site", sort = FALSE)

print(ggplot() +
  geom_polygon(data = world_sub, aes(x = long, y = lat, group = group), fill = "grey94", color = "grey70", linewidth = 0.25) +
  geom_point(data = plot_df, aes(x = lon, y = lat), size = 2.8) +
  geom_text(data = plot_df, aes(x = lon, y = lat, label = site), nudge_y = 0.06, size = 3) +
  coord_quickmap(xlim = range(plot_df$lon) + c(-1.0, 1.0), ylim = range(plot_df$lat) + c(-0.8, 0.8)) +
  labs(x = "Longitude", y = "Latitude", title = "ETBL-3 sampling sites") + theme_bw())

print(ggplot(ibd_df, aes(x = dist_km, y = fst)) + geom_point(size = 2) + geom_smooth(method = "lm", se = FALSE) + theme_bw() + labs(x = "Geographic distance (km)", y = expression(F[ST]), title = "ETBL-3 IBD"))

save(ETBL_3_fst, ETBL_3_coords, file = file.path("data", "ETBL-3.RData"))
