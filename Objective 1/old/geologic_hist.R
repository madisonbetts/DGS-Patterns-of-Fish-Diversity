# compute geologic regions following Gallager et al 2025
# cgh


# load libs
packages <- c("tidyverse", "elevatr", 
              "sf", "terra", "tigris", 
              "parallel", "ggridges")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
# install new packages
if (length(new_packages)) {
  install.packages(new_packages, dependencies = TRUE)
}
# load all packages
invisible(lapply(packages, library, character.only = TRUE))



states <- tigris::states(cb = TRUE, progress_bar = FALSE) |>
  filter(!STUSPS %in% c('HI', 'PR', 'AK', 'MP', 'GU', 'AS', 'VI')) |>
  st_transform(crs = 5070)

# read ice shp
ice <- st_read("Data/hist_shp/ice/ice18k.shp") |>
  st_transform(4326) |>
  st_crop(c(xmin = -125, ymin = 25, xmax = -65, ymax = 50)) |>
  st_transform(crs = 5070) |>
  filter(SYMB == "ICE")

# read ecoregion shp
ecor <- st_read("Data/hist_shp/Aggr_Ecoregions_2015/Aggr_Ecoregions_2015.shp") |>
  st_transform(crs = 5070)

ggplot(ice)+
  geom_sf(fill = "white")
ggplot(ecor)+
  geom_sf(aes(fill = WSA9_NAME)) +
  geom_sf(data = ice, fill = "white", alpha = 0.75)


# read in gendiv
gen_div <- read.csv("Data/filtered_data.csv") %>%
  select(x = "snapped_lon", y = "snapped_lat", Spec_id.x, Study_id, Spec_Latin_GBIF,Spec_Latin_GenDivRange, Genus_GBIF, Family_GBIF, Order_GBIF,
         Spec_common_EOL, Spec_id.y, Pop_id, Geog_1, Geog_2, Latitude, Longitude, 
         N, A, A_mean, A_tot, A_eff, A_private, N_genot, Ar, Ho, He, GD_Nei,
         D_clonal, F_is, F_is_sig, Ploidy) %>%
  filter(!is.na(x))


elev <- get_elev_point(gen_div, prj = 4326, src = "aws")
gen_div$elev <- elev$elevation

######################### nice way to calculate relief if needed
#gen_div$relief <- NA
#for (i in 1:nrow(gen_div)) {
#  x <- get_elev_raster(gen_div[i,1:2], z = 11, prj = 4326, clip = "bbox")
#  x <- rast(x)
#  elev_min <- minmax(x)[1]
#  elev_max <- minmax(x)[2]
#  gen_div$relief[i] <- elev_max - elev_min
#}

gen_div_sf <- st_as_sf(gen_div, coords = c("x", "y"), crs = 4326) |>
  st_transform(crs = 5070)

gen_div_sf <- gen_div_sf |>
  st_join(ice, join = st_intersects) |>
  st_join(ecor, join = st_intersects) |>
  mutate(
    Ice = if_else(SYMB == "ICE", "Glaciated", "Unglaciated", 
                  missing = "Unglaciated"),
    Tectonic = if_else(WSA9 %in% c("WMT", "XER"), 
                       "Active", "Stable"),
    SeaLevel = if_else((WSA9 == "CPL" & Latitude < 40.8), 
                       "Flooded", "Unflooded", missing = "Unflooded"),
    SeaShort = if_else(SeaLevel == "Unflooded", "", "- Flooded"),
    Topography = if_else(elev > 250 , 
                         "Highland", "Lowland"),
    History = paste(Tectonic, Ice, Topography, SeaShort, sep = " "),
    History = fct_relevel(History, 
                          c("Stable Unglaciated Lowland ", 
                            "Stable Unglaciated Highland ",
                            "Stable Unglaciated Lowland - Flooded",
                            "Stable Glaciated Lowland ",
                            "Stable Glaciated Highland ",
                            "Active Unglaciated Lowland ",
                            "Active Unglaciated Highland ",
                            "Active Glaciated Highland "))
  )

gen_div$Ice <- gen_div_sf$Ice
gen_div$Tectonic <- gen_div_sf$Tectonic
gen_div$SeaLevel <- gen_div_sf$SeaLevel
gen_div$Topography <- gen_div_sf$Topography
gen_div$History <- gen_div_sf$History




# save table
write.csv(gen_div, "Data/filtered_data.csv", row.names = FALSE)

# Create vector of colors to use in maps and figures
Colors <- c("dodgerblue", "mediumblue", "violet", "cadetblue1", "cadetblue", 
            "orange", "orange4", "wheat3", "darkorange")

# Map GenDivCoords and color by combined History
p <- ggplot() + 
  geom_sf(data = states, fill = 'white') +
  geom_sf(data = gen_div_sf, aes(color = History), alpha = 0.8) +
  labs(color = "Geologic history") +
  scale_color_manual(values = Colors) +
  theme_bw()

ggsave(p, filename = "Figures/geol_hist.pdf", width = 7, height = 5, dpi = 300)
