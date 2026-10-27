# initial filtering

setwd("~/grad_research/macrogen/DGS-Patterns-of-Fish-Diversity/Objective 1")

packages <- c("sf", "dplyr", 
              "readxl", "nhdplusTools", 
              "units", "rnaturalearth",
              "ape", "fishtree", "colorspace",
              "taxize", "phytools", "tigris", "ggplot2")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(new_packages)) {
  install.packages(new_packages, dependencies = TRUE)
}
invisible(lapply(packages, library, character.only = TRUE))

sf::sf_use_s2(FALSE)

sp  <- read_xlsx("Data/spec_tab_v2025-03-31.xlsx")
nrow(sp) # 1480

sp <- sp %>%
  filter(Life_form == "Fishes")
nrow(sp) #329

sp <- sp %>%
  filter(Habitat_breeding == "Freshwater") %>%
  filter(!(Spec_common_EOL == "Bull Shark"))
nrow(sp) #253

pop <- read_xlsx("Data/pop_tab_v2025-03-31.xlsx")

pop <- pop %>%
  filter(
    Geog_1 == "United States of America",
    Latitude < 50, Latitude > 25
  ) %>%
  filter(!(Latitude >= 46 & Longitude >= -80 & Longitude <= -72))

dat <- inner_join(sp, pop, by = "Study_id")
length(unique(dat$Study_id)) #110
length(unique(dat$Pop_id)) #1304

states <- tigris::states(cb = TRUE, progress_bar = FALSE) |>
  filter(!STUSPS %in% c('HI', 'PR', 'AK', 'MP', 'GU', 'AS', 'VI')) |>
  st_transform(crs = 5070)

dat_sf <- st_as_sf(dat, coords = c("Longitude", "Latitude"), crs = 4326) |>
  st_transform(crs = 5070)

# vis
ggplot() +
  geom_sf(data = states, fill = "white") +
  geom_sf(data = dat_sf, aes(color = He))+
  scale_color_viridis_c(
    option = "magma",
    limits = c(min(dat$He), max(dat$He)),
    breaks = seq(0.3, 0.7, by = 0.1)
  ) +
  theme_minimal()

# write to csv
filt_dat <- dat %>% 
  select(Data_source, Spec_id.x, Study_id, Spec_Latin_GenDivRange, DOI_study,
         DOI_data, Spec_Latin_publication, N_pops, N_loci, Spec_Latin_GBIF, Genus_GBIF,
         Family_GBIF, Order_GBIF, Class_GBIF, Phylum_GBIF, Spec_common_EOL, BIOME, Habitat_fishbase,
         Habitat_breeding, Habitat_adulthood, Pop_id, Geog_2, Latitude, Longitude,
         N, A, A_mean, A_tot, A_eff, A_private, N_genot, Ar, Ho, He)

write.csv(filt_dat, "Data/GenDiv_filtered.csv")



