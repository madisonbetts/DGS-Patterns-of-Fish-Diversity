######################################
# Full Dataset Construction Pipeline #
######################################

# authors: Cristian Hernandez, Parinaz

#########
# Setup #
#########

setwd("~/grad_research/macrogen/DGS-Patterns-of-Fish-Diversity/Objective 1")

packages <- c(
  "sf","dplyr","readxl","readr","stringr",
  "nhdplusTools","tigris","terra","elevatr",
  "FedData","units","tidyverse", "arrow", "bit64"
)

new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, dependencies = TRUE)

invisible(lapply(packages, library, character.only = TRUE))

sf::sf_use_s2(TRUE)

##############
# Get COMIDS #
##############

# read in filtered GenDivRange
dat <- read.csv("Data/GenDiv_filtered.csv")

dat$comid <- NA_real_
dat$raindrop_pathDist <- NA_real_

dat_sf <- st_as_sf(dat, coords = c("Longitude", "Latitude"), crs = 4326)


# loop through all pops to grab comid
# Totally throttles the API but this'll work for now
for (i in 1:nrow(dat_sf)){
  id <- tryCatch(discover_nhdplus_id(dat_sf[i,], raindrop = T)[1,], 
                 error = function(e) NULL)
  
  # skip failed requests
  if (is.null(id)) {
    message(paste("Failed at row", i))
    next
  }
  
  dat_sf$comid[i] <- id$comid
  dat_sf$raindrop_pathDist[i] <- id$raindrop_pathDist
  
  message(paste("Completed row", i))
  
}

# rm rows with NA in comid
dat_sf <- dat_sf %>%
  filter(!is.na(comid))
length(unique(dat_sf$Study_id)) # 103
length(unique(dat_sf$Pop_id)) # 1156

# rm rows with a path distance to comid of > 1000
dat_sf <- dat_sf %>%
  filter(raindrop_pathDist < 1000)
length(unique(dat_sf$Study_id)) # 101
length(unique(dat_sf$Pop_id)) # 1023

# quick plot
states <- tigris::states(cb = TRUE, progress_bar = FALSE) |>
  filter(!STUSPS %in% c('HI', 'PR', 'AK', 'MP', 'GU', 'AS', 'VI')) |>
  st_transform(crs = 5070)

ggplot() +
  geom_sf(data = states, fill = "white") +
  geom_sf(data = dat_sf, aes(color = He))

# checkpoint
saveRDS(dat_sf, file = "Data/GDR_w_comids.rds")


#################
# Pull GIS Data #
#################

dat_sf <- readRDS("Data/GDR_w_comids.rds")

# pull landscape data from NHD V2 using comids in dat_sf
file_path <- "../../../../../Downloads/nhdplusv2_landscape_variables_01302026_0952am.csv"

# get column names
col_names <- names(read_csv(file_path, n_max = 0, show_col_types = FALSE))

# default everything to float64
schema_list <- setNames(rep(list(float64()), length(col_names)),
                        col_names)

# override character columns
schema_list[["gnis_name"]] <- string()
schema_list[["ftype"]] <- string()
schema_list[["summer_temp_mean_class"]] <- string()

# build schema
custom_schema <- schema(schema_list)

# open csv
landsc <- open_csv_dataset(file_path, 
                           schema = custom_schema, 
                           skip = 1)

result <- landsc %>%
  filter(comid %in% dat_sf$comid) %>%
  collect()

# filter for relevant predictors
gis_dat <- result %>%
  select(comid, gnis_name, lengthkm, ftype, streamorder,
         matches("lc|nc"), fl_slope, summer_flow_mean_00_15,
         summer_temp_mean_raw, summer_temp_mean_class) %>%
  # calculate elevation and relief
  mutate(lc_meanelev = (lc_minelev + lc_maxelev) / 2,
         lc_relief = lc_maxelev - lc_minelev)

# checkpoint
saveRDS(gis_dat, file = "Data/filtered_gis_dat.rds")


dat_sf <- left_join(dat_sf, gis_dat, by = "comid")
dat_sf$Latitude <- st_coordinates(dat_sf$geometry)[,2]


ice <- st_read("Data/hist_shp/ice/ice18k.shp") %>%
  st_transform(4326) %>%
  filter(SYMB == "ICE")

ecor <- st_read("Data/hist_shp/Aggr_Ecoregions_2015/Aggr_Ecoregions_2015.shp") %>%
  st_transform(4326)

ice  <- st_make_valid(ice)
ecor <- st_make_valid(ecor)


dat_sf <- dat_sf |>
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
    Topography = if_else(lc_meanelev > 250 , 
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
                            "Active Glaciated Highland "
                            ))
  )






##############
# FishTraits #
##############

fish_traits <- read_excel("Data/FishTraits_14.3.xls") %>%
  mutate(
    ScientificName = str_to_lower(str_trim(paste(GENUS, SPECIES)))
  ) %>%
  select(ScientificName, AREAKM2, PATCHES, FECUNDITY, LONGEVITY, MAXTL)

dat_sf <- dat_sf %>%
  mutate(
    Spec_Latin_clean = str_to_lower(str_trim(Spec_Latin_GenDivRange))
  ) %>%
  left_join(fish_traits, by = c("Spec_Latin_clean" = "ScientificName"))



######################
# Remove weird stuff #
######################

taxa_to_remove <- c("Neogobius melanostomus", "Syngnathus scovelli")

dat_sf <- dat_sf %>%
  filter(!Spec_Latin_GenDivRange %in% taxa_to_remove) # 99 studies, 1016 sites

##########
# Export #
##########

final_df <- st_drop_geometry(dat_sf)

write_csv(final_df, "Data/final_fish_dataset.csv")

st_write(dat_sf,
         "Data/final_fish_dataset.gpkg",delete_dsn = TRUE)
