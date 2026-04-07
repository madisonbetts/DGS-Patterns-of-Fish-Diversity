######################################
# Full Dataset Construction Pipeline #
######################################

# authors: Cristian Hernandez, Parinaz
# Last modified: 4/7/2026 

#########
# Setup #
#########

packages <- c(
  "sf","dplyr","readxl","readr","stringr",
  "nhdplusTools","tigris","terra","elevatr",
  "FedData","units"
)

new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, dependencies = TRUE)

invisible(lapply(packages, library, character.only = TRUE))

sf::sf_use_s2(TRUE)


########################
# Load and filter Data #
########################

sp  <- read_xlsx("Data/spec_tab_v2025-03-31.xlsx")
pop <- read_xlsx("Data/pop_tab_v2025-03-31.xlsx")

pop <- pop %>%
  filter(
    Geog_1 == "United States of America",
    Latitude < 50, Latitude > 25,
    !grepl("lake|bay", Geog_2, ignore.case = TRUE)
  )

sp <- sp %>%
  filter(
    Life_form == "Fishes",
    Habitat_breeding == "Freshwater"
  )

dat <- inner_join(sp, pop, by = "Study_id")



######################################################
# if nearest feature is a not a coastline {         #
#   snap coords to that flowline                    #
#   extract stream order, catchment area, and slope #
# }                                                 #
#####################################################

# init vars to populate
dat$min_dist_to_flowline <- NA
dat$snapped_lon <- NA
dat$snapped_lat <- NA
dat$stream_order <- NA
dat$catch_area <- NA
dat$slope <- NA

# iterate through all rows
for (i in 1:nrow(dat)){
  lat  <- dat[i,]$Latitude
  long <- dat[i,]$Longitude
  
  # point in WGS84
  point <- st_sfc(st_point(c(long, lat)), crs = 4326)
  
  # extract huc
  huc <- tryCatch(get_huc(AOI = point), 
                  error = function(e) NULL)
  # get the flowlines
  flowlines <- tryCatch(get_nhdplus(AOI = huc[1,]), 
                        error = function(e) NULL)
  
  if (length(flowlines) > 0) {
    
    # project
    flowlines_proj <- st_transform(flowlines, 5070)
    point_proj <- st_transform(point, 5070)
    
    # compute and index distances
    dists <- st_distance(point_proj, flowlines_proj)
    min_idx <- which.min(dists)
    
    if (flowlines[min_idx,]$ftype != "Coastline"){
      
      # store min distance
      min_dist <- dists[min_idx]
      dat[i,]$min_dist_to_flowline <- as.numeric(min_dist)
      
      # only snap if within 500 m
      if (min_dist <= units::set_units(500, "m")) {
        
        nearest_line <- flowlines_proj[min_idx, ]
        
        # nearest point pair
        nearest_pair <- st_nearest_points(point_proj, nearest_line)
        
        # extract snapped point (second point)
        snapped_point_proj <- st_cast(nearest_pair, "POINT")[2]
        
        # back-transform
        snapped_point <- st_transform(snapped_point_proj, 4326)
        coords <- st_coordinates(snapped_point)
        
        # extract relevant data
        dat[i,]$snapped_lon <- coords[1]
        dat[i,]$snapped_lat <- coords[2]
        dat[i,]$stream_order <- nearest_line$streamorde
        dat[i,]$catch_area <- nearest_line$areasqkm
        dat[i,]$slope <- nearest_line$slope
        
      } else {
        # keep original if beyond threshold
        dat[i,]$snapped_lon <- NA
        dat[i,]$snapped_lat <- NA
      }
      
    } else {
      dat[i,]$min_dist_to_flowline <- NA
      dat[i,]$snapped_lon <- NA
      dat[i,]$snapped_lat <- NA
    }
    
  } else {
    dat[i,]$min_dist_to_flowline <- NA
    dat[i,]$snapped_lon <- NA
    dat[i,]$snapped_lat <- NA
  }
}

#################
# Convert to sf #
#################

gen_div <- dat %>%
  filter(!is.na(snapped_lon)) %>%
  rename(x = Longitude, y = Latitude)

gen_div_sf <- st_as_sf(gen_div, coords = c("x","y"), crs = 4326) %>%
  st_transform(5070)

##################### get rid of this
#library(ggplot2)
#ggplot()+
#  geom_sf(data = gen_div_sf)

###############################
# Geologic & Topographic Data #
###############################

ice <- st_read("Data/hist_shp/ice/ice18k.shp") %>%
  st_transform(5070) %>%
  filter(SYMB == "ICE")

ecor <- st_read("Data/hist_shp/Aggr_Ecoregions_2015/Aggr_Ecoregions_2015.shp") %>%
  st_transform(5070)

gen_div_sf <- gen_div_sf %>%
  st_join(ice) %>%
  st_join(ecor) %>%
  mutate(
    Ice = if_else(SYMB == "ICE","Glaciated","Unglaciated", missing = "Unglaciated"),
    Tectonic = if_else(WSA9 %in% c("WMT","XER"),"Active","Stable"),
    SeaLevel = if_else(WSA9 == "CPL" & snapped_lat < 40.8,"Flooded","Unflooded")
  )


elev <- get_elev_point(gen_div_sf, prj = 5070, src = "aws")
gen_div_sf$elev <- elev$elevation

gen_div_sf$relief <- NA
for (i in 1:nrow(gen_div_sf)) {
  x <- get_elev_raster(gen_div_sf[i,], z = 11, prj = 5070, clip = "bbox")
  x <- rast(x)
  elev_min <- minmax(x)[1]
  elev_max <- minmax(x)[2]
  gen_div_sf$relief[i] <- elev_max - elev_min
}

gen_div_sf$Topography <- NA
gen_div_sf$Topography <- ifelse(gen_div_sf$elev > 250, "Highland", "Lowland")


gen_div_sf <- gen_div_sf %>%
  mutate(
    SeaShort = if_else(SeaLevel == "Unflooded","", "- Flooded"),
    History = paste(Tectonic, Ice, Topography, SeaShort)
  )

##############
# FishTraits #
##############

fish_traits <- read_excel("Data/FishTraits_14.3.xls") %>%
  mutate(
    ScientificName = str_to_lower(str_trim(paste(GENUS, SPECIES)))
  ) %>%
  select(ScientificName, AREAKM2, PATCHES, FECUNDITY, LONGEVITY, MAXTL)

gen_div_sf <- gen_div_sf %>%
  mutate(
    Spec_Latin_clean = str_to_lower(str_trim(Spec_Latin_GenDivRange))
  ) %>%
  left_join(fish_traits, by = c("Spec_Latin_clean" = "ScientificName"))



#####################
# Remove weird taxa #
#####################

taxa_to_remove <- c("Neogobius melanostomus", "Syngnathus scovelli")

gen_div_sf <- gen_div_sf %>%
  filter(!Spec_Latin_GenDivRange %in% taxa_to_remove)


##########
# Export #
##########

final_df <- st_drop_geometry(gen_div_sf)

write_csv(final_df, "Data/final_fish_dataset.csv")

st_write(
  gen_div_sf,
  "Data/final_fish_dataset.gpkg",
  delete_dsn = TRUE
)

