# filtering for obj 1
# cgh

# install/load pacakges
packages <- c("sf", "dplyr", 
              "readxl", "nhdplusTools", 
              "units")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
# install new packages
if (length(new_packages)) {
  install.packages(new_packages, dependencies = TRUE)
}
# load all packages
invisible(lapply(packages, library, character.only = TRUE))

sf::sf_use_s2(TRUE)


# read in and filter data
sp <- read_xlsx("Data/spec_tab_v2025-03-31.xlsx")
pop <- read_excel("Data/pop_tab_v2025-03-31.xlsx")
# return only populations in the conterminous US
pop <- pop %>%
  filter(Geog_1 == "United States of America", Latitude < 50, Latitude > 25,
         !grepl("lake|bay", Geog_2, ignore.case = TRUE))
# return only freshwater fishes
sp <- sp %>% 
  filter(Life_form == "Fishes", 
         Habitat_breeding == "Freshwater", 
         Habitat_adulthood == "Freshwater")
# merge pop and dat by dat$Study_id
dat <- inner_join(sp, pop, by = "Study_id")



## get min distance to flowline for loop

# init
dat$min_dist_to_flowline <- NA
dat$snapped_lon <- NA
dat$snapped_lat <- NA

# iterate through all rows
for (i in 1:nrow(dat)){
  lat  <- dat[i,]$Latitude
  long <- dat[i,]$Longitude
  
  # point in WGS84
  point <- st_sfc(st_point(c(long, lat)), crs = 4326)
  
  # get flowlines
  huc <- tryCatch(get_huc(AOI = point), 
                  error = function(e) NULL)
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
      
      # only snap if within 200 m
      if (min_dist <= units::set_units(800, "m")) {
        
        nearest_line <- flowlines_proj[min_idx, ]
        
        # nearest point pair
        nearest_pair <- st_nearest_points(point_proj, nearest_line)
        
        # extract snapped point (second point)
        snapped_point_proj <- st_cast(nearest_pair, "POINT")[2]
        
        # back-transform
        snapped_point <- st_transform(snapped_point_proj, 4326)
        coords <- st_coordinates(snapped_point)
        
        dat[i,]$snapped_lon <- coords[1]
        dat[i,]$snapped_lat <- coords[2]
        
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



# write it to csv
write.csv(dat, "Data/filtered_data.csv", row.names = FALSE)

