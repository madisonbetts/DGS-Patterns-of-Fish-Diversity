# get comids

setwd("~/grad_research/macrogen/DGS-Patterns-of-Fish-Diversity/Objective 1")

packages <- c(
  "sf","dplyr","nhdplusTools", "tigris"
)

new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, dependencies = TRUE)

invisible(lapply(packages, library, character.only = TRUE))

sf::sf_use_s2(TRUE)

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

dat_sf <- dat_sf %>%
  filter(raindrop_pathDist < 1000)
length(unique(dat_sf$Study_id)) # 101
length(unique(dat_sf$Pop_id)) # 1023


states <- tigris::states(cb = TRUE, progress_bar = FALSE) |>
  filter(!STUSPS %in% c('HI', 'PR', 'AK', 'MP', 'GU', 'AS', 'VI')) |>
  st_transform(crs = 5070)

ggplot() +
  geom_sf(data = states, fill = "white") +
  geom_sf(data = dat_sf, aes(color = He))

saveRDS(dat_sf, file = "Data/GDR_w_comids.rds")




