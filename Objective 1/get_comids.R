# get comids


packages <- c(
  "sf","dplyr","readxl","readr","stringr",
  "nhdplusTools","tigris","terra","elevatr",
  "FedData","units","tidyverse", "arrow", "bit64"
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
for (i in 277:420){
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


dat_sf_subset <- dat_sf[277:420, ]

file_path <- "../../../../../Downloads/nhdplusv2_landscape_variables_01302026_0952am.csv"

# get column names
col_names <- names(
  read_csv(file_path, n_max = 0, show_col_types = FALSE)
)

# default everything to float64
schema_list <- setNames(
  rep(list(float64()), length(col_names)),
  col_names
)

# override character columns
schema_list[["gnis_name"]] <- string()
schema_list[["ftype"]] <- string()
schema_list[["summer_temp_mean_class"]] <- string()

# build schema
custom_schema <- schema(schema_list)

# open csv
landsc2 <- open_csv_dataset(
  file_path,
  schema = custom_schema,
  skip = 1
)
result <- landsc2 %>%
  filter(comid %in% dat_sf_subset$comid) %>%
  collect()




