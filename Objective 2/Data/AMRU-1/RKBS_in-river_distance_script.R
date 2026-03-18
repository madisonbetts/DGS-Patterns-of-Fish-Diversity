#--------------------------
# RKBS in-river distances
#--------------------------
library(terra)
library(sf)
library(dplyr)
library(FedData)
library(magrittr)

# read in rockbass data
load("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/datasets/RKBS-1/RKBS-1_genetic_distances.RData")

## -------------------------
## 1) lat/lon -> SpatVector
## -------------------------
rkbs_pts <- vect(
  rkbs_sites,
  geom = c("lon", "lat"),
  crs  = "EPSG:4326"
)

## -------------------------
## 2) download WBD for points
##    (this returns HUC12s)
## -------------------------
rkbs_wbd_huc12 <- get_wbd(
  template       = rkbs_pts,
  label          = "rkbs_points",
  extraction.dir = "data/FedData/wbd_rkbs",
  force.redo     = TRUE
) %>%
  vect()

# project pts to a UTM (meters0)
rkbs_pts <- rkbs_pts %>%
  project(crs(rkbs_wbd_huc12))

## reproject HUC12s to match point CRS
rkbs_wbd_huc12 <- project(rkbs_wbd_huc12, crs(rkbs_pts))



## inspect names to find HUC12 field if needed
names(rkbs_wbd_huc12)

## -------------------------
## 3) derive HUC2 from HUC12
## -------------------------
## FedData/WBD usually has a HUC12 column, but field names can vary.
## this grabs the first matching HUC column in a reasonably safe way.

huc_col <- names(rkbs_wbd_huc12)[grepl("^HUC", names(rkbs_wbd_huc12), ignore.case = TRUE)][1]

if (is.na(huc_col)) {
  stop("No HUC column found in rkbs_wbd_huc12")
}

rkbs_huc2 <- rkbs_wbd_huc12 |>
  mutate(HUC2 = substr(.data[[huc_col]], 1, 2)) |>
  group_by(HUC2) |>
  summarise(do_union = TRUE, .groups = "drop")

## optional check
print(rkbs_huc2)

## -------------------------
## 4) download NHD for HUC2 extent
## -------------------------
rkbs_nhd <- get_nhd(
  template       = rkbs_huc2,
  label          = "rkbs_huc2_extent",
  nhdplus        = FALSE,
  extraction.dir = "data/FedData/nhd_rkbs",
  force.redo     = FALSE
)

## -------------------------
## 5) inspect outputs
## -------------------------
names(rkbs_nhd)

## often includes things like:
## rkbs_nhd$NHDFlowline
## rkbs_nhd$NHDArea
## rkbs_nhd$NHDWaterbody
## rkbs_nhd$NHDLine

## quick map of points + HUC2 + flowlines if present
plot(st_geometry(rkbs_huc2), border = "red", lwd = 2)
plot(rkbs_sv, add = TRUE, pch = 16)

if ("NHDFlowline" %in% names(rkbs_nhd)) {
  plot(st_geometry(rkbs_nhd$NHDFlowline), add = TRUE, col = "blue")
}
