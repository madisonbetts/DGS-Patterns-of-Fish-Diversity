# -----------------------------
# SAVI-4
# Interior Highlands walleye (Sander vitreus)
# Berkman et al. 2023
# -----------------------------

library(ggplot2)
library(geosphere)
library(maps)

study_code <- "SAVI-4"
out_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SAVI-4/data"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

site_key <- data.frame(
  site = 1:25,
  site_name = c(
    "Clearwater Dam",
    "Markham Spring",
    "Deer Leap",
    "Doniphan",
    "Eleven Point River",
    "Beaver Lake",
    "Upper White",
    "Beaver Tailwater",
    "Powersite Dam",
    "Bull Shoals Dam",
    "James River",
    "Tecumseh",
    "Norfork Lake Dam",
    "South Fork",
    "Middle Fork",
    "Upper Greers Ferry",
    "Lower Greers Ferry",
    "Lake Hamilton",
    "Lake Catherine",
    "Smithville Reservoir",
    "Stockton Lake",
    "Mozingo Lake",
    "Niangua River",
    "Marais Des Cygnes River",
    "Fellows Lake"
  ),
  stringsAsFactors = FALSE
)

# Fig. 2 is displayed in a custom site order, not Table 1 order
fig_order <- c(
  1, 2, 3, 4, 5,
  6, 7, 8, 10, 9, 13, 12, 11,
  19, 18,
  16, 17, 15, 14,
  25, 24, 22, 23, 20, 21
)

fill_row <- function(mat, row, vals) {
  stopifnot(length(vals) == (row - 1))
  mat[row, seq_along(vals)] <- vals
  mat
}

fst_fig <- matrix(0, nrow = 25, ncol = 25)

fst_fig <- fill_row(fst_fig,  2, c(0.033))
fst_fig <- fill_row(fst_fig,  3, c(0.006, 0.025))
fst_fig <- fill_row(fst_fig,  4, c(0.017, 0.021, 0.010))
fst_fig <- fill_row(fst_fig,  5, c(0.004, 0.017, 0.007, 0.002))
fst_fig <- fill_row(fst_fig,  6, c(0.177, 0.166, 0.168, 0.157, 0.181))
fst_fig <- fill_row(fst_fig,  7, c(0.260, 0.228, 0.234, 0.221, 0.243, 0.000))
fst_fig <- fill_row(fst_fig,  8, c(0.164, 0.149, 0.157, 0.145, 0.166, 0.022, 0.058))
fst_fig <- fill_row(fst_fig,  9, c(0.167, 0.157, 0.160, 0.148, 0.170, 0.022, 0.052, 0.000))
fst_fig <- fill_row(fst_fig, 10, c(0.168, 0.164, 0.160, 0.150, 0.168, 0.023, 0.049, 0.003, 0.000))
fst_fig <- fill_row(fst_fig, 11, c(0.174, 0.164, 0.169, 0.158, 0.178, 0.023, 0.052, 0.002, 0.000, 0.001))
fst_fig <- fill_row(fst_fig, 12, c(0.141, 0.137, 0.133, 0.120, 0.144, 0.020, 0.064, 0.008, 0.005, 0.002, 0.003))
fst_fig <- fill_row(fst_fig, 13, c(0.257, 0.245, 0.240, 0.224, 0.249, 0.063, 0.083, 0.043, 0.038, 0.030, 0.036, 0.045))
fst_fig <- fill_row(fst_fig, 14, c(0.127, 0.129, 0.120, 0.109, 0.129, 0.035, 0.063, 0.021, 0.024, 0.027, 0.022, 0.023, 0.071))
fst_fig <- fill_row(fst_fig, 15, c(0.145, 0.138, 0.134, 0.126, 0.147, 0.032, 0.064, 0.017, 0.018, 0.017, 0.018, 0.017, 0.043, 0.004))
fst_fig <- fill_row(fst_fig, 16, c(0.132, 0.124, 0.130, 0.117, 0.139, 0.030, 0.062, 0.011, 0.006, 0.012, 0.006, 0.005, 0.051, 0.011, 0.009))
fst_fig <- fill_row(fst_fig, 17, c(0.167, 0.147, 0.156, 0.135, 0.164, 0.020, 0.058, 0.009, 0.004, 0.008, 0.002, 0.006, 0.034, 0.016, 0.003, 0.000))
fst_fig <- fill_row(fst_fig, 18, c(0.142, 0.135, 0.139, 0.126, 0.139, 0.027, 0.062, 0.020, 0.012, 0.008, 0.024, 0.015, 0.059, 0.023, 0.013, 0.005, 0.009))
fst_fig <- fill_row(fst_fig, 19, c(0.139, 0.108, 0.133, 0.109, 0.133, 0.022, 0.062, 0.006, 0.020, 0.014, 0.010, 0.010, 0.055, 0.009, 0.013, 0.005, 0.000, 0.004))
fst_fig <- fill_row(fst_fig, 20, c(0.251, 0.232, 0.232, 0.228, 0.241, 0.056, 0.073, 0.028, 0.046, 0.025, 0.034, 0.046, 0.072, 0.058, 0.045, 0.050, 0.060, 0.057, 0.052))
fst_fig <- fill_row(fst_fig, 21, c(0.228, 0.192, 0.205, 0.198, 0.220, 0.031, 0.066, 0.022, 0.018, 0.015, 0.012, 0.017, 0.004, 0.047, 0.028, 0.018, 0.009, 0.037, 0.030, 0.066))
fst_fig <- fill_row(fst_fig, 22, c(0.230, 0.209, 0.213, 0.200, 0.224, 0.033, 0.061, 0.024, 0.022, 0.018, 0.017, 0.023, 0.021, 0.056, 0.042, 0.041, 0.033, 0.048, 0.037, 0.052, 0.005))
fst_fig <- fill_row(fst_fig, 23, c(0.181, 0.171, 0.170, 0.164, 0.186, 0.027, 0.059, 0.014, 0.012, 0.007, 0.010, 0.007, 0.024, 0.041, 0.019, 0.014, 0.007, 0.029, 0.029, 0.041, 0.000, 0.011))
fst_fig <- fill_row(fst_fig, 24, c(0.235, 0.232, 0.225, 0.220, 0.235, 0.033, 0.084, 0.033, 0.035, 0.022, 0.026, 0.024, 0.025, 0.054, 0.040, 0.031, 0.033, 0.050, 0.050, 0.047, 0.000, 0.014, 0.006))
fst_fig <- fill_row(fst_fig, 25, c(0.213, 0.203, 0.198, 0.195, 0.215, 0.036, 0.074, 0.033, 0.021, 0.018, 0.021, 0.018, 0.018, 0.053, 0.028, 0.024, 0.020, 0.042, 0.047, 0.069, 0.000, 0.014, 0.002, 0.000))

fst_fig[upper.tri(fst_fig)] <- t(fst_fig)[upper.tri(fst_fig)]
diag(fst_fig) <- 0

perm_idx <- match(1:25, fig_order)
SAVI_4_fst <- fst_fig[perm_idx, perm_idx]
rownames(SAVI_4_fst) <- as.character(1:25)
colnames(SAVI_4_fst) <- as.character(1:25)

# workflow rules
SAVI_4_fst[SAVI_4_fst < 0] <- 0
diag(SAVI_4_fst) <- 0

stopifnot(identical(dim(SAVI_4_fst), c(25L, 25L)))
stopifnot(isTRUE(all.equal(SAVI_4_fst, t(SAVI_4_fst))))

SAVI_4_coords <- data.frame(
  site = 1:25,
  lat = c(
    37.1367166, 36.9812790, 36.6743450, 36.6208000, 36.6486000,
    36.4166700, 36.0485556, 36.4306314, 36.6594000, 36.3647000,
    36.8323000, 36.6231000, 36.2492000, 35.5667000, 35.6469000,
    35.5600000, 35.5206000, 34.4162000, 34.4471100, 39.4493690,
    37.6060523, 40.3542000, 37.6843056, 38.9188000, 37.3192540
  ),
  lon = c(
    -90.7731765, -90.6050070, -90.8858030, -90.8230000, -91.2008000,
    -93.8483300, -93.9742222, -93.8450372, -93.1247000, -92.5781000,
    -93.4466000, -92.2481000, -92.2378000, -92.3834899, -92.3081000,
    -92.1600000, -91.9997000, -93.0653000, -92.9199100, -94.5194078,
    -93.7994327, -94.7691000, -92.9246389, -94.7300000, -93.2147580
  )
)

SAVI_4_coords <- SAVI_4_coords[, c("site", "lat", "lon")]
stopifnot(identical(as.character(SAVI_4_coords$site), rownames(SAVI_4_fst)))

plot(NA,
     xlim = range(SAVI_4_coords$lon) + c(-2, 2),
     ylim = range(SAVI_4_coords$lat) + c(-2, 2),
     xlab = "Longitude",
     ylab = "Latitude",
     asp = 1,
     main = "SAVI-4 sampling sites")

map("world", regions = "Canada", add = TRUE, col = "grey95", fill = TRUE, border = "grey70")
map("state",
    regions = c("missouri", "arkansas", "kansas", "oklahoma", "tennessee", "kentucky", "illinois"),
    add = TRUE, col = "grey90", fill = TRUE, border = "grey60")
points(SAVI_4_coords$lon, SAVI_4_coords$lat, pch = 21, bg = "dodgerblue3", cex = 1.1)
text(SAVI_4_coords$lon, SAVI_4_coords$lat, labels = SAVI_4_coords$site, pos = 3, cex = 0.8)

dist_mat <- geosphere::distm(SAVI_4_coords[, c("lon", "lat")], fun = geosphere::distHaversine) / 1000
lower_idx <- lower.tri(SAVI_4_fst)
ibd_df <- data.frame(distance_km = dist_mat[lower_idx], fst = SAVI_4_fst[lower_idx])

plot(ibd_df$distance_km, ibd_df$fst,
     xlab = "Geographic distance (km)",
     ylab = "Pairwise FST",
     pch = 16,
     main = "SAVI-4 IBD plot")
ibd_lm <- lm(fst ~ distance_km, data = ibd_df)
abline(ibd_lm, col = "red", lwd = 2)

fst_name <- "SAVI_4_fst"
coords_name <- "SAVI_4_coords"
assign(fst_name, SAVI_4_fst, envir = .GlobalEnv)
assign(coords_name, SAVI_4_coords, envir = .GlobalEnv)

save(
  list = c(fst_name, coords_name),
  file = file.path(out_dir, "SAVI-4.RData")
)
