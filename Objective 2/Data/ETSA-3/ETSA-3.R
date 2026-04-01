# ========================================
# ETSA-3
# Etheostoma sagitta
# Kentucky Arrow Darter
#
# Builds:
#   ETSA_3_coords
#   ETSA_3_fst
#   ETSA_3_rivdists
#
# Saves to:
#   ETSA-3/data/ETSA-3.RData
#
# Plots in RStudio:
#   - site map
#   - IBD using in-river distance only
# ========================================

library(HWxtest)
library(adegenet)
library(hierfstat)
library(ggplot2)
library(geosphere)
library(maps)

study_code <- "ETSA-3"

base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETSA-3"
river_dir <- file.path(base_dir, "river_data")
data_dir  <- file.path(base_dir, "data")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

gen_file     <- file.path(river_dir, "KAD_Watson_GenepopFile_12popRearranged.txt")
coords_file  <- file.path(river_dir, "SiteCoords.csv")
rivdist_file <- file.path(river_dir, "riverdist_correctsitenumbers.csv")

stopifnot(file.exists(gen_file))
stopifnot(file.exists(coords_file))
stopifnot(file.exists(rivdist_file))

# -----------------------------
# 1) read exact site coordinates
# -----------------------------
site_info <- read.csv(coords_file, stringsAsFactors = FALSE)
names(site_info) <- trimws(names(site_info))

site_col   <- names(site_info)[tolower(names(site_info)) %in% c("site", "site_id", "id")]
stream_col <- names(site_info)[tolower(names(site_info)) %in% c("stream", "site name", "site_name")]
lat_col    <- names(site_info)[grepl("^lat$|latitude", tolower(names(site_info)))]
lon_col    <- names(site_info)[grepl("^long$|^lon$|longitude", tolower(names(site_info)))]

stopifnot(length(site_col) >= 1, length(stream_col) >= 1, length(lat_col) >= 1, length(lon_col) >= 1)

site_col   <- site_col[1]
stream_col <- stream_col[1]
lat_col    <- lat_col[1]
lon_col    <- lon_col[1]

site_info <- data.frame(
  site_orig = as.integer(site_info[[site_col]]),
  stream = as.character(site_info[[stream_col]]),
  lat = as.numeric(site_info[[lat_col]]),
  lon = as.numeric(site_info[[lon_col]]),
  stringsAsFactors = FALSE
)

expected_orig_order <- as.integer(c(1, 2, 3, 4, 8, 7, 9, 10, 11, 12, 13, 14))
stopifnot(length(site_info$site_orig) == length(expected_orig_order))
stopifnot(all(site_info$site_orig == expected_orig_order))

# -----------------------------
# 2) read Genepop and compute pairwise FST
# -----------------------------
et3_gen <- HWxtest::genepop.to.genind(
  name = gen_file,
  quiet = TRUE,
  ncode = 2
)

pop_sizes <- table(adegenet::pop(et3_gen))
stopifnot(length(pop_sizes) == nrow(site_info))

et3_hf <- hierfstat::genind2hierfstat(et3_gen)
fst_named <- as.matrix(hierfstat::pairwise.WCfst(et3_hf))

stopifnot(nrow(fst_named) == nrow(site_info))
stopifnot(ncol(fst_named) == nrow(site_info))

if (!is.null(rownames(fst_named)) && !is.null(colnames(fst_named))) {
  rn_num <- suppressWarnings(as.integer(gsub("[^0-9]", "", rownames(fst_named))))
  cn_num <- suppressWarnings(as.integer(gsub("[^0-9]", "", colnames(fst_named))))
  if (all(!is.na(rn_num)) && all(!is.na(cn_num)) &&
      setequal(rn_num, seq_len(nrow(site_info))) &&
      setequal(cn_num, seq_len(nrow(site_info)))) {
    ord <- order(rn_num)
    fst_named <- fst_named[ord, ord, drop = FALSE]
  }
}

fst_named[is.na(fst_named)] <- 0
fst_named[fst_named < 0] <- 0
diag(fst_named) <- 0
fst_named[lower.tri(fst_named)] <- t(fst_named)[lower.tri(fst_named)]

# -----------------------------
# 3) read in-river distance matrix
# -----------------------------
riv_raw <- read.csv(rivdist_file, stringsAsFactors = FALSE, check.names = FALSE)
names(riv_raw)[1] <- "row_id"

riv_row_ids <- suppressWarnings(as.integer(riv_raw$row_id))
riv_col_ids <- suppressWarnings(as.integer(names(riv_raw)[-1]))

keep_cols <- !is.na(riv_col_ids)
riv_col_ids <- riv_col_ids[keep_cols]

riv_vals <- as.matrix(riv_raw[, c(TRUE, keep_cols), drop = FALSE][, -1, drop = FALSE])
mode(riv_vals) <- "numeric"

riv_vals <- riv_vals[!is.na(riv_row_ids), , drop = FALSE]
riv_row_ids <- riv_row_ids[!is.na(riv_row_ids)]

rownames(riv_vals) <- riv_row_ids
colnames(riv_vals) <- riv_col_ids

stopifnot(setequal(as.integer(rownames(riv_vals)), site_info$site_orig))
stopifnot(setequal(as.integer(colnames(riv_vals)), site_info$site_orig))

riv_vals <- riv_vals[as.character(site_info$site_orig), as.character(site_info$site_orig), drop = FALSE]
riv_vals[lower.tri(riv_vals)] <- t(riv_vals)[lower.tri(riv_vals)]
diag(riv_vals) <- 0

# -----------------------------
# 4) build final objects
# -----------------------------
ETSA_3_coords <- data.frame(
  site = as.character(seq_len(nrow(site_info))),
  lat = site_info$lat,
  lon = site_info$lon,
  stringsAsFactors = FALSE
)

site_lookup <- data.frame(
  site = ETSA_3_coords$site,
  site_orig = site_info$site_orig,
  stream = site_info$stream,
  stringsAsFactors = FALSE
)

ETSA_3_fst <- fst_named
ETSA_3_rivdists <- riv_vals

rownames(ETSA_3_fst) <- ETSA_3_coords$site
colnames(ETSA_3_fst) <- ETSA_3_coords$site
rownames(ETSA_3_rivdists) <- ETSA_3_coords$site
colnames(ETSA_3_rivdists) <- ETSA_3_coords$site

# -----------------------------
# 5) pairwise dataframe for in-river IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1 = rownames(ETSA_3_fst)[row(ETSA_3_fst)[upper.tri(ETSA_3_fst)]],
  site2 = colnames(ETSA_3_fst)[col(ETSA_3_fst)[upper.tri(ETSA_3_fst)]],
  fst = ETSA_3_fst[upper.tri(ETSA_3_fst)],
  rivdist_km = ETSA_3_rivdists[upper.tri(ETSA_3_rivdists)] / 1000,
  stringsAsFactors = FALSE
)

ibd_df$site1_name <- site_lookup$stream[match(ibd_df$site1, site_lookup$site)]
ibd_df$site2_name <- site_lookup$stream[match(ibd_df$site2, site_lookup$site)]

# -----------------------------
# 6) map
# -----------------------------
state_map <- map_data("state")
plot_sites <- merge(ETSA_3_coords, site_lookup, by = "site", sort = FALSE)

lon_rng <- range(plot_sites$lon, na.rm = TRUE)
lat_rng <- range(plot_sites$lat, na.rm = TRUE)

x_pad <- max(0.20, diff(lon_rng) * 0.12)
y_pad <- max(0.15, diff(lat_rng) * 0.12)

xlim_use <- c(lon_rng[1] - x_pad, lon_rng[2] + x_pad)
ylim_use <- c(lat_rng[1] - y_pad, lat_rng[2] + y_pad)

map_plot <- ggplot() +
  geom_polygon(
    data = state_map,
    aes(x = long, y = lat, group = group),
    fill = "grey96",
    color = "grey55",
    linewidth = 0.25
  ) +
  geom_point(
    data = plot_sites,
    aes(x = lon, y = lat),
    size = 3
  ) +
  geom_text(
    data = plot_sites,
    aes(x = lon, y = lat, label = paste0(site, ". ", stream)),
    nudge_y = 0.012,
    size = 3.0
  ) +
  coord_fixed(xlim = xlim_use, ylim = ylim_use) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude", title = "ETSA-3 sampling locations")

print(map_plot)

ggsave(
  filename = file.path(base_dir, "ETSA-3_map.png"),
  plot = map_plot,
  width = 7.5,
  height = 5.75,
  dpi = 300
)

# -----------------------------
# 7) IBD plot using in-river distance only
# plotted in RStudio; not saved
# -----------------------------
ibd_plot_river <- ggplot(ibd_df, aes(x = rivdist_km, y = fst)) +
  geom_point(size = 2.7, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "In-river distance (km)",
    y = expression(F[ST]),
    title = "ETSA-3 isolation by distance (in-river)"
  )

print(ibd_plot_river)

# -----------------------------
# 8) checks + save
# -----------------------------
stopifnot(identical(rownames(ETSA_3_fst), ETSA_3_coords$site))
stopifnot(identical(colnames(ETSA_3_fst), ETSA_3_coords$site))
stopifnot(identical(rownames(ETSA_3_rivdists), ETSA_3_coords$site))
stopifnot(identical(colnames(ETSA_3_rivdists), ETSA_3_coords$site))
stopifnot(isTRUE(all.equal(ETSA_3_fst, t(ETSA_3_fst))))
stopifnot(isTRUE(all.equal(ETSA_3_rivdists, t(ETSA_3_rivdists))))

ETSA_3_rivdists = (ETSA_3_rivdists / 1000)


# export
save(
  ETSA_3_fst,
  ETSA_3_coords,
  ETSA_3_rivdists,
  file = file.path(data_dir, "ETSA-3.RData")
)

cat("\nETSA-3 completed successfully.\n")
print(site_lookup)
print(pop_sizes)
