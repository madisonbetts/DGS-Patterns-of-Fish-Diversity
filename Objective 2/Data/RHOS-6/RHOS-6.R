# RHOS-6 genotype-based workflow
# Rhinichthys osculus - Mussmann et al. ddRAD VCF
#
# Goal:
# 1) read the VCF as a genlight
# 2) assign each individual to a population from its sample name
# 3) attach best-available site coordinates to each individual / population
# 4) calculate pairwise population FST with dartR::gl.fst.pop
# 5) force diagonal = 0 and any negative FST = 0
# 6) save RHOS_6_fst and RHOS_6_coords
# 7) plot map + IBD in RStudio only

suppressPackageStartupMessages({
  library(vcfR)
  library(dartR)
  library(adegenet)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(geosphere)
})

# -----------------------------
# 0) paths
# -----------------------------
study_code <- "RHOS-6"
base_dir <- path.expand("~/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/RHOS-6")
out_dir   <- file.path(base_dir, "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

vcf_path <- file.path(base_dir, "spd_dv_input_files", "spd_dv.filt2.vcf")
stopifnot(file.exists(vcf_path))

obj_prefix  <- gsub("-", "_", study_code)
fst_name    <- paste0(obj_prefix, "_fst")
coords_name <- paste0(obj_prefix, "_coords")

# -----------------------------
# 1) read VCF
# -----------------------------
RHOS_6 <- vcfR2genlight(read.vcfR(vcf_path))

# -----------------------------
# 2) site/pop lookup
# Best-available approximate centroids for the named localities represented in Mussmann et al.
# Codes in sample names map to localities as follows:
#   AMA / damr = Upper Amargosa River at Beatty
#   AMC        = Amargosa Canyon
#   ASH        = Bradford Spring
#   ASR        = Rogers Spring
#   LVD / dwhs = Long Valley / Whitmore Hot Springs
#   PRC        = Poore Creek (Walker basin)
#   RFO        = Roberts Field
#   RUP        = Russi Pond
#   dhar       = Harris Ranch (Benton Valley)
#   dorb       = Owens River - Bishop
#   rhat       = Rogue River, MI (R. atratulus outgroup)
# -----------------------------
site_key <- tibble::tribble(
  ~pop_code, ~site_name,                            ~lat,      ~lon,
  "AMA",    "Upper Amargosa River at Beatty",    36.91398, -116.75551,
  "AMC",    "Amargosa Canyon",                   35.80020, -116.19520,
  "ASH",    "Bradford Spring",                   36.40214, -116.30236,
  "ASR",    "Rogers Spring",                     36.47920, -116.32530,
  "LVD",    "Long Valley / Whitmore Hot Springs",37.63333, -118.81667,
  "PRC",    "Poore Creek",                       38.58600, -119.63400,
  "RFO",    "Roberts Field",                     36.94600, -116.74400,
  "RUP",    "Russi Pond",                        37.36200, -118.40600,
  "HAR",    "Harris Ranch (Benton Valley)",      37.78700, -118.55300,
  "ORB",    "Owens River - Bishop",              37.37300, -118.39500,
  "RHAT",   "Rogue River, Michigan",             42.95600,  -85.60900
)

assign_pop <- function(x) {
  case_when(
    str_detect(x, "AMA|damr") ~ "AMA",
    str_detect(x, "AMC")      ~ "AMC",
    str_detect(x, "ASH")      ~ "ASH",
    str_detect(x, "ASR")      ~ "ASR",
    str_detect(x, "LVD|dwhs") ~ "LVD",
    str_detect(x, "PRC")      ~ "PRC",
    str_detect(x, "RFO")      ~ "RFO",
    str_detect(x, "RUP")      ~ "RUP",
    str_detect(x, "dhar")     ~ "HAR",
    str_detect(x, "dorb")     ~ "ORB",
    str_detect(x, "rhat")     ~ "RHAT",
    TRUE ~ NA_character_
  )
}

ind_meta <- data.frame(ind = indNames(RHOS_6), stringsAsFactors = FALSE) %>%
  mutate(pop_code = assign_pop(ind)) %>%
  left_join(site_key, by = "pop_code")

if (any(is.na(ind_meta$pop_code))) {
  stop("Some individuals were not assigned to a population: ",
       paste(ind_meta$ind[is.na(ind_meta$pop_code)], collapse = ", "))
}

# attach population factor to genlight
pop(RHOS_6) <- factor(ind_meta$pop_code, levels = unique(site_key$pop_code))

# -----------------------------
# 3) population/site coords object for project output
# Default: keep all RHOS populations except RHAT for FST output.
# If you want Death Valley ingroup only, also drop PRC below.
# -----------------------------
keep_for_fst <- ind_meta$pop_code != "RHAT"
# keep_for_fst <- !(ind_meta$pop_code %in% c("RHAT", "PRC"))  # optional stricter ingroup-only version

RHOS_6_fst_gl <- RHOS_6[keep_for_fst]
RHOS_6_fst_meta <- ind_meta[keep_for_fst, ]
pop(RHOS_6_fst_gl) <- droplevels(factor(RHOS_6_fst_meta$pop_code))

RHOS_6_coords <- RHOS_6_fst_meta %>%
  distinct(pop_code, site_name, lat, lon) %>%
  mutate(site_id = seq_len(n())) %>%
  select(site_id, pop_code, site_name, lat, lon)

# enforce project order to match the FST matrix
pop_levels <- RHOS_6_coords$pop_code
pop(RHOS_6_fst_gl) <- factor(pop(RHOS_6_fst_gl), levels = pop_levels)

# -----------------------------
# 4) pairwise FST with gl.fst.pop
# -----------------------------
fst_raw <- dartR::gl.fst.pop(RHOS_6_fst_gl, verbose = 0)

extract_fst_matrix <- function(x, target_levels) {
  # most likely: matrix / data.frame already
  if (is.matrix(x) || is.data.frame(x)) {
    m <- as.matrix(x)
  } else if (is.list(x)) {
    nms <- names(x)
    pick <- NULL
    for (cand in c("Fsts", "fst", "Fst", "fst.mat", "Fst.mat", "pairwise.fst", "PwFST")) {
      if (cand %in% nms) {
        pick <- x[[cand]]
        break
      }
    }
    if (is.null(pick)) {
      # try first matrix-like element
      idx <- which(vapply(x, function(z) is.matrix(z) || is.data.frame(z), logical(1)))[1]
      if (length(idx) == 0 || is.na(idx)) {
        stop("Could not find an FST matrix inside gl.fst.pop output. Inspect fst_raw with str().")
      }
      pick <- x[[idx]]
    }
    m <- as.matrix(pick)
  } else {
    stop("Unexpected gl.fst.pop output class: ", paste(class(x), collapse = ", "))
  }

  # if row/col names are missing but dimensions match, impose target levels
  if (is.null(rownames(m)) && nrow(m) == length(target_levels)) rownames(m) <- target_levels
  if (is.null(colnames(m)) && ncol(m) == length(target_levels)) colnames(m) <- target_levels

  # reorder if names are present
  if (!is.null(rownames(m)) && !is.null(colnames(m))) {
    if (all(target_levels %in% rownames(m)) && all(target_levels %in% colnames(m))) {
      m <- m[target_levels, target_levels, drop = FALSE]
    }
  }

  return(m)
}

fst_mat <- extract_fst_matrix(fst_raw, pop_levels)
mode(fst_mat) <- "numeric"

# project rules
if (nrow(fst_mat) != nrow(RHOS_6_coords)) {
  stop("FST matrix size does not match RHOS_6_coords rows.")
}

fst_mat[is.na(fst_mat)] <- 0
fst_mat <- (fst_mat + t(fst_mat)) / 2
fst_mat[fst_mat < 0] <- 0

diag(fst_mat) <- 0
rownames(fst_mat) <- RHOS_6_coords$site_id
colnames(fst_mat) <- RHOS_6_coords$site_id

RHOS_6_fst <- fst_mat
RHOS_6_coords <- RHOS_6_coords %>% select(site_id, lat, lon)

# -----------------------------
# 5) QC plots in RStudio only
# -----------------------------
# map
pts_sf <- st_as_sf(
  RHOS_6_fst_meta %>%
    distinct(pop_code, site_name, lat, lon),
  coords = c("lon", "lat"), crs = 4326
)

na <- rnaturalearth::ne_countries(country = c("United States of America", "Canada"), returnclass = "sf")

bb <- st_bbox(pts_sf)
x_pad <- diff(bb[c("xmin", "xmax")]) * 0.25
if (x_pad == 0) x_pad <- 0.5
y_pad <- diff(bb[c("ymin", "ymax")]) * 0.25
if (y_pad == 0) y_pad <- 0.5

plot_map <- ggplot() +
  geom_sf(data = na, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(data = pts_sf, size = 2) +
  ggrepel::geom_text_repel(
    data = cbind(st_drop_geometry(pts_sf), st_coordinates(pts_sf)),
    aes(X, Y, label = pop_code),
    size = 3
  ) +
  coord_sf(
    xlim = c(bb["xmin"] - x_pad, bb["xmax"] + x_pad),
    ylim = c(bb["ymin"] - y_pad, bb["ymax"] + y_pad),
    expand = FALSE
  ) +
  theme_bw() +
  labs(x = NULL, y = NULL, title = study_code)
print(plot_map)

# IBD using straight-line geographic distance as a QC plot
geo_mat <- geosphere::distm(
  x = RHOS_6_coords[, c("lon", "lat")],
  fun = geosphere::distHaversine
) / 1000
rownames(geo_mat) <- RHOS_6_coords$site_id
colnames(geo_mat) <- RHOS_6_coords$site_id

lower_idx <- lower.tri(RHOS_6_fst)
ibd_df <- data.frame(
  dist_km = geo_mat[lower_idx],
  fst = RHOS_6_fst[lower_idx]
)

plot_ibd <- ggplot(ibd_df, aes(dist_km, fst)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(x = "Geographic distance (km)", y = "Pairwise FST", title = paste0(study_code, " IBD QC"))
print(plot_ibd)

# -----------------------------
# 6) save outputs
# -----------------------------
assign(fst_name, RHOS_6_fst, envir = .GlobalEnv)
assign(coords_name, RHOS_6_coords, envir = .GlobalEnv)

colnames(RHOS_6_coords)[1] <- "site"

save(
  list = c(fst_name, coords_name),
  file = file.path(out_dir, paste0(study_code, ".RData"))
)

# -----------------------------
# 7) optional inspection
# -----------------------------
print(table(ind_meta$pop_code, useNA = "ifany"))
print(RHOS_6_coords)
print(round(RHOS_6_fst, 4))
