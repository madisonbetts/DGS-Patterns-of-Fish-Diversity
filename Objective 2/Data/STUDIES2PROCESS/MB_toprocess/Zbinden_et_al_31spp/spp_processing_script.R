library(terra)
library(vcfR)
library(dartR)
library(ggplot2)

# -------------------------------------------------------------------
# STEP 0: paths
# -------------------------------------------------------------------
src_base_dir <- path.expand(
  "~/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/STUDIES2PROCESS/MB_toprocess/Zbinden_et_al_31spp/osfstorage-archive/species"
)

out_base_dir <- path.expand(
  "~/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"
)

sites_path <- path.expand(
  "~/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/STUDIES2PROCESS/MB_toprocess/Zbinden_et_al_31spp/osfstorage-archive/spatial/sites/sites.shp"
)

streams_path <- path.expand(
  "~/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/STUDIES2PROCESS/MB_toprocess/Zbinden_et_al_31spp/osfstorage-archive/spatial/streams/streams.shp"
)

huc08_path <- path.expand(
  "~/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/STUDIES2PROCESS/MB_toprocess/Zbinden_et_al_31spp/osfstorage-archive/spatial/basins/HUC08/PROWBDHU08.shp"
)

# -------------------------------------------------------------------
# STEP 1: spatial data in decimal degrees
# -------------------------------------------------------------------
sites_raw <- vect(sites_path)
sites_ll <- project(sites_raw, "EPSG:4326")

streams_raw <- vect(streams_path)
streams_ll <- project(streams_raw, "EPSG:4326")

huc08_raw <- vect(huc08_path)
huc08_ll <- project(huc08_raw, "EPSG:4326")

# site lookup table used to join popmap site IDs to lon/lat
site_lookup <- data.frame(
  site_name = as.character(values(sites_ll)$Field1),
  lon = crds(sites_ll)[, 1],
  lat = crds(sites_ll)[, 2],
  stringsAsFactors = FALSE
)

stopifnot(!anyDuplicated(site_lookup$site_name))

# -------------------------------------------------------------------
# STEP 2: species key
# -------------------------------------------------------------------

# need to fix LEMA

spp_key <- data.frame(
  obj_name = c(
    "LASI_1",
    "LEMA_2", "LEME_1", "MIDO_3", "MISA_1",
    "COCA_1", "COHY_1",
    "FUCA_1", "FUOL_3",
    "CAAN_2", "CAOL_1", "CHER_1", "CYGA_1", "CYWH_1",
    "LUCH_1", "LUPI_1", "LUZO_1", "LYUM_1",
    "NOBO_1", "NONU_1", "NOPE_1", "NOTE_1",
    "PINO_1", "SEAT_3",
    "ETBL_1", "ETCA_5", "ETFL_1", "ETJU_1", "ETSP_1", "ETZO_1",
    "GAAF_1"
  ),
  study_code = c(
    "LASI-1",
    "LEMA-2", "LEME-1", "MIDO-3", "MISA-1",
    "COCA-1", "COHY-1",
    "FUCA-1", "FUOL-3",
    "CAAN-2", "CAOL-1", "CHER-1", "CYGA-1", "CYWH-1",
    "LUCH-1", "LUPI-1", "LUZO-1", "LYUM-1",
    "NOBO-1", "NONU-1", "NOPE-1", "NOTE-1",
    "PINO-1", "SEAT-3",
    "ETBL-1", "ETCA-5", "ETFL-1", "ETJU-1", "ETSP-1", "ETZO-1",
    "GAAF-1"
  ),
  family = c(
    "atherinopsidae",
    "centrarchidae", "centrarchidae", "centrarchidae", "centrarchidae",
    "cottidae", "cottidae",
    "fundulidae", "fundulidae",
    "leuciscidae", "leuciscidae", "leuciscidae", "leuciscidae", "leuciscidae",
    "leuciscidae", "leuciscidae", "leuciscidae", "leuciscidae",
    "leuciscidae", "leuciscidae", "leuciscidae", "leuciscidae",
    "leuciscidae", "leuciscidae",
    "percidae", "percidae", "percidae", "percidae", "percidae", "percidae",
    "poeciliidae"
  ),
  subfolder = c(
    "labsic",
    "lepmac", "lepmeg", "micdol", "micsal",
    "cotcar", "cothyp",
    "funcat", "funoli",
    "camano", "camoli", "chrery", "cypgal", "cypwhi",
    "luxchr", "luxpil", "luxzon", "lytumb",
    "notboo", "notnub", "notper", "nottel",
    "pimnot", "sematr",
    "ethble", "ethcae", "ethfla", "ethjul", "ethspe", "ethzon",
    "gamaff"
  ),
  stringsAsFactors = FALSE
)

# -------------------------------------------------------------------
# STEP 3: helper functions
# -------------------------------------------------------------------

# Read one VCF and convert to genlight
read_gl_from_vcf <- function(vcf_path) {
  vcf_obj <- read.vcfR(vcf_path, verbose = FALSE)
  gl_obj <- vcfR2genlight(vcf_obj)
  gl_obj
}

# Read one popmap.tsv
read_popmap <- function(popmap_path) {
  popmap_df <- read.delim(popmap_path, header = FALSE, stringsAsFactors = FALSE)
  colnames(popmap_df) <- c("individual", "site_name")
  popmap_df
}

# Attach pop slot and lat/lon to genlight
attach_pop_and_latlon <- function(gl_obj, popmap_df, site_lookup_df) {
  if (nrow(popmap_df) != nInd(gl_obj)) {
    stop("popmap row count does not match nInd(genlight)")
  }
  
  if (!all(popmap_df$individual == indNames(gl_obj))) {
    stop("popmap individual order does not match indNames(genlight)")
  }
  
  pop(gl_obj) <- as.factor(popmap_df$site_name)
  
  ind_site_df <- merge(
    popmap_df,
    site_lookup_df,
    by = "site_name",
    all.x = TRUE,
    sort = FALSE
  )
  
  # restore original individual order after merge
  ind_site_df <- ind_site_df[match(popmap_df$individual, ind_site_df$individual), ]
  
  if (anyNA(ind_site_df$lon) || anyNA(ind_site_df$lat)) {
    bad_sites <- unique(ind_site_df$site_name[is.na(ind_site_df$lon) | is.na(ind_site_df$lat)])
    stop("Unmatched site IDs in popmap: ", paste(bad_sites, collapse = ", "))
  }
  
  latlon_df <- ind_site_df[, c("lat", "lon")]
  rownames(latlon_df) <- ind_site_df$individual
  
  gl_obj@other$latlon <- latlon_df
  gl_obj@other$ind_sites <- ind_site_df
  
  gl_obj
}

# Filter individuals with too much missing data
filter_individual_missing <- function(gl_obj, threshold = 0.5) {
  gl_mat <- as.matrix(gl_obj)
  ind_missing_prop <- rowMeans(is.na(gl_mat))
  keep_ind <- ind_missing_prop <= threshold
  gl_obj[keep_ind, ]
}

# Filter loci with too much missing data
filter_locus_missing <- function(gl_obj, threshold = 0.3) {
  gl_mat <- as.matrix(gl_obj)
  loc_missing_prop <- colMeans(is.na(gl_mat))
  keep_loc <- loc_missing_prop <= threshold
  gl_obj[, keep_loc]
}

# Drop populations with < 2 individuals
drop_singleton_pops <- function(gl_obj) {
  pop_counts <- table(pop(gl_obj))
  keep_pops <- names(pop_counts)[pop_counts >= 2]
  gl_obj[as.character(pop(gl_obj)) %in% keep_pops, ]
}

# Subset internal @other objects after filtering
sync_other_slots <- function(gl_obj) {
  keep_ids <- indNames(gl_obj)
  
  if (!is.null(gl_obj@other$ind_sites)) {
    gl_obj@other$ind_sites <- gl_obj@other$ind_sites[
      match(keep_ids, gl_obj@other$ind_sites$individual),
    ]
  }
  
  if (!is.null(gl_obj@other$latlon)) {
    gl_obj@other$latlon <- gl_obj@other$latlon[keep_ids, , drop = FALSE]
  }
  
  gl_obj
}

# Coerce gl.fst.pop output to a matrix
as_fst_matrix <- function(fst_obj) {
  if (is.matrix(fst_obj)) {
    return(fst_obj)
  }
  
  if (is.data.frame(fst_obj)) {
    if (nrow(fst_obj) == ncol(fst_obj)) {
      return(as.matrix(fst_obj))
    }
    
    tmp <- fst_obj
    rownames(tmp) <- tmp[[1]]
    tmp[[1]] <- NULL
    
    if (nrow(tmp) == ncol(tmp)) {
      return(as.matrix(tmp))
    }
  }
  
  stop("Could not coerce gl.fst.pop output to square matrix")
}

# Build final coords df in exact FST order
build_final_coords_df <- function(fst_mat, site_lookup_df) {
  fst_site_names <- rownames(fst_mat)
  
  matched_sites <- site_lookup_df[match(fst_site_names, site_lookup_df$site_name), , drop = FALSE]
  
  if (anyNA(matched_sites$site_name)) {
    stop("Some FST site names could not be matched to the site lookup")
  }
  
  final_coords_df <- data.frame(
    site = seq_len(nrow(matched_sites)),
    lon = matched_sites$lon,
    lat = matched_sites$lat,
    stringsAsFactors = FALSE
  )
  
  rownames(fst_mat) <- final_coords_df$site
  colnames(fst_mat) <- final_coords_df$site
  
  list(fst = fst_mat, coords = final_coords_df)
}

# QC plot
make_qc_plot <- function(ind_sites_df, study_code, streams_v, huc08_v, out_png) {
  set.seed(1)
  
  plot_df <- ind_sites_df
  plot_df$lon_jit <- jitter(plot_df$lon, amount = 0.015)
  plot_df$lat_jit <- jitter(plot_df$lat, amount = 0.015)
  
  site_centroids_df <- unique(plot_df[, c("site_name", "lon", "lat")])
  
  p <- ggplot() +
    geom_spatvector(data = huc08_v, fill = NA, color = "grey40", linewidth = 0.3) +
    geom_spatvector(data = streams_v, color = "steelblue", linewidth = 0.15) +
    geom_point(
      data = plot_df,
      aes(x = lon_jit, y = lat_jit),
      color = "red",
      alpha = 0.6,
      size = 0.8
    ) +
    geom_point(
      data = site_centroids_df,
      aes(x = lon, y = lat),
      color = "black",
      size = 1.5
    ) +
    coord_sf() +
    theme_bw() +
    labs(
      title = paste0(study_code, " individual QC map"),
      x = "Longitude",
      y = "Latitude"
    )
  
  ggsave(out_png, plot = p, width = 7, height = 5, dpi = 300)
}

# -------------------------------------------------------------------
# STEP 4: main workflow loop
# -------------------------------------------------------------------
for (i in seq_len(nrow(spp_key))) {
  
  # ---------------------------------------------------------------
  # STEP 4.1: species metadata
  # ---------------------------------------------------------------
  obj_name_i <- spp_key$obj_name[i]
  study_code_i <- spp_key$study_code[i]
  family_i <- spp_key$family[i]
  subfolder_i <- spp_key$subfolder[i]
  
  cat("\n-----------------------------\n")
  cat("Processing", obj_name_i, "->", study_code_i, "\n")
  cat("-----------------------------\n")
  
  # ---------------------------------------------------------------
  # STEP 4.2: file paths for this species
  # ---------------------------------------------------------------
  vcf_path_i <- file.path(src_base_dir, family_i, subfolder_i, "data", "formatted.vcf")
  popmap_path_i <- file.path(src_base_dir, family_i, subfolder_i, "data", "popmap.tsv")
  out_dir_i <- file.path(out_base_dir, study_code_i, "data")
  
  if (!file.exists(vcf_path_i)) {
    stop("Missing VCF: ", vcf_path_i)
  }
  
  if (!file.exists(popmap_path_i)) {
    stop("Missing popmap: ", popmap_path_i)
  }
  
  dir.create(out_dir_i, recursive = TRUE, showWarnings = FALSE)
  
  # ---------------------------------------------------------------
  # STEP 4.3: read raw genotype and popmap
  # ---------------------------------------------------------------
  step_vcf_gl <- read_gl_from_vcf(vcf_path_i)
  step_popmap <- read_popmap(popmap_path_i)
  
  # debug objects if needed
  # step_vcf_gl
  # head(step_popmap)
  
  # ---------------------------------------------------------------
  # STEP 4.4: attach site IDs and lat/lon
  # ---------------------------------------------------------------
  step_gl_linked <- attach_pop_and_latlon(
    gl_obj = step_vcf_gl,
    popmap_df = step_popmap,
    site_lookup_df = site_lookup
  )
  
  # ---------------------------------------------------------------
  # STEP 4.5: filtering
  # ---------------------------------------------------------------
  # 1) individuals >50% missing removed
  step_gl_filt1 <- filter_individual_missing(step_gl_linked, threshold = 0.5)
  step_gl_filt1 <- sync_other_slots(step_gl_filt1)
  
  # 2) loci >30% missing removed
  step_gl_filt2 <- filter_locus_missing(step_gl_filt1, threshold = 0.3)
  step_gl_filt2 <- sync_other_slots(step_gl_filt2)
  
  # 3) singleton pops removed
  step_gl_filt3 <- drop_singleton_pops(step_gl_filt2)
  step_gl_filt3 <- sync_other_slots(step_gl_filt3)
  
  # 4) monomorphic loci removed
  step_gl_filt4 <- dartR::gl.filter.monomorphs(step_gl_filt3, verbose = 0)
  step_gl_filt4 <- sync_other_slots(step_gl_filt4)
  
  # 5) singleton pops checked again after locus filtering
  step_gl_final <- drop_singleton_pops(step_gl_filt4)
  step_gl_final <- sync_other_slots(step_gl_final)
  
  # quick debug summaries
  cat("nInd after filtering =", nInd(step_gl_final), "\n")
  cat("nPop after filtering =", nPop(step_gl_final), "\n")
  cat("nLoc after filtering =", nLoc(step_gl_final), "\n")
  
  # ---------------------------------------------------------------
  # STEP 4.6: FST
  # ---------------------------------------------------------------
  step_fst_raw <- dartR::gl.fst.pop(step_gl_final)
  step_fst_mat <- as_fst_matrix(step_fst_raw)
  
  step_fst_mat[is.na(step_fst_mat)] <- 0
  step_fst_mat <- (step_fst_mat + t(step_fst_mat)) / 2
  diag(step_fst_mat) <- 0
  step_fst_mat[step_fst_mat < 0] <- 0
  
  # ---------------------------------------------------------------
  # STEP 4.7: final coordinate dataframe and numeric site recode
  # ---------------------------------------------------------------
  step_final_objects <- build_final_coords_df(
    fst_mat = step_fst_mat,
    site_lookup_df = site_lookup
  )
  
  final_fst <- step_final_objects$fst
  final_coords <- step_final_objects$coords
  
  # ---------------------------------------------------------------
  # STEP 4.8: save RData
  # ---------------------------------------------------------------
  obj_prefix_i <- gsub("-", "_", study_code_i)
  fst_name_i <- paste0(obj_prefix_i, "_fst")
  coords_name_i <- paste0(obj_prefix_i, "_coords")
  
  assign(fst_name_i, final_fst, envir = .GlobalEnv)
  assign(coords_name_i, final_coords, envir = .GlobalEnv)
  
  save(
    list = c(fst_name_i, coords_name_i),
    file = file.path(out_dir_i, paste0(study_code_i, ".RData"))
  )
  
  # ---------------------------------------------------------------
  # STEP 4.9: QC plot
  # ---------------------------------------------------------------
  qc_plot_path_i <- file.path(out_dir_i, paste0(study_code_i, "_QC_individuals.png"))
  
  make_qc_plot(
    ind_sites_df = step_gl_final@other$ind_sites,
    study_code = study_code_i,
    streams_v = streams_ll,
    huc08_v = huc08_ll,
    out_png = qc_plot_path_i
  )
  
  cat("Saved:", file.path(out_dir_i, paste0(study_code_i, ".RData")), "\n")
  cat("QC plot:", qc_plot_path_i, "\n")
}