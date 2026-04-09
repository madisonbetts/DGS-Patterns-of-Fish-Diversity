
############################
## MAVE-1
## Gunn et al. 2023
## Hard-filtered diagnostic SNP workflow from raw genotype + metadata files
##
## Focal lineage: Neosho Bass (Micropterus velox)
##
## This script:
##   1) reads genotype_data.xlsx and metadata.xlsx
##   2) applies the repository filtering logic:
##        - drop loci with >20% missing
##        - drop samples with >20% missing
##        - remove likely duplicate focal-lineage samples (>95% identity)
##   3) reproduces a hard-filtered lineage dataset by removing the exact
##      number of hybrids reported per river in Gunn et al. (2023),
##      ranking candidate hybrids by diagnostic-SNP ancestry scores
##   4) computes retained-sample river centroids from exact metadata coords
##   5) calculates pairwise Weir & Cockerham FST among retained rivers
##   6) sets negative FST to 0, plots the sites and IBD, and saves RData
############################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(hierfstat)
  library(geosphere)
  library(ggplot2)
  library(maps)
  library(igraph)
})

# -----------------------------
# paths
# -----------------------------
out_dir  <- file.path("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data", "MAVE-1")
data_dir <- file.path(out_dir, "doi_10_5061_dryad_k6djh9wbs__v20240206/data")

geno_xlsx <- file.path(out_dir, "doi_10_5061_dryad_k6djh9wbs__v20240206/genotype_data.xlsx")
meta_xlsx <- file.path(out_dir, "doi_10_5061_dryad_k6djh9wbs__v20240206/metadata.xlsx")

if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

# -----------------------------
# focal settings
# -----------------------------
dataset_code   <- "MAVE-1"
focal_taxon    <- "NB"
river_order    <- c("Honey_Creek", "Spavinaw_Creek", "Baron_Fork", "Caney_Creek", "Lee_Creek")
site_name_map  <- c("Honey_Creek" = "Honey Creek", "Spavinaw_Creek" = "Spavinaw Creek", "Baron_Fork" = "Baron Fork", "Caney_Creek" = "Caney Creek", "Lee_Creek" = "Lee Creek")
drop_stage1_n  <- c("Honey_Creek" = 1, "Lee_Creek" = 1)
drop_stage2_n  <- c("Honey_Creek" = 1, "Spavinaw_Creek" = 3, "Baron_Fork" = 1)
drop_stage3_n  <- c()

# -----------------------------
# helpers
# -----------------------------
to_na_unknown <- function(df) {
  df[] <- lapply(df, function(x) {
    x <- as.character(x)
    x[x %in% c("unknown", "", "NA", "NaN", "<NA>")] <- NA_character_
    x
  })
  df
}

prop_missing_vec <- function(x) {
  mean(is.na(x))
}

sample_missing_prop <- function(df) {
  apply(df, 1, function(x) mean(is.na(x)))
}

pair_identity <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  n_ok <- sum(ok)
  if (n_ok == 0) {
    return(c(n_overlap = 0, p_ident = NA_real_))
  }
  c(
    n_overlap = n_ok,
    p_ident   = mean(x[ok] == y[ok])
  )
}

drop_duplicate_components <- function(meta_df, geno_df, focal_taxon,
                                      min_overlap = 150, ident_threshold = 0.95) {
  focal_idx <- which(meta_df$putative_taxon == focal_taxon)
  if (length(focal_idx) < 2) {
    return(integer(0))
  }

  edges <- list()

  for (ii in seq_len(length(focal_idx) - 1)) {
    i <- focal_idx[ii]
    xi <- unlist(geno_df[i, ], use.names = FALSE)

    for (jj in (ii + 1):length(focal_idx)) {
      j <- focal_idx[jj]
      yj <- unlist(geno_df[j, ], use.names = FALSE)

      pid <- pair_identity(xi, yj)
      same_taxon <- identical(meta_df$putative_taxon[i], meta_df$putative_taxon[j])

      if (same_taxon &&
          !is.na(pid["p_ident"]) &&
          pid["n_overlap"] >= min_overlap &&
          pid["p_ident"] > ident_threshold) {
        edges[[length(edges) + 1]] <- c(i, j)
      }
    }
  }

  if (length(edges) == 0) {
    return(integer(0))
  }

  edges_mat <- do.call(rbind, edges)
  g <- graph_from_edgelist(edges_mat, directed = FALSE)
  comps <- components(g)$membership
  drop_idx <- integer(0)

  for (cc in sort(unique(comps))) {
    members <- as.integer(names(comps)[comps == cc])
    if (length(members) < 2) next

    miss_prop <- sample_missing_prop(geno_df[members, , drop = FALSE])
    keep_one  <- members[which.min(miss_prop)]
    drop_idx  <- c(drop_idx, setdiff(members, keep_one))
  }

  sort(unique(drop_idx))
}

major_allele <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)

  alleles <- unlist(strsplit(x, ":", fixed = TRUE))
  tab <- sort(table(alleles), decreasing = TRUE)
  names(tab)[1]
}

ancestry_score <- function(meta_df, geno_df, loci, parent1_idx, parent2_idx) {
  loci <- intersect(loci, names(geno_df))
  score_list <- vector("list", length(loci))
  keep_loci   <- character(0)

  for (ll in seq_along(loci)) {
    loc <- loci[ll]
    p1_allele <- major_allele(geno_df[[loc]][parent1_idx])
    p2_allele <- major_allele(geno_df[[loc]][parent2_idx])

    if (is.na(p1_allele) || is.na(p2_allele) || identical(p1_allele, p2_allele)) {
      next
    }

    keep_loci <- c(keep_loci, loc)
    score_list[[ll]] <- vapply(geno_df[[loc]], function(gt) {
      if (is.na(gt)) return(NA_real_)
      parts <- strsplit(gt, ":", fixed = TRUE)[[1]]
      mean(parts == p2_allele)
    }, numeric(1))
  }

  keep_cols <- keep_loci
  if (length(keep_cols) == 0) {
    stop("No informative loci were retained for ancestry scoring.")
  }

  score_mat <- do.call(cbind, score_list[match(keep_cols, loci)])
  if (!is.matrix(score_mat)) score_mat <- matrix(score_mat, ncol = 1)

  rowMeans(score_mat, na.rm = TRUE)
}

fit_internal_centroids <- function(meta_df, geno_df, loci, species_levels = c("NB", "LRB", "OB")) {
  loci <- intersect(loci, names(geno_df))
  geno_chr <- geno_df[, loci, drop = FALSE]
  geno_chr[is.na(geno_chr)] <- "missing"

  mm <- model.matrix(~ . - 1, data = geno_chr)
  centroids <- lapply(species_levels, function(sp) {
    idx <- which(meta_df$putative_taxon == sp)
    if (length(idx) == 0) stop(paste("No rows for species", sp))
    colMeans(mm[idx, , drop = FALSE])
  })
  names(centroids) <- species_levels

  out <- matrix(NA_real_, nrow = nrow(mm), ncol = length(species_levels))
  colnames(out) <- species_levels

  for (jj in seq_along(species_levels)) {
    cc <- centroids[[jj]]
    diff_sq <- sweep(mm, 2, cc, FUN = "-")^2
    out[, jj] <- sqrt(rowSums(diff_sq))
  }

  as.data.frame(out)
}

choose_hybrids_by_count <- function(meta_df, score_col, focal_taxon, named_counts) {
  out <- integer(0)
  if (length(named_counts) == 0) return(out)

  for (rv in names(named_counts)) {
    n_drop <- unname(named_counts[[rv]])
    idx <- which(meta_df$putative_taxon == focal_taxon & meta_df$river == rv)
    if (length(idx) == 0 || n_drop == 0) next
    ord <- idx[order(meta_df[[score_col]][idx], decreasing = FALSE)]
    out <- c(out, head(ord, n_drop))
  }

  sort(unique(out))
}

choose_stage3_internal <- function(meta_df, dist_df, focal_taxon, named_counts) {
  out <- integer(0)
  if (length(named_counts) == 0) return(out)

  focal_col <- focal_taxon
  other_cols <- setdiff(colnames(dist_df), focal_col)

  # low margin to the nearest alternative lineage = most admixture-like
  margin <- vapply(seq_len(nrow(dist_df)), function(i) {
    min(dist_df[i, other_cols]) - dist_df[i, focal_col]
  }, numeric(1))

  for (rv in names(named_counts)) {
    n_drop <- unname(named_counts[[rv]])
    idx <- which(meta_df$putative_taxon == focal_taxon & meta_df$river == rv)
    if (length(idx) == 0 || n_drop == 0) next
    ord <- idx[order(margin[idx], decreasing = FALSE)]
    out <- c(out, head(ord, n_drop))
  }

  sort(unique(out))
}

make_hf_locus <- function(gt_vec) {
  gt_vec <- as.character(gt_vec)
  gt_vec[gt_vec %in% c(NA, "NA", "<NA>")] <- NA_character_

  non_missing <- gt_vec[!is.na(gt_vec)]
  alleles <- sort(unique(unlist(strsplit(non_missing, ":", fixed = TRUE))))

  allele_map <- setNames(seq_along(alleles), alleles)

  out <- vapply(gt_vec, function(gt) {
    if (is.na(gt)) return(NA_integer_)
    pp <- sort(strsplit(gt, ":", fixed = TRUE)[[1]])
    a1 <- allele_map[[pp[1]]]
    a2 <- allele_map[[pp[2]]]
    as.integer(sprintf("%03d%03d", a1, a2))
  }, integer(1))

  out
}

pretty_site_name <- function(river_code) {
  unname(site_name_map[river_code])
}

# -----------------------------
# read raw files
# -----------------------------
meta_raw <- readxl::read_excel(meta_xlsx) %>% as.data.frame(stringsAsFactors = FALSE)
geno_raw <- readxl::read_excel(geno_xlsx) %>% as.data.frame(stringsAsFactors = FALSE)

stopifnot(identical(meta_raw$sample_id, geno_raw$sample_id))

geno_calls <- geno_raw[, -1, drop = FALSE]
geno_calls <- to_na_unknown(geno_calls)

# -----------------------------
# repository-style filtering
# -----------------------------
locus_miss <- vapply(geno_calls, prop_missing_vec, numeric(1))
keep_loci  <- names(locus_miss)[locus_miss <= 0.20]

geno_locus <- geno_calls[, keep_loci, drop = FALSE]
sample_miss <- sample_missing_prop(geno_locus)
keep_samples <- sample_miss <= 0.20

meta_filt <- meta_raw[keep_samples, , drop = FALSE]
geno_filt <- geno_locus[keep_samples, , drop = FALSE]

dup_drop <- drop_duplicate_components(meta_filt, geno_filt, focal_taxon = focal_taxon)

if (length(dup_drop) > 0) {
  keep_dup <- setdiff(seq_len(nrow(meta_filt)), dup_drop)
  meta_filt <- meta_filt[keep_dup, , drop = FALSE]
  geno_filt <- geno_filt[keep_dup, , drop = FALSE]
}

rownames(meta_filt) <- NULL
rownames(geno_filt) <- NULL

# -----------------------------
# stage 1: SPB vs SMB-C hard filter
# paper threshold: remove reported per-river hybrid counts
# -----------------------------
spb_loci <- grep("^SM_SP_", names(geno_filt), value = TRUE)

score_smsp <- ancestry_score(
  meta_df    = meta_filt,
  geno_df    = geno_filt,
  loci       = spb_loci,
  parent1_idx = which(meta_filt$putative_taxon == "SPB"),
  parent2_idx = which(meta_filt$putative_taxon %in% c("SMB", "NB", "LRB", "OB"))
)

meta_filt$score_smsp <- score_smsp

drop1_idx <- choose_hybrids_by_count(
  meta_df      = meta_filt,
  score_col    = "score_smsp",
  focal_taxon  = focal_taxon,
  named_counts = drop_stage1_n
)

if (length(drop1_idx) > 0) {
  keep1 <- setdiff(seq_len(nrow(meta_filt)), drop1_idx)
  meta_s1 <- meta_filt[keep1, , drop = FALSE]
  geno_s1 <- geno_filt[keep1, , drop = FALSE]
} else {
  meta_s1 <- meta_filt
  geno_s1 <- geno_filt
}

rownames(meta_s1) <- NULL
rownames(geno_s1) <- NULL

# -----------------------------
# stage 2: SMB vs CIH hard filter
# -----------------------------
smb_loci <- grep("^N_SM_", names(geno_s1), value = TRUE)

score_nsm <- ancestry_score(
  meta_df     = meta_s1,
  geno_df     = geno_s1,
  loci        = smb_loci,
  parent1_idx = which(meta_s1$putative_taxon == "SMB"),
  parent2_idx = which(meta_s1$putative_taxon %in% c("NB", "LRB", "OB"))
)

meta_s1$score_nsm <- score_nsm

drop2_idx <- choose_hybrids_by_count(
  meta_df      = meta_s1,
  score_col    = "score_nsm",
  focal_taxon  = focal_taxon,
  named_counts = drop_stage2_n
)

if (length(drop2_idx) > 0) {
  keep2 <- setdiff(seq_len(nrow(meta_s1)), drop2_idx)
  meta_s2 <- meta_s1[keep2, , drop = FALSE]
  geno_s2 <- geno_s1[keep2, , drop = FALSE]
} else {
  meta_s2 <- meta_s1
  geno_s2 <- geno_s1
}

rownames(meta_s2) <- NULL
rownames(geno_s2) <- NULL

# -----------------------------
# stage 3: within-CIH hard filter, if needed
# only LRB has a reported CIH-stage hybrid to remove
# -----------------------------
cih_loci <- grep("^(NEO_|OULR_|OUOU_)", names(geno_s2), value = TRUE)

if (length(drop_stage3_n) > 0) {
  idx_internal <- which(meta_s2$putative_taxon %in% c("NB", "LRB", "OB"))
  meta_internal <- meta_s2[idx_internal, , drop = FALSE]
  geno_internal <- geno_s2[idx_internal, cih_loci, drop = FALSE]

  internal_dist <- fit_internal_centroids(
    meta_df = meta_internal,
    geno_df = geno_internal,
    loci    = cih_loci,
    species_levels = c("NB", "LRB", "OB")
  )

  drop3_local <- choose_stage3_internal(
    meta_df      = meta_internal,
    dist_df      = internal_dist,
    focal_taxon  = focal_taxon,
    named_counts = drop_stage3_n
  )

  drop3_idx <- which(meta_s2$sample_id %in% meta_internal$sample_id[drop3_local])
} else {
  drop3_idx <- integer(0)
}

if (length(drop3_idx) > 0) {
  keep3 <- setdiff(seq_len(nrow(meta_s2)), drop3_idx)
  meta_final_all <- meta_s2[keep3, , drop = FALSE]
  geno_final_all <- geno_s2[keep3, , drop = FALSE]
} else {
  meta_final_all <- meta_s2
  geno_final_all <- geno_s2
}

rownames(meta_final_all) <- NULL
rownames(geno_final_all) <- NULL

# -----------------------------
# keep focal lineage only
# -----------------------------
focal_idx <- which(meta_final_all$putative_taxon == focal_taxon)
meta_final <- meta_final_all[focal_idx, , drop = FALSE]
geno_final <- geno_final_all[focal_idx, cih_loci, drop = FALSE]

# enforce paper river order
meta_final$river <- factor(meta_final$river, levels = river_order)
ord_final <- order(meta_final$river, meta_final$sample_id)
meta_final <- meta_final[ord_final, , drop = FALSE]
geno_final <- geno_final[ord_final, , drop = FALSE]

# -----------------------------
# site lookup + retained-sample centroids
# exact coordinates from metadata, then averaged within river
# -----------------------------
site_lookup <- meta_final %>%
  mutate(site_name = pretty_site_name(as.character(river))) %>%
  group_by(river, site_name) %>%
  summarise(
    n_ind = n(),
    lat   = mean(latitude, na.rm = TRUE),
    lon   = mean(longitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(factor(river, levels = river_order)) %>%
  mutate(site_id = seq_len(n())) %>%
  select(site_id, river, site_name, n_ind, lat, lon)

# -----------------------------
# convert to hierfstat format
# -----------------------------
meta_final$site_id <- site_lookup$site_id[match(meta_final$river, site_lookup$river)]

hf_df <- data.frame(pop = meta_final$site_id, stringsAsFactors = FALSE)

for (loc in names(geno_final)) {
  hf_df[[loc]] <- make_hf_locus(geno_final[[loc]])
}

keep_hf_loci <- colSums(!is.na(hf_df[, -1, drop = FALSE])) > 0
hf_df <- hf_df[, c(TRUE, keep_hf_loci), drop = FALSE]

# -----------------------------
# pairwise FST
# -----------------------------
fst_mat <- hierfstat::pairwise.WCfst(hf_df)
fst_mat <- as.matrix(fst_mat)

fst_mat[lower.tri(fst_mat)] <- t(fst_mat)[lower.tri(fst_mat)]
diag(fst_mat) <- 0
fst_mat[fst_mat < 0] <- 0

rownames(fst_mat) <- site_lookup$site_id
colnames(fst_mat) <- site_lookup$site_id

coords_df <- site_lookup %>%
  select(site_id, lat, lon)

# rename objects to workflow-standard names
MAVE_1_fst <- fst_mat
MAVE_1_coords <- coords_df

# -----------------------------
# IBD dataframe
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(coords_df[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- coords_df$site_id
colnames(geo_dist_km) <- coords_df$site_id

ibd_df <- data.frame(
  site1   = rownames(fst_mat)[row(fst_mat)[upper.tri(fst_mat)]],
  site2   = colnames(fst_mat)[col(fst_mat)[upper.tri(fst_mat)]],
  fst     = fst_mat[upper.tri(fst_mat)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  stringsAsFactors = FALSE
)

# -----------------------------
# plots
# -----------------------------
usa <- map_data("state")
canada <- map_data("world", region = "Canada")

p_map <- ggplot() +
  geom_polygon(
    data = usa,
    aes(x = long, y = lat, group = group),
    fill = "gray95",
    color = "gray70",
    linewidth = 0.2
  ) +
  geom_polygon(
    data = canada,
    aes(x = long, y = lat, group = group),
    fill = "gray95",
    color = "gray70",
    linewidth = 0.2
  ) +
  geom_point(
    data = site_lookup,
    aes(x = lon, y = lat),
    size = 2.8
  ) +
  geom_text(
    data = site_lookup,
    aes(x = lon, y = lat, label = site_name),
    size = 3.0,
    nudge_y = 0.08
  ) +
  coord_quickmap(
    xlim = range(site_lookup$lon, na.rm = TRUE) + c(-4, 4),
    ylim = range(site_lookup$lat, na.rm = TRUE) + c(-2.5, 2.5)
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = paste0(dataset_code, " sampling sites")
  )

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.2, alpha = 0.75) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Straight-line distance (km)",
    y = expression(F[ST]),
    title = paste0(dataset_code, " isolation by distance")
  )

print(p_map)
print(p_ibd)

# -----------------------------
# quick checks
# -----------------------------
stopifnot(identical(rownames(fst_mat), as.character(coords_df$site_id)))
stopifnot(identical(colnames(fst_mat), as.character(coords_df$site_id)))
stopifnot(isTRUE(all.equal(fst_mat, t(fst_mat))))

cat("\nDataset:", dataset_code, "\n")
cat("Focal taxon:", focal_taxon, "\n")
cat("Retained sites:", nrow(site_lookup), "\n")
cat("Retained individuals:", nrow(meta_final), "\n")
cat("Retained loci:", ncol(hf_df) - 1, "\n\n")

print(site_lookup)

# -----------------------------
# save RData
# -----------------------------
save(
  MAVE_1_fst,
  MAVE_1_coords,
  file = file.path("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/MAVE-1/data/MAVE-1.RData")
)
