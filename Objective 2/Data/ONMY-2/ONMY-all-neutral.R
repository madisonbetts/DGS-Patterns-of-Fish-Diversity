
############################
## ONMY-all steelhead
## Hand et al. 2016 (Molecular Ecology)
## Neutral SNPs only
##
## This script:
##   1) reads all five neutral genotype sheets from SH_Raw_Genetic.xlsx
##   2) converts the two-column-per-locus GenAlEx-style layout into
##      hierfstat-ready diploid genotype codes
##   3) combines all individuals across all 79 populations
##   4) calculates pairwise Weir & Cockerham FST among all populations
##   5) forces symmetry and sets negative FST to 0
##   6) renames rows/cols with numeric site IDs
##   7) builds one site dataframe containing all populations
##
## Notes:
## - The five neutral sheets are:
##     Lo_Col_genos_neutral
##     Salmon_MPG_genos_neutral
##     Clearwater_MPG_genos_neutral
##     Cascades_MPG_genos_neutral
##     JohnDay_genos_neutral
## - Site names are mapped to pop codes in the same within-sheet order
##   as the supplement Table S1 / workbook population order.
## - This script does NOT assign coordinates because the workbook does not
##   contain site XYs. Add them later if you want mapping / IBD plots.
############################

suppressPackageStartupMessages({
  library(readxl)
  library(hierfstat)
})

# -----------------------------
# paths
# -----------------------------
xlsx_file <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ONMY-2/doi_10_5061_dryad_pr065__v20151211/DataforDryad/SH_Raw_Genetic.xlsx"

# -----------------------------
# sheets to combine
# -----------------------------
neutral_sheets <- c(
  "Lo_Col_genos_neutral",
  "Salmon_MPG_genos_neutral",
  "Clearwater_MPG_genos_neutral",
  "Cascades_MPG_genos_neutral",
  "JohnDay_genos_neutral"
)

# -----------------------------
# helper functions
# -----------------------------
fmt_allele <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "<NA>", "NULL")] <- NA_character_
  x
}

make_hf_genotype <- function(a1, a2) {
  a1 <- fmt_allele(a1)
  a2 <- fmt_allele(a2)

  bad1 <- !is.na(a1) & !grepl("^[0-9]+$", a1)
  bad2 <- !is.na(a2) & !grepl("^[0-9]+$", a2)

  if (any(bad1 | bad2)) {
    stop("Non-numeric allele codes detected in the raw genotype workbook.")
  }

  ifelse(
    is.na(a1) | is.na(a2),
    NA_integer_,
    as.integer(sprintf("%03d%03d", as.integer(a1), as.integer(a2)))
  )
}

sanitize_locus_names <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  make.unique(x, sep = "_")
}

read_neutral_sheet <- function(xlsx_file, sheet_name) {
  raw_xlsx <- readxl::read_excel(
    path = xlsx_file,
    sheet = sheet_name,
    col_names = FALSE
  )
  raw_xlsx <- as.data.frame(raw_xlsx, stringsAsFactors = FALSE)

  header_row <- raw_xlsx[3, ]
  ind_dat <- raw_xlsx[-c(1, 2, 3), , drop = FALSE]
  rownames(ind_dat) <- NULL
  names(ind_dat)[1:2] <- c("sample_id", "pop_code")

  locus_names <- as.character(unlist(header_row[3:ncol(raw_xlsx)]))
  locus_names <- trimws(locus_names)

  if (length(locus_names) %% 2 != 0) {
    stop(paste0("Expected an even number of genotype columns in sheet ", sheet_name, "."))
  }

  locus_pairs <- split(seq_along(locus_names), ceiling(seq_along(locus_names) / 2))

  # verify each locus appears as a clean allele pair
  pair_ok <- vapply(
    locus_pairs,
    function(idx) length(unique(locus_names[idx])) == 1,
    logical(1)
  )

  if (!all(pair_ok)) {
    stop(paste0("Not all genotype columns occur as clean allele pairs in sheet ", sheet_name, "."))
  }

  paired_names <- sanitize_locus_names(vapply(locus_pairs, function(idx) locus_names[idx[1]], character(1)))

  pop_codes <- unique(ind_dat$pop_code)
  pop_n <- as.integer(table(factor(ind_dat$pop_code, levels = pop_codes)))

  list(
    sheet_name = sheet_name,
    ind_dat = ind_dat,
    pop_codes = pop_codes,
    pop_n = pop_n,
    locus_pairs = locus_pairs,
    locus_names = paired_names
  )
}

sheet_to_hf <- function(sheet_obj, global_site_lookup) {
  ind_dat <- sheet_obj$ind_dat
  locus_pairs <- sheet_obj$locus_pairs
  locus_names <- sheet_obj$locus_names

  hf <- data.frame(
    pop = match(ind_dat$pop_code, global_site_lookup$pop_code),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(locus_pairs)) {
    cols <- locus_pairs[[i]] + 2L
    hf[[locus_names[i]]] <- make_hf_genotype(ind_dat[[cols[1]]], ind_dat[[cols[2]]])
  }

  hf
}

# -----------------------------
# site metadata in within-sheet order
# -----------------------------
site_meta <- list(
  Lo_Col_genos_neutral = data.frame(
    site_name = c(
      "Cowlitz R.", "Kalama R.", "E. F. Lewis R.", "N. F. Lewis R.",
      "Clackamas R.", "N. F. Eagle Cr.", "Still Cr.", "E. Fork Hood R.",
      "W. Fork Hood R.", "Coweeman R."
    ),
    metapop = "Lower Columbia",
    stringsAsFactors = FALSE
  ),
  Salmon_MPG_genos_neutral = data.frame(
    site_name = c(
      "Hazard Cr.", "Sawtooth Weir", "W. F. Yankee Fork R.", "Morgan Cr.",
      "Pahsimeroi Weir", "Hayden Cr.", "N. F. Salmon R.", "Marsh Cr.",
      "Sulphur Cr.", "Rapid R.", "Pistol Cr.", "Camas Cr.", "Loon Cr.",
      "Upper Big Cr.", "Chamberlain Cr.", "Bargamin Cr.",
      "East Fork S. F. Salmon R.", "Secesh R.", "Lick Cr.",
      "Stolle Meadows", "Boulder Cr.", "Slate Cr.", "Whitebird Cr."
    ),
    metapop = "Salmon",
    stringsAsFactors = FALSE
  ),
  Clearwater_MPG_genos_neutral = data.frame(
    site_name = c(
      "Mission Cr.", "Colt/ Lochsa R.", "Storm/ Lochsa R.", "Crooked Fork/Lochsa R.",
      "Lake/ Lochsa R.", "Fish/ Lochsa R.", "Canyon/ Lochsa R.",
      "Little Clearwater/ Selway R.", "Whitecap/ Selway R.", "Bear/ Selway R.",
      "N. F. Moose/ Selway R.", "Three Links/ Selway R.", "Gedney/ Selway R.",
      "OHara/ Selway R.", "Clear Cr.", "Crooked R.", "Tenmile Cr.",
      "Johns Cr.", "W. F. Potlatch R.", "E. F. Potlatch R.",
      "Big Bear Cr.", "Little Bear Cr."
    ),
    metapop = "Clearwater",
    stringsAsFactors = FALSE
  ),
  Cascades_MPG_genos_neutral = data.frame(
    site_name = c(
      "Shitike Cr.", "Buckhollow Cr.", "Trout Cr.", "Fifteen Cr.",
      "Lower Little Klickitat R.", "Lower Summit Cr.", "Upper Trout Cr.",
      "Deadcanyon Cr.", "Lower White Cr.", "Snyder Cr.", "Surveyor Cr.",
      "Swale Cr.", "Rock Cr.", "Squaw Cr."
    ),
    metapop = "Eastern Cascades",
    stringsAsFactors = FALSE
  ),
  JohnDay_genos_neutral = data.frame(
    site_name = c(
      "Beech Cr.", "Upper Mainstem J.D.", "Baldy Cr.", "Lower Mainstem J.D.",
      "Upper M. F. J.D.", "Granite Cr.", "Middle N. F. J.D.",
      "Big Wall Cr.", "Deer Cr.", "Murderers Cr."
    ),
    metapop = "John Day",
    stringsAsFactors = FALSE
  )
)

# -----------------------------
# read all sheets
# -----------------------------
sheet_objs <- lapply(neutral_sheets, function(sh) read_neutral_sheet(xlsx_file, sh))
names(sheet_objs) <- neutral_sheets

# -----------------------------
# verify same neutral locus panel across sheets
# -----------------------------
ref_loci <- sheet_objs[[1]]$locus_names
same_loci <- vapply(
  sheet_objs,
  function(x) identical(x$locus_names, ref_loci),
  logical(1)
)

if (!all(same_loci)) {
  stop("Neutral locus names are not identical across all sheets.")
}

# -----------------------------
# build one site lookup across all 79 populations
# -----------------------------
site_lookup_list <- vector("list", length(neutral_sheets))

for (i in seq_along(neutral_sheets)) {
  sh <- neutral_sheets[i]
  sh_obj <- sheet_objs[[sh]]
  sh_meta <- site_meta[[sh]]

  if (nrow(sh_meta) != length(sh_obj$pop_codes)) {
    stop(paste0("Site metadata length does not match pop code count for sheet ", sh, "."))
  }

  site_lookup_list[[i]] <- data.frame(
    sheet_name = sh,
    metapop = sh_meta$metapop,
    pop_code = sh_obj$pop_codes,
    site_name = sh_meta$site_name,
    n_ind = sh_obj$pop_n,
    stringsAsFactors = FALSE
  )
}

site_lookup <- do.call(rbind, site_lookup_list)
site_lookup$site_id <- seq_len(nrow(site_lookup))

# reorder columns
site_lookup <- site_lookup[, c("site_id", "sheet_name", "metapop", "pop_code", "site_name", "n_ind")]

# -----------------------------
# convert each sheet to hierfstat format and combine individuals
# -----------------------------
hf_list <- lapply(sheet_objs, function(x) sheet_to_hf(x, site_lookup))
ONMY_all_hf <- do.call(rbind, hf_list)

# drop loci that are entirely missing after combining
keep_loci <- colSums(!is.na(ONMY_all_hf[, -1, drop = FALSE])) > 0
ONMY_all_hf <- ONMY_all_hf[, c(TRUE, keep_loci), drop = FALSE]

# -----------------------------
# pairwise FST across all 79 populations
# -----------------------------
ONMY_all_fst <- hierfstat::pairwise.WCfst(ONMY_all_hf)
ONMY_all_fst <- as.matrix(ONMY_all_fst)

# force symmetry and cleanup
ONMY_all_fst[lower.tri(ONMY_all_fst)] <- t(ONMY_all_fst)[lower.tri(ONMY_all_fst)]
diag(ONMY_all_fst) <- 0
ONMY_all_fst[ONMY_all_fst < 0] <- 0

# replace row/col names with numeric site IDs
rownames(ONMY_all_fst) <- site_lookup$site_id
colnames(ONMY_all_fst) <- site_lookup$site_id

# -----------------------------
# final site dataframe
# -----------------------------
ONMY_all_sites <- site_lookup

# -----------------------------
# save objects
# -----------------------------

save(
  ONMY_2_fst,
  ONMY_2_fst,
  file = "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ONMY-2/data/ONMY-2.RData"
)

# -----------------------------
# quick output
# -----------------------------
cat("\nFinished combined neutral extraction across all five sheets.\n")
cat("Total populations:", nrow(ONMY_all_sites), "\n")
cat("Total individuals:", sum(ONMY_all_sites$n_ind), "\n")
cat("Total loci:", ncol(ONMY_all_hf) - 1, "\n\n")

print(head(ONMY_all_sites, 10))
print(dim(ONMY_all_fst))
print(round(ONMY_all_fst[1:min(10, nrow(ONMY_all_fst)), 1:min(10, ncol(ONMY_all_fst))], 5))



# plotting
ONMY_all_sites <- site_lookup


# -----------------------------
# add coordinates
# -----------------------------
# Fill these with your vetted lat/lon values.
# site_id must match ONMY_all_sites$site_id exactly.

site_coords <- data.frame(
  site_id = ONMY_all_sites$site_id,
  lat = NA_real_,
  lon = NA_real_
)

ONMY_all_sites <- merge(
  ONMY_all_sites,
  site_coords,
  by = "site_id",
  all.x = TRUE,
  sort = FALSE
)

# restore original order
ONMY_all_sites <- ONMY_all_sites[match(site_lookup$site_id, ONMY_all_sites$site_id), ]
rownames(ONMY_all_sites) <- NULL

# quick check
if (any(is.na(ONMY_all_sites$lat)) || any(is.na(ONMY_all_sites$lon))) {
  warning("Some site coordinates are missing. Map and IBD plots will fail until lat/lon are filled.")
}



# -----------------------------
# plot sites
# -----------------------------
library(ggplot2)
library(maps)

usa <- map_data("state")
canada <- map_data("world", region = "Canada")

p_sites <- ggplot() +
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
    data = ONMY_all_sites,
    aes(x = lon, y = lat, shape = metapop),
    size = 2.8
  ) +
  geom_text(
    data = ONMY_all_sites,
    aes(x = lon, y = lat, label = site_id),
    size = 2.7,
    nudge_y = 0.18
  ) +
  coord_quickmap(
    xlim = range(ONMY_all_sites$lon, na.rm = TRUE) + c(-2, 2),
    ylim = range(ONMY_all_sites$lat, na.rm = TRUE) + c(-2, 2)
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    shape = "Metapopulation",
    title = "ONMY all neutral sites"
  )

print(p_sites)

ggsave(
  filename = "ONMY_all_sites.png",
  plot = p_sites,
  width = 9,
  height = 7,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "ONMY_all_sites.pdf",
  plot = p_sites,
  width = 9,
  height = 7,
  bg = "white"
)