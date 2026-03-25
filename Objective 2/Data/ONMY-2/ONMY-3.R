############################
## ONMY-3 Eastern Cascades steelhead
## Hand et al. 2016 (Molecular Ecology)
## Neutral SNPs only
##
## This script:
##   1) reads the uploaded genepop file for the Eastern Cascades metapopulation
##   2) calculates pairwise Weir & Cockerham FST among populations
##   3) forces symmetry and sets negative FST to 0
##   4) renames rows/cols with numeric site IDs
##   5) builds a matching site dataframe
##
## NOTE:
## The genepop file contains 14 populations. I map those 14 POP codes to the
## 14 Eastern Cascades populations listed in Table S1, in the same order they
## appear in the supplement. That ordering assumption should be checked against
## the Dryad environmental file if you want a fully audited final version.
##
## Coordinates are left as NA here because they are not contained in the paper
## or the genepop file. Fill them after site vetting / georeferencing.
############################

suppressPackageStartupMessages({
  library(adegenet)
  library(hierfstat)
})

# -----------------------------
# file path
# -----------------------------
gen_file <- "/mnt/data/SH_MPG_Cascades_neutral_23pops_Feb26th_2014_GPOP"

# -----------------------------
# read genepop
# -----------------------------
ONMY_3_gen <- adegenet::read.genepop(gen_file, ncode = 3)

# population codes present in the file
pop_codes <- levels(ONMY_3_gen@pop)
print(pop_codes)

# -----------------------------
# site lookup
# -----------------------------
# Assumed to be in the same order as Eastern Cascades in Table S1.
site_lookup <- data.frame(
  pop_code = pop_codes,
  site_id  = seq_along(pop_codes),
  site_name = c(
    "Shitike Cr.",
    "Buckhollow Cr.",
    "Trout Cr.",
    "Fifteen Cr.",
    "Lower Little Klickitat R.",
    "Lower Summit Cr.",
    "Upper Trout Cr.",
    "Deadcanyon Cr.",
    "Lower White Cr.",
    "Snyder Cr.",
    "Surveyor Cr.",
    "Swale Cr.",
    "Rock Cr.",
    "Squaw Cr."
  ),
  tributary_region = c(
    "Deschutes",
    "Deschutes",
    "Deschutes",
    "Fifteen",
    "Little Klickitat",
    "Klickitat",
    "Klickitat",
    "Klickitat",
    "Klickitat",
    "Klickitat",
    "Klickitat",
    "Klickitat",
    "Rock",
    "Rock"
  ),
  lat = NA_real_,
  lon = NA_real_,
  stringsAsFactors = FALSE
)

# sanity check
if (nrow(site_lookup) != nPop(ONMY_3_gen)) {
  stop("site_lookup row count does not match the number of populations in the genepop file.")
}

# -----------------------------
# pairwise FST
# -----------------------------
ONMY_3_hf <- hierfstat::genind2hierfstat(ONMY_3_gen)
ONMY_3_fst <- hierfstat::pairwise.WCfst(ONMY_3_hf)

# force symmetry and clean negatives
ONMY_3_fst <- as.matrix(ONMY_3_fst)
ONMY_3_fst[lower.tri(ONMY_3_fst)] <- t(ONMY_3_fst)[lower.tri(ONMY_3_fst)]
diag(ONMY_3_fst) <- 0
ONMY_3_fst[ONMY_3_fst < 0] <- 0

# reorder to match site lookup / numeric site IDs
ONMY_3_fst <- ONMY_3_fst[site_lookup$pop_code, site_lookup$pop_code]
rownames(ONMY_3_fst) <- site_lookup$site_id
colnames(ONMY_3_fst) <- site_lookup$site_id

# -----------------------------
# coords dataframe
# -----------------------------
ONMY_3_coords <- site_lookup[, c("site_id", "site_name", "tributary_region", "lat", "lon")]

# -----------------------------
# optional pairwise dataframe for later IBD plotting
# -----------------------------
ONMY_3_ibd <- data.frame(
  site1 = rownames(ONMY_3_fst)[row(ONMY_3_fst)[upper.tri(ONMY_3_fst)]],
  site2 = colnames(ONMY_3_fst)[col(ONMY_3_fst)[upper.tri(ONMY_3_fst)]],
  fst   = ONMY_3_fst[upper.tri(ONMY_3_fst)],
  stringsAsFactors = FALSE
)

# -----------------------------
# save objects
# -----------------------------
save(
  ONMY_3_fst,
  ONMY_3_coords,
  ONMY_3_ibd,
  file = "/mnt/data/ONMY-3.RData"
)

# -----------------------------
# quick output
# -----------------------------
cat("\nFinished ONMY-3 extraction from genepop file.\n")
cat("Populations in genepop:", nPop(ONMY_3_gen), "\n")
cat("Individuals:", nInd(ONMY_3_gen), "\n\n")

print(ONMY_3_coords)
print(round(ONMY_3_fst, 5))
