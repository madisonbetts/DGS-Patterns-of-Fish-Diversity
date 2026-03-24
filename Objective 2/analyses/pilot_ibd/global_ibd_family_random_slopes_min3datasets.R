# -------------------------------
# Global IBD workflow across all datasets
# random intercepts + slopes by Family,
# restricted to families represented by >= 3 datasets
# -------------------------------

library(terra)
library(readxl)
library(dplyr)
library(ggplot2)
library(lme4)
library(stringr)
library(tibble)

base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"
metadata_file <- file.path(base_dir, "Study_metadata.xlsx")

# -------------------------------
# metadata
# -------------------------------
meta <- read_excel(metadata_file)
names(meta) <- trimws(names(meta))

if (!"Study_ID" %in% names(meta)) stop("Column 'Study_ID' not found in Study_metadata.xlsx")
if (!"Family" %in% names(meta)) stop("Column 'Family' not found in Study_metadata.xlsx")

meta <- meta %>%
  mutate(study_code = as.character(Study_ID)) %>%
  select(study_code, Family, Common_name, Spp)

# -------------------------------
# list valid study folders
# -------------------------------
study_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
study_dirs <- study_dirs[grepl("^[A-Za-z0-9]+-[0-9]+$", basename(study_dirs))]

# -------------------------------
# helper: identify coordinate df
# -------------------------------
find_lon_lat_cols <- function(df) {
  nms <- names(df)
  nms_lower <- tolower(nms)

  lon_idx <- grep("(^lon$|^long$|^longitude$|lon|long|longitude)", nms_lower)
  lat_idx <- grep("(^lat$|^latitude$|lat|latitude)", nms_lower)

  if (length(lon_idx) == 0 || length(lat_idx) == 0) return(NULL)

  out <- df[, c(nms[lon_idx[1]], nms[lat_idx[1]]), drop = FALSE]
  names(out) <- c("lon", "lat")

  site_idx <- grep("(^site$|site|pop|locality|location|id)", nms_lower)
  if (length(site_idx) > 0) {
    out$site <- as.character(df[[nms[site_idx[1]]]])
  } else {
    out$site <- rownames(df)
  }

  out <- out %>%
    mutate(
      lon = suppressWarnings(as.numeric(lon)),
      lat = suppressWarnings(as.numeric(lat)),
      site = as.character(site)
    ) %>%
    filter(!is.na(lon), !is.na(lat))

  if (nrow(out) == 0) return(NULL)
  out
}

# -------------------------------
# helper: identify fst matrix + coords df
# -------------------------------
extract_fst_and_coords <- function(rdata_file) {
  e <- new.env()
  load(rdata_file, envir = e)
  obj_names <- ls(e)

  fst_candidates <- obj_names[grepl("_fst$|fst", tolower(obj_names))]
  coord_candidates <- c(
    obj_names[grepl("_coords$|coords|sites", tolower(obj_names))],
    obj_names[!grepl("_fst$|fst", tolower(obj_names))]
  ) %>% unique()

  fst <- NULL
  coords <- NULL
  fst_name <- NA_character_
  coords_name <- NA_character_

  for (nm in fst_candidates) {
    obj <- get(nm, envir = e)
    if (is.matrix(obj) || is.data.frame(obj)) {
      m <- as.matrix(obj)
      if (nrow(m) == ncol(m) && !is.null(rownames(m)) && !is.null(colnames(m))) {
        fst <- m
        fst_name <- nm
        break
      }
    }
  }

  for (nm in coord_candidates) {
    obj <- get(nm, envir = e)

    if (is.data.frame(obj)) {
      tmp <- find_lon_lat_cols(obj)
      if (!is.null(tmp)) {
        coords <- tmp
        coords_name <- nm
        break
      }
    }

    if (is.list(obj) && !is.data.frame(obj)) {
      sub_nms <- names(obj)
      if (length(sub_nms) > 0) {
        for (sub_nm in sub_nms) {
          sub_obj <- obj[[sub_nm]]
          if (is.data.frame(sub_obj)) {
            tmp <- find_lon_lat_cols(sub_obj)
            if (!is.null(tmp)) {
              coords <- tmp
              coords_name <- paste0(nm, "$", sub_nm)
              break
            }
          }
        }
      }
      if (!is.null(coords)) break
    }
  }

  list(
    fst = fst,
    coords = coords,
    fst_name = fst_name,
    coords_name = coords_name
  )
}

# -------------------------------
# build pooled pairwise dataframe
# -------------------------------
all_ibd_list <- list()
skipped <- character(0)

for (d in study_dirs) {

  study_code <- basename(d)
  rdata_file <- list.files(file.path(d, "data"), pattern = "\\.RData$", full.names = TRUE)

  if (length(rdata_file) == 0) {
    skipped <- c(skipped, study_code)
    next
  }

  objs <- extract_fst_and_coords(rdata_file[1])
  fst <- objs$fst
  coords <- objs$coords

  if (is.null(fst) || is.null(coords)) {
    skipped <- c(skipped, study_code)
    next
  }

  fst <- as.matrix(fst)
  mode(fst) <- "numeric"

  if (is.null(rownames(fst)) || is.null(colnames(fst))) {
    skipped <- c(skipped, study_code)
    next
  }

  if ("site" %in% names(coords)) {
    coords <- coords[match(rownames(fst), coords$site), , drop = FALSE]
  } else {
    coords <- coords[match(rownames(fst), rownames(coords)), , drop = FALSE]
  }

  if (nrow(coords) != nrow(fst) || any(is.na(coords$lon)) || any(is.na(coords$lat))) {
    skipped <- c(skipped, study_code)
    next
  }

  dmat <- terra::distance(
    terra::vect(coords[, c("lon", "lat")], crs = "EPSG:4326")
  ) / 1000
  dmat <- as.matrix(dmat)

  pair_df <- data.frame(
    study_code = study_code,
    site1 = rownames(fst)[row(fst)[upper.tri(fst)]],
    site2 = colnames(fst)[col(fst)[upper.tri(fst)]],
    dist_km = dmat[upper.tri(dmat)],
    fst = fst[upper.tri(fst)],
    stringsAsFactors = FALSE
  )

  pair_df <- pair_df %>%
    filter(!is.na(dist_km), !is.na(fst)) %>%
    mutate(species_code = sub("-[0-9]+$", "", study_code))

  all_ibd_list[[study_code]] <- pair_df
}

all_ibd <- bind_rows(all_ibd_list) %>%
  left_join(meta, by = "study_code") %>%
  filter(!is.na(Family)) %>%
  mutate(fst = ifelse(fst < 0, 0, fst))

# -------------------------------
# keep only families with >= 3 datasets
# -------------------------------
family_dataset_counts <- all_ibd %>%
  distinct(study_code, Family) %>%
  count(Family, name = "n_datasets") %>%
  arrange(desc(n_datasets), Family)

eligible_families <- family_dataset_counts %>%
  filter(n_datasets >= 3) %>%
  pull(Family)

all_ibd_mod <- all_ibd %>%
  filter(Family %in% eligible_families) %>%
  mutate(
    Family = factor(Family, levels = family_dataset_counts$Family[family_dataset_counts$Family %in% eligible_families]),
    study_code = factor(study_code)
  )

# -------------------------------
# summaries
# -------------------------------
cat("\n================ IBD DATA SUMMARY ================\n")
cat("All datasets with usable IBD data: ", n_distinct(all_ibd$study_code), "\n", sep = "")
cat("All families with usable IBD data: ", n_distinct(all_ibd$Family), "\n", sep = "")
cat("Families eligible for random slopes (>= 3 datasets): ", length(eligible_families), "\n", sep = "")
cat("Datasets retained for model: ", n_distinct(all_ibd_mod$study_code), "\n", sep = "")
cat("Pairwise comparisons retained for model: ", nrow(all_ibd_mod), "\n", sep = "")

cat("\nEligible families:\n")
print(family_dataset_counts %>% filter(Family %in% eligible_families))

if (length(skipped) > 0) {
  cat("\nSkipped datasets:\n")
  cat(paste0("  ", sort(unique(skipped))), sep = "\n")
  cat("\n")
}
cat("==================================================\n\n")

# -------------------------------
# scale distance for modeling
# -------------------------------
dist_center <- mean(all_ibd_mod$dist_km, na.rm = TRUE)
dist_scale  <- sd(all_ibd_mod$dist_km, na.rm = TRUE)

all_ibd_mod <- all_ibd_mod %>%
  mutate(dist_sc = (dist_km - dist_center) / dist_scale)

# -------------------------------
# mixed model:
# fixed global slope
# random intercept + slope by family
# random intercept by dataset
# -------------------------------
ibd_mod <- lmer(
  fst ~ dist_sc + (dist_sc | Family) + (1 | study_code),
  data = all_ibd_mod,
  REML = TRUE
)

cat("\n================ MODEL SUMMARY ===================\n")
print(summary(ibd_mod))
cat("==================================================\n\n")

# -------------------------------
# extract family random effects
# -------------------------------
re_family <- ranef(ibd_mod)$Family %>%
  rownames_to_column("Family") %>%
  rename(
    rand_intercept = `(Intercept)`,
    rand_slope = dist_sc
  )

fixef_vals <- fixef(ibd_mod)

family_coef <- re_family %>%
  mutate(
    intercept = fixef_vals["(Intercept)"] + rand_intercept,
    slope = fixef_vals["dist_sc"] + rand_slope
  ) %>%
  left_join(
    all_ibd_mod %>% distinct(study_code, Family) %>% count(Family, name = "n_datasets"),
    by = "Family"
  ) %>%
  left_join(
    all_ibd_mod %>% count(Family, name = "n_pairs"),
    by = "Family"
  ) %>%
  arrange(desc(slope))

cat("\n=========== FAMILY RANDOM EFFECTS ================\n")
print(family_coef)
cat("==================================================\n\n")

write.csv(
  family_coef,
  file.path(base_dir, "family_random_effects_min3datasets.csv"),
  row.names = FALSE
)

# -------------------------------
# prediction grid for family lines
# -------------------------------
dist_grid <- seq(min(all_ibd_mod$dist_km), max(all_ibd_mod$dist_km), length.out = 200)

pred_lines <- expand.grid(
  dist_km = dist_grid,
  Family = levels(all_ibd_mod$Family),
  stringsAsFactors = FALSE
) %>%
  mutate(dist_sc = (dist_km - dist_center) / dist_scale) %>%
  left_join(family_coef %>% select(Family, intercept, slope), by = "Family") %>%
  mutate(pred_fst = intercept + slope * dist_sc)

# -------------------------------
# plot 1: pooled scatter + family fitted lines
# -------------------------------
p1 <- ggplot(all_ibd_mod, aes(x = dist_km, y = fst, color = Family)) +
  geom_point(alpha = 0.22, size = 1.1) +
  geom_line(
    data = pred_lines,
    aes(x = dist_km, y = pred_fst, color = Family),
    linewidth = 1.1,
    alpha = 0.95
  ) +
  theme_classic() +
  labs(
    x = "Geographic distance (km)",
    y = expression(F[ST]),
    title = "Global IBD with family-specific random-slope fits",
    subtitle = "Only families with >= 3 datasets retained"
  )

print(p1)

# -------------------------------
# plot 2: family random intercepts
# -------------------------------
p2_dat <- family_coef %>%
  arrange(intercept) %>%
  mutate(Family = factor(Family, levels = Family))

p2 <- ggplot(p2_dat, aes(x = intercept, y = Family)) +
  geom_vline(xintercept = fixef_vals["(Intercept)"], linetype = 2, color = "grey50") +
  geom_point(size = 3) +
  theme_classic() +
  labs(
    x = "Family-specific intercept",
    y = NULL,
    title = "Random intercepts by family",
    subtitle = "Only families with >= 3 datasets retained"
  )

print(p2)

# -------------------------------
# plot 3: family random slopes
# -------------------------------
p3_dat <- family_coef %>%
  arrange(slope) %>%
  mutate(Family = factor(Family, levels = Family))

p3 <- ggplot(p3_dat, aes(x = slope, y = Family)) +
  geom_vline(xintercept = fixef_vals["dist_sc"], linetype = 2, color = "grey50") +
  geom_point(size = 3) +
  theme_classic() +
  labs(
    x = "Family-specific slope on scaled distance",
    y = NULL,
    title = "Random slopes by family",
    subtitle = "Only families with >= 3 datasets retained"
  )

print(p3)

# -------------------------------
# save outputs
# -------------------------------
ggsave(
  file.path(base_dir, "global_ibd_family_lines_min3datasets.png"),
  p1,
  width = 11,
  height = 7,
  dpi = 400
)

ggsave(
  file.path(base_dir, "family_random_intercepts_min3datasets.png"),
  p2,
  width = 8,
  height = 6,
  dpi = 400
)

ggsave(
  file.path(base_dir, "family_random_slopes_min3datasets.png"),
  p3,
  width = 8,
  height = 6,
  dpi = 400
)
