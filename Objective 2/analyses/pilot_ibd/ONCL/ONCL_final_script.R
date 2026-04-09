# ONCL (e.g., Oncorhynchus clarkii) studies — same workflow as SAFO

library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(geosphere)

# -----------------------------
# paths
# -----------------------------
data_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"

out_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/analyses/pilot_ibd/ONCL/"

study_dirs <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)
study_dirs <- study_dirs[grepl("^ONCL-", study_dirs)]

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------
# helper functions
# -----------------------------
read_study <- function(study_code, data_dir) {
  rdata_path <- file.path(data_dir, study_code, "data", paste0(study_code, ".RData"))

  if (!file.exists(rdata_path)) {
    stop("Missing file: ", rdata_path)
  }

  e <- new.env(parent = emptyenv())
  load(rdata_path, envir = e)
  obj_names <- ls(e)

  fst <- get(obj_names[sapply(obj_names, function(x) {
    obj <- get(x, envir = e)
    is.matrix(obj) && is.numeric(obj) && nrow(obj) == ncol(obj)
  })][1], envir = e)

  coords <- get(obj_names[sapply(obj_names, function(x) {
    obj <- get(x, envir = e)
    is.data.frame(obj) && all(c("lat", "lon") %in% names(obj))
  })][1], envir = e)

  if (!"site_id" %in% names(coords)) coords$site_id <- seq_len(nrow(coords))

  if (nrow(coords) != nrow(fst)) {
    stop("Mismatch in ", study_code)
  }

  list(study_code = study_code, fst = fst, coords = coords)
}

study_to_pairs <- function(study) {
  idx <- which(lower.tri(study$fst), arr.ind = TRUE)

  tibble(
    study_code = study$study_code,
    fst = study$fst[idx],
    lon1 = study$coords$lon[idx[,1]],
    lat1 = study$coords$lat[idx[,1]],
    lon2 = study$coords$lon[idx[,2]],
    lat2 = study$coords$lat[idx[,2]],
    site_id1 = study$coords$site_id[idx[,1]],
    site_id2 = study$coords$site_id[idx[,2]]
  ) %>%
    mutate(
      distance_km = geosphere::distHaversine(
        matrix(c(lon1, lat1), ncol = 2),
        matrix(c(lon2, lat2), ncol = 2)
      ) / 1000,
      fst = ifelse(fst < 0, 0, fst)
    )
}

# -----------------------------
# read all studies
# -----------------------------
all_pairs <- map_dfr(study_dirs, ~ study_to_pairs(read_study(.x, data_dir))) %>%
  filter(is.finite(fst), is.finite(distance_km))

# order by dataset size (largest first)
study_order <- all_pairs %>%
  count(study_code, name = "n") %>%
  arrange(desc(n)) %>%
  pull(study_code)

all_pairs <- all_pairs %>%
  mutate(study_code = factor(study_code, levels = study_order)) %>%
  arrange(study_code)

# stats
n_unique_sites <- all_pairs %>%
  select(study_code, site_id1, site_id2) %>%
  pivot_longer(cols = c(site_id1, site_id2), values_to = "site_id") %>%
  distinct(study_code, site_id) %>%
  nrow()

n_pairwise <- nrow(all_pairs)

annot_br <- paste0(
  "Unique sites: ", format(n_unique_sites, big.mark=","), "\n",
  "Pairwise comparisons: ", format(n_pairwise, big.mark=",")
)

x_min <- min(all_pairs$distance_km)
x_max <- max(all_pairs$distance_km)
y_min <- min(all_pairs$fst)
y_max <- max(all_pairs$fst)

# plot
p <- ggplot(all_pairs, aes(distance_km, fst, color = study_code)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(
    data = all_pairs,
    aes(distance_km, fst),
    method = "lm",
    se = FALSE,
    inherit.aes = FALSE,
    color = "blue",
    linewidth = 1
  ) +
  annotate("text",
           x = x_max*0.98,
           y = y_min + 0.04*(y_max - y_min),
           label = annot_br,
           hjust = 1,
           vjust = 0,
           size = 4) +
  annotate("text",
           x = x_min,
           y = y_max,
           label = "Cutthroat Trout (ONCL)",
           hjust = 0,
           vjust = 1,
           size = 5) +
  theme_bw(base_size = 12) +
  labs(
    x = "Straight-line distance among sites (km)",
    y = expression(pairwise~F[ST]),
    color = "Study"
  )

print(p)

# save
ggsave(file.path(out_dir, "ONCL_all_studies_fst_vs_distance.png"),
       p, width = 9, height = 6, dpi = 300)

write.csv(all_pairs,
          file.path(out_dir, "ONCL_all_pairs.csv"),
          row.names = FALSE)
