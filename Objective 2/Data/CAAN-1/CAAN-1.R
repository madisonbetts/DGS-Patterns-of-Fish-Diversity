############################
## CAAN_1 workflow
## 1) site coordinates
## 2) pairwise FST matrix
## 3) pairwise river-distance matrix (km)
## 4) long-format site-pair data frame
## 5) IBD plot using river distance
############################

library(ggplot2)
library(ggrepel)

# -------------------------
# 1) Site coordinates
# -------------------------
# These are the sampling coordinates for the 20 sites.
# Object name requested: CAAN_1_coords

CAAN_1_coords <- data.frame(
  site = c(
    "BRUS","CLCR","COWS","HOGC","HUZZ","LIBE","SUGA","SULO","SWAN","SYLA",
    "CHER","FLAT","LONG","MAST","MINK","OSAG","SPID","TOWN","WEFO","WHAR"
  ),
  lat = c(
    36.131211, 36.122016, 36.965145, 36.151714, 36.231834,
    36.799673, 36.291167, 36.387254, 36.700713, 35.994070,
    35.989033, 36.075829, 36.292399, 35.853857, 35.893545,
    36.163822, 36.441026, 36.044794, 35.848494, 36.021875
  ),
  lon = c(
    -93.948083, -92.903423, -92.726148, -92.960925, -92.990015,
    -92.908690, -92.919110, -92.970669, -93.097324, -92.210800,
    -93.827983, -93.103605, -93.282247, -94.004442, -93.582264,
    -93.353172, -93.841951, -94.176363, -94.105501, -93.631034
  ),
  stringsAsFactors = FALSE
)

CAAN_1_coords

# -------------------------
# Quick map of sampling sites
# -------------------------
ggplot(CAAN_1_coords, aes(lon, lat)) +
  geom_point(size = 2.5) +
  geom_text_repel(aes(label = site), size = 3.5) +
  coord_fixed() +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "CAAN_1 sampling localities"
  )

# -------------------------
# 2) Pairwise FST matrix
# -------------------------
# Enter the lower triangle of the pairwise FST matrix row-wise.
# Object name requested: CAAN_1_fst

# duplicate sites
sites <- CAAN_1_coords$site

# fst matrices
fst_vals <- c(
  0.00,
  0.01, 0.00,
  0.13, 0.00, 0.08,
  0.00, 0.00, 0.00, 0.04,
  0.00, 0.05, 0.18, 0.10, 0.05,
  0.63, 0.00, 0.10, 0.00, 0.00, 0.02,
  0.00, 0.10, 0.12, 0.25, 0.00, 0.08, 0.15,
  0.11, 0.17, 0.13, 0.33, 0.20, 0.40, 0.35, 0.21,
  0.00, 0.00, 0.00, 0.02, 0.00, 0.15, 0.01, 0.12, 0.05,
  0.00, 0.13, 0.14, 0.23, 0.05, 0.06, 0.16, 0.00, 0.19, 0.12,
  0.05, 0.00, 0.00, 0.00, 0.02, 0.11, 0.00, 0.20, 0.20, 0.00, 0.18,
  0.00, 0.00, 0.00, 0.07, 0.00, 0.00, 0.04, 0.00, 0.07, 0.00, 0.00, 0.00,
  0.00, 0.03, 0.02, 0.14, 0.02, 0.14, 0.13, 0.00, 0.00, 0.00, 0.00, 0.05, 0.00,
  0.42, 0.42, 0.44, 0.58, 0.48, 0.67, 0.59, 0.53, 0.07, 0.33, 0.44, 0.49, 0.40, 0.22,
  0.08, 0.01, 0.01, 0.14, 0.10, 0.28, 0.19, 0.20, 0.00, 0.00, 0.19, 0.01, 0.01, 0.00, 0.24,
  0.24, 0.13, 0.29, 0.12, 0.23, 0.15, 0.02, 0.34, 0.47, 0.16, 0.30, 0.12, 0.20, 0.27, 0.67, 0.30,
  0.03, 0.00, 0.00, 0.00, 0.03, 0.14, 0.05, 0.16, 0.10, 0.00, 0.15, 0.00, 0.00, 0.00, 0.35, 0.00, 0.18,
  0.24, 0.31, 0.28, 0.46, 0.33, 0.52, 0.46, 0.33, 0.00, 0.20, 0.27, 0.35, 0.23, 0.05, 0.10, 0.14, 0.57, 0.22,
  0.00, 0.18, 0.19, 0.31, 0.10, 0.18, 0.22, 0.00, 0.18, 0.14, 0.00, 0.25, 0.00, 0.00, 0.42, 0.21, 0.38, 0.18, 0.23
)

# make raw matrix
CAAN_1_fst <- matrix(
  NA_real_,
  nrow = length(sites),
  ncol = length(sites),
  dimnames = list(sites, sites)
)

# fill in opposite side 
CAAN_1_fst[lower.tri(CAAN_1_fst)] <- fst_vals
CAAN_1_fst[upper.tri(CAAN_1_fst)] <- t(CAAN_1_fst)[upper.tri(CAAN_1_fst)]
diag(CAAN_1_fst) <- 0

CAAN_1_fst # QC

# -------------------------
# 3) Pairwise river-distance matrix (km)
# -------------------------

# river distances values
rivdist_vals <- c(
  448,
  293, 289,
  447, 13, 288,
  446, 40, 287, 42,
  263, 259, 30, 258, 257,
  430, 39, 270, 39, 19, 240,
  296, 176, 137, 181, 187, 107, 184,
  226, 222, 67, 221, 220, 37, 203, 70,
  392, 148, 307, 157, 172, 277, 151, 204, 166,
  41, 468, 313, 467, 466, 283, 450, 316, 246, 412,
  504, 245, 344, 232, 228, 314, 228, 260, 277, 208, 523,
  233, 312, 157, 311, 310, 126, 293, 160, 90, 255, 253, 367,
  63, 492, 337, 491, 490, 306, 473, 340, 270, 435, 57, 547, 277,
  273, 446, 291, 445, 444, 261, 428, 294, 224, 390, 293, 501, 225, 317,
  427, 409, 254, 408, 407, 223, 390, 257, 187, 352, 248, 464, 192, 272, 114,
  74, 374, 219, 373, 372, 189, 355, 222, 152, 318, 94, 429, 159, 118, 199, 154,
  47, 443, 287, 442, 441, 258, 425, 291, 221, 387, 41, 498, 234, 41, 243, 214, 93,
  71, 497, 343, 497, 496, 312, 479, 346, 276, 441, 66, 553, 283, 67, 323, 278, 124, 34,
  97, 517, 362, 516, 515, 332, 499, 365, 295, 461, 120, 572, 302, 143, 340, 297, 143, 121, 147
)

# make empty matrix
CAAN_1_rivdists <- matrix(
  NA_real_,
  nrow = length(sites),
  ncol = length(sites),
  dimnames = list(sites, sites)
)

# fill in opposite side of matrix
CAAN_1_rivdists[lower.tri(CAAN_1_rivdists)] <- rivdist_vals
CAAN_1_rivdists[upper.tri(CAAN_1_rivdists)] <- t(CAAN_1_rivdists)[upper.tri(CAAN_1_rivdists)]
diag(CAAN_1_rivdists) <- 0

CAAN_1_rivdists

# -------------------------
# 4) Convert matrices to long format for IBD
# -------------------------

# Keep only one copy of each pair (lower triangle).
pair_idx <- lower.tri(CAAN_1_fst)

# 
CAAN_1_ibd_df <- data.frame(
  site1 = rownames(CAAN_1_fst)[row(CAAN_1_fst)[pair_idx]],
  site2 = colnames(CAAN_1_fst)[col(CAAN_1_fst)[pair_idx]],
  fst = CAAN_1_fst[pair_idx],
  river_distance_km = CAAN_1_rivdists[pair_idx],
  stringsAsFactors = FALSE
)

head(CAAN_1_ibd_df)

# -------------------------
# 5) IBD plot using river distance
# -------------------------
ggplot(CAAN_1_ibd_df, aes(river_distance_km, fst)) +
  geom_point(size = 2.3) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  labs(
    x = "River distance (km)",
    y = expression(pairwise~F[ST]),
    title = "Isolation by distance for CAAN_1"
  )