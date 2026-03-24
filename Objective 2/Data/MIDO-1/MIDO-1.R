# -----------------------------
# MIDO-1 Smallmouth Bass
# Euclide et al. 2020
# Table 2 FST + Fig 1 sites
# -----------------------------

library(dplyr)
library(geosphere)
library(ggplot2)
library(maps)

# -----------------------------
# paths
# -----------------------------
setwd("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/MIDO-1")

# -----------------------------
# site key in FST / paper order
# -----------------------------
site_key <- data.frame(
  site = 1:32,
  code = c(
    "APS1", "APS2", "LBDN", "BBDN", "BOAR", "BPS1", "BPS2", "CEDR",
    "CHCF", "EPFC", "FOXR", "GALR", "LKEG", "LSTB", "LWAU", "LWIR",
    "MANI", "MENA", "MENB", "MILR", "MINE", "NBLM", "NEBL", "OCON",
    "PALL", "PESH", "ROWB", "SAWH", "SHEB", "USHF", "WASH", "WTWI"
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# sites (site, lat, lon ONLY)
# best-available named waterbody anchors matched to paper sites
# site numbering preserved from current script
# -----------------------------
MIDO_1_sites <- data.frame(
  site = 1:32,
  lat = c(
    43.3789,  # 1  APS1  Lake Wisconsin
    43.3100,  # 2  APS2  Lower Wisconsin River
    45.7666,  # 3  LBDN  Little Bay de Noc
    45.7833,  # 4  BBDN  Big Bay de Noc
    44.6753,  # 5  BOAR  Boardman River
    43.1840,  # 6  BPS1  Lower Wisconsin River
    43.1650,  # 7  BPS2  Lower Wisconsin River
    45.4111,  # 8  MENA  Upper Menominee River
    44.8650,  # 9  CHCF  Chippewa River
    45.1569,  # 10 MENB  Lower Menominee River
    44.5200,  # 11 FOXR  Fox River
    42.4170,  # 12 GALR  Galena River
    42.9617,  # 13 LKEG  Lake Kegonsa
    44.8275,  # 14 PESH  Peshtigo River
    43.0019,  # 15 LWAU  Lake Waubesa
    43.1160,  # 16 LWIR  Lower Wisconsin River
    45.9578,  # 17 MANI  Manistique
    45.7500,  # 18 CEDR  Cedar River
    45.0960,  # 19 OCON  Oconto River
    43.1160,  # 20 MILR  Milwaukee River
    42.8580,  # 21 MINE  Pecatonica River
    45.1800,  # 22 WASH  Washington Island
    46.0520,  # 23 NEBL  Nebish Lake
    44.8850,  # 24 EPFC  Ephraim
    46.0500,  # 25 PALL  Pallette Lake
    45.0500,  # 26 LSTB  Sturgeon Bay
    45.1900,  # 27 ROWB  Rowleys Bay
    44.8600,  # 28 SAWH  Sawyer Harbor
    43.7430,  # 29 SHEB  Sheboygan River
    42.5050,  # 30 USHF  Illinois–Fox River
    45.3990,  # 31 NBLM  North Bay
    44.1820   # 32 WTWI  West Twin River
  ),
  lon = c(
    -89.7569,  # 1  APS1
    -89.9600,  # 2  APS2
    -87.0126,  # 3  LBDN
    -86.7001,  # 4  BBDN
    -85.6309,  # 5  BOAR
    -89.7050,  # 6  BPS1
    -89.6670,  # 7  BPS2
    -87.3546,  # 8  MENA
    -91.9300,  # 9  CHCF
    -87.1709,  # 10 MENB
    -88.0150,  # 11 FOXR
    -90.4290,  # 12 GALR
    -89.3157,  # 13 LKEG
    -87.3770,  # 14 PESH
    -89.3172,  # 15 LWAU
    -90.7050,  # 16 LWIR
    -86.2463,  # 17 MANI
    -88.0400,  # 18 CEDR
    -87.6100,  # 19 OCON
    -87.9100,  # 20 MILR
    -90.1840,  # 21 MINE
    -87.0580,  # 22 WASH
    -89.6760,  # 23 NEBL
    -87.8650,  # 24 EPFC
    -89.6400,  # 25 PALL
    -87.7450,  # 26 LSTB
    -86.9900,  # 27 ROWB
    -87.3400,  # 28 SAWH
    -87.7140,  # 29 SHEB
    -88.9450,  # 30 USHF
    -86.9100,  # 31 NBLM
    -87.5880   # 32 WTWI
  )
)

# -----------------------------
# FST matrix (Table 2)
# -----------------------------
MIDO_1_fst <- matrix(0, 32, 32)
rownames(MIDO_1_fst) <- colnames(MIDO_1_fst) <- site_key$site

MIDO_1_fst[2,1] <- 0.04
MIDO_1_fst[3,1] <- 0.17; MIDO_1_fst[3,2] <- 0.12
MIDO_1_fst[4,1] <- 0.20; MIDO_1_fst[4,2] <- 0.17; MIDO_1_fst[4,3] <- 0.01
MIDO_1_fst[5,1] <- 0.27; MIDO_1_fst[5,2] <- 0.22; MIDO_1_fst[5,3] <- 0.08; MIDO_1_fst[5,4] <- 0.08
MIDO_1_fst[6,1] <- 0.01; MIDO_1_fst[6,2] <- 0.03; MIDO_1_fst[6,3] <- 0.15; MIDO_1_fst[6,4] <- 0.19; MIDO_1_fst[6,5] <- 0.25
MIDO_1_fst[7,1] <- 0.00; MIDO_1_fst[7,2] <- 0.03; MIDO_1_fst[7,3] <- 0.15; MIDO_1_fst[7,4] <- 0.19; MIDO_1_fst[7,5] <- 0.25; MIDO_1_fst[7,6] <- 0.00
MIDO_1_fst[8,1] <- 0.18; MIDO_1_fst[8,2] <- 0.14; MIDO_1_fst[8,3] <- 0.03; MIDO_1_fst[8,4] <- 0.05; MIDO_1_fst[8,5] <- 0.12; MIDO_1_fst[8,6] <- 0.17; MIDO_1_fst[8,7] <- 0.17
MIDO_1_fst[9,1] <- 0.07; MIDO_1_fst[9,2] <- 0.07; MIDO_1_fst[9,3] <- 0.17; MIDO_1_fst[9,4] <- 0.21; MIDO_1_fst[9,5] <- 0.26; MIDO_1_fst[9,6] <- 0.06; MIDO_1_fst[9,7] <- 0.05; MIDO_1_fst[9,8] <- 0.18
MIDO_1_fst[10,1] <- 0.19; MIDO_1_fst[10,2] <- 0.16; MIDO_1_fst[10,3] <- 0.04; MIDO_1_fst[10,4] <- 0.02; MIDO_1_fst[10,5] <- 0.08; MIDO_1_fst[10,6] <- 0.17; MIDO_1_fst[10,7] <- 0.17; MIDO_1_fst[10,8] <- 0.05; MIDO_1_fst[10,9] <- 0.19
MIDO_1_fst[11,1] <- 0.12; MIDO_1_fst[11,2] <- 0.10; MIDO_1_fst[11,3] <- 0.10; MIDO_1_fst[11,4] <- 0.12; MIDO_1_fst[11,5] <- 0.19; MIDO_1_fst[11,6] <- 0.10; MIDO_1_fst[11,7] <- 0.10; MIDO_1_fst[11,8] <- 0.08; MIDO_1_fst[11,9] <- 0.13; MIDO_1_fst[11,10] <- 0.09
MIDO_1_fst[12,1] <- 0.06; MIDO_1_fst[12,2] <- 0.04; MIDO_1_fst[12,3] <- 0.10; MIDO_1_fst[12,4] <- 0.15; MIDO_1_fst[12,5] <- 0.22; MIDO_1_fst[12,6] <- 0.04; MIDO_1_fst[12,7] <- 0.05; MIDO_1_fst[12,8] <- 0.13; MIDO_1_fst[12,9] <- 0.07; MIDO_1_fst[12,10] <- 0.15; MIDO_1_fst[12,11] <- 0.07
MIDO_1_fst[13,1] <- 0.19; MIDO_1_fst[13,2] <- 0.13; MIDO_1_fst[13,3] <- 0.23; MIDO_1_fst[13,4] <- 0.27; MIDO_1_fst[13,5] <- 0.34; MIDO_1_fst[13,6] <- 0.13; MIDO_1_fst[13,7] <- 0.16; MIDO_1_fst[13,8] <- 0.22; MIDO_1_fst[13,9] <- 0.17; MIDO_1_fst[13,10] <- 0.26; MIDO_1_fst[13,11] <- 0.19; MIDO_1_fst[13,12] <- 0.14
MIDO_1_fst[14,1] <- 0.16; MIDO_1_fst[14,2] <- 0.12; MIDO_1_fst[14,3] <- 0.05; MIDO_1_fst[14,4] <- 0.06; MIDO_1_fst[14,5] <- 0.09; MIDO_1_fst[14,6] <- 0.15; MIDO_1_fst[14,7] <- 0.15; MIDO_1_fst[14,8] <- 0.06; MIDO_1_fst[14,9] <- 0.16; MIDO_1_fst[14,10] <- 0.02; MIDO_1_fst[14,11] <- 0.07; MIDO_1_fst[14,12] <- 0.11; MIDO_1_fst[14,13] <- 0.25
MIDO_1_fst[15,1] <- 0.11; MIDO_1_fst[15,2] <- 0.09; MIDO_1_fst[15,3] <- 0.20; MIDO_1_fst[15,4] <- 0.23; MIDO_1_fst[15,5] <- 0.30; MIDO_1_fst[15,6] <- 0.08; MIDO_1_fst[15,7] <- 0.09; MIDO_1_fst[15,8] <- 0.21; MIDO_1_fst[15,9] <- 0.11; MIDO_1_fst[15,10] <- 0.22; MIDO_1_fst[15,11] <- 0.14; MIDO_1_fst[15,12] <- 0.10; MIDO_1_fst[15,13] <- 0.11; MIDO_1_fst[15,14] <- 0.20
MIDO_1_fst[16,1] <- 0.04; MIDO_1_fst[16,2] <- 0.02; MIDO_1_fst[16,3] <- 0.13; MIDO_1_fst[16,4] <- 0.18; MIDO_1_fst[16,5] <- 0.24; MIDO_1_fst[16,6] <- 0.01; MIDO_1_fst[16,7] <- 0.02; MIDO_1_fst[16,8] <- 0.16; MIDO_1_fst[16,9] <- 0.05; MIDO_1_fst[16,10] <- 0.17; MIDO_1_fst[16,11] <- 0.10; MIDO_1_fst[16,12] <- 0.03; MIDO_1_fst[16,13] <- 0.10; MIDO_1_fst[16,14] <- 0.15; MIDO_1_fst[16,15] <- 0.06
MIDO_1_fst[17,1] <- 0.18; MIDO_1_fst[17,2] <- 0.13; MIDO_1_fst[17,3] <- 0.04; MIDO_1_fst[17,4] <- 0.05; MIDO_1_fst[17,5] <- 0.10; MIDO_1_fst[17,6] <- 0.16; MIDO_1_fst[17,7] <- 0.16; MIDO_1_fst[17,8] <- 0.06; MIDO_1_fst[17,9] <- 0.16; MIDO_1_fst[17,10] <- 0.06; MIDO_1_fst[17,11] <- 0.11; MIDO_1_fst[17,12] <- 0.11; MIDO_1_fst[17,13] <- 0.24; MIDO_1_fst[17,14] <- 0.05; MIDO_1_fst[17,15] <- 0.20; MIDO_1_fst[17,16] <- 0.14
MIDO_1_fst[18,1] <- 0.15; MIDO_1_fst[18,2] <- 0.11; MIDO_1_fst[18,3] <- 0.03; MIDO_1_fst[18,4] <- 0.06; MIDO_1_fst[18,5] <- 0.07; MIDO_1_fst[18,6] <- 0.13; MIDO_1_fst[18,7] <- 0.14; MIDO_1_fst[18,8] <- 0.07; MIDO_1_fst[18,9] <- 0.14; MIDO_1_fst[18,10] <- 0.06; MIDO_1_fst[18,11] <- 0.11; MIDO_1_fst[18,12] <- 0.10; MIDO_1_fst[18,13] <- 0.19; MIDO_1_fst[18,14] <- 0.05; MIDO_1_fst[18,15] <- 0.17; MIDO_1_fst[18,16] <- 0.12; MIDO_1_fst[18,17] <- 0.04
MIDO_1_fst[19,1] <- 0.12; MIDO_1_fst[19,2] <- 0.11; MIDO_1_fst[19,3] <- 0.04; MIDO_1_fst[19,4] <- 0.06; MIDO_1_fst[19,5] <- 0.11; MIDO_1_fst[19,6] <- 0.11; MIDO_1_fst[19,7] <- 0.10; MIDO_1_fst[19,8] <- 0.05; MIDO_1_fst[19,9] <- 0.11; MIDO_1_fst[19,10] <- 0.04; MIDO_1_fst[19,11] <- 0.06; MIDO_1_fst[19,12] <- 0.08; MIDO_1_fst[19,13] <- 0.19; MIDO_1_fst[19,14] <- 0.06; MIDO_1_fst[19,15] <- 0.16; MIDO_1_fst[19,16] <- 0.10; MIDO_1_fst[19,17] <- 0.06; MIDO_1_fst[19,18] <- 0.04
MIDO_1_fst[20,1] <- 0.12; MIDO_1_fst[20,2] <- 0.10; MIDO_1_fst[20,3] <- 0.07; MIDO_1_fst[20,4] <- 0.09; MIDO_1_fst[20,5] <- 0.11; MIDO_1_fst[20,6] <- 0.11; MIDO_1_fst[20,7] <- 0.11; MIDO_1_fst[20,8] <- 0.08; MIDO_1_fst[20,9] <- 0.12; MIDO_1_fst[20,10] <- 0.07; MIDO_1_fst[20,11] <- 0.08; MIDO_1_fst[20,12] <- 0.10; MIDO_1_fst[20,13] <- 0.21; MIDO_1_fst[20,14] <- 0.05; MIDO_1_fst[20,15] <- 0.16; MIDO_1_fst[20,16] <- 0.11; MIDO_1_fst[20,17] <- 0.06; MIDO_1_fst[20,18] <- 0.06; MIDO_1_fst[20,19] <- 0.05
MIDO_1_fst[21,1] <- 0.06; MIDO_1_fst[21,2] <- 0.04; MIDO_1_fst[21,3] <- 0.17; MIDO_1_fst[21,4] <- 0.21; MIDO_1_fst[21,5] <- 0.29; MIDO_1_fst[21,6] <- 0.03; MIDO_1_fst[21,7] <- 0.04; MIDO_1_fst[21,8] <- 0.20; MIDO_1_fst[21,9] <- 0.07; MIDO_1_fst[21,10] <- 0.20; MIDO_1_fst[21,11] <- 0.12; MIDO_1_fst[21,12] <- 0.05; MIDO_1_fst[21,13] <- 0.14; MIDO_1_fst[21,14] <- 0.17; MIDO_1_fst[21,15] <- 0.09; MIDO_1_fst[21,16] <- 0.02; MIDO_1_fst[21,17] <- 0.17; MIDO_1_fst[21,18] <- 0.16; MIDO_1_fst[21,19] <- 0.14; MIDO_1_fst[21,20] <- 0.14
MIDO_1_fst[22,1] <- 0.20; MIDO_1_fst[22,2] <- 0.18; MIDO_1_fst[22,3] <- 0.05; MIDO_1_fst[22,4] <- 0.04; MIDO_1_fst[22,5] <- 0.08; MIDO_1_fst[22,6] <- 0.20; MIDO_1_fst[22,7] <- 0.20; MIDO_1_fst[22,8] <- 0.04; MIDO_1_fst[22,9] <- 0.22; MIDO_1_fst[22,10] <- 0.04; MIDO_1_fst[22,11] <- 0.12; MIDO_1_fst[22,12] <- 0.18; MIDO_1_fst[22,13] <- 0.29; MIDO_1_fst[22,14] <- 0.07; MIDO_1_fst[22,15] <- 0.25; MIDO_1_fst[22,16] <- 0.20; MIDO_1_fst[22,17] <- 0.06; MIDO_1_fst[22,18] <- 0.07; MIDO_1_fst[22,19] <- 0.06; MIDO_1_fst[22,20] <- 0.09; MIDO_1_fst[22,21] <- 0.23
MIDO_1_fst[23,1] <- 0.13; MIDO_1_fst[23,2] <- 0.14; MIDO_1_fst[23,3] <- 0.25; MIDO_1_fst[23,4] <- 0.28; MIDO_1_fst[23,5] <- 0.35; MIDO_1_fst[23,6] <- 0.11; MIDO_1_fst[23,7] <- 0.12; MIDO_1_fst[23,8] <- 0.26; MIDO_1_fst[23,9] <- 0.12; MIDO_1_fst[23,10] <- 0.25; MIDO_1_fst[23,11] <- 0.19; MIDO_1_fst[23,12] <- 0.14; MIDO_1_fst[23,13] <- 0.24; MIDO_1_fst[23,14] <- 0.23; MIDO_1_fst[23,15] <- 0.20; MIDO_1_fst[23,16] <- 0.13; MIDO_1_fst[23,17] <- 0.26; MIDO_1_fst[23,18] <- 0.21; MIDO_1_fst[23,19] <- 0.19; MIDO_1_fst[23,20] <- 0.15; MIDO_1_fst[23,21] <- 0.16; MIDO_1_fst[23,22] <- 0.30
MIDO_1_fst[24,1] <- 0.07; MIDO_1_fst[24,2] <- 0.07; MIDO_1_fst[24,3] <- 0.06; MIDO_1_fst[24,4] <- 0.08; MIDO_1_fst[24,5] <- 0.14; MIDO_1_fst[24,6] <- 0.06; MIDO_1_fst[24,7] <- 0.05; MIDO_1_fst[24,8] <- 0.07; MIDO_1_fst[24,9] <- 0.11; MIDO_1_fst[24,10] <- 0.06; MIDO_1_fst[24,11] <- 0.07; MIDO_1_fst[24,12] <- 0.07; MIDO_1_fst[24,13] <- 0.19; MIDO_1_fst[24,14] <- 0.08; MIDO_1_fst[24,15] <- 0.14; MIDO_1_fst[24,16] <- 0.07; MIDO_1_fst[24,17] <- 0.08; MIDO_1_fst[24,18] <- 0.07; MIDO_1_fst[24,19] <- 0.03; MIDO_1_fst[24,20] <- 0.06; MIDO_1_fst[24,21] <- 0.10; MIDO_1_fst[24,22] <- 0.09; MIDO_1_fst[24,23] <- 0.18
MIDO_1_fst[25,1] <- 0.36; MIDO_1_fst[25,2] <- 0.33; MIDO_1_fst[25,3] <- 0.50; MIDO_1_fst[25,4] <- 0.49; MIDO_1_fst[25,5] <- 0.55; MIDO_1_fst[25,6] <- 0.35; MIDO_1_fst[25,7] <- 0.36; MIDO_1_fst[25,8] <- 0.43; MIDO_1_fst[25,9] <- 0.34; MIDO_1_fst[25,10] <- 0.47; MIDO_1_fst[25,11] <- 0.42; MIDO_1_fst[25,12] <- 0.30; MIDO_1_fst[25,13] <- 0.48; MIDO_1_fst[25,14] <- 0.44; MIDO_1_fst[25,15] <- 0.45; MIDO_1_fst[25,16] <- 0.36; MIDO_1_fst[25,17] <- 0.51; MIDO_1_fst[25,18] <- 0.36; MIDO_1_fst[25,19] <- 0.40; MIDO_1_fst[25,20] <- 0.40; MIDO_1_fst[25,21] <- 0.37; MIDO_1_fst[25,22] <- 0.52; MIDO_1_fst[25,23] <- 0.34; MIDO_1_fst[25,24] <- 0.46
MIDO_1_fst[26,1] <- 0.09; MIDO_1_fst[26,2] <- 0.09; MIDO_1_fst[26,3] <- 0.06; MIDO_1_fst[26,4] <- 0.09; MIDO_1_fst[26,5] <- 0.15; MIDO_1_fst[26,6] <- 0.08; MIDO_1_fst[26,7] <- 0.08; MIDO_1_fst[26,8] <- 0.08; MIDO_1_fst[26,9] <- 0.13; MIDO_1_fst[26,10] <- 0.06; MIDO_1_fst[26,11] <- 0.05; MIDO_1_fst[26,12] <- 0.07; MIDO_1_fst[26,13] <- 0.22; MIDO_1_fst[26,14] <- 0.07; MIDO_1_fst[26,15] <- 0.17; MIDO_1_fst[26,16] <- 0.09; MIDO_1_fst[26,17] <- 0.10; MIDO_1_fst[26,18] <- 0.06; MIDO_1_fst[26,19] <- 0.02; MIDO_1_fst[26,20] <- 0.05; MIDO_1_fst[26,21] <- 0.12; MIDO_1_fst[26,22] <- 0.10; MIDO_1_fst[26,23] <- 0.18; MIDO_1_fst[26,24] <- 0.01; MIDO_1_fst[26,25] <- 0.45
MIDO_1_fst[27,1] <- 0.19; MIDO_1_fst[27,2] <- 0.17; MIDO_1_fst[27,3] <- 0.05; MIDO_1_fst[27,4] <- 0.03; MIDO_1_fst[27,5] <- 0.09; MIDO_1_fst[27,6] <- 0.19; MIDO_1_fst[27,7] <- 0.19; MIDO_1_fst[27,8] <- 0.05; MIDO_1_fst[27,9] <- 0.22; MIDO_1_fst[27,10] <- 0.03; MIDO_1_fst[27,11] <- 0.11; MIDO_1_fst[27,12] <- 0.17; MIDO_1_fst[27,13] <- 0.28; MIDO_1_fst[27,14] <- 0.07; MIDO_1_fst[27,15] <- 0.24; MIDO_1_fst[27,16] <- 0.19; MIDO_1_fst[27,17] <- 0.07; MIDO_1_fst[27,18] <- 0.09; MIDO_1_fst[27,19] <- 0.07; MIDO_1_fst[27,20] <- 0.10; MIDO_1_fst[27,21] <- 0.22; MIDO_1_fst[27,22] <- 0.00; MIDO_1_fst[27,23] <- 0.29; MIDO_1_fst[27,24] <- 0.08; MIDO_1_fst[27,25] <- 0.49; MIDO_1_fst[27,26] <- 0.09
MIDO_1_fst[28,1] <- 0.14; MIDO_1_fst[28,2] <- 0.12; MIDO_1_fst[28,3] <- 0.03; MIDO_1_fst[28,4] <- 0.04; MIDO_1_fst[28,5] <- 0.09; MIDO_1_fst[28,6] <- 0.13; MIDO_1_fst[28,7] <- 0.13; MIDO_1_fst[28,8] <- 0.04; MIDO_1_fst[28,9] <- 0.16; MIDO_1_fst[28,10] <- 0.01; MIDO_1_fst[28,11] <- 0.06; MIDO_1_fst[28,12] <- 0.11; MIDO_1_fst[28,13] <- 0.23; MIDO_1_fst[28,14] <- 0.01; MIDO_1_fst[28,15] <- 0.19; MIDO_1_fst[28,16] <- 0.13; MIDO_1_fst[28,17] <- 0.05; MIDO_1_fst[28,18] <- 0.06; MIDO_1_fst[28,19] <- 0.03; MIDO_1_fst[28,20] <- 0.05; MIDO_1_fst[28,21] <- 0.17; MIDO_1_fst[28,22] <- 0.04; MIDO_1_fst[28,23] <- 0.21; MIDO_1_fst[28,24] <- 0.04; MIDO_1_fst[28,25] <- 0.44; MIDO_1_fst[28,26] <- 0.05; MIDO_1_fst[28,27] <- 0.03
MIDO_1_fst[29,1] <- 0.10; MIDO_1_fst[29,2] <- 0.07; MIDO_1_fst[29,3] <- 0.06; MIDO_1_fst[29,4] <- 0.10; MIDO_1_fst[29,5] <- 0.13; MIDO_1_fst[29,6] <- 0.07; MIDO_1_fst[29,7] <- 0.08; MIDO_1_fst[29,8] <- 0.10; MIDO_1_fst[29,9] <- 0.11; MIDO_1_fst[29,10] <- 0.09; MIDO_1_fst[29,11] <- 0.07; MIDO_1_fst[29,12] <- 0.07; MIDO_1_fst[29,13] <- 0.17; MIDO_1_fst[29,14] <- 0.07; MIDO_1_fst[29,15] <- 0.14; MIDO_1_fst[29,16] <- 0.07; MIDO_1_fst[29,17] <- 0.08; MIDO_1_fst[29,18] <- 0.06; MIDO_1_fst[29,19] <- 0.06; MIDO_1_fst[29,20] <- 0.05; MIDO_1_fst[29,21] <- 0.11; MIDO_1_fst[29,22] <- 0.11; MIDO_1_fst[29,23] <- 0.17; MIDO_1_fst[29,24] <- 0.03; MIDO_1_fst[29,25] <- 0.37; MIDO_1_fst[29,26] <- 0.04; MIDO_1_fst[29,27] <- 0.12; MIDO_1_fst[29,28] <- 0.06
MIDO_1_fst[30,1] <- 0.10; MIDO_1_fst[30,2] <- 0.09; MIDO_1_fst[30,3] <- 0.18; MIDO_1_fst[30,4] <- 0.21; MIDO_1_fst[30,5] <- 0.27; MIDO_1_fst[30,6] <- 0.08; MIDO_1_fst[30,7] <- 0.08; MIDO_1_fst[30,8] <- 0.19; MIDO_1_fst[30,9] <- 0.09; MIDO_1_fst[30,10] <- 0.19; MIDO_1_fst[30,11] <- 0.12; MIDO_1_fst[30,12] <- 0.08; MIDO_1_fst[30,13] <- 0.16; MIDO_1_fst[30,14] <- 0.17; MIDO_1_fst[30,15] <- 0.09; MIDO_1_fst[30,16] <- 0.06; MIDO_1_fst[30,17] <- 0.17; MIDO_1_fst[30,18] <- 0.16; MIDO_1_fst[30,19] <- 0.12; MIDO_1_fst[30,20] <- 0.14; MIDO_1_fst[30,21] <- 0.07; MIDO_1_fst[30,22] <- 0.23; MIDO_1_fst[30,23] <- 0.19; MIDO_1_fst[30,24] <- 0.12; MIDO_1_fst[30,25] <- 0.42; MIDO_1_fst[30,26] <- 0.15; MIDO_1_fst[30,27] <- 0.22; MIDO_1_fst[30,28] <- 0.17; MIDO_1_fst[30,29] <- 0.11
MIDO_1_fst[31,1] <- 0.23; MIDO_1_fst[31,2] <- 0.20; MIDO_1_fst[31,3] <- 0.04; MIDO_1_fst[31,4] <- 0.02; MIDO_1_fst[31,5] <- 0.07; MIDO_1_fst[31,6] <- 0.21; MIDO_1_fst[31,7] <- 0.22; MIDO_1_fst[31,8] <- 0.07; MIDO_1_fst[31,9] <- 0.24; MIDO_1_fst[31,10] <- 0.03; MIDO_1_fst[31,11] <- 0.13; MIDO_1_fst[31,12] <- 0.18; MIDO_1_fst[31,13] <- 0.29; MIDO_1_fst[31,14] <- 0.07; MIDO_1_fst[31,15] <- 0.26; MIDO_1_fst[31,16] <- 0.21; MIDO_1_fst[31,17] <- 0.08; MIDO_1_fst[31,18] <- 0.08; MIDO_1_fst[31,19] <- 0.08; MIDO_1_fst[31,20] <- 0.10; MIDO_1_fst[31,21] <- 0.24; MIDO_1_fst[31,22] <- 0.02; MIDO_1_fst[31,23] <- 0.30; MIDO_1_fst[31,24] <- 0.10; MIDO_1_fst[31,25] <- 0.48; MIDO_1_fst[31,26] <- 0.09; MIDO_1_fst[31,27] <- 0.03; MIDO_1_fst[31,28] <- 0.05; MIDO_1_fst[31,29] <- 0.11; MIDO_1_fst[31,30] <- 0.24
MIDO_1_fst[32,1] <- 0.15; MIDO_1_fst[32,2] <- 0.15; MIDO_1_fst[32,3] <- 0.16; MIDO_1_fst[32,4] <- 0.19; MIDO_1_fst[32,5] <- 0.19; MIDO_1_fst[32,6] <- 0.14; MIDO_1_fst[32,7] <- 0.14; MIDO_1_fst[32,8] <- 0.17; MIDO_1_fst[32,9] <- 0.14; MIDO_1_fst[32,10] <- 0.14; MIDO_1_fst[32,11] <- 0.15; MIDO_1_fst[32,12] <- 0.14; MIDO_1_fst[32,13] <- 0.28; MIDO_1_fst[32,14] <- 0.12; MIDO_1_fst[32,15] <- 0.21; MIDO_1_fst[32,16] <- 0.15; MIDO_1_fst[32,17] <- 0.14; MIDO_1_fst[32,18] <- 0.12; MIDO_1_fst[32,19] <- 0.12; MIDO_1_fst[32,20] <- 0.09; MIDO_1_fst[32,21] <- 0.19; MIDO_1_fst[32,22] <- 0.19; MIDO_1_fst[32,23] <- 0.21; MIDO_1_fst[32,24] <- 0.11; MIDO_1_fst[32,25] <- 0.52; MIDO_1_fst[32,26] <- 0.12; MIDO_1_fst[32,27] <- 0.18; MIDO_1_fst[32,28] <- 0.12; MIDO_1_fst[32,29] <- 0.08; MIDO_1_fst[32,30] <- 0.18; MIDO_1_fst[32,31] <- 0.19

MIDO_1_fst[upper.tri(MIDO_1_fst)] <- t(MIDO_1_fst)[upper.tri(MIDO_1_fst)]
MIDO_1_fst[MIDO_1_fst < 0] <- 0
diag(MIDO_1_fst) <- 0

# -----------------------------
# pairwise geographic distances for IBD
# -----------------------------
coords <- MIDO_1_sites %>% select(lon, lat)
geo_dist_km <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
rownames(geo_dist_km) <- colnames(geo_dist_km) <- MIDO_1_sites$site

ibd_df <- data.frame(
  site1   = rownames(MIDO_1_fst)[row(MIDO_1_fst)[upper.tri(MIDO_1_fst)]],
  site2   = colnames(MIDO_1_fst)[col(MIDO_1_fst)[upper.tri(MIDO_1_fst)]],
  fst     = MIDO_1_fst[upper.tri(MIDO_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# map of sampling sites
# -----------------------------
world_map <- map_data("world")
usa_map   <- subset(world_map, region %in% c("USA", "Canada"))

map_df <- MIDO_1_sites %>%
  mutate(label = as.character(site))

p_map <- ggplot() +
  geom_polygon(
    data = usa_map,
    aes(x = long, y = lat, group = group),
    fill = "grey92",
    color = "grey55",
    linewidth = 0.2
  ) +
  geom_point(
    data = map_df,
    aes(x = lon, y = lat),
    color = "red3",
    size = 3
  ) +
  geom_text(
    data = map_df,
    aes(x = lon, y = lat, label = label),
    nudge_y = 0.12,
    size = 3.0
  ) +
  coord_fixed(
    ratio = 1.3,
    xlim = range(map_df$lon) + c(-0.6, 0.6),
    ylim = range(map_df$lat) + c(-0.4, 0.4)
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "MIDO-1 sampling sites"
  )

# -----------------------------
# IBD plot
# -----------------------------
p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Geographic distance (km)",
    y = expression(F[ST]),
    title = "MIDO-1 isolation by distance"
  )

# -----------------------------
# save RData
# -----------------------------
save(
  MIDO_1_fst,
  MIDO_1_sites,
  file = "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/MIDO-1/data/MIDO-1.RData"
)

# -----------------------------
# print outputs
# -----------------------------
print(MIDO_1_sites)
print(round(MIDO_1_fst, 4))
print(head(ibd_df))
print(p_map)
print(p_ibd)
