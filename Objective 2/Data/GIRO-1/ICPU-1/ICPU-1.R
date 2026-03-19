#========================================
# Ictalurus punctatus
# Wabash + Ohio rivers
# pairwise FST + IBD workflow
#========================================

library(ggplot2)
library(vegan)

#------------------------------
# 1) site coordinates
#------------------------------
sites <- data.frame(
  site = c("WA_19", "WA_136", "WA_189", "WA_268",
           "OH_1478", "OH_1511", "OH_1548"),
  lat = c(37.90736, 38.29545, 38.64305, 39.12247,
          37.23654, 37.12979, 37.22777),
  lon = c(-88.09838, -87.83276, -87.6168, -87.64498,
          -88.47449, -88.42943, -88.95386),
  stringsAsFactors = FALSE
)

sites

#------------------------------
# 2) pairwise FST matrix
# values above diagonal from the table
#------------------------------
site_order <- c("WA_268", "WA_189", "WA_136", "WA_19",
                "OH_1478", "OH_1511", "OH_1548")

fst_matrix <- matrix(NA, nrow = 7, ncol = 7,
                     dimnames = list(site_order, site_order))

diag(fst_matrix) <- 0

fst_matrix["WA_268", "WA_189"]   <- 0.009
fst_matrix["WA_268", "WA_136"]   <- 0.072
fst_matrix["WA_268", "WA_19"]    <- 0.123
fst_matrix["WA_268", "OH_1478"]  <- 0.173
fst_matrix["WA_268", "OH_1511"]  <- 0.157
fst_matrix["WA_268", "OH_1548"]  <- 0.184

fst_matrix["WA_189", "WA_136"]   <- 0.029
fst_matrix["WA_189", "WA_19"]    <- 0.080
fst_matrix["WA_189", "OH_1478"]  <- 0.142
fst_matrix["WA_189", "OH_1511"]  <- 0.102
fst_matrix["WA_189", "OH_1548"]  <- 0.123

fst_matrix["WA_136", "WA_19"]    <- 0.018
fst_matrix["WA_136", "OH_1478"]  <- 0.069
fst_matrix["WA_136", "OH_1511"]  <- 0.049
fst_matrix["WA_136", "OH_1548"]  <- 0.074

fst_matrix["WA_19", "OH_1478"]   <- 0.066
fst_matrix["WA_19", "OH_1511"]   <- 0.050
fst_matrix["WA_19", "OH_1548"]   <- 0.087

fst_matrix["OH_1478", "OH_1511"] <- 0.021
fst_matrix["OH_1478", "OH_1548"] <- 0.172

fst_matrix["OH_1511", "OH_1548"] <- 0.113

# mirror upper triangle to lower triangle
fst_matrix[lower.tri(fst_matrix)] <- t(fst_matrix)[lower.tri(fst_matrix)]

fst_matrix <- as.matrix(fst_matrix)
round(fst_matrix, 3)

#------------------------------
# 3) geographic distance matrix (km)
# values below diagonal from the table
#------------------------------
geo_matrix <- matrix(NA, nrow = 7, ncol = 7,
                     dimnames = list(site_order, site_order))

diag(geo_matrix) <- 0

geo_matrix["WA_189", "WA_268"]   <- 79

geo_matrix["WA_136", "WA_268"]   <- 132
geo_matrix["WA_136", "WA_189"]   <- 53

geo_matrix["WA_19", "WA_268"]    <- 249
geo_matrix["WA_19", "WA_189"]    <- 170
geo_matrix["WA_19", "WA_136"]    <- 117

geo_matrix["OH_1478", "WA_268"]  <- 379
geo_matrix["OH_1478", "WA_189"]  <- 299
geo_matrix["OH_1478", "WA_136"]  <- 246
geo_matrix["OH_1478", "WA_19"]   <- 129

geo_matrix["OH_1511", "WA_268"]  <- 412
geo_matrix["OH_1511", "WA_189"]  <- 332
geo_matrix["OH_1511", "WA_136"]  <- 279
geo_matrix["OH_1511", "WA_19"]   <- 162
geo_matrix["OH_1511", "OH_1478"] <- 33

geo_matrix["OH_1548", "WA_268"]  <- 449
geo_matrix["OH_1548", "WA_189"]  <- 369
geo_matrix["OH_1548", "WA_136"]  <- 316
geo_matrix["OH_1548", "WA_19"]   <- 199
geo_matrix["OH_1548", "OH_1478"] <- 70
geo_matrix["OH_1548", "OH_1511"] <- 37

# mirror lower triangle to upper triangle
geo_matrix[upper.tri(geo_matrix)] <- t(geo_matrix)[upper.tri(geo_matrix)]

geo_matrix <- as.matrix(geo_matrix)
geo_matrix

#------------------------------
# 4) pairwise dataframe for plotting
#------------------------------
ibd_df <- data.frame(
  site1 = rownames(fst_matrix)[row(fst_matrix)[lower.tri(fst_matrix)]],
  site2 = colnames(fst_matrix)[col(fst_matrix)[lower.tri(fst_matrix)]],
  distance_km = geo_matrix[lower.tri(geo_matrix)],
  fst = fst_matrix[lower.tri(fst_matrix)]
)

ibd_df

#------------------------------
# 5) IBD plot
#------------------------------
p <- ggplot(ibd_df, aes(x = distance_km, y = fst)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  labs(
    x = "Geographic distance (km)",
    y = expression(pairwise~F[ST]),
    title = "Ictalurus punctatus: isolation by distance",
    subtitle = "Wabash and Ohio river sites"
  )

print(p)

#------------------------------
# 6) simple correlation on pairwise values
#------------------------------
cor_test <- cor.test(ibd_df$distance_km, ibd_df$fst, method = "pearson")
cor_test

#------------------------------
# 7) Mantel test
#------------------------------
mantel_res <- mantel(
  as.dist(geo_matrix),
  as.dist(fst_matrix),
  method = "pearson",
  permutations = 9999
)

mantel_res

#------------------------------
# 8) optional: write outputs
#------------------------------
write.csv(sites,
          "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/datasets/ICTP-1/ictalurus_punctatus_sites.csv",
          row.names = FALSE)

write.csv(fst_matrix,
          "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/datasets/ICTP-1/ictalurus_punctatus_fst_matrix.csv",
          row.names = TRUE)

write.csv(geo_matrix,
          "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/datasets/ICTP-1/ictalurus_punctatus_geo_matrix_km.csv",
          row.names = TRUE)

save(sites, fst_matrix, geo_matrix, ibd_df, mantel_res,
     file = "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/datasets/ICTP-1/ictalurus_punctatus_ibd.RData")