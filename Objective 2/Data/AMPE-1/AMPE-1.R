sand_darter_sites <- data.frame(
  site = c("ER","EF","BC","DC","RED","LK","LM","HR","SC","MA","Syd","TH","GR","RAS","RR","CC"),
  name = c(
    "Eel River","East Fork White River","Big Creek","Deer Creek","Red River","Licking River","Little Muskingum River",
    "Hocking River main","Salt Creek","Maumee River main","Sydenham River","Thames River upper","Grand River upper",
    "Riviere au Saumon","Richelieu River","Champlain Canal"
  ),
  lat = c(
    40.828056,
    39.138611,
    38.809167,
    39.500556,
    37.819722,
    38.208333,
    39.411667,
    39.300833,
    39.433333,
    41.084167,
    42.646944,
    42.931944,
    43.127778,
    44.999167,
    45.635000,
    43.352500
  ),
  lon = c(
    -86.113889,
    -85.893889,
    -85.643889,
    -86.930278,
    -83.575833,
    -83.680278,
    -81.358611,
    -81.963889,
    -82.680000,
    -85.019722,
    -82.009722,
    -81.426389,
    -80.199167,
    -74.510556,
    -73.190556,
    -73.495556
  )
)


#---
# AI
#

library(ggplot2)
library(ggrepel)

# -----------------------------
# 1) site list in decimal degrees
# -----------------------------
sand_darter_sites <- data.frame(
  site = c("ER","EF","BC","DC","LK","RED","LM","HR","SC","MA","Syd","TH","GR","RAS","RR","CC"),
  full_name = c(
    "Eel River",
    "East Fork White River",
    "Big Creek",
    "Deer Creek",
    "Licking River",
    "Red River",
    "Little Muskingum River (site 1)",
    "Hocking River main (site 1)",
    "Salt Creek (site 1)",
    "Maumee River main (site 1)",
    "Sydenham River",
    "Thames River upper1",
    "Grand River upper1",
    "Riviere au Saumon",
    "Richelieu River1",
    "Champlain Canal"
  ),
  lat = c(
    40.828056,
    39.138611,
    38.809167,
    39.500556,
    38.208333,
    37.819722,
    39.411667,
    39.300833,
    39.433333,
    41.084167,
    42.646944,
    42.931944,
    43.127778,
    44.999167,
    45.635000,
    43.352500
  ),
  lon = c(
    -86.113889,
    -85.893889,
    -85.643889,
    -86.930278,
    -83.680278,
    -83.575833,
    -81.358611,
    -81.963889,
    -82.680000,
    -85.019722,
    -82.009722,
    -81.426389,
    -80.199167,
    -74.510556,
    -73.190556,
    -73.495556
  )
)

sand_darter_sites$label <- paste0(sand_darter_sites$site, " — ", sand_darter_sites$full_name)

sand_darter_sites

# -----------------------------
# 2) plot the sites
# -----------------------------
ggplot(sand_darter_sites, aes(x = lon, y = lat)) +
  geom_point(size = 2.8) +
  geom_text_repel(
    aes(label = label),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.2
  ) +
  coord_fixed(1.15) +
  labs(
    x = "Longitude",
    y = "Latitude",
    #title = "Eastern sand darter populations in the FST matrix"
  ) +
  theme_bw()

# -----------------------------
# 3) reconstruct the FST matrix
#    (values below diagonal from screenshot)
# -----------------------------
pops <- c("ER","EF","BC","DC","LK","RED","LM","HR","SC","MA","Syd","TH","GR","RAS","RR","CC")

fst_mat <- matrix(0, nrow = length(pops), ncol = length(pops),
                  dimnames = list(pops, pops))

# fill lower triangle row by row
fst_mat["EF","ER"]   <- 0.075

fst_mat["BC","ER"]   <- 0.085
fst_mat["BC","EF"]   <- 0.011

fst_mat["DC","ER"]   <- 0.076
fst_mat["DC","EF"]   <- 0.024
fst_mat["DC","BC"]   <- 0.009

fst_mat["LK","ER"]   <- 0.160
fst_mat["LK","EF"]   <- 0.103
fst_mat["LK","BC"]   <- 0.081
fst_mat["LK","DC"]   <- 0.078

fst_mat["RED","ER"]  <- 0.144
fst_mat["RED","EF"]  <- 0.089
fst_mat["RED","BC"]  <- 0.069
fst_mat["RED","DC"]  <- 0.063
fst_mat["RED","LK"]  <- 0.032

fst_mat["LM","ER"]   <- 0.103
fst_mat["LM","EF"]   <- 0.072
fst_mat["LM","BC"]   <- 0.063
fst_mat["LM","DC"]   <- 0.042
fst_mat["LM","LK"]   <- 0.075
fst_mat["LM","RED"]  <- 0.049

fst_mat["HR","ER"]   <- 0.164
fst_mat["HR","EF"]   <- 0.119
fst_mat["HR","BC"]   <- 0.085
fst_mat["HR","DC"]   <- 0.073
fst_mat["HR","LK"]   <- 0.080
fst_mat["HR","RED"]  <- 0.046
fst_mat["HR","LM"]   <- 0.053

fst_mat["SC","ER"]   <- 0.153
fst_mat["SC","EF"]   <- 0.139
fst_mat["SC","BC"]   <- 0.123
fst_mat["SC","DC"]   <- 0.112
fst_mat["SC","LK"]   <- 0.069
fst_mat["SC","RED"]  <- 0.060
fst_mat["SC","LM"]   <- 0.075
fst_mat["SC","HR"]   <- 0.081

fst_mat["MA","ER"]   <- 0.081
fst_mat["MA","EF"]   <- 0.047
fst_mat["MA","BC"]   <- 0.058
fst_mat["MA","DC"]   <- 0.077
fst_mat["MA","LK"]   <- 0.148
fst_mat["MA","RED"]  <- 0.145
fst_mat["MA","LM"]   <- 0.120
fst_mat["MA","HR"]   <- 0.165
fst_mat["MA","SC"]   <- 0.162

fst_mat["Syd","ER"]  <- 0.062
fst_mat["Syd","EF"]  <- 0.071
fst_mat["Syd","BC"]  <- 0.084
fst_mat["Syd","DC"]  <- 0.084
fst_mat["Syd","LK"]  <- 0.172
fst_mat["Syd","RED"] <- 0.159
fst_mat["Syd","LM"]  <- 0.121
fst_mat["Syd","HR"]  <- 0.175
fst_mat["Syd","SC"]  <- 0.154
fst_mat["Syd","MA"]  <- 0.054

fst_mat["TH","ER"]   <- 0.053
fst_mat["TH","EF"]   <- 0.047
fst_mat["TH","BC"]   <- 0.054
fst_mat["TH","DC"]   <- 0.053
fst_mat["TH","LK"]   <- 0.123
fst_mat["TH","RED"]  <- 0.110
fst_mat["TH","LM"]   <- 0.083
fst_mat["TH","HR"]   <- 0.126
fst_mat["TH","SC"]   <- 0.134
fst_mat["TH","MA"]   <- 0.050
fst_mat["TH","Syd"]  <- 0.021

fst_mat["GR","ER"]   <- 0.099
fst_mat["GR","EF"]   <- 0.077
fst_mat["GR","BC"]   <- 0.090
fst_mat["GR","DC"]   <- 0.088
fst_mat["GR","LK"]   <- 0.156
fst_mat["GR","RED"]  <- 0.149
fst_mat["GR","LM"]   <- 0.109
fst_mat["GR","HR"]   <- 0.168
fst_mat["GR","SC"]   <- 0.165
fst_mat["GR","MA"]   <- 0.090
fst_mat["GR","Syd"]  <- 0.044
fst_mat["GR","TH"]   <- 0.055

fst_mat["RAS","ER"]  <- 0.114
fst_mat["RAS","EF"]  <- 0.070
fst_mat["RAS","BC"]  <- 0.056
fst_mat["RAS","DC"]  <- 0.060
fst_mat["RAS","LK"]  <- 0.159
fst_mat["RAS","RED"] <- 0.147
fst_mat["RAS","LM"]  <- 0.115
fst_mat["RAS","HR"]  <- 0.130
fst_mat["RAS","SC"]  <- 0.171
fst_mat["RAS","MA"]  <- 0.096
fst_mat["RAS","Syd"] <- 0.116
fst_mat["RAS","TH"]  <- 0.081
fst_mat["RAS","GR"]  <- 0.105

fst_mat["RR","ER"]   <- 0.148
fst_mat["RR","EF"]   <- 0.096
fst_mat["RR","BC"]   <- 0.098
fst_mat["RR","DC"]   <- 0.086
fst_mat["RR","LK"]   <- 0.184
fst_mat["RR","RED"]  <- 0.170
fst_mat["RR","LM"]   <- 0.118
fst_mat["RR","HR"]   <- 0.146
fst_mat["RR","SC"]   <- 0.190
fst_mat["RR","MA"]   <- 0.125
fst_mat["RR","Syd"]  <- 0.143
fst_mat["RR","TH"]   <- 0.098
fst_mat["RR","GR"]   <- 0.093
fst_mat["RR","RAS"]  <- 0.060

fst_mat["CC","ER"]   <- 0.259
fst_mat["CC","EF"]   <- 0.170
fst_mat["CC","BC"]   <- 0.184
fst_mat["CC","DC"]   <- 0.193
fst_mat["CC","LK"]   <- 0.279
fst_mat["CC","RED"]  <- 0.267
fst_mat["CC","LM"]   <- 0.224
fst_mat["CC","HR"]   <- 0.237
fst_mat["CC","SC"]   <- 0.281
fst_mat["CC","MA"]   <- 0.243
fst_mat["CC","Syd"]  <- 0.289
fst_mat["CC","TH"]   <- 0.204
fst_mat["CC","GR"]   <- 0.205
fst_mat["CC","RAS"]  <- 0.155
fst_mat["CC","RR"]   <- 0.175

# mirror lower triangle to upper triangle
fst_mat[upper.tri(fst_mat)] <- t(fst_mat)[upper.tri(fst_mat)]

# set diagonal to NA to match typical pairwise FST display
diag(fst_mat) <- NA

fst_mat

# optional rounded print
round(fst_mat, 3)

# -----------------------------
# 4) optional save objects
# -----------------------------
# save(sand_darter_sites, fst_mat, file = "sand_darter_sites_fst.RData")
# write.csv(sand_darter_sites, "sand_darter_sites.csv", row.names = FALSE)
# write.csv(fst_mat, "sand_darter_fst_matrix.csv")

library(geosphere)
library(ggplot2)

# -----------------------------
# pairwise straight-line distance (km)
# -----------------------------
dist_mat <- geosphere::distm(
  x = sand_darter_sites[, c("lon", "lat")],
  fun = geosphere::distGeo
) / 1000

rownames(dist_mat) <- sand_darter_sites$site
colnames(dist_mat) <- sand_darter_sites$site

# reorder to match fst_mat if needed
dist_mat <- dist_mat[rownames(fst_mat), colnames(fst_mat)]

# -----------------------------
# pull lower triangle into a df
# -----------------------------
idx <- lower.tri(fst_mat)

ibd_df <- data.frame(
  site1 = rownames(fst_mat)[row(fst_mat)[idx]],
  site2 = colnames(fst_mat)[col(fst_mat)[idx]],
  fst = fst_mat[idx],
  distance_km = dist_mat[idx]
)

ibd_df <- ibd_df[complete.cases(ibd_df), ]

# inspect
ibd_df

# -----------------------------
# IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = distance_km, y = fst)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Straight-line distance (km)",
    y = expression(pairwise~F[ST]),
    title = "Isolation by distance"
  ) +
  theme_bw()