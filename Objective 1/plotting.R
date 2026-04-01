#################  need to write

# install/load pacakges
packages <- c("sf", "dplyr", 
              "readxl", "nhdplusTools", 
              "units", "rnaturalearth",
              "ape", "fishtree", "colorspace",
              "taxize", "phytools", "tigris", "ggplot2")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
# install new packages
if (length(new_packages)) {
  install.packages(new_packages, dependencies = TRUE)
}
# load all packages
invisible(lapply(packages, library, character.only = TRUE))

sf::sf_use_s2(FALSE)

dat <- read.csv("Data/filtered_data.csv")

# some plots quick
He_lat <- ggplot(dat, aes(y, He)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Latitude", y = "Expected Heterozygosity (He)") +
  theme_classic()
ggsave(He_lat, filename = "Figures/He_lat.pdf", width = 7, height = 5, dpi = 300)

Ho_lat <- ggplot(dat, aes(y, Ho)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Latitude", y = "Observed Heterozygosity (Ho)") +
  theme_classic()
ggsave(Ho_lat, filename = "Figures/Ho_lat.pdf", width = 7, height = 5, dpi = 300)


# shape files for maps
states <- tigris::states(cb = TRUE, progress_bar = FALSE) |>
  filter(!STUSPS %in% c('HI', 'PR', 'AK', 'MP', 'GU', 'AS', 'VI')) |>
  st_transform(crs = 5070)
# plot


coords <- dat |>
  st_as_sf(coords = c("x", "y"), crs = 4326) |>
  st_transform(5070) |>
  st_coordinates()

dat$x_5070 <- coords[,1]
dat$y_5070 <- coords[,2]



He_gdr <- ggplot() +
  geom_sf(data = states, fill = "white") +
  geom_point(data = dat, aes(x = x_5070, y = y_5070, color = He))+
  scale_color_viridis_c(
    option = "magma",
    limits = c(min(dat$He), max(dat$He)),
    breaks = seq(0.3, 0.7, by = 0.1)
  ) +
  theme_void()
ggsave(He_gdr, filename = "Figures/He_gdr.pdf", width = 7, height = 5, dpi = 300)

# but dont include populations who do not have an Ho value
Ho_gdr <- ggplot() +
  geom_sf(data = states, fill = "white") +
  geom_point(data = dat %>% filter(!is.na(Ho)), aes(x = x_5070, y = y_5070, color = Ho))+
  scale_color_viridis_c(
    option = "magma",
    breaks = seq(0.3, 0.7, by = 0.1)) +
  theme_void()
ggsave(Ho_gdr, filename = "Figures/Ho_gdr.pdf", width = 7, height = 5, dpi = 300)

####################### pull phylogeny and summarize He/Ho ################################################


dat$Spec_Latin_GenDivRange <- gsub(" ", "_", dat$Spec_Latin_GenDivRange)

phy <- fishtree_phylogeny(species = unique(dat$Spec_Latin_GenDivRange))

#plot.phylo(phy, cex = 0.5)


########## add missing tips to tree following taxonomy ###########
# need to write



summary <-  dat %>% 
  group_by(Spec_Latin_GenDivRange) %>% 
  summarise(mean_He = mean(He), sd_He = sd(He), mean_Ho = mean(Ho), sd_Ho = sd(Ho)) %>% 
  filter(Spec_Latin_GenDivRange %in% phy$tip.label)

# get orders from taxize
summary$tax <- NA
for (i in 1:nrow(summary)){
  summary$tax[i] <- tax_name(summary$Spec_Latin_GenDivRange[i], get = "order", db = "ncbi")$order
}

# ensure the order of the summary matches the order of the tips in the phylogeny
summary <- summary[match(phy$tip.label, summary$Spec_Latin_GenDivRange),]

# pull taxonomy and species names
order_vec <- summary$tax
names(order_vec) <- phy$tip.label
all(names(order_vec) == phy$tip.label)


########################## plotting the phylogeny with He ####################################################

# prep data for plot
phy <- paintSubTree(phy, node = length(phy$tip.label) + 1, state = "other")

# setup colors - ensure "other" is black
orders <- unique(order_vec)
cols <- setNames(grDevices::hcl.colors(length(orders), "Dark3"), orders)
cols["other"] <- "black" 

# paint the tree
for (ord in orders) {
  tips_in_order <- names(order_vec)[order_vec == ord]
  
  if (length(tips_in_order) > 1) {
    mrca_node <- getMRCA(phy, tips_in_order)
    phy <- paintSubTree(phy, mrca_node, ord)
  } else {
    tip_node <- which(phy$tip.label == tips_in_order)
    phy <- paintBranches(phy, tip_node, ord)
  }
}

box_cols <- setNames(
  grDevices::hcl.colors(length(orders), "Dark3"),
  orders
)
box_cols <- box_cols[order_vec]
# named box_cols based on tip labels
names(box_cols) <- phy$tip.label


# boxplot
# create named he_vec for all obs in dat
he_vec_all <- setNames(dat$He, dat$Spec_Latin_GenDivRange)

# open a png for this plot
png("Figures/phylo_boxplot.png", width = 8, height = 6, units = "in", res = 600)

par(oma = c(0, 0, 0, 0)) 
plotTree.boxplot(phy, he_vec_all,
                 args.plotTree = list(ftype = "off", mar = c(5, 0.5, 4, 0)),
                 args.boxplot = list(col = box_cols,
                                     mar  = c(5, 0, 4, 1),  
                                     xlab = "He",
                                     ylim = c(0, 1)))
par(mfg=c(1,1))
plot(phy,cols,ftype = "off",mar=c(5.2, 0.5, 4.25, 0))



leg_cols <- cols[names(cols) != "other"]
# reverse order
leg_cols <- rev(leg_cols)
legend(x = "bottomleft", legend = names(leg_cols),
       pch = 22, pt.cex = 2, pt.bg = leg_cols,
       box.col = "transparent")

# switch to boxplot panel
par(mfg = c(1, 2))  # row 1, column 2 of mfrow


# get tip y positions
pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
yy <- pp$yy[1:length(phy$tip.label)]

#abline(h = yy, col = rgb(0,0,0,0.1), lwd = 1)

usr <- par("usr")  # plotting limits
xrange <- diff(usr[1:2])

for(i in seq_along(yy)) {
  if(i %% 2 == 0) {
    rect(
      xleft   = usr[1],
      xright  = usr[1] + 0.932 * xrange,
      ybottom = yy[i] - 0.5,
      ytop    = yy[i] + 0.5,
      col     = rgb(0,0,0,0.1),
      border  = NA
    )
  }
}

mtext(
  "Heterozygosity of Select North American \n Freshwater Fishes"
)

dev.off()

