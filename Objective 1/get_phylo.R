# Get the phylogeny for the PGLMM
# cgh


#install/load pacakges
packages <- c("dplyr", "ape", "fishtree",
              "taxize", "phytools")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
# install new packages
if (length(new_packages)) {
  install.packages(new_packages, dependencies = TRUE)
}
# load all packages
invisible(lapply(packages, library, character.only = TRUE))



# load filtered data
dat <- read.csv("Data/final_fish_dataset.csv")

# pull FishTreee phylo for species in filtered dataset
phy <- fishtree_phylogeny(species = unique(dat$Spec_Latin_GenDivRange))

# get names of missing fish
not_in <- unique(dat$Spec_Latin_GenDivRange)[!(gsub(" ", "_", unique(dat$Spec_Latin_GenDivRange)) %in% phy$tip.label)]


###########################################################
# graft in species that are missing in FishTree phylogeny #
###########################################################

# Catostomus discobolus
# Divergence time from https://doi.org/10.7717/peerj.5168
phy <- bind.tip(phy,
                  tip.label = "Catostomus_discobolus",
                  where = which(phy$tip.label == "Catostomus_warnerensis"),
                edge.length = 15, position = 15)
plot(phy)

# Cyprinodon nevadensis
# Divergence time from https://doi.org/10.1643/CG-03-093R3
phy <- bind.tip(phy,
                tip.label = "Cyprinodon_nevadensis",
                where = which(phy$tip.label == "Cyprinodon_eremus"),
                edge.length = 5, position = 5)
plot(phy)

# Erimyzon claviformes
# div time from https://doi.org/10.7717/peerj.5168
phy <- bind.tip(phy,
                tip.label = "Erimyzon_claviformis",
                where = getMRCA(phy, c("Catostomus_discobolus", "Catostomus_warnerensis")),
                edge.length = 30+15, position = 30)

plot(phy)

# Fundulus heteroclitus
# Divergence time from https://doi.org/10.11646/zootaxa.4250.6.5
phy <- bind.tip(phy,
                tip.label = "Fundulus_heteroclitus",
                where = which(phy$tip.label == "Fundulus_grandis"),
                edge.length = 6, position = 6)
plot(phy)

# Salevelinus alpinus
# Divergence time from https://doi.org/10.1111/bij.12559
phy <- bind.tip(phy,
                tip.label = "Salvelinus_alpinus",
                where = which(phy$tip.label == "Salvelinus_confluentus"),
                edge.length = 3, position = 3)
plot(phy)


# write tree
write.tree(phy, file = "Data/complete_tree.nwk")


