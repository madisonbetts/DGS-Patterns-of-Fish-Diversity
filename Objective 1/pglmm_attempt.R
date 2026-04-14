##### PGLMM #####

# first shot
library(dplyr)
library(ape)
library(glmmTMB) 
library(bbmle) 
library(ggplot2)

# load in dataset, remove unnecessary columns
dat <- read.csv("Data/final_fish_dataset.csv")

dat_lm <- dat %>%
  select(
    Study_id, Spec_Latin_GenDivRange, N_pops, N_loci, 
         Order_GBIF, Family_GBIF, Habitat_adulthood, Habitat_breeding, 
         Habitat_fishbase, Pop_id, N, A, A_mean, A_tot, A_eff, A_private, 
         N_genot, Ar, Ho, He, GD_Nei, D_clonal, F_is, F_is_sig, Ploidy,
         snapped_lon, snapped_lat, stream_order, catch_area, slope, elev,
         WSA9, Ice, Tectonic, SeaLevel, SeaShort, Topography, History,
         AREAKM2, PATCHES, FECUNDITY, LONGEVITY, MAXTL
    ) %>%
  mutate(
    Anadromy = ifelse(Habitat_adulthood == "Marine", 1, 0)
  )

dat_lm$Spec_Latin_GenDivRange <- gsub(" ", "_", dat_lm$Spec_Latin_GenDivRange)

# testing for multicolinearity
lm_vif <- lm(
  He ~ N + N_loci + N_pops +
    slope + elev + catch_area + stream_order +
    Ice + Tectonic + Topography + SeaLevel +
    AREAKM2 + PATCHES +
    FECUNDITY + LONGEVITY + MAXTL + Anadromy +
    snapped_lat + snapped_lon,
  data = dat_lm
)

car::vif(lm_vif) #snapped_lon & tectonism -> removed snapped_lon

# load phylogeny
phylo <- read.tree("Data/complete_tree.nwk")

# phylogenetic variance-covariance matrix
phylo_varcov <- vcv(phylo[[1]])

dat_lm$Spec_Latin_GenDivRange <- factor(
  dat_lm$Spec_Latin_GenDivRange,
  levels = rownames(phylo_varcov)
)

dat_lm$phylo_all <- factor("all")


#################### NULL ########################
null_model <- glmmTMB(
  He ~ 1 +
    (1 | Study_id) +
    (1 | Spec_Latin_GenDivRange) +
    propto(0 + Spec_Latin_GenDivRange | phylo_all, phylo_varcov),
  family = beta_family(link = "logit"),
  data = dat_lm
)




# He glmm + full predictor set
lm1 <- glmmTMB(
  He ~ 
    # local factors
    slope + elev + catch_area + stream_order +
    # historical factors
    History +
    # life history
    #AREAKM2 + PATCHES +
    #FECUNDITY + LONGEVITY + MAXTL + Anadromy +
    snapped_lat +
    
    # random effects
    (1 | Study_id) +
    (1 | Spec_Latin_GenDivRange) +
    propto(0 + Spec_Latin_GenDivRange | phylo_all, phylo_varcov),
  family = beta_family(link = "logit"),
  data = dat_lm
)
summary(lm1)
################### trying phyr###############################

library(phyr)
library(INLA)
library(MASS)


phyr_pglmm_full <- pglmm(
  He ~ 
    # local factors
    slope + elev + catch_area + stream_order +
    
    # historical factors
    History + snapped_lat +
    
    # life history
    AREAKM2 + PATCHES +
    FECUNDITY + LONGEVITY + MAXTL + Anadromy
    
    # random effects
    (1 | Study_id) +
    (1 | Spec_Latin_GenDivRange__), 
  data = dat_lm, 
  family = "gaussian", 
  cov_ranef = list(Spec_Latin_GenDivRange = phylo[[1]])
)
summary(phyr_pglmm_full)




phyr_pglmm_loc <- pglmm(
  He ~ 
    # local factors
    slope + elev + catch_area + stream_order +
    
    # random effects
    (1 | Study_id) +
    (1 | Spec_Latin_GenDivRange__), 
  data = dat_lm, 
  family = "gaussian", 
  cov_ranef = list(Spec_Latin_GenDivRange = phylo[[1]])
)
summary(phyr_pglmm_loc)


phyr_pglmm_hist <- pglmm(
  He ~ 
    # historical factors
    History + snapped_lat +
    
    # random effects
    (1 | Study_id) +
    (1 | Spec_Latin_GenDivRange__), 
  data = dat_lm, 
  family = "gaussian", 
  cov_ranef = list(Spec_Latin_GenDivRange = phylo[[1]])
)
summary(phyr_pglmm_hist)


phyr_pglmm_lh <- pglmm(
  He ~ 
    # life history
    FECUNDITY + LONGEVITY + MAXTL + Anadromy +
    
    # random effects
    (1 | Study_id) +
    (1 | Spec_Latin_GenDivRange__), 
  data = dat_lm, 
  family = "gaussian", 
  cov_ranef = list(Spec_Latin_GenDivRange = phylo[[1]])
)
summary(phyr_pglmm_lh)


  
phyr_pglmm_bg <- pglmm(
  He ~ 
    # biogeography
    AREAKM2 + PATCHES +
      
    # random effects
    (1 | Study_id) +
    (1 | Spec_Latin_GenDivRange__), 
  data = dat_lm, 
  family = "gaussian", 
  cov_ranef = list(Spec_Latin_GenDivRange = phylo[[1]])
)
summary(phyr_pglmm_bg)
