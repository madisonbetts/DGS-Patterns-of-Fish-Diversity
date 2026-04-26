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
    Anadromy = ifelse(Habitat_adulthood == "Marine", 1, 0),
    AREAKM2 = ifelse(AREAKM2 %in% c(-555, -1, -999), NA_real_, AREAKM2),
    PATCHES = ifelse(PATCHES %in% c(-555, -1, -999), NA_real_, PATCHES),
    FECUNDITY = ifelse(FECUNDITY %in% c(-555, -1, -999), NA_real_, FECUNDITY),
    LONGEVITY = ifelse(LONGEVITY %in% c(-555, -1, -999), NA_real_, LONGEVITY),
    MAXTL = ifelse(MAXTL %in% c(-555, -1, -999), NA_real_, MAXTL)
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

lm_vif <- lm(
  He ~ N + N_loci + N_pops +
    slope + elev + catch_area + stream_order +
    Ice + Tectonic + Topography + SeaLevel +
    AREAKM2 + PATCHES +
    FECUNDITY + LONGEVITY + MAXTL + Anadromy +
    snapped_lat,
  data = dat_lm
)
car::vif(lm_vif)

# load phylogeny
phylo <- read.tree("Data/complete_tree.nwk")

# phylogenetic variance-covariance matrix
phylo_varcov <- vcv(phylo)

dat_lm$Spec_Latin_GenDivRange <- factor(
  dat_lm$Spec_Latin_GenDivRange,
  levels = rownames(phylo_varcov)
)

################### trying phyr###############################

library(phyr)
library(INLA)


#########################
# null
phyr_pglmm_null <- pglmm(
  He ~ 1 +
    # random effects
    (1 | Study_id) +
    (1 | Spec_Latin_GenDivRange__), 
  data = dat_lm, 
  family = "gaussian", 
  cov_ranef = list(Spec_Latin_GenDivRange = phylo), 
  REML = F
)
summary(phyr_pglmm_null)
#########################

#########################
m1 <- pglmm(
  He ~ 
    # local factors
    slope + elev + catch_area + stream_order +
    
    # historical factors
    History + snapped_lat +
    
    # Species Biogeography
    AREAKM2 + PATCHES +
    
    # Species Life History
    FECUNDITY + LONGEVITY + MAXTL + Anadromy
    
    # random effects
    (1 | Study_id) +
    (1 | Spec_Latin_GenDivRange__), 
  data = dat_lm, 
  family = "gaussian", 
  cov_ranef = list(Spec_Latin_GenDivRange = phylo),
  REML = F
)
summary(m1)
#########################

#########################
m2 <- pglmm(
  He ~ 
    # local factors
    slope + elev + catch_area + stream_order +
    
    # historical factors
    History + snapped_lat +
    
    # Species Biogeography
    AREAKM2 + PATCHES +
  
    # random effects
    (1 | Study_id) +
    (1 | Spec_Latin_GenDivRange__), 
  data = dat_lm, 
  family = "gaussian", 
  cov_ranef = list(Spec_Latin_GenDivRange = phylo),
  REML = F
)
summary(m2)
#########################

#########################
m3 <- pglmm(
  He ~ 
    # historical factors
    History + snapped_lat +
    
    # Species Biogeography
    AREAKM2 + PATCHES +
    
    # Species Life History
    FECUNDITY + LONGEVITY + MAXTL + Anadromy
  
    # random effects
    (1 | Study_id) +
    (1 | Spec_Latin_GenDivRange__), 
  data = dat_lm, 
  family = "gaussian", 
  cov_ranef = list(Spec_Latin_GenDivRange = phylo),
  REML = F
)
summary(m3)
#########################

#########################
m4 <- pglmm(
  He ~ 
    # local factors
    slope + elev + catch_area + stream_order +
        
    # Species Biogeography
    AREAKM2 + PATCHES +
    
    # Species Life History
    FECUNDITY + LONGEVITY + MAXTL + Anadromy
  
    # random effects
    (1 | Study_id) +
    (1 | Spec_Latin_GenDivRange__), 
  data = dat_lm, 
  family = "gaussian", 
  cov_ranef = list(Spec_Latin_GenDivRange = phylo), 
  REML = F
)
summary(m4)
#########################

#########################
m5 <- pglmm(
  He ~ 
    # local factors
    slope + elev + catch_area + stream_order +
    
    # historical factors
    History + snapped_lat +
    
    # Species Life History
    FECUNDITY + LONGEVITY + MAXTL + Anadromy
  
    # random effects
    (1 | Study_id) +
    (1 | Spec_Latin_GenDivRange__), 
  data = dat_lm, 
  family = "gaussian", 
  cov_ranef = list(Spec_Latin_GenDivRange = phylo), 
  REML = F
)
summary(m5)
#########################

#########################
m6 <- pglmm(
  He ~ 
    # local factors
    slope + elev + catch_area + stream_order +
    
    # historical factors
    History + snapped_lat +

    # random effects
    (1 | Study_id) +
    (1 | Spec_Latin_GenDivRange__), 
  data = dat_lm, 
  family = "gaussian", 
  cov_ranef = list(Spec_Latin_GenDivRange = phylo), 
  REML = F
)
summary(m6)
#########################

#########################
m7 <- pglmm(
  He ~ 
    # Species Biogeography
    AREAKM2 + PATCHES +
    
    # Species Life History
    FECUNDITY + LONGEVITY + MAXTL + Anadromy
  
    # random effects
    (1 | Study_id) +
    (1 | Spec_Latin_GenDivRange__), 
  data = dat_lm, 
  family = "gaussian", 
  cov_ranef = list(Spec_Latin_GenDivRange = phylo),
  REML = F
)
summary(m7)
#########################

#########################
m8 <- pglmm(
  He ~ 
    # local factors
    slope + elev + catch_area + stream_order +
    
    # random effects
    (1 | Study_id) +
    (1 | Spec_Latin_GenDivRange__), 
  data = dat_lm, 
  family = "gaussian", 
  cov_ranef = list(Spec_Latin_GenDivRange = phylo),
  REML = F
)
summary(m8)
#########################

#########################
m9 <- pglmm(
  He ~ 
    # historical factors
    History + snapped_lat +
    
    # random effects
    (1 | Study_id) +
    (1 | Spec_Latin_GenDivRange__), 
  data = dat_lm, 
  family = "gaussian", 
  cov_ranef = list(Spec_Latin_GenDivRange = phylo), 
  REML = F
)
summary(m9)

m10 <- pglmm(
  He ~ 
    # biogeography
    AREAKM2 + PATCHES +
    
    # random effects
    (1 | Study_id) +
    (1 | Spec_Latin_GenDivRange__), 
  data = dat_lm, 
  family = "gaussian", 
  cov_ranef = list(Spec_Latin_GenDivRange = phylo), 
  REML = F
)
summary(m10)
#########################

#########################
m11 <- pglmm(
  He ~ 
    # life history
    FECUNDITY + LONGEVITY + MAXTL + Anadromy +
    
    # random effects
    (1 | Study_id) +
    (1 | Spec_Latin_GenDivRange__), 
  data = dat_lm, 
  family = "gaussian", 
  cov_ranef = list(Spec_Latin_GenDivRange = phylo),
  REML = F
)
summary(m11)


phyr_pglmm_null$AIC
m1$AIC
m2$AIC
m3$AIC
m4$AIC
m5$AIC
m6$AIC
m7$AIC
m8$AIC
m9$AIC
m10$AIC
m11$AIC
  
# Best model with REML for fixed effect terms
m6_REML <- pglmm(
  He ~ 
    # local factors
    slope + elev + catch_area + stream_order +
    
    # historical factors
    History + snapped_lat +
    
    # random effects
    (1 | Study_id) +
    (1 | Spec_Latin_GenDivRange__), 
  data = dat_lm, 
  family = "gaussian", 
  cov_ranef = list(Spec_Latin_GenDivRange = phylo), 
  REML = T
)
summary(m6_REML)
summary(m6)


fixef_df <- data.frame(
  term = rownames(m6_REML$B),
  estimate = as.numeric(m6_REML$B),
  se = as.numeric(m6_REML$B.se),
  pval = as.numeric(m6_REML$B.pvalue)
) %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se
  )

fixef_df <- fixef_df %>%
  mutate(
    sig = case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    term_clean = gsub("^History", "", term),
    term_clean = gsub("\\(Intercept\\)", "Intercept", term_clean)
  )


fixef_df <- fixef_df %>%
  arrange(estimate) %>%
  mutate(term_clean = factor(term_clean, levels = term_clean))

fixef_df <- fixef_df %>%
  filter(term != "(Intercept)")

fixef_df$group <- ifelse(grepl("History", fixef_df$term),
                         "History",
                         "Local / Spatial")

ggplot(filter(fixef_df, group == "History"), aes(x = estimate, y = term_clean)) +
  # zero reference line
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  # confidence intervals
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  # points colored by significance
  geom_point(size = 3) +
  geom_text(
    aes(label = sig),
    nudge_y = 0.25,
    size = 4,
    vjust = 0
  ) +
  # labels
  labs(
    x = "Effect size (estimate ± 95% CI)",
    y = "",
    shape = "Significance",
    title = "Historical Fixed Effects"
  ) +
  # clean theme
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank()
  )

ggplot(filter(fixef_df, group == "Local / Spatial"), aes(x = estimate, y = term_clean)) +
  # zero reference line
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  # confidence intervals
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  # points colored by significance
  geom_point(size = 3) +
  geom_text(
    aes(label = sig),
    nudge_y = 0.25,
    size = 4,
    vjust = 0
  ) +
  # labels
  labs(
    x = "Effect size (estimate ± 95% CI)",
    y = "",
    shape = "Significance",
    title = "Local Fixed Effects"
  ) +
  # clean theme
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank()
  )


################################# an attempt at TMBglmm
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
    #FECUNDITY + LONGEVITY + MAXTL + 
    Anadromy +
    snapped_lat +
    
    # random effects
    (1 | Study_id) +
    (1 | Spec_Latin_GenDivRange) +
    propto(0 + Spec_Latin_GenDivRange | phylo_all, phylo_varcov),
  family = beta_family(link = "logit"),
  data = dat_lm
)
summary(lm1)
