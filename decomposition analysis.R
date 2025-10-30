# clearing workspace
rm(list = ls())

# ------------------------------------------------------------------------------

## Header ##
# Title: Animal carrion decomposition in a subtropical forest biodiversity experiment
# Author: Finn Rehling 
# Date: 30th October 2025

# ------------------------------------------------------------------------------

## Comments ##
# For simplicity, I leave out all diagnostic and plotting procedures

# ------------------------------------------------------------------------------
# 1. Load libraries and dataset
# ------------------------------------------------------------------------------

# Directory
setwd() # set working directory

# load
require(glmmTMB)
require(ggeffects)
require(car)
require(ggplot2)
require(tidyr)
require(dplyr)
require(DHARMa)
require(psych)
require(irr)
require(lme4)
require(vegan)
require(permute)
require(RVAideMemoire)
require(grid)
require(stringr)
require(tibble)
require(gridExtra)

#carrion <- read.csv()
str(carrion)

# ------------------------------------------------------------------------------
# 2. Data preparation
# ------------------------------------------------------------------------------

# factor / numeric conversion
carrion$lgtsr       = log2(carrion$tree.richness)
carrion$location    = as.character(carrion$location)
carrion$site        = as.factor(carrion$site)
carrion$plot        = as.factor(carrion$plot)
carrion$year        = as.factor(carrion$year)
carrion$days.fac    = as.factor(carrion$days.fac)
carrion$tree.species  = as.factor(carrion$tree.species)
carrion$slope.steepness = as.numeric(carrion$slope.steepness)
carrion$prime.deco  = as.factor(carrion$prime.deco)
carrion$burried = as.factor(carrion$burried)
carrion$prime.deco = as.factor(carrion$prime.deco)
carrion$sc.initial = scale(carrion$initial.mass)
carrion$sc.cc = scale(carrion$canopy.cover)
carrion$sc.sl = scale(carrion$slope.steepness)
carrion$sc.tsr = scale(carrion$lgtsr)
carrion$sc.initial = scale(carrion$initial.mass)

# ------------------------------------------------------------------------------
# 3. Inter-rater reliability
# ------------------------------------------------------------------------------
alpha = psych::alpha(carrion[,c("deco.EN", "deco.MN", "deco.FR", "deco.JC")])
alpha  # n = 2622

# ------------------------------------------------------------------------------
# 4. Adjust decomposition scores for different photo days
# ------------------------------------------------------------------------------

# please note: because photos in 2023 after four days contained too many fully
# decomposed mice this led to a right censoring and a loss of variation
# we therefore excluded photos taken after 4 days from the analysis and only
# continued with those after 2 days. 

# However, a few photos were also taken after 1 or 3 days, affecting their 
# decomposition. Therefore, we adjusted the values of deco.raw -> deco.adj
# and corrected all values such that they would represent the decomposition
# stage after two days

#carrion$days.ord = factor(carrion$days, levels = 1:7, ordered = TRUE)
#m.cor1 = glmmTMB(deco~days.ord+(1|site/plot), 
 #                carrion[carrion$year == "2023" & carrion$days.fac == "early",], 
  #               family=tweedie)

# now adjust values for the differences between photos taken on different days
#pred <- ggpredict(m.cor1, terms = "days.ord")
#day2_value <- pred$predicted[pred$x == 2]

# Calculate adjustment factors for each day relative to day 2
#pred$adjustment <- day2_value - pred$predicted

# Step 4: Apply the adjustment to deco in your dataset
#carrion$deco.adj = carrion$deco
#carrion$deco.adj[carrion$year == "2023" & carrion$days.fac == "early"] <- 
 # carrion$deco[carrion$year == "2023" & carrion$days.fac == "early"] + 
  #pred$adjustment[match(carrion$days[carrion$year == "2023" & carrion$days.fac == "early"], 
   #                     pred$x)]

# Step 5: There is one plot in 2024, of which photos were taken only after four days. We will
# divide the scores by two to adjust the scores, with a min value of 1.
#carrion$deco.adj[carrion$days.fac == "early" & carrion$year == "2024" & carrion$days == 4] = 
 # ifelse(carrion$deco[carrion$days.fac == "early" & carrion$year == "2024" & carrion$days == 4]/2 >= 1,
  #       carrion$deco[carrion$days.fac == "early" & carrion$year == "2024" & carrion$days == 4]/2, 1)

# raw decomposition scores vs adjusted decomposition scores
#p.adj = ggplot()+
 # geom_jitter(carrion[carrion$days.fac == "early",], mapping=aes(days, deco.adj), shape=21, fill="goldenrod3", size=2, width=0.1)+
  #geom_jitter(carrion[carrion$days.fac == "early",], mapping=aes(days, deco), shape=21, fill="grey70", size=2, width=0.1)+
  #labs(y="Decomposition score", x="Number of days")+
  #annotate("text", x = 4.5, y = 6.5, 
   #        label = paste0("raw scores"),
    #       hjust = 1, size = 4, colour="grey70")+
  #annotate("text", x = 4.5, y = 6.2, 
   #        label = paste0("adjusted scores"),
    #       hjust = 1, size = 4, colour="goldenrod3")+
#  theme_classic()
#p.adj

# ------------------------------------------------------------------------------
# 5. Insect counts & scaling
# ------------------------------------------------------------------------------
# scale the number of decomposer for each year,
# i.e. flies for 2023
# and ants and other for 2024
# in addition, we will scale ants and no.other from 2024, and exchange values 
# from 2023 with those from 2024, as we did not sample them in 2023.

## 1.  Build a one-row-per plot table with the 2024 insect counts
##    if every 2024 replicate already carries the same value you can
##      use `distinct()`; otherwise take a summary (mean, max, etc.).

insects_24 <- carrion %>%                       # your full data frame
  filter(year == 2024) %>%                      # keep 2024 rows only
  group_by(site, plot) %>%                      # << keys that define a plot
  summarise(ants_24  = first(no.ants),          # or `mean(.., na.rm = TRUE)`
            other_24 = first(no.other),         #        ???
            .groups  = "drop")

carrion <- carrion %>%
  mutate(across(c(site, plot), as.factor))
insects_24 <- insects_24 %>%
  mutate(across(c(site, plot), as.factor))

carrion$no.ants[carrion$year == 2023] = NA
carrion$no.other[carrion$year == 2023] = NA

carrion <- carrion %>%
  left_join(insects_24, by = c("site", "plot")) %>%
  mutate(
    ants_24  = if_else(year == 2024, no.ants, ants_24),
    other_24 = if_else(year == 2024, no.other, other_24)
  )

carrion = carrion %>% 
  group_by(year) %>%
  mutate(
    lgants24      = log10(ants_24+1),
    lgflies = log10(no.flies+1),
    sc.lgants24   = as.numeric(scale(lgants24)),       # centre & sd???scale within yr
    sc.other24    = as.numeric(scale(other_24)),        # ditto for ???other??? insects
    sc.lgflies = as.numeric(scale(lgflies))
  ) %>% 
  ungroup()

str(carrion)

# extracting scaling parameters
scaling_params <- carrion %>%
  group_by(year) %>%
  summarise(
    mean_lgants24 = mean(log10(ants_24 + 1), na.rm = TRUE),
    sd_lgants24   = sd(log10(ants_24 + 1), na.rm = TRUE),
    mean_lgflies  = mean(lgflies, na.rm = TRUE),
    sd_lgflies    = sd(lgflies, na.rm = TRUE)
  )

scaling_params

b_ants_23  <- function(x) {
  lg <- x * scaling_params$sd_lgants[scaling_params$year == "2023"] + scaling_params$mean_lgants[scaling_params$year == "2023"]
  10^lg+1 
}
b_ants_24  <- function(x) {
  lg <- x * scaling_params$sd_lgants[scaling_params$year == "2024"] + scaling_params$mean_lgants[scaling_params$year == "2024"]
  10^lg+1 
}

b_flies_23 <- function(x) {
  lg <- x * scaling_params$sd_lgflies[scaling_params$year == "2023"] + scaling_params$mean_lgflies[scaling_params$year == "2023"]
  10^lg 
}

b_flies_24 <- function(x) {
  lg <- x * scaling_params$sd_lgflies[scaling_params$year == "2024"] + scaling_params$mean_lgflies[scaling_params$year == "2024"]
  10^lg 
}

#This fills ants_24/other_24 for both years but leaves the original 2023 no.ants, no.other intact.

# ------------------------------------------------------------------------------
# 6. Probability of main decomposer
# ------------------------------------------------------------------------------

### probability that flies were main decomposer in 2023 ~ environment
### probability that ants were main decomposer in 2024 ~ environment

#flies from 2023
mflies23 = glmmTMB(fly.prob ~ 
                     sc.initial+
                     sc.cc+
                     sc.sl+sc.tsr+
                     sc.lgflies+
                     sc.lgants24+
                     (1|site) 
                   , carrion[carrion$days.fac == "early" & carrion$year == "2023",], family=binomial)
Anova(mflies23, type="3")

mants24 = glmmTMB(ant.prob ~ sc.initial+
                    sc.cc+
                    sc.sl+
                    sc.tsr+
                    sc.lgants24+
                    sc.lgflies+
                    (1|site) 
                  , carrion[carrion$days.fac == "early" & carrion$year == "2024",], family=binomial)
Anova(mants24, type=3)

# ------------------------------------------------------------------------------
# 9. GLMM: decomposition (remain/loss)
# ------------------------------------------------------------------------------

# decomposition
mloss1 = glmmTMB(cbind(remain, loss)~sc.initial*year+
                   sc.cc*year+
                   sc.sl*year+
                   sc.tsr*year+
                   prime.deco*year+
                   (1|site/plot)+(1|tree.species)
                 , carrion[carrion$days.fac == "early",], family=betabinomial, 
                 dispformula = ~year+sc.cc+prime.deco)

Anova(mloss1, type="3")

