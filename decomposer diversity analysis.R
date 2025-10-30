
# clearing workspace
rm(list = ls())

# ------------------------------------------------------------------------------
## Header ##
# Title: Carrion diversity in a subtropical forest biodiversity experiment
# Author: Finn Rehling
# Date: 30th October 2025

# ------------------------------------------------------------------------------

## Comments ##
# For simplicity, I leave out all diagnostic and plotting procedures

# ------------------------------------------------------------------------------
# 1. Load libraries and set working directory
# ------------------------------------------------------------------------------

# Directory
setwd() # set working directory

# load
library(car)
library(dplyr)
library(DHARMa)
library(ggplot2)
library(ggeffects)
library(grid)
library(gridExtra)
library(glmmTMB)
library(iNEXT)
library(irr)
library(lme4)
library(permute)
library(psych)
library(RVAideMemoire)
library(stringr)
library(tibble)
library(tidyr)
library(vegan)

# ------------------------------------------------------------------------------
# 2. Read and prepare data
# ------------------------------------------------------------------------------

#decom <- read.csv()
str(decom)

sum(table(unique(decom$plot[decom$site == "A"]))) # 48 plots
sum(table(unique(decom$plot[decom$site == "B"]))) # 48 plots

decom$tree.richness[decom$tree.richness>16] = 16  # the 24-tree species richness plots will be reduced to 16 tree species, 
# because the 16 tree species are closer to the realized tree species richness of those plots after 15 years of growth.
decom$lgtsr = log2(decom$tree.richness) #log-transformation
decom$year = as.character(decom$year) # year as categorical factor

# ------------------------------------------------------------------------------
# 3. Prepare decomposer matrices (flies, ants, other)
# ------------------------------------------------------------------------------

# long-format
decom_long = decom[,c(1:3, 7:83)] %>%
  pivot_longer(cols = fly.1:formicidae.k, names_to = "species", values_to = "count") %>%
  filter(count > 0)
decom_long

#dataset on flies in 2023
# take only those traps when you have at least two species and five individuals
flies.23 <- decom_long %>%
  filter(year == 2023, str_detect(species, "^fly")) %>%
  pivot_wider(names_from = species, values_from = count, values_fill = 0) %>%
  unite("site.plot", site, plot, remove = FALSE) %>%
  select(-site, -plot, -year) %>%
  column_to_rownames("site.plot") %>%
  as.data.frame() %>%
  filter(
    rowSums(. != 0) >= 2,   # at least two non???zero species counts
    rowSums(.)       >= 5    # at least five individuals in total
  )

#dataset on flies in 2024
# take only those traps when you have at least two species and five individuals
flies.24 <- decom_long %>%
  filter(year == 2024, str_detect(species, "^fly")) %>%
  pivot_wider(names_from = species, values_from = count, values_fill = 0) %>%
  unite("site.plot", site, plot, remove = FALSE) %>%
  select(-site, -plot, -year) %>%
  column_to_rownames("site.plot") %>%
  as.data.frame() %>%
  filter(
    rowSums(. != 0) >= 2,   # at least two non???zero species counts
    rowSums(.)       >= 5    # at least five individuals in total
  )

#transpone
flies.23.t = t(flies.23)
flies.24.t = t(flies.24)

# rarefaction and extrapolation
require(iNEXT)

### --- Helper: check if a plot can realistically reach the coverage level ---
is_coverage_reachable <- function(vec, coverage_level) {
  vec <- vec[vec > 0]
  if (length(vec) < 2) return(FALSE)
  
  current_cov <- DataInfo(vec, datatype = "abundance")$SC
  if (coverage_level <= current_cov) return(TRUE)
  
  size_needed <- tryCatch(
    coverage_to_size(vec, datatype = "abundance", level = coverage_level),
    error = function(e) Inf
  )
  size_needed <= 2 * sum(vec)  # no more than double observed individuals
}

### --- get_inext_coverage_summary ---
get_inext_coverage_summary <- function(comm_matrix, coverage_level, year_label, prefix = "flies_") {
  valid_plots <- colnames(comm_matrix)[
    sapply(colnames(comm_matrix), function(name) {
      is_coverage_reachable(comm_matrix[, name], coverage_level)
    })
  ]
  
  if (length(valid_plots) == 0) return(NULL)
  
  filtered_matrix <- comm_matrix[, valid_plots, drop = FALSE]
  
  cov_out <- estimateD(filtered_matrix, q = 0, datatype = "abundance",
                       base = "coverage", level = coverage_level, conf = 0.95)
  
  summary <- cov_out %>%
    filter(Order.q == 0) %>%
    select(Assemblage, Estimator = qD, SC) %>%
    mutate(Coverage_Level = coverage_level, year = year_label) %>%
    separate(Assemblage, into = c("site", "plot"), sep = "_") %>%
    rename_with(~ paste0(prefix, .), c("Estimator", "SC"))
  
  return(summary)
}

# a helper function for iNEXT, which gives up estimates for species richness 
# (q = 0) depending on year, and the species group
flies.23.sum <- get_inext_coverage_summary(flies.23.t, 0.8, "2023", prefix = "flies_")
flies.24.sum <- get_inext_coverage_summary(flies.24.t, 0.8, "2024", prefix = "flies_")

flies_summary <- bind_rows(flies.23.sum, flies.24.sum)

# flies_summary contains the species richness (q = 0) for flies in the traps
# at a coverage level of 0.8

#merge
decom <- left_join(decom, flies_summary,
                   by = c("site", "plot", "year"))

#reducing dataset on ants in 2024
ants.24 <- decom_long %>%
  filter(year == 2024, str_detect(species, "^formicidae")) %>%
  pivot_wider(names_from = species, values_from = count, values_fill = 0) %>%
  unite("site.plot", site, plot, remove = FALSE) %>%
  select(-site, -plot, -year) %>%
  column_to_rownames("site.plot") %>%
  as.data.frame() %>%
  filter(
    rowSums(. != 0) >= 2,   # at least two non???zero species counts
    rowSums(.)       >= 5    # at least five individuals in total
  )

ants.24.t = t(ants.24)
# species coverage standardized to 0.8
ants.24.sum <- get_inext_coverage_summary(ants.24.t, 0.8, "2024", prefix = "ants_")

# Merge ants summary (only fills 2024 rows; 2023 rows get NA)
decom <- decom %>%
  left_join(ants.24.sum, by = c("site", "plot", "year"))

# reducing dataset for other decomposer in 2024
other.24 <- decom_long %>%
  filter(year == 2024) %>%
  filter(!str_detect(species, "^fly") & 
           !str_detect(species, "^formicidae")) %>%
  pivot_wider(names_from = species, values_from = count, values_fill = 0) %>%
  unite("site.plot", site, plot, remove = FALSE) %>%
  select(-site, -plot, -year) %>%
  column_to_rownames("site.plot") %>%
  as.data.frame() %>%
  filter(
    rowSums(. != 0) >= 2,   # at least two non???zero species counts
    rowSums(.)       >= 5    # at least five individuals in total
  )

other.24.t <- t(other.24)

other.24.sum <- get_inext_coverage_summary(
  comm_matrix = other.24.t, 
  coverage_level = 0.6, 
  year_label = "2024", 
  prefix = "other_"
)

# Merge other decomposer
decom <- decom %>%
  left_join(other.24.sum, by = c("site", "plot", "year"))

# Species coverage
mean(decom$flies_SC, na.rm=T)
mean(decom$ants_SC, na.rm=T)
mean(decom$other_SC, na.rm=T)

# ------------------------------------------------------------------------------
# 4. GLMM analyses
# ------------------------------------------------------------------------------

decom$sc.cc = scale(decom$canopy.cover)
decom$sc.sl = scale(decom$slope.steepness)
decom$sc.tsr = scale(decom$lgtsr)

# Identify the fly columns (assuming names start with "fly.")
fly_cols <- grep("^fly\\.", names(decom), value = TRUE)

# Sum flies per row and log10-transform
decom$no.flies <- rowSums(decom[, fly_cols], na.rm = TRUE)
decom$lgflies <- log10(decom$no.flies + 1)  # +1 avoids log(0) for some traps 

sum(table(unique(decom$plot)))

#fly abundance
m.flies = glmmTMB(no.flies ~ sc.tsr*year+
                    sc.cc*year+
                    sc.sl*year+
                    (1|site),
                  decom, family=nbinom1, dispformula = ~year)
Anova(m.flies, type=3)

plot(ggpredict(m.flies, terms=c("sc.cc", "year")), show_data=T)
plot(ggpredict(m.flies, terms=c("sc.sl", "year")), show_data=T)
plot(ggpredict(m.flies, terms=c("sc.tsr", "year")), show_data=T)

#fly species richness
hist(decom$flies_Estimator)

m.spflies = glmmTMB(flies_Estimator ~ 
                      lgflies*year+
                      sc.tsr*year+
                      sc.cc*year+
                      sc.sl*year+
                      (1|site),family= Gamma(link = "log"),
                    decom)
Anova(m.spflies, type=3)
plot(ggpredict(m.spflies, terms=c("sc.cc", "year")), show_data=T)

# ant analyses
# Identify the ant columns (assuming names start with "ant.")
ant_cols <- grep("^formicidae\\.", names(decom), value = TRUE)

# Sum ants per row and log10-transform
decom$no.ants <- rowSums(decom[, ant_cols], na.rm = TRUE)
decom$lgants <- log10(decom$no.ants + 1)  # +1 avoids log(0) for some traps 

#ant abundance
m.ants = glmmTMB(no.ants ~ sc.tsr+
                   sc.cc+
                   sc.sl+
                   (1|site),
                 decom[decom$year=="2024",], family=nbinom1)
Anova(m.ants, type=3)

plot(ggpredict(m.ants, terms=c("sc.cc")), show_data=T)
plot(ggpredict(m.ants, terms=c("sc.sl")), show_data=T)
plot(ggpredict(m.ants, terms=c("sc.tsr")), show_data=T)

#species richness
m.spants = glmmTMB(ants_Estimator ~ lgants+sc.tsr+
                     sc.cc+
                     sc.sl+
                     (1|site), family= Gamma(link = "log"),
                   decom[decom$year=="2024",])
Anova(m.spants, type=3)
plot(ggpredict(m.spants, terms=c("sc.cc")), show_data=T)

# other analyses
other_cols <- names(decom)[
  grepl("^(coleoptera|lepidoptera|orthoptera|blattoidea|dermaptera|hemiptera|mecoptera|scutigeromorpha)\\.", names(decom))
]
other_cols

# Sum other per row and log10-transform
decom$no.other <- rowSums(decom[, other_cols], na.rm = TRUE)
decom$lgother <- log10(decom$no.other + 1)  # +1 avoids log(0) for some traps 

#other abundance
m.other = glmmTMB(no.other ~ sc.tsr+
                    sc.cc+
                    sc.sl+
                    (1|site),
                  decom[decom$year=="2024",], family=nbinom1)
Anova(m.other, type=3)

plot(ggpredict(m.other, terms=c("sc.cc")), show_data=T)
plot(ggpredict(m.other, terms=c("sc.sl")), show_data=T)
plot(ggpredict(m.other, terms=c("sc.tsr")), show_data=T)

#gamma distribution
m.spother = glmmTMB(other_Estimator ~ 
                      lgother+
                      sc.cc+sc.tsr+
                      sc.sl, family= Gamma(link = "log"),
                    decom[decom$year=="2024",])
Anova(m.spother, type=3)

# ------------------------------------------------------------------------------
# 6. Community composition: NMDS + adonis
# ------------------------------------------------------------------------------

## exclude plots with partially knocked over traps
# i.e I28, F25, Q5, U20
# and exclude plots with 0 flies
fly.com = decom %>% filter(decom$no.flies>5)
fly.com.23 = decom %>% filter(decom$no.flies>5 & year == "2023")
fly.com.24 = decom %>% filter(decom$no.flies>5 & year == "2024")

# stata = decom.com$site to account for nestedness of data structure
perm.23 <- how(plots = Plots(strata = fly.com.23$site), nperm=9999)
perm.24 <- how(plots = Plots(strata = fly.com.24$site), nperm=9999)

# hellinger transformation to weight common species
fly.comp.23 = fly.com.23[,c(7:22)]
fly.comp.24 = fly.com.24[,c(7:22)]

all.mds.23 <- metaMDS(fly.comp.23)
all.mds.24 <- metaMDS(fly.comp.24)
fly.com.23$NMDS1 = as.data.frame(scores(all.mds.23)$sites)[,1]
fly.com.23$NMDS2 = as.data.frame(scores(all.mds.23)$sites)[,2]
fly.com.24$NMDS1 = as.data.frame(scores(all.mds.24)$sites)[,1]
fly.com.24$NMDS2 = as.data.frame(scores(all.mds.24)$sites)[,2]

set.seed(1234)
fly.bray.23 <- adonis2(fly.comp.23~ sc.cc+sc.tsr+sc.sl, 
                       data = fly.com.23,
                       permutations = perm.23,
                       method = "bray",
                       by="margin")

fly.bray.24 <- adonis2(fly.comp.24~ sc.cc+sc.tsr+sc.sl, 
                       data = fly.com.24,
                       permutations = perm.24,
                       method = "bray",
                       by="margin")
fly.bray.23
fly.bray.24

cor.test(fly.com.23$NMDS1, fly.com.23$sc.cc)
cor.test(fly.com.23$NMDS2, fly.com.23$sc.cc)

#### other carrion decomposer
other.com = decom %>% filter(year == "2024" & no.other > 5)
perm1 <- how(plots = Plots(strata = other.com$site), nperm=9999)
other.comp = other.com[,c(23:70)]

all.mds <- metaMDS(other.comp) 
other.com$NMDS1 = as.data.frame(scores(all.mds)$sites)[,1]
other.com$NMDS2 = as.data.frame(scores(all.mds)$sites)[,2]

set.seed(1234)
other.bray <- adonis2(other.comp~sc.cc+sc.tsr+sc.sl, 
                      data = other.com,
                      permutations = perm1,
                      method = "bray",
                      by="margin")
other.bray

# Correlation between NMDS1 and the environmental variables
cor.test(other.com$canopy.cover, other.com$NMDS1, method="spearman")

# ant community composition analysis
ants.com2 = decom %>% filter(year == "2024" & no.ants > 5)

# stata = decom.com$site to account for nestedness of data structure
perm1 <- how(plots = Plots(strata = ants.com2$site), nperm=9999)
# hellinger transformation to weight common species
# Note: No other insects were found in plot M7 in 2023 and 2024
ants.comp = ants.com2[,c(71:83)]

set.seed(1234)
all.mds <- metaMDS(ants.comp) 
ants.com2$NMDS1 = as.data.frame(scores(all.mds)$sites)[,1]
ants.com2$NMDS2 = as.data.frame(scores(all.mds)$sites)[,2]
ants.bray <- adonis2(ants.comp~sc.cc+sc.tsr+sc.sl, 
                     data = ants.com2,
                     permutations = perm1,
                     method = "bray",
                     by="margin")
ants.bray

