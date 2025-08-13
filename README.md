# Realistic-Null-Models
# R programmed functions to the realistic null models by Lopez-Angulo et al. (2021; Ecography) to detect assembly processes on the structure of plant communities.


######################################################################################################
# Manuscript: "A dimmer shade of pale: revealing the faint signature of local assembly processes 
#              on the structure of strongly filtered plant communities"
######################################################################################################

# ---------------------------------------------------------------------------------------------------
# 0. LOAD REQUIRED PACKAGES
# ---------------------------------------------------------------------------------------------------
# Each package serves a specific purpose in the workflow:
# - psych: PCA with varimax rotation for environmental variables
# - vegan: Beals smoothing for co-occurrence probabilities
# - MASS: stepAIC for GLM model selection
# - picante: functional diversity metrics (adapted)
# - lme4, lmerTest: mixed-effects models with significance testing
# - emmeans: post-hoc comparisons
# - ggplot2, ggpubr: plotting and figure combination
library(psych)
library(vegan)
library(MASS)
library(picante)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(ggpubr)

# ---------------------------------------------------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------------------------------------------------
# This script requires the file 'DataSets.Rdata' available from:
# López-Angulo et al. 2019. Alpine vegetation dataset from three contrasting mountain ranges.
# DOI: https://doi.org/10.1016/j.dib.2019.104816
#
# The file contains:
# - plant_cover: % cover matrix of species (sites x species)
# - funct_traits: data.frame of functional traits (scaled, species x traits)
# - predictors: standardized environmental variables (11 climatic + 3 topographic)
#
# NOTE: 5 plots with <90% cumulative cover of species with measured traits were removed beforehand.
load("DataSets.Rdata")

# Ensure reproducibility of randomization steps
set.seed(123)

# ---------------------------------------------------------------------------------------------------
# 2. AUXILIARY FUNCTIONS
# ---------------------------------------------------------------------------------------------------

## 2.1. Generate null community matrices while keeping local species richness constant
# cover: observed species cover matrix
# prob: species occurrence probabilities (can be uniform, environmental, or co-occurrence based)
# nsim: number of null matrices to generate
simu.com <- function(cover, prob, nsim) {
  richness <- rowSums(cover > 0)
  nsp <- ncol(cover)
  result <- vector("list", nsim)
  for(i in seq_len(nsim)) {
    null.assemblage <- t(mapply(function(r, p) {
      sp_sel <- sort(sample(seq_len(nsp), size = r, prob = p))
      seq_len(nsp) %in% sp_sel
    }, richness, data.frame(t(prob))))
    com.sim <- null.assemblage
    for(sp in seq_len(ncol(com.sim))) {
      com.sim[null.assemblage[,sp] > 0, sp] <- sample(
        cover[cover[,sp] > 0, sp],
        size = sum(null.assemblage[,sp] > 0),
        replace = TRUE
      )
    }
    colnames(com.sim) <- colnames(cover)
    rownames(com.sim) <- rownames(cover)
    result[[i]] <- com.sim
  }
  return(result)
}

## 2.2. "Melodic" function (Bello et al. 2016) to compute Rao's quadratic entropy
# Calculates MPD and RaoQ for both abundance and presence-absence data
melodic <- function(samp, dis, type="both") {
  if(!is.matrix(samp)) samp <- as.matrix(samp)
  if(!is.matrix(dis))  dis  <- as.matrix(dis)
  if(is.null(colnames(samp)) | is.null(colnames(dis))) stop("Missing 'colnames' in 'samp' and/or 'dis'")
  
  N <- nrow(samp)
  melodic <- list()
  if(type %in% c("both","abundance")) {
    melodic$abundance <- list(mpd=numeric(N), rao=numeric(N), simpson=numeric(N))
  }
  if(type %in% c("both","presence")) {
    melodic$presence <- list(mpd=numeric(N), rao=numeric(N), simpson=numeric(N))
  }
  for(i in seq_len(N)) {
    spp <- names(samp[i, samp[i,] > 0])
    if(length(spp) > 1) {
      d <- dis[spp, spp]
      if(type %in% c("both","abundance")) {
        w <- samp[i, spp] / sum(samp[i, spp])
        W <- outer(w, w)
        melodic$abundance$mpd[i] <- weighted.mean(d[lower.tri(d)], W[lower.tri(W)])
        melodic$abundance$rao[i] <- sum(W * d)
      }
      if(type %in% c("both","presence")) {
        w <- rep(1, length(spp)) / length(spp)
        W <- outer(w, w)
        melodic$presence$mpd[i] <- weighted.mean(d[lower.tri(d)], W[lower.tri(W)])
        melodic$presence$rao[i] <- sum(W * d)
      }
    } else {
      if(type %in% c("both","abundance")) melodic$abundance$rao[i] <- 0
      if(type %in% c("both","presence"))  melodic$presence$rao[i]  <- 0
    }
  }
  return(melodic)
}

## 2.3. Compute standardized effect size (SES) of RaoQ given a null model
# nullcom: list of null community matrices
# cover: observed species cover matrix
# traits: functional trait matrix
# trait_index: column index of the focal trait in 'traits'
ses.rao <- function(nullcom, cover, traits, trait_index) {
  runs <- length(nullcom)
  dis <- dist(as.data.frame(traits)[, trait_index])
  rao.obs <- melodic(as.matrix(cover), as.matrix(dis), type="abundance")$abundance$rao
  rand <- simplify2array(lapply(nullcom, function(x) {
    melodic(as.matrix(x[, rownames(traits)]), as.matrix(dis), type="abundance")$abundance$rao
  }))
  rao.mean <- apply(rand, 1, mean, na.rm=TRUE)
  rao.sd   <- apply(rand, 1, sd, na.rm=TRUE)
  rao.z    <- (rao.obs - rao.mean) / rao.sd
  rao.rank <- apply(cbind(rao.obs, rand), 1, rank)[1,]
  rao.perc <- rao.rank / (runs + 1)
  data.frame(ntaxa=specnumber(cover), runs=runs,
             rao.obs, rao.rand.mean=rao.mean, rao.rand.sd=rao.sd,
             rao.obs.z=rao.z, rao.obs.perc=rao.perc,
             row.names=row.names(cover))
}

# ---------------------------------------------------------------------------------------------------
# 3. BUILD SPECIES OCCURRENCE PROBABILITY MODELS
# ---------------------------------------------------------------------------------------------------

## 3.1. Stochastic model: uniform probability of occurrence for all species
stochastic_prob <- matrix(sum(plant_cover > 0) / (nrow(plant_cover)*ncol(plant_cover)),
                          nrow=nrow(plant_cover), ncol=ncol(plant_cover))

## 3.2. Independence model: GLMs based on environmental predictors
# PCA of environmental variables (excluding first 3 variables)
pca_env1 <- principal(scale(predictors[,-c(1:3)]), nfactors=4, rotate="varimax")
pca_env <- data.frame(c1=pca_env1$scores[,1], c2=pca_env1$scores[,2],
                      c3=pca_env1$scores[,3], c4=pca_env1$scores[,4])
pca2 <- pca_env^2
pca_predictors <- cbind(pca_env, pca2)
colnames(pca_predictors) <- c("c1","c2","c3","c4","c12","c22","c32","c42")

# Fit a binomial GLM for each species, select best model via AIC
list_binomial_glms <- vector("list", ncol(plant_cover))
for(i in seq_len(ncol(plant_cover))) {
  glm_i <- glm((plant_cover[,i] > 0) ~ ., pca_predictors, family=binomial)
  list_binomial_glms[[i]] <- stepAIC(glm_i, trace=FALSE)
}
glm_prob <- as.data.frame(lapply(list_binomial_glms, predict, type="response"))
colnames(glm_prob) <- colnames(plant_cover)

## 3.3. Co-occurrence model: Beals smoothing
cooccurr_prob <- beals(plant_cover, type=0, include=FALSE)

# ---------------------------------------------------------------------------------------------------
# 4. GENERATE NULL COMMUNITIES
# ---------------------------------------------------------------------------------------------------
stochastic_matrices   <- simu.com(plant_cover, stochastic_prob, nsim=999)
independence_matrices <- simu.com(plant_cover, glm_prob, nsim=999)
cooccurrence_matrices <- simu.com(plant_cover, cooccurr_prob, nsim=999)

# ---------------------------------------------------------------------------------------------------
# 5. CALCULATE SES FOR EACH TRAIT UNDER EACH NULL MODEL
# ---------------------------------------------------------------------------------------------------
# Trait indices in 'funct_traits'
traits_names <- c("seed_mass"=1, "leaf_thickness"=2, "LDMC"=3, "height"=5)
null_models <- list(stochastic=stochastic_matrices,
                    independence=independence_matrices,
                    cooccurrence=cooccurrence_matrices)

SES_results <- list()
for(trait in names(traits_names)) {
  SES_results[[trait]] <- lapply(null_models, function(nm) {
    ses.rao(nm, plant_cover, funct_traits, traits_names[trait])
  })
}

# ---------------------------------------------------------------------------------------------------
# 6. STATISTICAL TESTING EXAMPLE (LMM for seed mass SES)
# ---------------------------------------------------------------------------------------------------
# Prepare data for mixed-effects modeling
null_data <- do.call(rbind, replicate(length(null_models), predictors[,1:2], simplify=FALSE))
null_data$null <- factor(rep(names(null_models), each=nrow(plant_cover)))
null_data$sample <- factor(rep(1:nrow(plant_cover), times=length(null_models)))

SES_seed <- do.call(rbind, SES_results$seed_mass)
lm_seed <- lmer(rao.obs.z ~ null + (1|sample) + (1|plot), data=cbind(SES_seed, null_data))
summary(lm_seed)
emmeans(lm_seed, pairwise~null, adjust="tukey")
