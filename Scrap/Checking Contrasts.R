
##### CALCULATING TREATMENT EFFECT ESTIMATES

library(tidyverse)

# spit out helpful summary stats from samples of posterior distribution
good_info <- function(x){
  x2 <- x[which(x != 0)]
  c(median(x2, na.rm = TRUE), # median
  as.numeric(quantile(x2, 0.025, na.rm = TRUE)), # lower 95th CI
  as.numeric(quantile(x2, 0.975, na.rm = TRUE)), # upper 95th CI
  sum(x2 > 0, na.rm = TRUE)/sum(!is.na(x2), na.rm = TRUE), # probability of positive value
  sum(x2 < 0, na.rm = TRUE)/sum(!is.na(x2), na.rm = TRUE)) # probability of negative value
}

# bring in stanfit objects
BH_fit <- readRDS("BH_fit.rds")
germ_ests <- readRDS("germ_ests.rds")


# Calculate Metrics -------------------------------------------------------

# extract posterior samples from stanfit objects & put 'em in arrays

lams_df <- as.data.frame(BH_fit)[,grepl("lam", names(as.data.frame(BH_fit)))]
lams_array <- array(as.numeric(unlist(lams_df)), dim = c(nrow(lams_df), 3, 2, 2))

alphas_df <- as.data.frame(BH_fit)[,grepl("alpha", names(as.data.frame(BH_fit)))]
alphas_array <- array(as.numeric(unlist(alphas_df)), dim = c(nrow(lams_df), 3, 3, 2, 2))

germ_df <- as.data.frame(germ_ests)[,1:12]
germ_array <- array(as.numeric(unlist(germ_df)), dim = c(nrow(germ_df), 3, 2, 2))

# calculate competition metrics
ND <- CR1 <- CR2 <- CR <- DR1 <- DR2 <- DR <- FI <- array(NA, dim = c(nrow(lams_df), 3, 3, 2, 2))
for(i in 1:3){
  for(j in 1:3){
    ND[,i,j,,] <- 1 - sqrt((alphas_array[,i,j,,]*alphas_array[,j,i,,])/(alphas_array[,i,i,,]*alphas_array[,j,j,,]))
    CR1[,i,j,,] <- sqrt((alphas_array[,j,i,,]*alphas_array[,j,j,,])/(alphas_array[,i,i,,]*alphas_array[,i,j,,]))
    CR2[,i,j,,] <- sqrt((alphas_array[,i,i,,]*alphas_array[,i,j,,])/(alphas_array[,j,i,,]*alphas_array[,j,j,,]))
    DR1[,i,j,,] <- (germ_array[,i,,]*lams_array[,i,,])/(germ_array[,j,,]*lams_array[,j,,])
    DR2[,i,j,,] <- (germ_array[,j,,]*lams_array[,j,,])/(germ_array[,i,,]*lams_array[,i,,])
    FI[,i,j,,] <- ifelse(DR1[,i,j,,]*CR1[,i,j,,] > DR2[,i,j,,]*CR2[,i,j,,], 
                         DR1[,i,j,,]*CR1[,i,j,,], DR2[,i,j,,]*CR2[,i,j,,])
    CR[,i,j,,] <- ifelse(DR1[,i,j,,]*CR1[,i,j,,] > DR2[,i,j,,]*CR2[,i,j,,], 
                         CR1[,i,j,,], CR2[,i,j,,])
    DR[,i,j,,] <- ifelse(DR1[,i,j,,]*CR1[,i,j,,] > DR2[,i,j,,]*CR2[,i,j,,], 
                         DR1[,i,j,,], DR2[,i,j,,])
  }
}

# Calculate Contrasts -----------------------------------------------------

# subtracting experimental treatment from control

# germination

Grhizo_cont <- germ_array[,,2,] - germ_array[,,1,]
Gnitro_cont <- germ_array[,,,2] - germ_array[,,,1]
Ginter_cont <- Grhizo_cont[,,2] - Grhizo_cont[,,1]

# intrinsic growth rates

Lrhizo_cont <- lams_array[,,2,] - lams_array[,,1,]
Lnitro_cont <- lams_array[,,,2] - lams_array[,,,1]
Linter_cont <- Lrhizo_cont[,,2] - Lrhizo_cont[,,1]

# competition coefficients

Arhizo_cont <- alphas_array[,,,2,] - alphas_array[,,,1,]
Anitro_cont <- alphas_array[,,,,2] - alphas_array[,,,,1]
Ainter_cont <- Arhizo_cont[,,,2] - Arhizo_cont[,,,1]

# niche differences

NDrhizo_cont <- ND[,,,2,] - ND[,,,1,]
NDnitro_cont <- ND[,,,,2] - ND[,,,,1]
NDinter_cont <- NDrhizo_cont[,,,2] - NDrhizo_cont[,,,1]

# competitive ratio

CRrhizo_cont <- CR[,,,2,] - CR[,,,1,]
CRnitro_cont <- CR[,,,,2] - CR[,,,,1]
CRinter_cont <- CRrhizo_cont[,,,2] - CRrhizo_cont[,,,1]

# demographic ratio

DRrhizo_cont <- DR[,,,2,] - DR[,,,1,]
DRnitro_cont <- DR[,,,,2] - DR[,,,,1]
DRinter_cont <- DRrhizo_cont[,,,2] - DRrhizo_cont[,,,1]

# fitness inequalities

FIrhizo_cont <- FI[,,,2,] - FI[,,,1,]
FInitro_cont <- FI[,,,,2] - FI[,,,,1]
FIinter_cont <- FIrhizo_cont[,,,2] - FIrhizo_cont[,,,1]


# Compile into dataframes -------------------------------------------------

## Overall Effects --------------------------------------------------------

effect.df <- cbind.data.frame(rbind(good_info2(Grhizo_cont), good_info2(Gnitro_cont), 
                                    good_info2(Ginter_cont), good_info2(Lrhizo_cont), 
                                    good_info2(Lnitro_cont), good_info2(Linter_cont), 
                                    good_info2(Arhizo_cont), good_info2(Anitro_cont), 
                                    good_info2(Ainter_cont), good_info2(NDrhizo_cont), 
                                    good_info2(NDnitro_cont), good_info2(NDinter_cont), 
                                    good_info2(CRrhizo_cont), good_info2(CRnitro_cont), 
                                    good_info2(CRinter_cont), good_info2(DRrhizo_cont), 
                                    good_info2(DRnitro_cont), good_info2(DRinter_cont), 
                                    good_info2(FIrhizo_cont), good_info2(FInitro_cont), 
                                    good_info2(FIinter_cont)),
                              pred = rep(c("rhizo", "nitro", "inter"), 7), 
                              resp = rep(c("germ", "lam", "alph", 
                                           "ND", "CR", "DR", "FI"), each = 3))

names(effect.df) <- c("median", "lower", "upper", 
                      "pos_pd", "neg_pd", "pred", "resp")

# which treatment effects are we rather confident in? 

confident_effects <- function(x, p = .95){
  # x is dataframe of effects; p is probability cutoff
  x[which(x$pos_pd > p | x$neg_pd > p),]
}

confident_effects(effect.df, .8)

## Species Specific Effects -----------------------------------------------

# functions to get good_info2 output into binded up data

get_sp <- function(x) t(apply(x, 2, good_info2))
get_sp_pair <- function(x) apply(apply(x, c(2,3), good_info2), 1, identity)

sp.effect.df <- cbind.data.frame(rbind(get_sp(Grhizo_cont), get_sp(Gnitro_cont), 
                                       get_sp(Ginter_cont), get_sp(Lrhizo_cont), 
                                       get_sp(Lnitro_cont), get_sp(Linter_cont), 
                                       get_sp_pair(Arhizo_cont), 
                                       get_sp_pair(Anitro_cont), 
                                       get_sp_pair(Ainter_cont), 
                                       get_sp_pair(NDrhizo_cont), 
                                       get_sp_pair(NDnitro_cont), 
                                       get_sp_pair(NDinter_cont), 
                                       get_sp_pair(CRrhizo_cont), 
                                       get_sp_pair(CRnitro_cont), 
                                       get_sp_pair(CRinter_cont), 
                                       get_sp_pair(DRrhizo_cont), 
                                       get_sp_pair(DRnitro_cont), 
                                       get_sp_pair(DRinter_cont), 
                                       get_sp_pair(FIrhizo_cont), 
                                       get_sp_pair(FInitro_cont), 
                                       get_sp_pair(FIinter_cont)),
                              pred = c(rep(rep(c("rhizo", "nitro", "inter"), each = 3), 2), 
                                       rep(rep(c("rhizo", "nitro", "inter"), each = 9), 5)), 
                              resp = c(rep(c("germ", "lam"), each = 9), 
                                       rep(c("alph", "ND", "CR", "DR", "FI"), each = 27)), 
                              foc = rep(LETTERS[1:3], 51), 
                              comp = c(rep(NA, 18), rep(rep(LETTERS[1:3], each = 3), 15)))

names(sp.effect.df) <- c("median", "lower", "upper", 
                      "pos_pd", "neg_pd", "pred", "resp", "foc", "comp")

sp.effect.df <- sp.effect.df[which(sp.effect.df$median != "NA"),]

# which treatment effects are we rather confident in? 

confident_effects(sp.effect.df, .8)
