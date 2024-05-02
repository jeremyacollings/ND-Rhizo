
good_info <- function(x){
  print(paste("Median =", median(x, na.rm = TRUE), ""))
  print(paste("Lower 95% CI =", as.numeric(quantile(x, 0.025, na.rm = TRUE)), ""))
  print(paste("Upper 95% CI =", as.numeric(quantile(x, 0.975, na.rm = TRUE)), ""))
  print(paste("Prop of samples > 0 =", sum(x > 0, na.rm = TRUE)/sum(!is.na(x), na.rm = TRUE), ""))
  ggplot() + geom_histogram(aes(x))
}

# calculate ND, DR, CR, and FI

BH_fit <- readRDS("BH_fit.rds")
germ_ests <- readRDS("germ_ests.rds")

lams_df <- as.data.frame(BH_fit)[,1:12]
lams_array <- array(as.numeric(unlist(lams_df)), dim = c(nrow(lams_df), 3, 2, 2))

alphas_df <- as.data.frame(BH_fit)[,25:60]
alphas_array <- array(as.numeric(unlist(alphas_df)), dim = c(nrow(lams_df), 3, 3, 2, 2))

germ_df <- as.data.frame(germ_ests)[,1:12]
germ_array <- array(as.numeric(unlist(germ_df)), dim = c(nrow(germ_df), 3, 2, 2))

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

# some contrasts

NDrhizo_cont_zn <- ND[,,,1,1] - ND[,,,2,1] # sterile - rhizo; no nitrogen
NDrhizo_cont_n <- ND[,,,1,2] - ND[,,,2,2] # sterile - rhizo; nitrogen
NDnitro_cont_s <- ND[,,,1,1] - ND[,,,1,2] # no nitrogen - nitrogen; sterile
NDnitr_cont_r <- ND[,,,2,1] - ND[,,,2,2] # no nitrogen - nitrogen; rhizobia
NDinter_cont <- NDrhizo_cont_zn - NDrhizo_cont_n

CRrhizo_cont_zn <- CR[,,,1,1] - CR[,,,2,1] # sterile - rhizo; no nitrogen
CRrhizo_cont_n <- CR[,,,1,2] - CR[,,,2,2] # sterile - rhizo; nitrogen
CRnitro_cont_s <- CR[,,,1,1] - CR[,,,1,2] # no nitrogen - nitrogen; sterile
CRnitr_cont_r <- CR[,,,2,1] - CR[,,,2,2] # no nitrogen - nitrogen; rhizobia
CRinter_cont <- CRrhizo_cont_zn - CRrhizo_cont_n

DRrhizo_cont_zn <- DR[,,,1,1] - DR[,,,2,1] # sterile - rhizo; no nitrogen
DRrhizo_cont_n <- DR[,,,1,2] - DR[,,,2,2] # sterile - rhizo; nitrogen
DRnitro_cont_s <- DR[,,,1,1] - DR[,,,1,2] # no nitrogen - nitrogen; sterile
DRnitr_cont_r <- DR[,,,2,1] - DR[,,,2,2] # no nitrogen - nitrogen; rhizobia
DRinter_cont <- DRrhizo_cont_zn - DRrhizo_cont_n

FIrhizo_cont_zn <- FI[,,,1,1] - FI[,,,2,1] # sterile - rhizo; no nitrogen
FIrhizo_cont_n <- FI[,,,1,2] - FI[,,,2,2] # sterile - rhizo; nitrogen
FInitro_cont_s <- FI[,,,1,1] - FI[,,,1,2] # no nitrogen - nitrogen; sterile
FInitr_cont_r <- FI[,,,2,1] - FI[,,,2,2] # no nitrogen - nitrogen; rhizobia
FIinter_cont <- FIrhizo_cont_zn - FIrhizo_cont_n

# is there an effect of rhizobia (when no nitrogen was added)
# on the competitive ratio between Acmispon and Collinsia
good_info(CRrhizo_cont_zn[,1,2]) # probably not

# search for the interesting ones

