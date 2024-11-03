
########## Calculating Metrics ##########

set.seed(6)

library(tidyverse)
library(tidybayes)
library(reticulate)
np <- import("numpy",convert=F)
py_available() 
py_numpy_available()

# loads the relevant python code
source_python("numerical_NFD.py")

RickUB <- readRDS(file.path("RDS_Files", "fit_RickUB.RDS"))
germ <- readRDS(file.path("RDS_Files", "germ.RDS"))

# Calculate ND and FD -----------------------------------------------------

# extract posterior samples from stanfit objects & put 'em in arrays

lams_df <- as.data.frame(RickUB)[,grepl("lam", names(as.data.frame(RickUB)))]
lams_array <- array(as.numeric(unlist(lams_df)), dim = c(nrow(lams_df), 3, 2, 2))

alphas_df <- as.data.frame(RickUB)[,grepl("alpha", names(as.data.frame(RickUB)))]
alphas_array <- array(as.numeric(unlist(alphas_df)), dim = c(nrow(lams_df), 3, 3, 2, 2))
# Godoy et al. (2014) ND & FI is undefined for any facilitation
# so replace facilitation with very slight competition
alphas_array2 <- alphas_array; alphas_array2[which(alphas_array2 <= 0)] <- 1*10^-5

germ_df <- as.data.frame(germ)[,grepl("g", names(as.data.frame(germ)))]
germ_array <- array(as.numeric(unlist(germ_df)), dim = c(nrow(germ_df), 3))

# growth model for Spaak & De Laender (2021) calculations
growth_mod <- function(N,g, lam, A){
  return(log((1 - g) + (g*lam*exp(-(A%*%(g*N))))))
  }

# initialize arrays to store metrics 
NDg <- CR1 <- CR2 <- CR <- DR1 <- DR2 <- DR <- FIg <- array(NA, dim = c(nrow(lams_df), 3, 3, 2, 2))
NDs <- FIs <- array(NA, dim = c(nrow(lams_df), 3, 3, 2, 2))
for(i in 1:3){
  for(j in 1:3){
    # Calculate Godoy et al. 2014 ND & FI
    CR1[,i,j,,] <- sqrt((alphas_array2[,j,i,,]*alphas_array2[,j,j,,])/(alphas_array2[,i,i,,]*alphas_array2[,i,j,,]))
    CR2[,i,j,,] <- sqrt((alphas_array2[,i,i,,]*alphas_array2[,i,j,,])/(alphas_array2[,j,i,,]*alphas_array2[,j,j,,]))
    DR1[,i,j,,] <- (germ_array[,i]*lams_array[,i,,])/(germ_array[,j]*lams_array[,j,,])
    DR2[,i,j,,] <- (germ_array[,j]*lams_array[,j,,])/(germ_array[,i]*lams_array[,i,,])
    FIg[,i,j,,] <- ifelse(DR1[,i,j,,]*CR1[,i,j,,] > DR2[,i,j,,]*CR2[,i,j,,], 
                          DR1[,i,j,,]*CR1[,i,j,,], DR2[,i,j,,]*CR2[,i,j,,])
    CR[,i,j,,] <- ifelse(DR1[,i,j,,]*CR1[,i,j,,] > DR2[,i,j,,]*CR2[,i,j,,],
                         CR1[,i,j,,], CR2[,i,j,,])
    DR[,i,j,,] <- ifelse(DR1[,i,j,,]*CR1[,i,j,,] > DR2[,i,j,,]*CR2[,i,j,,],
                         DR1[,i,j,,], DR2[,i,j,,])
    NDg[,i,j,,] <- 1 - sqrt((alphas_array2[,i,j,,]*alphas_array2[,j,i,,])/
                              (alphas_array2[,i,i,,]*alphas_array2[,j,j,,]))
    for(r in 1:2){
      for(n in 1:2){
        for(iter in 1:dim(alphas_array)[1]){
          # prep data for Spaak & De Laender (2021) calculations
          temp = list(c(germ_array[iter,i],
                        germ_array[iter,j]),
                      c(lams_array[iter,i,r,n],
                        lams_array[iter,j,r,n]),
                      matrix(c(alphas_array[iter,i,i,r,n], 
                               alphas_array[iter,i,j,r,n],
                               alphas_array[iter,j,j,r,n],
                               alphas_array[iter,j,i,r,n]), 
                             ncol = 2, nrow = 2, byrow = TRUE))
          
          # check if self facilitation occurs and if so, replace with
          # very slight competition
          #if(temp[[3]][1,1] < 0) temp[[3]][1,1] <- 1*10^-5
          #if(temp[[3]][2,2] < 0) temp[[3]][2,2] <- 1*10^-5
          
          # calculate Spaak & De Laender (2021) ND & FI
          mod <- try(NFD_model(growth_mod, n_spec = 2, 
                           args = temp, from_R = TRUE))
          if(inherits(mod, "try-error")) {}
          else{
            NDs[iter,i,j,r,n] <- mod$ND[1]
            FIs[iter,i,j,r,n] <- mod$FD[1]
          }
        }
      }
    }
  }
  }

output_list <- list(CR = CR, DR = DR, 
                    NDg = NDg, FIg = FIg, 
                    NDs = NDs, FIs = FIs)
write_rds(output_list, file.path("RDS_Files", "output_list.RDS"))
