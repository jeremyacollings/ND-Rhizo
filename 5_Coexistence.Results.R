
########## COEXISTENCE RESULTS ##########

library(tidyverse)

output_list <- read_rds(file.path("RDS_Files", "output_list.RDS"))

# Coexistence Graphs ------------------------------------------------------

# Godoy et al.

samps <- sample(1:dim(output_list[['NDg']])[1], 100)

god.samps <- cbind.data.frame(ND = c(output_list[['NDg']][samps,,,,]), 
                 FI = c(output_list[['FIg']][samps,,,,]), 
                 foc = rep(rep(c("AA", "CG", "LP"), each = 100), 3*2*2), 
                 comp = rep(rep(c("AA", "CG", "LP"), each = 100*3), 2*2), 
                 rhizo = rep(rep(c("cont", "add"), each = 100*3*3), 2), 
                 nitro = rep(c("uninoc", "inoc"), each = 100*3*3*2))
god.samps$pair <- interaction(do.call(pmin, god.samps[,c("foc", "comp")]), 
                              do.call(pmax, god.samps[,c("foc", "comp")]))
god.samps$treat <- paste(god.samps$rhizo, god.samps$nitro)

god.meds <- cbind.data.frame(ND = c(apply(output_list[['NDg']], c(2,3,4,5), median, na.rm = TRUE)), 
                             FI = c(apply(output_list[['FIg']], c(2,3,4,5), median, na.rm = TRUE)), 
                             foc = rep(c("AA", "CG", "LP"), 3*2*2), 
                             comp = rep(rep(c("AA", "CG", "LP"), each = 3), 2*2), 
                             rhizo = rep(rep(c("cont", "add"), each = 3*3), 2), 
                             nitro = rep(c("uninoc", "inoc"), each = 3*3*2))
god.meds$pair <- interaction(do.call(pmin, god.meds[,c("foc", "comp")]), 
                              do.call(pmax, god.meds[,c("foc", "comp")]))
god.meds$treat <- paste(god.meds$rhizo, god.meds$nitro)

god.samps %>%
  filter(foc != comp) %>%
  ggplot(aes(x = ND, y = FI, color = pair, shape = treat)) + 
  geom_point(alpha = .25) + 
  xlim(-50, 1) + ylim(1, 50) + 
  stat_function(fun = function(x) -x+1, color = "black") + 
  stat_function(fun = function(x) 1/(-x+1), color = "black") + 
  theme_classic(base_size = 15) + 
  geom_point(data = god.meds[which(god.meds$foc != god.meds$comp),], 
             size = 5)
  
# Spaak & De Laender

spaak.samps <- cbind.data.frame(ND = c(output_list[['NDs']][samps,,,,]), 
                              FI = c(output_list[['FIs']][samps,,,,]), 
                              foc = rep(rep(c("AA", "CG", "LP"), each = 100), 3*2*2), 
                              comp = rep(rep(c("AA", "CG", "LP"), each = 100*3), 2*2), 
                              rhizo = rep(rep(c("uninoc", "inoc"), each = 100*3*3), 2), 
                              nitro = rep(c("cont", "add"), each = 100*3*3*2))
spaak.samps$pair <- interaction(do.call(pmin, spaak.samps[,c("foc", "comp")]), 
                               do.call(pmax, spaak.samps[,c("foc", "comp")]))
spaak.samps$treat <- paste(spaak.samps$rhizo, spaak.samps$nitro)

spaak.meds <- cbind.data.frame(ND = c(apply(output_list[['NDs']], c(2,3,4,5), median, na.rm = TRUE)), 
                             FI = c(apply(output_list[['FIs']], c(2,3,4,5), median, na.rm = TRUE)), 
                             foc = rep(c("AA", "CG", "LP"), 3*2*2), 
                             comp = rep(rep(c("AA", "CG", "LP"), each = 3), 2*2), 
                             rhizo = rep(rep(c("uninoc", "inoc"), each = 3*3), 2), 
                             nitro = rep(c("cont", "add"), each = 3*3*2))
spaak.meds$pair <- interaction(do.call(pmin, spaak.meds[,c("foc", "comp")]), 
                             do.call(pmax, spaak.meds[,c("foc", "comp")]))
spaak.meds$treat <- paste(spaak.meds$rhizo, spaak.meds$nitro)

spaak.samps %>%
  filter(foc != comp) %>%
  ggplot(aes(x = ND, y = -FI, color = pair, shape = treat)) + 
  geom_point(alpha = .1) + 
  xlim(-.15, 1.5) + ylim(-1, 1) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 1, linetype = "dashed") + 
  stat_function(fun = function(x) x/(1-x), color = "black") +
  theme_classic(base_size = 15) + 
  geom_point(data = spaak.meds[which(spaak.meds$foc != spaak.meds$comp),], 
             size = 3) + 
  scale_shape_manual(values = c(0, 15, 1, 16), 
                     name = "Treatment", 
                     labels = c("Add.Inoc", "Add.Uninoc", 
                                "Cont.Inoc", "Cont.Uninoc"))

# Coexistence Probabilities -----------------------------------------------

# get coexistence probabilities per pair per treatment

persist <- -output_list[["FIs"]] < output_list[["NDs"]]/
  (1 - output_list[["NDs"]])

outcome <- array(NA, dim = c(3, dim(persist)[c(1,4,5)]))
counter = 1
pairs = c()
for(i in 1:2){
  for(j in (i+1):3){
    coex[counter,,,] <- ifelse(persist[,i,j,,] & persist[,j,i,,],  # coexistence = 1
                               1, ifelse(persist[,i,j,,] & !persist[,j,i,,], # i excludes j = 2
                                         2, ifelse(!persist[,i,j,,] & persist[,j,i,,], # j excludes i = 3
                                                   3, ifelse(!persist[,i,j,,] & !persist[,j,i,,], 4, NA))))  # pe = 4
    counter <- counter + 1
    pairs = c(pairs, paste(i,j, sep = "."))
  }
}

# what proportion of samples are predicted to coexiste?
apply(coex, c(1,3,4), function(x) sum(x == 1, na.rm = TRUE)/
        length(x))

# what proportion of samples are predicted to exclude?
apply(coex, c(1,3,4), function(x) sum(x %in% c(2,3), na.rm = TRUE)/
        length(x))

# what proportion of samples are predicted to experience priority effects?
apply(coex, c(1,3,4), function(x) sum(x == 4, na.rm = TRUE)/
        length(x))

# what proportion of samples are NA?
apply(coex, c(1,3,4), function(x) sum(is.na(x))/
        length(x))

# get probability of nitrogen switching coexistence outcome

sum(coex[,,,1] != coex[,,,2] & !is.na(coex[,,,1]) & !is.na(coex[,,,2]))/
  sum(!is.na(coex[,,,1]) & !is.na(coex[,,,2]))

apply(coex, 1, FUN = function(x) sum(x[,,1] != x[,,2] & !is.na(x[,,1]) & !is.na(x[,,2]))/
        sum(!is.na(x[,,1]) & !is.na(x[,,2])))

# get probability of rhizobia switching coexistence outcome

sum(coex[,,1,] != coex[,,2,] & !is.na(coex[,,1,]) & !is.na(coex[,,2,]))/
  sum(!is.na(coex[,,1,]) & !is.na(coex[,,2,]))

apply(coex, 1, FUN = function(x) sum(x[,1,] != x[,2,] & !is.na(x[,1,]) & !is.na(x[,2,]))/
        sum(!is.na(x[,1,]) & !is.na(x[,2,])))

# get probability of interaction effect on coexistence outcome

sum(coex[,,1,2] != coex[,,2,2] & coex[,,1,1] == coex[,,2,1] &
      !is.na(coex[,,1,1]) & !is.na(coex[,,2,1]) &
      !is.na(coex[,,1,2]) & !is.na(coex[,,2,2]))/
  sum(!is.na(coex[,,1,1]) & !is.na(coex[,,2,1]) &
        !is.na(coex[,,1,2]) & !is.na(coex[,,2,2]))

apply(coex, 1, FUN = function(x) sum(x[,1,2] != x[,2,2] & x[,1,1] == x[,2,1] &
                                       !is.na(x[,1,1]) & !is.na(x[,2,1]) &
                                       !is.na(x[,1,2]) & !is.na(x[,2,2]))/
        sum(!is.na(x[,1,1]) & !is.na(x[,2,1]) &
              !is.na(x[,1,2]) & !is.na(x[,2,2])))
