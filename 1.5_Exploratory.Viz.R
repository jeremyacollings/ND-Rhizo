
########## EXPLORATORY DATA VISUALIZATION ##########

set.seed(6)

library(tidyverse)

source("1_Data.Cleaning.R")

ggplot(data = fec_dat) + 
  geom_histogram(aes(x = fit, fill = sp), position = "identity", 
                 alpha = .6) + 
  theme_classic()

ggplot(data = fec_dat, aes(x = nitrogen, y = fit, color = rhizobia)) + 
  geom_jitter(position = position_jitterdodge()) + geom_boxplot() + facet_grid(~ sp) + 
  theme_classic()

ggplot(data = fec_dat, aes(x = rhizobia, y = fit, color = nitrogen)) + 
  geom_jitter(position = position_jitterdodge()) + geom_boxplot() + facet_grid(~ sp) + 
  theme_classic()

# these figs are simplified looks into the data and should not
# perfectly reflect the models we'll fit

fec_plots <- list()
for(foc in c("A", "B", "C")){
  fec_plots[[foc]][["A"]] <- ggplot(data = fec_dat[which(fec_dat$sp == foc & fec_dat$competitor %in% c("A", "ALONE")),], 
                                    aes(x = A_comp, y = fit)) + 
    geom_jitter(width = .25) + 
    stat_smooth(method = "nls",
                formula = y ~ a/(1+b*x),
                method.args = list(start = list(a = mean(fec_dat$fit[1:nrow(fec_dat)], na.rm = TRUE), b = .1)),
                se = FALSE, na.rm = TRUE) + theme_classic() +
    xlab("Competitor Density") + 
    facet_grid(nitrogen ~ rhizobia) + 
    ggtitle(paste(foc, "A", sep = "-"))
  
  fec_plots[[foc]][["B"]] <- ggplot(data = fec_dat[which(fec_dat$sp == foc & fec_dat$competitor %in% c("B", "ALONE")),], 
                                    aes(x = B_comp, y = fit)) + 
    geom_jitter(width = .25) + 
    stat_smooth(method = "nls",
                formula = y ~ a/(1+b*x),
                method.args = list(start = list(a = mean(fec_dat$fit[1:nrow(fec_dat)], na.rm = TRUE), b = .1)),
                se = FALSE, na.rm = TRUE) + theme_classic() +
    xlab("Competitor Density") + 
    facet_grid(nitrogen ~ rhizobia) + 
    ggtitle(paste(foc, "B", sep = "-"))
  
  fec_plots[[foc]][["C"]] <- ggplot(data = fec_dat[which(fec_dat$sp == foc & fec_dat$competitor %in% c("C", "ALONE")),], 
                                    aes(x = C_comp, y = fit)) + 
    geom_jitter(width = .25) + 
    stat_smooth(method = "nls",
                formula = y ~ a/(1+b*x),
                method.args = list(start = list(a = mean(fec_dat$fit[1:nrow(fec_dat)], na.rm = TRUE), b = .1)),
                se = FALSE, na.rm = TRUE) + theme_classic() +
    xlab("Competitor Density") + 
    facet_grid(nitrogen ~ rhizobia) + 
    ggtitle(paste(foc, "C", sep = "-"))
}

fec_plots
