
########## DATA CLEANING ##########

set.seed(6)

library(tidyverse)
library(readxl)

# Load in data
fec_dat <- read_excel(file.path("Data", "NRdatMP.xlsx"), 
                      sheet = "fecundity")
germ_dat <- read_excel(file.path("Data", "NRdatMP.xlsx"), 
                       sheet = "germdat")
seeds <- read_excel(file.path("Data", "NRdatMP.xlsx"), 
                    sheet = "seedavgerages", 
                    col_names = FALSE)


# Fix case inconsistency 
fec_dat$competitor <- toupper(fec_dat$competitor)

# Convert pot-level germ data to seed-level germ data

germ_dat2 <- pmap_dfr(germ_dat, 
                      function(pot_ID, nitrogen, rhizobia, 
                               focal, competitor, density, 
                               sp, baseline,
                               n.germ, n.seeded) {
                        data.frame(pot_ID = pot_ID, 
                                   nitrogen = nitrogen,
                                   rhizobia = rhizobia,
                                   focal = focal, 
                                   competitor = competitor, 
                                   density = density, 
                                   sp = sp,
                                   baseline = baseline, 
                                   germ = c( rep(1, n.germ),
                                             rep(0, n.seeded - n.germ) ) )
                      }
)

# Get total seed output using mean seed counts
# Note: not very reproducible... fix seeds data file and code at some point

fec_dat$fit <- NA

fec_dat$fit <- ifelse(fec_dat$sp == "A" & fec_dat$nitrogen == "nitrogen", 
                      fec_dat$flower_count * as.numeric(seeds[1, 2]), 
                      ifelse(fec_dat$sp == "A" & fec_dat$nitrogen == "zero nitrogen", 
                             fec_dat$flower_count * as.numeric(seeds[2, 2]), fec_dat$fit ))

fec_dat$fit <- ifelse(fec_dat$sp == "B" & fec_dat$nitrogen == "nitrogen", 
                      fec_dat$flower_count * as.numeric(seeds[3, 2]), 
                      ifelse(fec_dat$sp == "B" & fec_dat$nitrogen == "zero nitrogen", 
                             fec_dat$flower_count * as.numeric(seeds[4, 2]), fec_dat$fit ))

fec_dat$fit <- ifelse(fec_dat$sp == "C" & fec_dat$nitrogen == "nitrogen", 
                      fec_dat$flower_count * as.numeric(seeds[5, 2]), 
                      ifelse(fec_dat$sp == "C" & fec_dat$nitrogen == "zero nitrogen", 
                             fec_dat$flower_count * as.numeric(seeds[6, 2]), fec_dat$fit ))
