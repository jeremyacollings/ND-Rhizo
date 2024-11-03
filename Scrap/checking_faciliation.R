
library(tidybayes)

par_dat <- spread_draws(RickUB, lam[foc, rhizo, nitro], 
                        alpha[foc, comp, rhizo, nitro])

tapply(par_dat$alpha[which(par_dat$rhizo == 2 & par_dat$nitro == 2)],
       list(par_dat$foc[which(par_dat$rhizo == 2 & par_dat$nitro == 2)],
            par_dat$comp[which(par_dat$rhizo == 2 & par_dat$nitro == 2)]),
       median)

exp(abs(tapply(par_dat$alpha[which(par_dat$rhizo == 2 & par_dat$nitro == 2)],
       list(par_dat$foc[which(par_dat$rhizo == 2 & par_dat$nitro == 2)],
            par_dat$comp[which(par_dat$rhizo == 2 & par_dat$nitro == 2)]),
       median)))

# Probability of Competition

tapply(par_dat$alpha, list(par_dat$foc, par_dat$comp, 
                           par_dat$rhizo, par_dat$nitro), 
       function(x) sum(x > 0)/length(x))

# Probability of Facilitation

tapply(par_dat$alpha, list(par_dat$foc, par_dat$comp, 
                           par_dat$rhizo, par_dat$nitro), 
       function(x) sum(x < 0)/length(x))

# Probability of Rhizobia switching interaction type

p.switch.interaction <- function(df, var, vals = c(1, 2),
                                 interaction = "alpha", 
                                 focal = "foc", competitor = "comp"){
  for(f in unique(df[,focal])){
    for(co in unique(df[,competitor])){
      df_temp <- df[which(df[,focal] == f & 
                            df[,competitor == co]),]
      sum(sign(df_temp[which(df_temp[,var] == vals[1]),interaction]) !=
            sign(df_temp[which(df_temp[,var] == vals[2]),interaction]))/
        nrow(df_temp[which(df_temp[,var] == vals[1]),interaction])
    }
  }
}

p.switch.interaction(par_dat, "rhizo")
