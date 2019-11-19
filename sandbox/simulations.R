###--------------------------------------------------------###
###--------------------------------------------------------###
### Initial simulations for amputing/imputing
###--------------------------------------------------------###
###--------------------------------------------------------###

library(mice)
library(metafor)
library(tidyverse)

source("./sandbox/helperfuns_simulation.R")
source("./sandbox/imputation_functions.R")

###---Get complete df
bd <- metafor::dat.bangertdrowns2004

bd <- bd %>% 
  select(y = yi, v = vi, x = feedback) %>%
  na.omit()

##-----------MCAR SIMS
nsims = 2
m = 5
props = c(0.2, 0.5, 0.75)
mcar_results = list()

for(i in seq(props)){
  mcar_results[[i]] = list()
  prop = props[i]
  
  for(j in 1:nsims){
    
    # Step 1: Ampute data
    amp_obj = ampute(bd, mech = "MCAR",
                     patterns = c(1, 1, 0),
                     prop = prop)
    
    data = amp_obj$amp
    
    # Step 2: Complete cases
    cc_mod = rma(y, v, mods = ~ x, method = "PM", data = data)
    
    # Step 3: Mice default
    eval(metafor:::.mice)
    mice_def <- mice(data)
    midef_mods = with(mice_def, 
                      rma(y, v,
                          mods = ~ x,
                          method = "PM"))
    mice_def = pool(midef_mods)
    
    # Step 4: Mice RF
    mice_rf <- mice(data, method = "rf")
    mirf_mods = with(mice_rf, 
                      rma(y, v,
                          mods = ~ x,
                          method = "PM"))
    mice_rf = pool(mirf_mods)
    
    # Step 5: Compatible imputations
    compat = impute_x_onevar(data)
    compat_mods = lapply(1:m, FUN=function(i)
                       rma(y, v,
                           mods = ~ x,
                           method = "PM", data = compat[[i]]))
    comp = pool(compat_mods)
    
    mcar_results[[i]][[j]] = list(prop = prop, 
                             m = m, 
                             cc = cc_mod, 
                             mi_def = mice_def, 
                             mi_rf = mice_rf, 
                             comp = comp)
    
  }
}




##-----------MCAR SIMS
mar_results = list()

for(i in seq(props)){
  mar_results[[i]] = list()
  prop = props[i]
  
  for(j in 1:nsims){
    
    # Step 1: Ampute data
    amp_obj = ampute(bd, mech = "MAR",
                     patterns = c(1, 1, 0),
                     prop = prop)
    
    data = amp_obj$amp
    
    # Step 2: Complete cases
    cc_mod = rma(y, v, mods = ~ x, method = "PM", data = data)
    
    # Step 3: Mice default
    eval(metafor:::.mice)
    mice_def <- mice(data)
    midef_mods = with(mice_def, 
                      rma(y, v,
                          mods = ~ x,
                          method = "PM"))
    mice_def = pool(midef_mods)
    
    # Step 4: Mice RF
    mice_rf <- mice(data, method = "rf")
    mirf_mods = with(mice_rf, 
                     rma(y, v,
                         mods = ~ x,
                         method = "PM"))
    mice_rf = pool(mirf_mods)
    
    # Step 5: Compatible imputations
    compat = impute_x_onevar(data)
    compat_mods = lapply(1:m, FUN=function(i)
      rma(y, v,
          mods = ~ x,
          method = "PM", data = compat[[i]]))
    comp = pool(compat_mods)
    
    mar_results[[i]][[j]] = list(prop = prop, 
                                 m = m, 
                                 cc = cc_mod, 
                                 mi_def = mice_def, 
                                 mi_rf = mice_rf, 
                                 comp = comp)
    
  }
}

