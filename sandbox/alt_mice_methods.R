###--------------------------------------------------------###
###--------------------------------------------------------###
### Practice simulations for amputing/imputing
###--------------------------------------------------------###
###--------------------------------------------------------###

### Load libraries
library(mice) # for multiple imputation 
library(metafor) # for meta-regression
library(tidyverse) # for data wrangling
library(randomForest)


### Get the dataset
bd <- metafor::dat.bangertdrowns2004

bd <- bd %>% 
  select(yi, vi, feedback, wic, imag) %>%
  na.omit()

### Fit MR model
full_model <- rma(yi, vi, 
                  mods = ~ feedback + wic + imag, 
                  method = "PM", 
                  data = bd)

### Ampute data
# This is data were some of the covariates are missing

# First we need to describe missingness patterns
# we will have six patterns, one where each covariate is missing
# and one where two are missing at a time
miss_pattern = matrix(c(rep(1, 6), 
                        rep(1, 6), 
                        c(1, 0, 0, 1, 1, 0), 
                        c(0, 1, 0, 1, 0, 1), 
                        c(0, 0, 1, 0, 1, 1)), 
                      ncol = 5)
miss_pattern

# Use the ampute() function to delete variabes
bd_mis_obj <- ampute(bd, # data to be amputed
                     prop = 0.5, # proportion of missingness we want
                     patterns = miss_pattern, # missingness patterns
                     mech = "MCAR") # data missing mcar. You can set this to "MAR" or "MNAR"

# Extract the data frame with missingness
bd_mis <- bd_mis_obj$amp

# Complete case model
cc_model <- rma(yi, vi, 
                mods = ~ feedback + wic + imag, 
                method = "PM", 
                data = bd_mis)

### Run MICE and analyze data
# Here we impute the missing values with mice
bd_mice <- mice(bd_mis, printFlag = FALSE)

bd_mice_rf <- mice(bd_mis, method = "rf", , printFlag = FALSE)

# fit models to imputed datasets
mi_mods <- with(bd_mice, 
                rma(yi, vi,
                    mods = ~ feedback + wic + imag,
                    method = "PM"))

# pool the estimates using Rubin's rules
eval(metafor:::.mice)
pooled_models <- pool(mi_mods) 
# UH OH! There's an error when we call pool()--a MICE function--mi_mods--a list of metafor objects
# Try updating metafor (since it's supposed to work with MICE)
# devtools::install_github("wviechtb/metafor")
# Github version doesn't work!
# Quick function for pooling models in  helperfuns file
 source("./sandbox/helperfuns_simulation.R")
 pooled_models <- pool_meta(mi_mods)


# loop --------------------------------------------------------------------

bd_mice_rf <- mice(bd_mis, method = "rf", , printFlag = FALSE)
 
# Received the following warning: 
# In randomForest.default(x = xobs, y = yobs, ntree = 1,  ... :
#                           The response has five or fewer unique values.  Are you sure you want to do regression?
#                           
# fit models to imputed datasets
 
 mi_mods_rf <- with(bd_mice_rf, 
                 rma(yi, vi,
                     mods = ~ feedback + wic + imag,
                     method = "PM"))
 
 # pool the estimates using Rubin's rules
 source("./sandbox/helperfuns_simulation.R")
 pooled_models_rf <- pool_meta(mi_mods_rf)
 

# end loop -----------------------------------------------------------------

  
### Compare MICE to MR model
full_model
pooled_models_rf
pooled_models
cc_model

# Beatrice's misc
# cc_model_df <- tibble(estimate = cc_model$beta, 
#                       se = cc_model$se,
#                       zval = cc_model$zval,
#                       pval = cc_model$pval,
#                       ci.lb = cc_model$ci.lb,
#                       ci.ub = cc_model$ci.ub)



