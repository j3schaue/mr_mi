library(metafor)
library(mice)
library(metami)
library(tidyverse)
library(doParallel)
# source("./metami/R/initialize.R")
# source("./metami/R/draw_parameters.R")
# source("./metami/R/impute_x.R")
# source("./metami/R/mi.R")

# cl <- parallel::makeCluster(3)
# doParallel::registerDoParallel(cl)

nsims <- 2
k <- 60
tau2 <- 0
beta0 <- -0.05
beta1 <- 0.2
beta2 <- 0.15
m <- 10
prop <- 0.1 * 1:2
fixed_effects <- ifelse(tau2 == 0, TRUE, FALSE)

Sigma <- matrix(c(1, .3, .3, 1), ncol = 2)

set.seed(5897)
dat <- tibble(n = sample(60:300, k, replace = TRUE)) %>%
  mutate(v = 4/n) %>%
  select(v)


simulations <- list()
for(j in seq_along(prop)) {#foreach(pp = prop) %dopar% {
  
  pp <- prop[j]
  
  res <- list()
  mi_pmm <- list()
  mi_rf <- list()
  cc_res <- list()
  for(i in 1:nsims){
    
    # ampute
    xtib <- as_tibble(matrix(MASS::mvrnorm(k, c(0, 0), Sigma), ncol = 2)) %>%
      mutate(x2 = ifelse(V2 > 0.07, 1, 0)) %>%
      select(x1 = V1, x2)
    
    dat_sim <- dat %>%
      bind_cols(xtib) %>%
      mutate(t = rnorm(k, beta0 + beta1 * x1 + beta2 * x2, sqrt(tau2 + v))) %>%
      select(t, x1, x2, v) 
    
    amp_obj <- mice::ampute(dat_sim, prop = pp, 
                            patterns = matrix(c(1, 0, 0, 1, 
                                                1, 0, 1, 1,
                                                1, 1, 0, 1), 
                                              ncol = 4, 
                                              byrow = TRUE), 
                            mech = "MAR", 
                            weights = matrix(c(.7, 0, 0, .3, 
                                               .6, 0, .2, .2,
                                               .6, .2, 0, .2), 
                                             ncol = 4, 
                                             byrow = TRUE)) 
    
    dat_amp <- amp_obj$amp
    
    miss_pcts <- apply(dat_amp, 2, FUN = function(x) mean(is.na(x)))
    
    while(any(miss_pcts == 1)){
      
      amp_obj <- mice::ampute(dat_sim, prop = pp, 
                              patterns = matrix(c(1, 0, 0, 1, 
                                                  1, 0, 1, 1,
                                                  1, 1, 0, 1), 
                                                ncol = 4, 
                                                byrow = TRUE), 
                              mech = "MAR", 
                              weights = matrix(c(.7, 0, 0, .3, 
                                                 .6, 0, .2, .2,
                                                 .6, .2, 0, .2), 
                                               ncol = 4, 
                                               byrow = TRUE)) 
      
      dat_amp <- amp_obj$amp
      
      miss_pcts <- apply(dat_amp, 2, FUN = function(x) mean(is.na(x)))
      
    } 
    
    print(paste0("Proportion: ", pp, ". Simulation: ", i))
    
    imp_mid <- suppressMessages(
      MetaMice(t = dat_amp %>% select(t), 
               x = dat_amp %>% select(x1, x2), 
               v = dat_amp %>% select(v), 
               m = m, 
               variables = c("x1", "x2"),
               types = c("continuous", "binary"), 
               fixed = fixed_effects, 
               burn_in = 2)
    )
    
    mid_def <- mice::mice(dat_amp, m = 10)
    mid_rf <- mice::mice(dat_amp, m = 10, method = "rf")
    mid_norm <- mice::mice(dat_amp, m = 10, method = c("norm", "logreg"))
    
    res[[i]] <- pool(with(imp_mid, rma(yi = t, vi = v, mods = ~ x1 + x2, method = method_opt)))
    mi_pmm[[i]] <- pool(with(mid_def, rma(yi = t, vi = v, mods = ~ x1 + x2, method = method_opt)))
    mi_rf[[i]] <- pool(with(mid_rf, rma(yi = t, vi = v, mods = ~ x1 + x2, method = method_opt)))
    cc_res[[i]] <- rma(yi = t, vi = v, mods = ~ x1 + x2, method = "FE", data = dat_amp)
    norm_res[[i]] <- pool(with(mid_norm, rma(yi = t, vi = v, mods = ~ x1 + x2, method = method_opt)))
    
  }

  # results_comp <- lapply(seq_along(res), 
  #                        FUN = function(i) res[[i]]$pooled %>% mutate(sim = i)) %>%
  #   bind_rows() %>%
  #   mutate(prop = pp, 
  #          imputation_type = "compatible")
  # 
  # results_pmm <- lapply(seq_along(mi_pmm), 
  #                       FUN = function(i) mi_pmm[[i]]$pooled %>% mutate(sim = i)) %>%
  #   bind_rows() %>%
  #   mutate(prop = pp, 
  #          imputation_type = "pmm")
  # 
  # results_rf <- lapply(seq_along(mi_rf), 
  #                      FUN = function(i) mi_rf[[i]]$pooled %>% mutate(sim = i)) %>%
  #   bind_rows() %>%
  #   mutate(prop = pp, 
  #          imputation_type = "rf")
  # 
  # results_cc <- lapply(seq_along(cc_res), 
  #                      FUN = function(i) cc_res[[i]]$pooled %>% mutate(sim = i)) %>%
  #   bind_rows() %>%
  #   mutate(prop = pp, 
  #          imputation_type = "cc")
  # 
  # results <- bind_rows(results_comp, results_rf, results_pmm, results_cc)
  simulations[[j]] <- list(compatible = res, 
                           complete_case = cc_res, 
                           mice_rf = mi_rf, 
                           mice_pmm = mi_pmm, 
                           mice_norm = norm_res,
                           prop = c(pp))
}

simulations[[j + 1]] <- dat

# saveRDS(simulations,
#         paste0("./results/simulations_mice_tau2_", tau2, "_normal.RDS")
# )

