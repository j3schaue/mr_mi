library(metafor)
library(mice)
library(tidyverse)
library(metami)
# source("./metami/R/initialize.R")
# source("./metami/R/draw_parameters.R")
# source("./metami/R/impute_x.R")
# source("./metami/R/mi.R")
# cl <- parallel::makeCluster(2)
# doParallel::registerDoParallel(cl)

nsims <- 5
k <- 60
tau2 <- 0 #1/50
beta0 <- -0.05
beta1 <- 0.2
beta2 <- 0.15
m <- 10

Sigma <- matrix(c(1, .3, .3, 1), ncol = 2)

# set.seed(5897)
dat <- tibble(n = sample(60:300, k, replace = TRUE)) %>%
  mutate(v = 4/n) %>%
  select(v)
  

res <- list()
for(i in 1:nsims){
  
  # ampute
  xtib <- as_tibble(matrix(MASS::mvrnorm(k, c(-0.2, 0.2), Sigma), ncol = 2)) %>%
    mutate(#x1 = ifelse(V1 > -.1, 1, 0),
           x2 = ifelse(V2 > 0.15, 1, 0)) %>%
    select(x1 = V1, x2)
  
  dat_sim <- dat %>%
    bind_cols(xtib) %>%
    mutate(t = rnorm(k, beta0 + beta1 * x1 + beta2 * x2, sqrt(tau2 + v))) %>%
    select(t, x1, x2, v) 
  
  amp_obj <- mice::ampute(dat_sim, prop = 0.3, 
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
  
  imp_mid <- MetaMice(t = dat_amp %>% select(t), 
                      x = dat_amp %>% select(x1, x2), 
                      v = dat_amp %>% select(v), 
                      m = m, 
                      variables = c("x1", "x2") ,
                      types = c("continuous", "binary"), 
                      fixed = TRUE)
  
  res[[i]] <- pool(with(imp_mid, rma(yi = t, vi = v, mods = ~ x1 + x2, method = "FE")))
  
}

results <- lapply(seq_along(res), 
                  FUN = function(i) res[[i]]$pooled %>% mutate(sim = i)) %>%
  bind_rows()

ggplot(results) + 
  stat_density(aes(estimate)) +
  facet_grid(~term)

results %>%
  group_by(term) %>%
  summarize(mn = mean(estimate))