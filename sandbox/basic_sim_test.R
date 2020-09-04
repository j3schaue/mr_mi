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

nsims <- 200
k <- 60
tau2 <- 1/50
beta0 <- -0.1
beta1 <- 0.2
gamma0 <- -1
gamma1 <- 40
m <- 10

# set.seed(5897)
dat <- tibble(n = sample(60:300, k, replace = TRUE)) %>%
  mutate(v = 4/n, 
         x = rbinom(k, 1, .6),#rbinom(k, 1, exp(gamma0 + gamma1 * v)/(1 + exp(gamma0 + gamma1 * v))),
         t = rnorm(k, beta0 + beta1 * x, sqrt(tau2 + v))) %>%
  select(v) 

res <- list()
for(i in 1:nsims){
  
  # ampute
  dat_sim <- dat %>%
  mutate(x = rnorm(k, v, 1),#rbinom(k, 1, exp(gamma0 + gamma1 * v)/(1 + exp(gamma0 + gamma1 * v))),
         t = rnorm(k, beta0 + beta1 * x, sqrt(tau2 + v))) %>%
    select(t, x, v)

  amp_obj <- mice::ampute(dat_sim, prop = 0.3, 
                          patterns = matrix(c(1, 0, 1), ncol = 3), 
                          mech = "MAR") 
  dat_amp <- amp_obj$amp

  imp_mid <- MetaMi(t = dat_amp %>% select(t), 
                    x = dat_amp %>% select(x), 
                    v = dat_amp %>% select(v), 
                    m = m, 
                    variable = "x", 
                    type = "continuous", 
                    fixed = FALSE)
  
  res[[i]] <- pool(with(imp_mid, rma(yi = t, vi = v, mods = ~ x, method = "FE")))
    
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
