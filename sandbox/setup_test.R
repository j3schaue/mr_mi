library(tidyverse)
# library(doParallel)

# registerDoParallel(cores = 4)

bd = metafor::dat.bangertdrowns2004 %>%
  select(t = yi, v = vi, feedback, wic)

t = bd %>% select(t)
v = bd %>% select(v)
x = bd %>% select(wic)
missing <- "wic"
variable <- "wic"
fixed <- FALSE
type = "binary"
m = 5
variables <- c("wic", "feedback")
types <- c("binary", "continuous")
burn_in = 3

# source("./R/initialize.R")
# source("./R/draw_parameters.R")
# source("./R/impute_x.R")
# source("./R/mi.R")
# # t = t %>% slice(to_impute)
# # v = v %>% slice(to_impute)
# # x = x %>% slice(to_impute)


mm1 <- MetaMi(t = t, x = x %>% select(wic), v = v, m = m, variable = "wic", type = "binary", fixed = TRUE)
mm2 <- MetaMi(t = t, x = x %>% select(wic), v = v, m = m, variable = "wic", type = "binary", fixed = FALSE)
mm3 <- MetaMi(t = t, x = x %>% select(wic), v = v, m = m, variable = "wic", type = "continuous", fixed = TRUE)
mm4 <- MetaMi(t = t, x = x %>% select(wic), v = v, m = m, variable = "wic", type = "continuous", fixed = FALSE)



mmi1 <- MetaMice(t = t, x = x, v = v, m = m,
                 variables = c("wic", "feedback"), types = c("binary", "continuous"), fixed = TRUE,
                 burn_in = burn_in)
mmi2 <- MetaMice(t = t, x = x, v = v, m = m,
                 variables = c("wic", "feedback"), types = c("binary", "binary"), fixed = FALSE,
                 burn_in = burn_in)

pool(with(mmi2, rma(yi = t, vi = v, mods = ~ wic + feedback)))
rma(yi = t, vi = v, mods = ~ wic + feedback, data = bd)
