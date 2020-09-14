
#' Draw parameters from partial imputation model
#'
#' @param x a tibble of covariates
#' @param v a tibble with one column (v) of variances
#' @param variable a string indicating which variable in \code(x) to model
#' @param m a number of draws from the posterior to make
#' @param type string indicating 'binary' or 'continuous'
#' @return tibble of draws from auxiliary model
GetPostPsi <- function(x, v, variable, m = 5, type = NULL){

  dat <- bind_cols(x, v)
  fmla <- paste(variable, "~ .")

  if(type == "continuous"){

    post <- suppressMessages(
             rstanarm::stan_lm(fmla,
                               data = dat,
                               na.action = na.omit,
                               prior = rstanarm::R2(location = .05, what = "mean"),
                               seed = 1234)
            )

    psi <- as_tibble(post) %>%
      rename(intrcpt = `(Intercept)`) %>%
      select(-R2, -`log-fit_ratio`) %>%
      sample_n(m)

  } else if(type == "binary"){

    t_prior <- rstanarm::student_t(df = 5, location = 0, scale = 10)

    post <- suppressMessages(
             rstanarm::stan_glm(fmla,
                                data = dat,
                                family = binomial(link = "logit"),
                                prior = t_prior, prior_intercept = t_prior,
                                seed = 1234)
    )

    psi <- as_tibble(post) %>%
      rename(intrcpt = `(Intercept)`) %>%
      sample_n(m)

  }

  return(psi)

}


#' Draw parameters from substantive model
#'
#' @param t a tibble with one column (t) of effect estimates
#' @param x a tibble of covariates
#' @param v a tibble with one column (v) of estimation variances
#' @param m the number of posterior draws to generate
#' @param fixed boolean indicates fixed/random effects model
#' @return tibble of posterior draws from a MR model
GetPostEta <- function(t, x, v, m = 5, fixed = TRUE){

  dat <- bind_cols(t, v, x) %>%
    mutate(study = 1:nrow(.))

  x_vars <- names(x)

  if(fixed){# For fixed-effects models, draw beta only.
    # under flat prior, beta is normal with mean hat(beta)
    # and covariance matrix vcov(beta)

    # fit regression mod
    mod <- metafor::rma(yi = t, vi = v,
                        mods = ~ . - t - v - study,
                        data = dat,
                        method = "FE")

    # extract coeficient estimates and covariance matrix
    mu_vec <- coef(mod)
    sigma_mat <- vcov(mod)

    # draw beta ~ N(beta_hat, VC(beta_hat))
    beta <- MASS::mvrnorm(n = m, mu = mu_vec, Sigma = sigma_mat)

    # set param draws
    if(m == 1){
      eta <- as_tibble_row(beta) %>%
        mutate(tau2 = 0)
    } else {
      eta <- as_tibble(beta)  %>%
        mutate(tau2 = 0)
    }

  } else {

    dat <- dat %>%
      mutate(se_t = sqrt(v)) %>%
      select(-v)

    mod <- brms::brm(t | se(se_t) ~ . + (1 | study),
                      data = dat)

    samps <- suppressMessages(
        brms::posterior_samples(mod) %>%
          sample_n(m) %>%
          mutate(intrcpt = b_Intercept,
                 tau2 = sd_study__Intercept^2)
    )

    beta <- samps %>% select("intrcpt", paste0("b_", x_vars))
    names(beta) <- gsub("b_", "", names(beta))
    tau2 <- samps %>% select(tau2)

    eta <- as_tibble(c(beta, tau2))

  }

  return(eta)
}
