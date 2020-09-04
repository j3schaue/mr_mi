
#' Draw new values of a single covariate (given other covariates)
#'
#' @param t a tibble of effect sizes (one column called 't')
#' @param x a tibble of covariates
#' @param v a tibble of estimation variances (one column called 'v')
#' @param eta tibble of draws from the posterior of the MR model
#' @param psi tibble of draws from the aux model
#' @param variable string indicating which covariate to impute
#' @param type string indicating covariate type ('binary', 'continuous')
#' @return tibble with imputed values for \code(variable) in \code(x)
ImputeX <- function(t, x, v, eta, psi, variable, type){

  # standardize regression parameters
  tau2 <- eta %>% select(tau2)
  beta <- eta %>% select(-tau2)

  dat <- bind_cols(t, x, v) %>%
    mutate(row = 1:nrow(.))

  if(type == "binary"){

    reg_dat <- dat %>%
      mutate(intrcpt = 1,
             !!variable := 1) %>%
      select(t, names(beta))

    xb1 <- tibble(xb1 = as.vector(as.matrix(reg_dat %>% select(-t)) %*% t(as.matrix(beta))))

    reg_dat <- dat %>%
      mutate(intrcpt = 1,
             !!variable := 0) %>%
      select(t, names(beta))

    xb0 <- tibble(xb0 = as.vector(as.matrix(reg_dat %>% select(-t)) %*% t(as.matrix(beta))))

    aux_dat <- dat %>%
      mutate(intrcpt = 1) %>%
      select(names(psi))

    xp <- tibble(xp = as.vector(as.matrix(aux_dat) %*% t(as.matrix(psi))))

    imputed_vals <- reg_dat %>%
      bind_cols(xb1, xb0, xp) %>%
      bind_cols(dat %>% select(v), tau2) %>%
      mutate(p1_reg = dnorm(t, xb1, sqrt(tau2 + v)),
             p0_reg = dnorm(t, xb0, sqrt(tau2 + v)),
             prob_reg = p1_reg/(p1_reg + p0_reg),
             p1_aux = exp(xp)/(1 + exp(xp)),
             p0_aux = 1 - p1_aux,
             p1 <- prob_reg * p1_aux,
             p0 <- (1 - prob_reg) * p0_aux,
             prob_1 = p1/(p1 + p0),
             # prob_1 = prob_reg,
             imputation = rbinom(nrow(.), 1, prob_1))

  } else if(type == "continuous"){

    reg_dat <- dat %>%
      mutate(intrcpt = 1) %>%
      select(names(beta), t, -all_of(variable))

    beta_j <- beta %>%
      select(beta_j = all_of(variable))

    beta_reg <- beta %>%
      mutate(t = -1) %>%
      select(-all_of(variable))

    xb <- tibble(xb = as.vector(as.matrix(reg_dat) %*% t(as.matrix(beta_reg))))

    psi_df <- psi %>%
      select(-sigma)

    aux_dat <- dat %>%
      mutate(intrcpt = 1) %>%
      select(names(psi_df))

    xp <- tibble(xp = as.vector(as.matrix(aux_dat) %*% t(as.matrix(psi_df))))

    imputed_vals <- dat %>%
      select(-all_of(variable)) %>%
      bind_cols(xb, xp,
                psi %>% select(sigma),
                tau2,
                beta_j) %>%
      mutate(mu_reg = xb/(-beta_j),
             v_reg = (tau2 + v)/(beta_j^2),
             mu_aux = xp,
             v_aux = sigma^2,
             muimp = (mu_reg/v_reg + mu_aux/v_aux)/(1/v_reg + 1/v_aux),
             vimp = 1/(1/v_reg + 1/v_aux),
             imputation = rnorm(nrow(.), muimp, sqrt(vimp)))

    # A <- tibble(A = as.vector(as.matrix(reg_dat) %*% t(as.matrix(beta_reg))))
    #
    # psi_df <- psi %>%
    #   select(-sigma)
    #
    # aux_dat <- dat %>%
    #   mutate(intrcpt = 1) %>%
    #   select(names(psi_df))
    #
    # B <- tibble(B = as.vector(as.matrix(aux_dat) %*% t(as.matrix(psi_df))))
    #
    # imputed_vals <- dat %>%
    #   select(-all_of(variable)) %>%
    #   bind_cols(A, B,
    #             psi %>% select(sigma),
    #             tau2,
    #             beta %>% select(all_of(variable))) %>%
    #   mutate(C = !!sym(variable) * sigma,
    #          D = sqrt(tau2 + v),
    #          lambda = sigma^2 * D^2,
    #          norm_var = sigma^2 * D^2,
    #          den = C^2 + D^2,
    #          M = (C * sigma * A + D^2 * B)/den,
    #          imputation = rnorm(nrow(.), M, sqrt(lambda/den)))

  }



  return(imputed_vals$imputation)

}

