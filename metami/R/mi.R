#' Generate imputations for a single covariate
#'
#' @param t a tibble of effect sizes (one column called 't')
#' @param v a tibble of estimation variances (one column called 'v')
#' @param x a tibble of covariates
#' @param m number of imputations to generate
#' @param variable a string indicating which variable in \code(x) to impute
#' @param type string indicating covariate type ('binary', 'continuous')
#' @param fixed boolean indicates fixed/random effects model
#' @param to_impute optional argument, list of row indices to impute
#' @return mids object of imputed datasets
MetaMi <- function(t, x, v,
                   m = 5,
                   variable = NULL,
                   type,
                   fixed = TRUE,
                   to_impute = NULL){

  # Sanitize inputs
  orig_data <- bind_cols(t, v, x) %>%
    mutate(.imp = 0,
           .id = 1:nrow(.))

  if(is.null(variable)){
    variable <- x %>%
      select_if(anyNA) %>%
      names()

    if(length(variable) > 1){
      print("You have missingness in multiple columns. Use MetaMice().")
      return(NULL)
    } else if(length(variable) == 0){
      print("There are no missing covariates.")
      return(NULL)
    }
  }

  # Get complete covariates
  xj_not <- x %>% select(-all_of(variable))

  # Get column to impute
  xj <- x %>% select(all_of(variable))

  if(is.null(to_impute)){
    to_impute <- which(is.na(xj))
  }

  params <- foreach(i = 1:2) %dopar% {
    if(i == 1){
      tmp_params <- GetPostPsi(x, v, variable = variable, type = type, m = m)
    } else {
      tmp_params <- GetPostEta(t, x, v, fixed = fixed, m = m)
    }
    tmp_params
  }

  post_psi <- params[[1]]
  post_eta <- params[[2]]

  if(m > 1){

    out <- foreach(i = 1:m) %dopar% {

      xj_tmp <- xj

      # Draw parameters
      eta <- post_eta %>%
        slice(i)

      psi <- post_psi %>%
        slice(i)

      # Draw imputations
      t_imp <- t %>% slice(to_impute)
      x_imp <- x %>% slice(to_impute)
      v_imp <- v %>% slice(to_impute)
      xj_new <- ImputeX(t = t_imp,
                        x = x_imp,
                        v = v_imp,
                        eta = eta,
                        psi = psi,
                        variable = variable,
                        type = type)

      xj_tmp[[variable]][to_impute] <- xj_new

      # out[[i]] <- bind_cols(t, v, xj_not, xj_tmp) %>%
      #   mutate(.imp = i,
      #          .id = 1:nrow(.))
      bind_cols(t, v, xj_not, xj_tmp) %>%
        mutate(.imp = i,
               .id = 1:nrow(.))

    }
  } else {

    xj_tmp <- xj

    # Draw parameters
    eta <- post_eta %>%
      slice(1)

    psi <- post_psi %>%
      slice(1)

    # Draw imputations
    t_imp <- t %>% slice(to_impute)
    x_imp <- x %>% slice(to_impute)
    v_imp <- v %>% slice(to_impute)
    xj_new <- ImputeX(t = t_imp,
                      x = x_imp,
                      v = v_imp,
                      eta = eta,
                      psi = psi,
                      variable = variable,
                      type = type)

    xj_tmp[[variable]][to_impute] <- xj_new


    out <- bind_cols(t, v, xj_not, xj_tmp) %>%
      mutate(.imp = 1,
             .id = 1:nrow(.))
  }

  imps <- mice::as.mids(bind_rows(orig_data, out))
  return(imps)
}


#' Generate imputations for a multiple covariates
#'
#' @param t a tibble of effect sizes (one column called 't')
#' @param v a tibble of estimation variances (one column called 'v')
#' @param x a tibble of covariates
#' @param m number of imputations to generate
#' @param variables list of strings indicating which variables in \code(x) to impute
#' @param types list of strings indicating covariate types ('binary', 'continuous')
#' @param fixed boolean indicates fixed/random effects model
#' @param burn_in number of burn-in iterations for chained equations
#' @return mids object of imputed datasets
MetaMice <- function(t, x, v, m = 5, variables, types, fixed = TRUE, burn_in = 10){

  orig_data <- bind_cols(t, v, x) %>%
    mutate(.imp = 0,
           .id = 1:nrow(.))

  if(is.null(names(types))){
    names(types) <- variables
  }

  x_vars <- names(x)

  missing_vals <- lapply(x_vars, function(xvar) which(is.na(x[[xvar]])))
  names(missing_vals) <- x_vars

  x_t <- Initialize(x)

  out <- list()

  for(i in 1:(burn_in + m)){
    for(var_name in variables){

      new_dat <- MetaMi(t = t, x = x_t, v = v,
                        variable = var_name,
                        type = types[[var_name]],
                        fixed = fixed,
                        m = 1,
                        to_impute = missing_vals[[var_name]]) %>%
        complete("long") %>%
        filter(.imp == 1)

      x_t <- new_dat %>%
        select(all_of(x_vars))

    }

    if(i > burn_in){

      print(paste("Imputed dataset", i - burn_in, "added to stack"))
      out[[i - burn_in]] <- new_dat %>%
        mutate(.imp = i - burn_in)

    }
  }

  imps <- mice::as.mids(bind_rows(orig_data, out))
  return(imps)

}
