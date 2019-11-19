###-------------------------------------------------------###
###-------------------------------------------------------###
###-------------------------------------------------------###
### Some potentially useful functions for generating 
### compatible imputations
###-------------------------------------------------------###
###-------------------------------------------------------###
###-------------------------------------------------------###

###---Load packages
require(fastDummies)
require(tidyverse)
require(rstan)


#-------------------------------------------------------#
#-------------------------------------------------------#
factor_to_dummy <- function(data, var){
  
  other_vars = data %>% select(-var)
  dummy = data %>% select(var) %>%
    fastDummies::dummy_cols(remove_first_dummy = TRUE) %>%
    select(-var) %>%
    cbind(other_vars, .)
  levs = levels(data[[var]])
  
  return(list(data=dummy, levels=levs))

}


dummy_to_factor <- function(data, var, levels){
  
  other_vars = data %>% select_if(!grepl(var, names(data)))
  
  dummy_codes = data %>% 
    select_if(grepl(var, names(data)))
  
  fct_var = factor(as.matrix(dummy_codes) %*% 1:I(length(levels) - 1), labels = levels)
  
  other_vars[[var]] = fct_var
  
  return(other_vars)
}



#-------------------------------------------------------#
#' @name mr_likelihood
#' @param t double, effect estimate
#' @param v double, estimation error variance
#' @param x vector, x value (covariate)
#'          does not include intercept (bias) column
#' @param beta vector of regression coefficients
#' @param tau2 double, random effect variance
#' @return value of the likelihood function
#' @note the vector 'x' does not include an intercept
#'       column for input. This is appended afterward.
#-------------------------------------------------------#
mr_likelihood = function(t, v, x, beta, tau2){
  
  # Set up design vector
  x = c(1, x)
  
  # Check that beta and x have the same dimension
  if(length(x) != length(beta)){
    
    # If dimensions don't match, kill the process
    print("You need x to have length of length(beta) - 1.")
    return(NULL)
    
  } else {
    
    # If dimensions match, get the mean for t
    mu <- beta %*% x
    
    # Return the value of the likelihood
    return(dnorm(t, mean = mu, sd = sqrt(tau2 + v)))
  }
  
}


#-------------------------------------------------------#
#' @name prop_x_cat
#' @param t, double, effect estimate
#' @param v, double, effect estimate variance
#' @param beta, vector, regression coefficients
#' @param tau2, double, variance component
#' @return vector of probabilities for each level that
#'         x could take. 
#'         Assumes x is categorical with p levels
#' @note DOES NOT INCLUDE PRIOR INFO ON X. 
#' @note ONLY WORKS FOR A SINGLE COVARIATE!
#-------------------------------------------------------#
prop_x_cat = function(t, v, beta, tau2){
  
  # Set up the possible values X could take
  p = length(beta) - 1
  xs = rbind(rep(0, p), diag(1, p))
  # Note that we do not add a unit vector here.
  # This gets added in the likelihood function.
  
  # Get the likelihood for each category
  prop = sapply(1:nrow(xs), 
                FUN = function(i){
                  mr_likelihood(t, v, xs[i,], beta, tau2)
                })
  
  return(prop/sum(prop))
}


#-------------------------------------------------------#
#' @name draw_x_cat
#' @param t, double, effect estimate
#' @param v, double, effect estimate variance
#' @param beta, vector, regression coefficients
#' @param tau2, double, variance component
#' @return vector drawn from the posterior predictive
#'         distribution of X.
#' @note This returns a vector W/O intercept. 
#-------------------------------------------------------#
draw_x_cat = function(t, v, beta, tau2){
  
  prop = prop_x_cat(t, v, beta, tau2)
  xvec = t(rmultinom(1, 1, prop))
  xvec[1,1] = 0 # Set first entry to zero as you
                # would for a design matrix 
                # w/dummy-coded variables
  
  return(xvec)
}


#-------------------------------------------------------#
#-------------------------------------------------------#
draw_x_cts <- function(t, v, beta, tau2){
  
  z = rnorm(1, beta[1] - t, sqrt(v + tau2))
  return(-z/beta[1])
  
}


#-------------------------------------------------------#
#' @name get_posterior_stan
#' @param t, vector of effect estimates (length k)
#' @param x, matrix of predictors (k x p)
#' @param v, vector of effect variances (length k)
#' @param fixed, boolean for whether we are doing a 
#'           fixed- or random-effects model
#'           default is FALSE.
#' @note x is NOT a design matrix! It is a (dummy-coded)
#'       matrix of predictors.
#-------------------------------------------------------#
get_posterior_stan = function(t, x, v, fixed = FALSE){
  
  # Get dimensions of the design matrix
  k = nrow(x); p = ncol(x)
  
  # Write out the MR model in STAN
  if(fixed){
    
    # Fixed-effects model in STAN
    stan_mod = "
      data {
        int k;
        int p;
        vector[k] t;
        vector<lower=0>[k] v;
        matrix[k, p] x;
      }
      
      parameters {
        real alpha; // intercept
        vector[p] beta; // regression coefficient
      }
      
      model {
        t ~ normal(x * beta + alpha, sqrt(v));
      }
  	"
    
  } else {
    
    # Random-effects model in STAN
    stan_mod = "
      data {
        int k;
        int p;
        vector[k] t;
        vector<lower=0>[k] v;
        matrix[k, p] x;
      }
      
      parameters {
        vector[k] theta; // effect parameter
        real alpha; // intercept
        vector[p] beta; // regression coefficient
        real<lower=0> tau; // between study variance
      }
      
      model {
        t ~ normal(theta, sqrt(v));
        theta ~ normal(x * beta + alpha, tau);
      }
  	"
  }
  
  # Run the MR model
  data = list(t = t, v = v, x = x, k = k, p = p)
  fit <- stan(model_code = stan_mod, data = data, 
              iter = 1000, chains = 4)
  
  # Extract the posterior
  posterior_dist = as.data.frame(rstan::extract(fit))
  
  return(posterior_dist)
}


#-------------------------------------------------------#
#-------------------------------------------------------#
prep_df <- function(data, categorical = TRUE){
  
  # Sort out the data that is missing vs. not missing
  comp_df = data %>% drop_na()
  mis_df = data %>% setdiff(comp_df)
  inds = which(is.na(data$x)) # get missing data indices
  
  # Select only the predictors
  comp_df_x = comp_df %>% 
    select(-y, -v)
  
  if(categorical){
    
    if(sort(unique(comp_df$x)) != c(0, 1)){
      ## STUFF FOR DUMMY CODING
      
      # Get the # of predictors
      nx = ncol(comp_df_x)
      
      # Create dummy codes
      comp_df_x = comp_df_x %>% 
        fastDummies::dummy_cols(remove_first_dummy = TRUE)
      
      # Get updated column count
      ncx = ncol(comp_df_x)
      
      # Select only the dummy-coded variables
      comp_df_x = comp_df_x[, (nx + 1):ncx] %>%
        as.matrix()
    }
    
  }
    
  return(list(comp_df = comp_df, mis_df = mis_df, comp_df_x = comp_df_x, inds = inds))
}


#-------------------------------------------------------#
#' @name impute_x_onevar
#' @param df, data frame
#' @param m, number of imputations
#' @param fixed, boolean, indicates fixed-effects model
#' @param true_beta, boolean, simulation check allows
#'             us to use the true tau2
#' @param beta, vector             
#-------------------------------------------------------#
impute_x_onevar = function(data, m = 5, fixed = FALSE,
                    categorical = TRUE, 
                    beta = NULL, tau2 = NULL){
  
  ###---Prep the DF
  prepped <- prep_df(data)
  comp_df_x = prepped$comp_df_x
  comp_df = prepped$comp_df
  mis_df = prepped$mis_df
  inds = prepped$inds

    
  ###---Get the posterior of the MR parameters
  post = get_posterior_stan(t = comp_df$y, 
                            x = comp_df_x, 
                            v = comp_df$v, 
                            fixed = fixed)
  
  if(fixed){
    
    # Extract parameters into a DF
    params = post %>%
      sample_n(m) %>%
      select_if(grepl("alpha|beta", names(.))) %>%
      mutate(tau2 = 0)
    
  } else {
    
    params = post %>%
      sample_n(m) %>%
      mutate(tau2 = tau^2) %>%
      select_if(grepl("tau2|beta|alpha", names(.)))
    
  }
    
  
  ###---Generate imputations
  dfs = list()
  
  # Posterior draws of coefficients as a matrix
  post_betas = params %>% 
    select_if(grepl("alpha|beta", names(.))) %>%
    as.matrix()
  
  # Generate m imputed datasets
  for(i in 1:m){
    
    # get the data that are missing covariates
    tmp = mis_df %>%
      select(y, v) 
    
    if(categorical){
      # Draw all of the xs as a matrix
      xMatrix = sapply(1:nrow(tmp), 
                       FUN = function(j){
                         draw_x_cat(t = tmp$y[j], 
                                    v = tmp$v[j], 
                                    beta = post_betas[i,],
                                    tau2 = params[i, "tau2"])}
                       ) %>%
        t()
      
      # Convert the matrix of imputed dummy variables into a factor
      Xvals = sapply(1:nrow(xMatrix), 
                     FUN = function(k){
                       xx = xMatrix[k,]
                       if(sum(xx) == 0){
                         return(0)
                       } else{
                         return(which(xx == 1) - 1)
                       }
                     })
    } else {
      
      Xvals = sapply(1:nrow(tmp), 
                     FUN = function(i){
                       draw_x_cts(t = tmp$y[j], 
                                  v = tmp$v[j], 
                                  beta = post_betas[i,], 
                                  tau2 = params[i, "tau2"])
                     })
      
    }
      
    imp_df = tmp %>%
      mutate(x = Xvals)
    
    dfs[[i]] = bind_rows(comp_df, imp_df)
  }
  
  return(dfs)
  
}
