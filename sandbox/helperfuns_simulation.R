###-------------------------------------------------------###
###-------------------------------------------------------###
###-------------------------------------------------------###
### Some potentially useful functions for imputation
### and simulations
###-------------------------------------------------------###
###-------------------------------------------------------###
###-------------------------------------------------------###

library(metafor)
library(mice)
library(tidyverse)

#-----------------------------------------------------------#
#' @name  ampute_df
#' @param data data_frame used in a meta-regression. 
#'             assume yi is in col 1 and vi is in col 2.
#' @param prob double indicating desired pct of missingness
#'             default 40%
#' @param pattern matrix/data_frame of missingness patterns
#'                default to complete yi, vi, but missing xi
#' @param freq vector of #pattern denoting proportion of 
#'             missingness that follows each pattern. 
#'             default to null
#' @param mech string missingness mechanism: 
#'             "MAR", "MCAR", "MNAR"
#' @param weights matrix/data_frame of variable weights for
#'                determining missingness. see MICE 
#' @param bycases logical, denotes if missingness pct is by 
#'                case or across all cells. Default TRUE
#' @param df_only logical for if only the amputed data should 
#'                be included in the output. Default TRUE
#' @return out list that includes amputed data
#'               may return a DF if exclude_orig = FALSE
#-----------------------------------------------------------#
ampute_df <- function(data, 
                      prob = 0.4, 
                      pattern = NULL, 
                      freq = NULL, 
                      mech = "MAR",
                      weights = NULL,
                      bycases = TRUE,
                      df_only = TRUE){
  
  # Set the missingness pattern if it is not specified
  nc <- ncol(data)
  if(is.null(pattern)){
    
    # If they specify a frequency but not a pattern, throw an error message.
    if(!is.null(freq)){
      print("You have not set a pattern for missing data, 
             but you have set a frequency for that pattern. 
             You must set either both to NULL or set freq to NULL.")
      return(NULL)
    }
    
    # Set all possible patterns where only cols 1-2 are NOT missing data
    pattern <- sapply(1:(nc - 1), FUN=function(i){
        pat = rep(0, nc)
        pat[c(1:2, i+1)] = 1
        return(pat)
      }) %>%
      matrix(nrow = nc - 1, byrow = TRUE)
    # This is a matrix with the first two cols as 1s
    # The first row for cols 3+ are all 0
    # Cols 3+ x Rows 2+ are the identify matrix
    
    # Set frequency as each pattern being equally probable
    freq = rep(1/nrow(pattern), nrow(pattern))
  
  }
  
  # Actually delete cells
  amp_df <- mice::ampute(data, prop = prob, patterns = pattern,
                        freq = freq, mech = mech, weights = weights, 
                        bycases = bycases)
  
  # Set output
  if(df_only){
    
    # If excluding original, only return amputed df
    out <- amp_df$amp
    
  } else {
    
    # If including original, stash original and amputed df in a list
    out <- amp_df
    
  }
  
  # Return a value
  return(out)
  
}


###---Example
# dat = metafor::dat.bangertdrowns2004 %>%
#   select(yi, vi, wic, feedback, pers) %>%
#   na.omit()
# 
# foo = ampute_df(dat)
# foo



#-----------------------------------------------------------#
#' @name  pool_meta
#' @param mir object with list of metafor models fit to MICE
#'            imputations
#' @return data frame with pooled parameter estimates
#' @note This is a quick and dirty pooling function that 
#'       only uses marginal variances and not the entire 
#'       covariance matrix! 
#'       You should re-code it later to use the cov. matrix!
#-----------------------------------------------------------#
pool_meta = function(mir){

  ### Create tibble with 3 columns:
  # Column 1 is the analyses run on each MI dataset
  if(is.mira(mir)){
    analysis_results = tibble(analyses = mir$analyses)
  } else {
    analysis_results = tibble(analyses = mir$analyses)
  }
  analysis_results %>% 
    # Column 2 is a dataframe of point estimates for each MI dataset
    mutate(ests = map(.x = analyses, 
                      .f = function(x){
                        data.frame(param = c(row.names(x$beta), "tau2"),
                                   value = c(x$beta, x$tau2), 
                                   stat = "estimate")
                        
                      }), 
           ### Column 3 is a dataframe of variances for each MI dataset
           vars = map(.x = analyses, 
                      .f = function(x){
                        data.frame(param = c(row.names(x$beta), "tau2"),
                                   value = c(x$se^2, x$se.tau2^2), 
                                   stat = "variance")
                      })) 
  
  
  # Create a data frame that takes the mean within-df variances
  within_vals = do.call(rbind, analysis_results$vars) %>%
    group_by(param) %>%
    summarize(Ubar = mean(value))
  
  # Create a data frame that gives the pooled point estimates
  # and takes the variance between datasets
  between_vals = do.call(rbind, analysis_results$ests) %>%
    group_by(param) %>% # for each parameter
    summarize(est = mean(value), # get the mean across dfs
              B = var(value)) # get the between-df variance
  
  # Join the within- and between- results  
  out = left_join(between_vals, within_vals) %>%
    mutate(var = B + B/length(analysis_results$analyses) + Ubar, # formula for pooling variances
           se = sqrt(var)) # get standard error
  
  return(out %>% dplyr::select(param, est, var, se, B, Ubar))
}
