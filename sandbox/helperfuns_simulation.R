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
