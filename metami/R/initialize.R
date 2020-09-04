#' Initialize dataset for chained imputations
#'
#' @param x a tibble of covariates
#' @return the same tibble \code(x) with missing values filled in
Initialize <- function(x){

  is_missing <- x %>%
    select_if(anyNA) %>%
    names()

  out <- x %>%
    tidyr::fill(all_of(is_missing), .direction = "updown")

  return(out)
}
