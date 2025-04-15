#' Helper function for getting confidence intervals on J fits
#'
#' @description
#' Used by fitJ to return a confidence interval on a fitted J value.
#' Similar to fitJ, but optimizes a function that has its maximum some distance
#' below the maximum likelihood. That distance is 1.92 for a 95% confidence interval.
#'
#'
#' @param log10J log of J (community size).
#' @param occs matrix of the number of observations in each species at each time.
#' @param ages vector containing the ages of each time slice, in years, from
#' oldest to youngest.
#' @param ML maximum likelihood, passed to the function from fitJ.
#' @param confidence value indicating degree of confidence desired. Default is 95%.
#' @param sampled boolean indicating whether occs represents a sampled
#' community (TRUE) or instead represents true species abundances (FALSE).
#' @param generationtime time between generations, in years.
#'
#' @returns the value of the function at the given logJ.
CIfunc <- function(log10J,occs,ages,ML,confidence=0.95,sampled=TRUE,generationtime=1){
  diff <- qchisq(p=confidence,d=1)/2
  out <- xxprob(log10J,occs,ages,sampled,generationtime)
  return(abs(ML-diff-out))}
