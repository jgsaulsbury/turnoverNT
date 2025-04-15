#' Find the maximum-likelihood value of J for a community change dataset
#'
#' @description
#' Uses optim() to fit the J parameter in Hubbell's neutral theory to a dataset.
#' Optionally, outputs a confidence interval.
#'
#' @param occs matrix of the number of observations in each species at each time.
#' One column for each species, one row for each time slice. Time goes from oldest
#' at the bottom to youngest at the top.
#' @param ages vector containing the ages of each time slice, in years, from
#' oldest to youngest.
#' @param sampled boolean indicating whether occs represents a sampled
#' community (TRUE) or instead represents true species abundances (FALSE).
#' @param generationtime time between generations, in years.
#' @param CI boolean indicating whether or not to calculate a confidence interval
#' on the maximum-likelihood estimate of J.
#' @param searchinterval consider values of log10J within this interval.
#'
#' @returns a list containing "loglik", "J", and (optionally) "CI"
#' @export
#'
#' @examples
#' mat <- matrix(data=c(52,12,160,109,30,401,93,31,355),nrow=3)
#' fitJ(occs=mat,ages=c(200,100,0),CI=TRUE)
fitJ <- function(occs,ages,sampled=TRUE,generationtime=1,CI=FALSE,searchinterval=c(1,9)){
  #error handling
  if(dim(occs)[1] != length(ages)){stop("'ages' must have length equal to the number of rows of 'occs'")}
  op <- stats::optimize(xxprob,interval=c(searchinterval[1],searchinterval[2]),occs=occs,ages=ages,sampled=sampled,generationtime=generationtime,maximum=TRUE)
  out <- list("loglik"=op$objective,"J"=10^op$maximum)
  if(CI){
    left <- stats::optimize(CIfunc,interval=c(searchinterval[1],log10(out$J)),occs=occs,ages=ages,ML=out$loglik,
                            sampled=sampled,generationtime=generationtime,maximum=FALSE)$minimum
    right <- stats::optimize(CIfunc,interval=c(log10(out$J),searchinterval[2]),occs=occs,ages=ages,ML=out$loglik,
                             sampled=sampled,generationtime=generationtime,maximum=FALSE)$minimum
    out$CI <- 10^c(left,right)}
  return(out)}
