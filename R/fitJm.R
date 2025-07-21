#' Find the maximum-likelihood values of J and m for a community change dataset
#'
#' @description
#' Uses optim() to fit the J and m parameters in Hubbell's neutral theory to
#' a dataset, where J is community size and m is the rate at which vacancies in
#' the community are filled by migrants from a static metacommunity.
#'
#' @details
#' optim() takes the midpoint of the J and m search bounds as the starting value.
#'
#' @param occs matrix of the number of observations in each species at each time.
#' One column for each species, one row for each time slice. Time goes from oldest
#' at the bottom to youngest at the top.
#' @param ages vector containing the ages of each time slice, in years, from
#' oldest to youngest.
#' @param sampled boolean indicating whether occs represents a sampled
#' community (TRUE) or instead represents true species abundances (FALSE).
#' @param generationtime time between generations, in years.
#' @param lowerbounds vector of length 2 giving the lower search bounds for log10J and log10m.
#' @param upperbounds vector of length 2 giving the upper search bounds for log10J and log10m.
#'
#' @returns a list containing "loglik", "J", and "m"
#' @export
#'
#' @examples
#' set.seed(10)
#' sim <- simNT(c(1000,1000,1000,1000),ts=seq(0,2000,50),m=0.001,ss=1000) #simulate under neutral theory with migration
#' fitJm(occs=sim$simulation,ages=sim$times)
fitJm <- function(occs,ages,sampled=TRUE,generationtime=1,lowerbounds=c(1,-30),upperbounds=c(10,-1)){
  if(dim(occs)[1] != length(ages)){stop("'ages' must have length equal to the number of rows of 'occs'")}
  op <- stats::optim(par=c(5,-3),fn=xxprobm,method="L-BFGS-B",lower=lowerbounds,
                     upper=upperbounds,control=list(fnscale=-1),occs=occs,ages=ages,
                     sampled=sampled,generationtime=generationtime)
  out <- list("loglik"=op$value,"J"=10^op$par[1],"m"=10^op$par[2])
  return(out)}
