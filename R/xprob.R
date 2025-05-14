#' Probability of a single transition under neutral theory
#'
#' @description
#' Gives the log-likelihood of a single observed transition in community compositions
#' under Hubbell's neutral theory, for a given ratio of J (community size) to
#' t (# generations).
#'
#' @details
#' Used by xxprob. Sums of m and n must each be in the interval (0,1]. Allows
#' for incomplete sampling; if nothing is given for ss, assumes m and n represent
#' true relative abundances.
#'
#' @param m vector of relative abundances at start of transition.
#' @param n vector of relative abundances at end of transition.
#' @param Jt ratio of J (community size) to t (# generations).
#' @param ss single value or vector of length 2 giving number of samples at
#' times for m and n. If a single value given, assumes it applies to both m and n.
#'
#' @importFrom cbinom dcbinom
#'
#' @returns loglik value.
#' @export
#'
#' @examples
#'xprob(m=c(0.1,0.1,0.8),n=c(0.2,0.3,0.5),Jt=10,ss=c(1000,1000)) #-0.01480662
#'
xprob <- function(m,n,Jt,ss=NA){
  #error handling
  tol <- 1E-7
  if(sum(m)<=0|sum(n)<=0|sum(m)>1+tol|sum(n)>1+tol){
    stop("Sums of m and n must each be greater than 0 and no greater than 1")}
  if(length(m)!=length(n)){stop("m and n must have the same length")}
  if(!Jt>0){stop("Jt must be positive")}
  #body
  if(!any(is.na(ss))){ #if ss given...
    if(length(ss)==1){ss <- rep(ss,2)} #duplicate if a single value given
    sizes <- c(ss,Jt)
    size <- min(sizes) / sum(min(sizes)/sizes)
  } else {
    size <- Jt}
  order <- rev(order(m)) #sort by abundance of m
  m <- m[order]
  n <- n[order]
  while(sum(m) > 1-tol | sum(n) > 1-tol){ #remove last species if sum is 1, not needed for calculation
    m <- m[-length(m)]
    n <- n[-length(n)]}
  m <- m[!n==0];n <- n[!n==0] #remove indices where n is zero (zeros in m removed in previous step)
  n <- n*size + 0.5 #moving onto scale of cbinom
  out <- ifelse(length(m)==0,0, #check to make sure there is at least one transition left
                cbinom::dcbinom(x=n[1],size=size,prob=m[1],log=T)+log(size)) #for taxon 1 (multiplied to make sense for prob densities)
  for(i in seq_len(length(m))[-1]){ #for every other taxon i
      out <- out + dcbinom(x=n[i],size=size-sum(n[1:i-1]-0.5),prob=m[i]/(1-sum(m[1:i-1])),log=T)+log(size)}
  return(out)}
