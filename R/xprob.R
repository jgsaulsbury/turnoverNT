#xprob gives the probability of a sampled community transition under some ratio of t and J
#takes two vectors of relative abundance m and n
#takes a value t/J for the ratio of t (# generations) to J (community size)
#ss is a vector of length 2 giving the number of samples at times for m & n
#conditional probability for multiple binomial draws is:


#' Probability of a single transition under neutral theory
#'
#' @description
#' Gives the probability of a single observed transition in community compositions
#' under Hubbell's neutral theory, for a given ratio of t (# generations) to J
#' (community size)
#'
#' @details
#' Used by xxprob. Sums of m and n must each be in the interval (0,1]. Allows
#' for incomplete sampling; if nothing is given for ss, assumes m and n represent
#' true relative abundances.
#'
#' @param m vector of relative abundances at start of transition.
#' @param n vector of relative abundances at end of transition.
#' @param tJ ratio of t (# generations) to J (community size).
#' @param ss vector of length 2 giving number of samples at times for m and n.
#'
#' @importFrom cbinom dcbinom
#'
#' @returns loglik value.
#' @export
#'
#' @examples
#'xprob(m=c(0.1,0.1,0.8),n=c(0.2,0.3,0.5),tJ=0.1,ss=c(1000,1000))
#'
xprob <- function(m,n,tJ,ss=NA){
  if(sum(m)<=0|sum(n)<=0|sum(m)>1|sum(n)>1){
    stop("Sums of m and n must each be greater than 0 and no greater than 1")}
  if(length(m)!=length(n)){stop("m and n must have the same length")}
  J <- 1/tJ
  if(!any(is.na(ss))){ #if ss given...
    sizes <- c(ss,1/tJ)
    size <- min(sizes) / sum(min(sizes)/sizes)
  } else {
    size <- 1/tJ}
  while(sum(m)==1 | sum(n)==1){ #remove last species if sum is 1, not needed for calculation
    m <- m[-length(m)]
    n <- n[-length(n)]}
  n <- n*size + 0.5 #moving onto scale of cbinom
  out <- ifelse(m[1]>0 & n[1]>0,cbinom::dcbinom(x=n[1],size=size,prob=m[1],log=T)+log(size+1),0) #for sp 1 (multiplied to make sense for prob densities)
  for(i in seq_len(length(m))[-1]){ #for every other sp i
    if(m[i]>0 & n[i]>0){ #if sp i not at abundance 0...
      out <- out + dcbinom(x=n[i],size=size-sum(n[1:i-1]-0.5),prob=m[i]/(1-sum(m[1:i-1])),log=T)+log(size+1)}}

  return(out)}
