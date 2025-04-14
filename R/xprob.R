#xprob gives the probability of a sampled community transition under some ratio of t and J
#takes two vectors of relative abundance m and n
#takes a value t/J for the ratio of t (# generations) to J (community size)
#ss is a vector of length 2 giving the number of samples at times for m & n
#conditional probability for multiple binomial draws is:
# P(n[i]) = B(N = J - m[1:i-1], p = m[i]/(J-m[1:i-1]))
xprob <- function(m,n,tJ,ss=NA){
  J <- 1/tJ
  if(!any(is.na(ss))){ #if ss given...
    sizes <- c(ss,1/tJ)
    size <- min(sizes) / sum(min(sizes)/sizes)
  } else {
    size <- 1/tJ}
  while(sum(m)==1 | sum(n)==1){ #remove last species if sum is 1, not needed for calculation
    m <- head(m,-1)
    n <- head(n,-1)}
  n <- n*size + 0.5 #moving onto scale of cbinom
  out <- ifelse(m[1]>0 & n[1]>0,dcbinom(x=n[1],size=size,prob=m[1],log=T)+log(size+1),0) #for sp 1 (multiplied to make sense for prob densities)
  for(i in seq_len(length(m))[-1]){ #for every other sp i
    if(m[i]>0 & n[i]>0){ #if sp i not at abundance 0...
      out <- out + dcbinom(x=n[i],size=size-sum(n[1:i-1]-0.5),prob=m[i]/(1-sum(m[1:i-1])),log=T)+log(size+1)}} #conditional binomial
  return(out)}
