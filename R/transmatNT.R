#transmat returns a transition matrix for a species in a neutral community with non-overlapping generations
#gives the true expectations for abundance of single species after some time
#useful for comparing with our approximation but not exported
transmat <- function(J){
  mat <- matrix(0,J+1,J+1)
  for(i in seq(0,J)){
    mat[i+1,] <- stats::dbinom(size=J,prob=i/J,x=seq(0,J))}
  return(mat)}
