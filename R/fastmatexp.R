#fastmatexp takes a matrix mat to a power n
#uses the exponentiation by squaring algorithm
#use with transmatNT() to get true transition probs for a single species
fastmatexp <- function(mat,n){
  if(n==0){
    return(diag(dim(mat)[1])) #identity matrix
  } else if(n==1) {
    return(mat)
  } else if(n %% 2 == 0){ #else if n is even
    return(fastexp(mat %*% mat,n/2))
  }else{ #else if n is odd
    return(mat %*% fastexp(mat %*% mat,(n-1)/2))
  }}
