#Helper function for getting confidence intervals on J fits
#Used by fitJ to return a confidence interval on a fitted J value.
#Similar to fitJ, but optimizes a function that has its maximum some distance
#below the maximum likelihood. That distance is 1.92 for a 95% confidence interval.
#confidence indicates degree of confidence desired. Default is 95% (0.95)
#returns the value of the function at the given logJ.
CIfunc <- function(log10J,occs,ages,ML,confidence=0.95,sampled=TRUE,generationtime=1){
  diff <- stats::qchisq(p=confidence,d=1)/2
  out <- xxprob(log10J,occs,ages,sampled,generationtime)
  return(abs(ML-diff-out))}
