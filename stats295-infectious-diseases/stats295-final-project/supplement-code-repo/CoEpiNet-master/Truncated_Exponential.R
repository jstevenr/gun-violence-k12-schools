# 05/08/2019
# Truncated Exponential Distribution

# 1. a function to generate samples from Exp(rate) truncated at (a,b)
sam_texp <- function(n, rate, a=0, b=Inf){
  quans = runif(n)
  sams = -(log(exp(-rate*a)-quans*(exp(-rate*a)-exp(-rate*b))))/rate
  return(sams)
}


# 2. a function to evaluate density of truncated exponential
dtexp <- function(x, rate, a=0, b=Inf, log=F){
  if(any(x < a | x > b)){
    stop("Truncated Exponential: input values not within the bounds!\n")
  }
  dx = rate * exp(-rate * x)/(exp(-rate * a) - exp(-rate * b))
  if(length(x)==1){
    res = dx
  }else{
    res = prod(dx)
  }
  if(log){
    return(log(res))
  }else{
    return(res)
  }
}

# # checking
# x = sam_texp(20,0.5,1,7)
# log(prod(TruncatedDistributions::dtexp(x,0.5,1,7)))
# dtexp(x,0.5,1,7,log = T)
# # seems to be working
