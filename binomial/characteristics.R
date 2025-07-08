
operating_characteristic <- function(test, th) {

  M = test[["info"]]$M
  N = length(M)
  vals = list()
  partial_sums_M = cumsum(M)
  sum_M = partial_sums_M[N]

  accept = test[["data"]]$accept
  cont = test[["data"]]$cont

  if (N == 0)
    return(vals)

  z0 = dbinom(0:sum_M, sum_M, th)
  vals[[N]] = z0 * ifelse(accept[[N]], 1, 0)

if(N>1)
  for(n in (N - 1):1){
    s = seq(0, partial_sums_M[n])
    z0 = dbinom(s, partial_sums_M[n], th)

    v = 0
    for(x in 0:M[n + 1]){

      # shift by 1 since arrays are indexed from 1
      t = vals[[n + 1]][s + 1 + x]

      h = exp(lchoose(partial_sums_M[n], s)-
          lchoose(partial_sums_M[n + 1], s + x) + lchoose(M[n + 1], x) )

      v = v + t * h
    }

    vals[[n]] = z0 * ifelse(accept[[n]], 1, 0) * ifelse(cont[[n]], 0, 1) + v * ifelse(cont[[n]], 1, 0)
  }

  v = 0
  for(x in 0:M[1]){

    # shift by 1 since arrays are indexed from 1
    t = vals[[1]][x + 1]
    v = v + t
  }

  v
}



average_sample_number <- function(test, th) {

  M = test[["info"]]$M
  N = length(M)
  vals = list()
  partial_sums_M = cumsum(M)
  sum_M = partial_sums_M[N]

  cont = test[["data"]]$cont

  if (N == 0)
    return(vals)

  vals[[N]] = rep(0, sum_M + 1)

if(N>1)
  for(n in (N - 1):1){
    s = seq(0, partial_sums_M[n])
    z = M[n + 1] * dbinom(s, partial_sums_M[n], th)

    v = 0
    for(x in 0:M[n + 1]){

      # shift by 1 since arrays are indexed from 1
      t = vals[[n + 1]][s + 1 + x]

      h = exp(lchoose(partial_sums_M[n], s)-
          lchoose(partial_sums_M[n + 1], s + x) + lchoose(M[n + 1], x) )

      v = v + t * h
    }

    vals[[n]] = (z + v) * ifelse(cont[[n]], 1, 0)
  }

  v = 0
  for(x in 0:M[1]){

    # shift by 1 since arrays are indexed from 1
    t = vals[[1]][1 + x]

    v = v + t
  }

  v + M[1]
}

average_sample_cost <- function(test, th, cost_fn=function(x, k) x) {

  M = test[["info"]]$M
  N = length(M)
  vals = list()
  partial_sums_M = cumsum(M)
  sum_M = partial_sums_M[N]

  cont = test[["data"]]$cont

  if (N == 0)
    return(vals)

  vals[[N]] = rep(0, sum_M + 1)

if(N>1)
  for(n in (N - 1):1){
    s = seq(0, partial_sums_M[n])
    z = cost_fn(M[n + 1], n+1) * dbinom(s, partial_sums_M[n], th)

    v = 0
    for(x in 0:M[n + 1]){

      # shift by 1 since arrays are indexed from 1
      t = vals[[n + 1]][s + 1 + x]

      h = exp(lchoose(partial_sums_M[n], s) -
          +lchoose(partial_sums_M[n + 1], s + x)+ lchoose(M[n + 1], x))

      v = v + t * h
    }

    vals[[n]] = (z + v) * ifelse(cont[[n]], 1, 0)
  }

  v = 0
  for(x in 0:M[1]){

    # shift by 1 since arrays are indexed from 1
    t = vals[[1]][1 + x]

    v = v + t
  }

  v + cost_fn(M[1], 1)
}
if(sys.nframe()==0){
  source("./optimal_plan/binomial/plan.R")
  cost<-function(n){
c0=0
c1=1
c1*n+c0}
# cost<-function(n){
# n}

l0=160
l1=160
th0=0.4
th1=0.6
# gsizes=seq(1,37,2)
  M=rep(30,5)
print("group sequential")
test=calculate_plan(M,l0,l1,th0,th1, gamma=0.5)
#   MC(test,th0,100000)
alpha=1-operating_characteristic(test,th0)
 asn0=average_sample_number(test,th0)
 beta=operating_characteristic(test,th1)
 asn1=average_sample_number(test,th1)

  print(paste(alpha,beta,asn0,asn1))
print(l0*alpha+l1*beta+asn0+asn1)
  #
  #
  #   MC(test,th1,100000)

# 5 groups
# test




}