
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
