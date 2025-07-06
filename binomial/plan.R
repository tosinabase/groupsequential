

calculate_plan <- function(M, l0, l1, th0, th1, gam, thgam, cost_fn=function(x,k) x) {

  N = length(M)
  w = list()
  accept = list()
  cont = list()
  partial_sums_M = cumsum(M)
  sum_M = partial_sums_M[N]

  if (N == 0)
    return(w)

  z0 = l0 * dbinom(0:sum_M, sum_M, th0)
  z1 = l1 * dbinom(0:sum_M, sum_M, th1)
  w[[N]] = pmin(z0, z1)

  accept[[N]] = z0 >= z1
  cont[[N]] = rep(FALSE, sum_M + 1)

  if(N > 1)
    for(n in (N - 1):1){
      s = seq(0, partial_sums_M[n])
      z0 = l0 * dbinom(s, partial_sums_M[n], th0)
      z1 = l1 * dbinom(s, partial_sums_M[n], th1)

      z = cost_fn(M[n + 1], n+1) * sapply(thgam, function(x) dbinom(s, partial_sums_M[n], x))%*%gam
      z = c(z)

      v = 0
      for(x in 0:M[n + 1]){

        # shift by 1 since arrays are indexed from 1
        t = w[[n + 1]][s + 1 + x]

        h = exp(lchoose(partial_sums_M[n], s)-lchoose(partial_sums_M[n + 1], s + x) + lchoose(M[n + 1], x))


        v = v + t * h
      }

      w[[n]] = pmin(z0, z1, z + v)

      accept[[n]] = z0 >= z1
      cont[[n]] = z + v < pmin(z0, z1)
    }

  list(
    info=list(M=M, l0=l0, l1=l1, th0=th0, th1=th1, gam=gam, thgam=thgam),
    data=list(w=w, cont=cont, accept=accept)
  )

}

calculate_simplified <- function(M, l0, l1, th0, th1, gamma) {
  N = length(M)

  accept = list()
  cont = list()
  partial_sums_M = cumsum(M)
  sum_M = partial_sums_M[N]


  z0 = l0 * dbinom(0:sum_M, sum_M, th0)
  z1 = l1 * dbinom(0:sum_M, sum_M, th1)

  accept[[N]] = z0 >= z1
  cont[[N]] = rep(FALSE, sum_M + 1)
   for (n in 1:(N-1)) {
    s = seq(0, partial_sums_M[n])
    z0 = l0 * dbinom(s, partial_sums_M[n], th0)
    z1 = l1 * dbinom(s, partial_sums_M[n], th1)
    z =M[n+1]*  (gamma * dbinom(s, partial_sums_M[n], th0) + (1 - gamma) * dbinom(s, partial_sums_M[n], th1))
     cont[[n]]=z<pmin(z0,z1)
    accept[[n]] = z0 >= z1
  }

  list(
    info=list(M=M, l0=l0, l1=l1, th0=th0, th1=th1, gam=gam, thgam=thgam),
    data=list(cont=cont, accept=accept)
  )

}




if (sys.nframe() == 0){
  source("./optimal_plan/binomial/characteristics.R")

  M <- rep(30,5)
  l0 <- 160
  l1 <- 160
  th0 <- 0.4
  th1 <- 0.6

  test = calculate_plan(M, l0, l1, th0, th1, gam=c(0.5, 0.5), thgam=c(0.5, 0.5))
  1 - operating_characteristic(test,th0)
  operating_characteristic(test,th1)
  test
}
