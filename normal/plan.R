
require("pracma")
# require("Bolstad2")
source("./normal/characteristics.R")


d <- function (s, thx, m) {
  # probability density function without common multiple
  if (length(s) != 1) stop("s should be a number here")

  exp(thx * s - (m * thx ^ 2) / 2 )
}


L <- function (s, l0, l1, th0, th1, M, N, sums_M){
  # first (at step N - 1) integral value (found analiticaly)
  if (length(s) != 1) stop("s should be a number here")

  B = tail(sums_M, 1) * (th0 + th1) / 2 + log(l0 / l1) / (th1 - th0)
  rp1 = l0 * d(s, th0, sums_M[N - 1]) * (1 - erf(((B - s)/ M[N] - th0) * sqrt(M[N] / 2)))
  rp2 = l1 * d(s, th1, sums_M[N - 1]) * (1 + erf(((B - s)/ M[N] - th1) * sqrt(M[N] / 2)))

  (rp1 + rp2) / 2
}


find_first_continuation_interval <- function (l0, l1, th0, th1, gam, thgam, M, N, sums_M, costs){

  current_loss <- function(s){ min(l0 * d(s, th0, sums_M[N - 1]), l1 * d(s, th1, sums_M[N - 1])) }
  future_loss <- function(s){
    tmp = sapply(thgam, function(x) d(s, x, sums_M[N - 1]))%*%gam
    costs[N]*tmp + L(s, l0, l1, th0, th1, M, N, sums_M)
  }

  to_solve <- function (s){ current_loss(s) - future_loss(s) }

  c = calc_center_bound(N - 1, sums_M[N - 1], l0, l1, th0, th1)  # center_bound
  interval_const = 10

  left_solution = uniroot(to_solve, c(c - interval_const, c), extendInt="upX")
  right_solution = uniroot(to_solve, c(c, c + interval_const), extendInt="downX")

  return(c(left_solution$root, right_solution$root))
}


find_continuation_interval <- function (stepdata, l0, l1, th0, th1, gam, thgam, M, sums_M, costs){
  # all calculatuions exactly for n (previous step minus one)
  n = stepdata$n - 1

  current_loss <- function(s){ min(l0 * d(s, th0, sums_M[n]), l1 * d(s, th1, sums_M[n])) }

  future_loss <- function(s){
    tmp = sapply(thgam, function(x) d(s, x, sums_M[n]))%*%gam
    costs[n + 1] * tmp + calculate_integral(s, stepdata, l0, l1, th0, th1, M, sums_M)
  }

  to_solve <- function (s){ current_loss(s) - future_loss(s) }

  c = calc_center_bound(n, sums_M[n], l0, l1, th0, th1)  # center_bound

  a = calc_left_bound(n, sums_M[n], l1, th0, th1)  # left_bound

  interval_const = 10
  left_solution = uniroot(to_solve, c(a, c))
  right_solution = uniroot(to_solve, c(c, c + interval_const), extendInt="downX")


  return(c(left_solution$root, right_solution$root))
}

calculate_integral <- function (s, stepdata, l0, l1, th0, th1, M, sums_M){
  # s should be a number

  if (length(s) != 1) stop("s should be a number here")

  n = stepdata$n

  # left bound of interval
  a = head(stepdata$grid, 1)
  # right bound of interval
  b = tail(stepdata$grid, 1)


  # numerical part of integral
  int3 = intgr(stepdata$grid, stepdata$val, s, M[n])
  # analytical parts of integral
  int1 = l0 * d(s, th0, sums_M[n - 1]) * (1 - erf(((b - s)/ M[n] - th0) * sqrt(M[n]/2 ))) / 2
  int2 = l1 * d(s, th1, sums_M[n - 1]) * (1 + erf(((a - s)/ M[n] - th1) * sqrt(M[n]/2 ))) / 2

  parts = c(int1, int2, int3)
  sum(parts)
}

calculate_plan <- function(M, l0, l1, th0, th1, gam, thgam,
                           costs=NULL, cost_fn=NULL,
                           precision=0.005, print_output=TRUE){
  # M is vector of groups sizes
  # It is possible set costs as a vector (costs) or as a function (cost_fn)
      # costs is the vector of group costs or cost function,
      # cost_fn is functon of costs


  if(!is.null(cost_fn)){
    costs = sapply(M, FUN=cost_fn)
  }
  if(is.null(costs)){
      costs = M
  }

  test = list()

  N = length(M)
  sums_M = cumsum(M)

  info_list = list(l0=l0, l1=l1, th0=th0, th1=th1, M=M, costs=costs, N=N, gam=gam, thgam=thgam,
            precision=precision)

  if (N == 0)
    return(list(info=info_list, data=test))

  n = N - 1
  acceptance_constant = log(l0 / l1) / (th1 - th0) + tail(sums_M, n=1) * (th0 + th1) / 2
  test[[n + 1]] = list(n=N, const=acceptance_constant, grid=acceptance_constant, val=NA)

  info_list[["acceptance_constant"]] = acceptance_constant

  if (print_output) {
    print(N)
    print(acceptance_constant)
  }

  if (N == 1){
    return(list(info=info_list, data=test))
  }

  continuation_interval = find_first_continuation_interval(l0, l1, th0, th1, gam, thgam, M, N, sums_M, costs)
  a = continuation_interval[1]
  b = continuation_interval[2]
  nint = ceiling((b - a) / precision)
  grid = seq(a, b, length=nint + 1)
  tmp = sapply(thgam, function(x) Vectorize(d, vectorize.args="s")(grid, x, sums_M[N - 1]))%*%gam
  val = costs[N] * ((tmp)) + Vectorize(L, vectorize.args="s")(grid, l0, l1, th0, th1, M, N, sums_M)
  stepdata = list(n=N - 1, grid=grid, val=val)
  test[[N - 1]] = stepdata

  if (print_output) {
    print(N - 1)
    print(continuation_interval)
  }

  if (N == 2){
    return(list(info=info_list, data=test))
  }

  repeat{
    n = stepdata$n - 1
    continuation_interval = find_continuation_interval(stepdata, l0, l1, th0, th1, gam, thgam, M, sums_M, costs)

    if (print_output) {
      print(n)
      print(continuation_interval)
    }

    a = continuation_interval[1]
    b = continuation_interval[2]
    nint = ceiling((b - a) / precision)
    grid = seq(a, b, length=nint + 1)
    tmp = sapply(thgam, function(x) Vectorize(d, vectorize.args="s")(grid, x, sums_M[n]))%*%gam
    val = costs[n+1] * tmp + Vectorize(calculate_integral, vectorize.args="s")(grid, stepdata, l0, l1, th0, th1, M, sums_M)

    stepdata = list(n=n, grid=grid, val=val)
    test[[n]] = stepdata

    if(n == 1)
      break
  }
  return(list(info=info_list, data=test))
}


calc_center_bound <- function (n, m, l0, l1, th0, th1){
  bound(l0/l1, th0, th1, n, m)
}

calc_left_bound <- function(n, m, l1, th0, th1){
  bound(1/l1, th0, th1, n, m)
}

bound <- function (l, tx0, tx1, n, m){
  log(l) / (tx1 - tx0) +  (b_func(tx1, m) - b_func(tx0, m)) / (tx1 - tx0)
}

b_func <- function (thx, m){
  # for normal distribution
  m * thx^2 / 2
}


monte_carlo_simulation <- function(M, test, hyp, costs=NULL, cost_fn=NULL, nMC=10000){
  # hyp = true parameter value
  # returns relative frequency of acceptations
  # the average run length and its standard deviationtest

  if(!is.null(cost_fn)){
    costs = sapply(M, FUN=cost_fn)
  }

  if(is.null(costs)){
    stop("One of the variables costs or cost_fn has to be set")
  }

  H = length(M)
  const = test$info$acceptance_constant
  ss = 0
  totaccepted = 0
  totn = 0
  tot_stages = 0
  tot_cost = 0

  for(i in 1:nMC){
    # accumulated sum for test run
    s = 0

    # number of steps in current run
    n = 0

    # accumulated cost in current run
    cost_accum = 0

    #number of accepting in the run (0 or 1)
    accepted = 0

    # run starts
    for (stage in 1:H){
      #generate
      # summand = sum(rnorm(M[stage], mean=hyp))
      summand = rnorm(1, mean=M[stage]*hyp, sd=sqrt(M[stage]))

      # accumulated sum
      s = s + summand

      # one step more
      n = n + M[stage]

      # cost accumulation
      cost_accum = cost_accum + costs[stage]

      # the last stage of the run
      if(stage == H){
        # if s is below of the decision constant then accept
        if(s <= const)
          accepted = accepted + 1
      }
      else {
        # stop by the optimal stopping rule
        if(s < head(test$data[[stage]]$grid, 1) | s > tail(test$data[[stage]]$grid, 1)){
          if (s <= head(test$data[[stage]]$grid, 1))
            accepted = accepted + 1
          # if s is in the acceptation
          # region then accept
          # accepted or rejected, stop the run; no more stages
          break
        }
      }
    }

    totaccepted = totaccepted + accepted
    totn = totn + n
    tot_stages = tot_stages + stage
    tot_cost = tot_cost + cost_accum
    ss = ss + n^2
  }
  OC = (totaccepted) / as.double(nMC)
  ESS = totn / as.double(nMC)
  EGN = tot_stages / as.double(nMC)
  ASC = tot_cost / as.double(nMC)
  # stdev = sqrt((ss - totn^2 / as.double(nMC)) / as.double(nMC - 1))

  return(c(OC, ESS, EGN, ASC))
}
