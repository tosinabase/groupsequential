

intgr <- function(spts, vls, t, M){
  # Backward induction integration
  # numerical integral of \int u(x)*dnorm(x-t,sd=sqrt(M)) dx
  # spts is grid of x values, vls corresponding values of u
  PhiVec = pnorm(spts - t, sd=sqrt(M))
  phiVec = dnorm(spts - t, sd=sqrt(M))
  diffPhi = diff(PhiVec)
  diffphi = diff(phiVec)
  tmp = sum(head(vls, -1) * diffPhi)
  vec = diffPhi * head(spts - t, -1) + M * diffphi
  tmp = tmp - sum((diff(vls)/diff(spts)) * vec)
  return(tmp)
}



operating_characteristic <- function(test, th) {


  int_d <- function(s, n, M, d_vals){
    # d values integration
    intgr(
      d_vals[[n + 1]]$grid - s,
      d_vals[[n + 1]]$val * exp((d_vals[[n + 1]]$grid - s) * th) * exp(-th^2 * M[n + 1] / 2),
      0,
      M[n + 1]
    )
  }

  F <- function(s, n, M, d_vals){
    pnorm(d_vals[[n + 1]]$grid[1] - s, mean=M[n + 1] * th, sd=sqrt(M[n + 1]))
  }

  first_F <- function(const, s, n, M){
    pnorm(const - s, mean=M[n + 1] * th, sd=sqrt(M[n + 1]))
  }

  M = test[["info"]]$M
  N = length(M)


  test_vals = test[["data"]]
  d_vals = list()

  acceptance_constant = test[["info"]]$acceptance_constant

  if(N == 1){
    return(
      pnorm(acceptance_constant, mean=M[1] * th, sd=sqrt(M[1]))
    )
  }

  grid = test_vals[[N - 1]]$grid
  val = Vectorize(first_F, vectorize.args="s")(acceptance_constant, grid, N - 1, M)

  stepdata = list(n=N - 1, grid=grid, val=val)
  d_vals[[N - 1]] = stepdata

  if(N > 2)
    for (n in (N - 2):1){
      p1 = Vectorize(F, vectorize.args="s")(test_vals[[n]]$grid, n, M, d_vals)
      p2 = Vectorize(int_d, vectorize.args="s")(test_vals[[n]]$grid, n, M, d_vals)

      grid = test_vals[[n]]$grid
      val = p1 + p2

      stepdata = list(n=n, grid=grid, val=val)
      d_vals[[n]] = stepdata
    }

  F(0 ,0, M, d_vals) + int_d(0, 0, M, d_vals)
}


average_sampling_cost <- function(test, th, costs=NULL, cost_fn=NULL){
  # costs is the vector of group costs
  # If costs=NULL the parameter test[["info"]]$costs will be used.

  if(!is.null(cost_fn)){
    costs = sapply(M, FUN=cost_fn)
  }

  if(is.null(costs)){
    stop("One of the variables costs or cost_fn has to be set")
  }

  int_q <- function(s, n, M, q_vals){
    # q values integration
    intgr(
      q_vals[[n + 1]]$grid - s,
      q_vals[[n + 1]]$val * exp((q_vals[[n + 1]]$grid - s) * th) * exp(-th^2 * M[n + 1] / 2),
      0,
      M[n + 1]
    )
  }

  M = test[["info"]]$M
  N = length(M)
  test_vals = test[["data"]]
  q_vals = list()

  if (N == 1){
    return(costs[1])
  }

  grid = test_vals[[N - 1]]$grid
  val = rep(costs[N], length(test_vals[[N - 1]]$grid))

  stepdata = list(n=N - 1, grid=grid, val=val)
  q_vals[[N - 1]] = stepdata

  if(N > 2)
    for (n in (N - 2):1){
      p = Vectorize(int_q, vectorize.args="s")(test_vals[[n]]$grid, n, M, q_vals)
      
      grid = test_vals[[n]]$grid
      val = costs[n + 1] + p

      stepdata = list(n=n, grid=grid, val=val)
      q_vals[[n]] = stepdata
    }

  costs[1] + int_q(0, 0, M, q_vals)
}


prob_to_stop_after <- function(k, test, th){
  # k < N

  int_t <- function(s, n, M, t_vals){
    # t values integration
    intgr(
      t_vals[[n + 1]]$grid - s,
      t_vals[[n + 1]]$val * exp((t_vals[[n + 1]]$grid - s) * th) * exp(-th^2 * M[n + 1] / 2),
      0,
      M[n + 1]
    )
  }

  M = test[["info"]]$M
  N = length(M)
  test_vals = test[["data"]]
  t_vals = list()

  if (N == 1){
    return(0)
  }


  grid = test_vals[[k]]$grid
  val = rep(1, length(test_vals[[k]]$grid))

  stepdata = list(n=k, grid=grid, val=val)
  t_vals[[k]] = stepdata

  if(k > 1)
    for (n in (k - 1):1){

      grid = test_vals[[n]]$grid
      val = Vectorize(int_t, vectorize.args="s")(test_vals[[n]]$grid, n, M, t_vals)

      stepdata = list(n=n, grid=grid, val=val)
      t_vals[[n]] = stepdata
    }

  int_t(0, 0, M, t_vals)

}


asc_using_prob_to_stop_after <- function(test, th, costs=NULL, cost_fn=NULL){
  # Average sampling cost, realized through prob_to_stop_after,
  # computationally not perfect but can be used for testing prob_to_stop_after or average_sampling_cost

  M = test[["info"]]$M
  N = test[["info"]]$N

  if(!is.null(cost_fn)){
    costs = sapply(M, FUN=cost_fn)
  }

  if(is.null(costs)){
    stop("One of the variables costs or cost_fn has to be set")
  }

  s = costs[1]
  for(k in 2:N){
    s = s + costs[k] * prob_to_stop_after(k - 1, test, th)
  }

  s
}

step_quantile <- function (test, th, q){

  n = 0
  p = 1
  if (p <= 1 - q) return(0)

  M = test[["info"]]$M
  N = length(M)

  repeat{
    n = n + 1
    if (n == N) return(N)

    p = prob_to_stop_after(n, test, th)
    if (p <= 1 - q) return(n)
  }
}


sample_number_quantile <- function (test, th, q){

  n = 0
  p = 1
  if (p <= 1 - q) return(0)

  M = test[["info"]]$M
  sums_M = cumsum(M)
  N = length(M)

  repeat{
    n = n + 1
    if (n == N) return(sums_M[N])

    p = prob_to_stop_after(n, test, th)
    if (p <= 1 - q) return(sums_M[n])
  }
}

