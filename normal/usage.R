# The usage example

source("./normal/plan.R")
source("./normal/characteristics.R")



plan_summary <- function(test, M, l0, l1, th0, th1, costs, cost_fn) {

  result_list = list(
    "lambda0" = l0,
    "lambda1" = l1,
    "theta0" = th0,
    "theta1" = th1,
    "M" = M,
    "costs" = test[["info"]]$costs,
    "N" = length(M),
    "alpha" = 1 - operating_characteristic(test, th0),
    "beta" = operating_characteristic(test, th1),
    "AGN0" = average_sampling_cost(test, th0, cost_fn=function (x) 1),
    "AGN1" = average_sampling_cost(test, th1, cost_fn=function (x) 1),
    "ASN0" = average_sampling_cost(test, th0, costs=M),
    "ASN1" = average_sampling_cost(test, th1, costs=M),
    "ASС0" = average_sampling_cost(test, th0, costs, cost_fn),
    "ASС1" = average_sampling_cost(test, th1, costs, cost_fn)
  )
  return (result_list)
}

if (sys.nframe() == 0){

  # input parameters can be set here

  l0 <- 530.3
  l1 <- 530.3
  th0 <- -0.1
  th1 <- 0.1
  precision = 0.01

  # thgam = c(th0, th1)
  # gam = c(0.5, 0.5)

  thgam = 0
  gam = 1

  M <- c(30, 20, 30, 25, 50)
  # M <- rep(1, 100)

  # costs = M
  costs = c(5, 5, 17, 32.55, 67.88)


  cost_fn = NULL
  # cost_fn <- function (x) 1 + 2*x


  test = calculate_plan(M, l0, l1, th0, th1, gam, thgam,
                        costs=costs,
                        cost_fn=cost_fn,
                        precision=precision, print_output=TRUE)

  results = plan_summary(test, M, l0, l1, th0, th1, costs, cost_fn)
  print(results)

  for (name in names(results)){
    print(paste(name, results[[name]]))

  }

  asc_0_with_p = asc_using_prob_to_stop_after(test, th0, costs=costs, cost_fn=cost_fn)
  print("asc_0 using prob_to_stop_after")
  print(asc_0_with_p)

  res = monte_carlo_simulation(M, test, th0, costs=costs, cost_fn=cost_fn, nMC=100000)
  print(res)



}
