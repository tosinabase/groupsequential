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
    "N" = length(M),
    "alpha" = 1 - operating_characteristic(test, th0),
    "beta" = operating_characteristic(test, th1),
    "AGN0" = average_sampling_cost(test, th0, cost_fn=function (x) 1),
    "AGN1" = average_sampling_cost(test, th1, cost_fn=function (x) 1),
    "ASN0" = average_sampling_cost(test, th0, costs),
    "ASN1" = average_sampling_cost(test, th1, costs)
  )
  return (result_list)
}

print("An example of optimal plan computation")
# input parameters:
# Lagrangian multipliers
  l0 <- 1526 # lambda0
  l1 <- 326  # lambda1
# hypothesized values:
  th0 <- 0    # theta0
  th1 <- 0.21 # theta1
  precision = 0.05   #  grid step
  thgam = c(th0,th1) #  points at which ASC is weighted
  gam = c(0.5, 0.5)  #  weights for weighting ASC
  M = c(100, 100)    #  sample size for the first and the second step
  costs = c(100, 100) # group costs
  cost_fn <- function (x) x # cost of a group group of x observations
  # main function call for optimal plan cumputation:
  test = calculate_plan(M, l0, l1, th0, th1, gam, thgam,
                        costs=NULL,
                        cost_fn=cost_fn,
                        precision=precision, print_output=FALSE)
# report on the computed plan
  results = plan_summary(test, M, l0, l1, th0, th1, costs, cost_fn)

# print the computed plan characteristics:
  for (name in names(results)){
    print(paste(name, results[[name]]))
  }

print("Continuation interval for the first step")
    print(paste(head(test$data[[1]]$grid, 1), tail(test$data[[1]]$grid, 1)))
print("Critical constant for the second step")
  print(test$data[[2]]$const)



