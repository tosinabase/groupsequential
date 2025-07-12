# The usage example

source("./binomial/plan.R")
source("./binomial/characteristics.R")



plan_summary <- function(test, M, l0, l1, th0, th1, gam, thgam, cost_fn) {

  result_list = list(
    "lambda0" = l0,
    "lambda1" = l1,
    "theta0" = th0,
    "theta1" = th1,
    "gam" = gam,
    "thgam" = thgam,
    "M" = M,
    "N" = length(M),
    "alpha" = 1 - operating_characteristic(test, th0),
    "beta" = operating_characteristic(test, th1),
    "ASN0" = average_sample_number(test, th0),
    "ASN1" = average_sample_number(test, th1),
    "ASC0" = expected_sampling_cost(test, th0, cost_fn=cost_fn),
    "ASC1" = expected_sampling_cost(test, th1, cost_fn=cost_fn),
    "Aver."=  (expected_sampling_cost(test, th0, cost_fn=cost_fn)+
               expected_sampling_cost(test, th1, cost_fn=cost_fn))/2
  )
  return (result_list)
}

if (sys.nframe() == 0){

  # input parameters can be set here

  cost_fn = function (M, k){M}

  M <- c(100)
  l0 <- 30
  l1 <- 45
  th0 <- 0.935
  th1 <- 0.64

  gam = c(0.5, 0.5)
  thgam = c(th0, th1)




  test = calculate_plan(M, l0, l1, th0, th1, gam, thgam)

  results = plan_summary(test, M, l0, l1, th0, th1, gam, thgam, cost_fn)
  print(results)

  for (name in names(results)){
    print(paste(name, results[[name]]))

  }

}
