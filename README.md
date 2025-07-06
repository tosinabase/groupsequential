# Group-sequential hypothesis testing: an R project for optimal test design


This repository contains accompanying R program code the article "Farkhshatov, F., Novikov, A. Optimal group sequential hypothesis testing for exponential families of distributions: a numerical approach".

## Contains

* [binomial](binomial) for Bernoulli model
  * [binomial/characteristics.R](binomial/characteristics.R)
  * [binomial/plan.R](binomial/plan.R)
  * [binomial/usage.R](binomial/usage.R)

* [normal](normal) for normal model
  * [normal/characteristics.R](normal/characteristics.R)
  * [normal/plan.R](normal/plan.R)
  * [normal/usage.R](normal/usage.R)
  * [normal/example.R](normal/example.R)

## Run
Usage examples can be found in the following files
* [binomial/usage.R](binomial/usage.R)
* [normal/usage.R](normal/usage.R)

### Example - Two-step test for normal model ([normal/example.R](normal/example.R)):
```r
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

cost_fn <- function (x) x # cost of a group group of x observations

# main function call for optimal plan cumputation:
test = calculate_plan(M, l0, l1, th0, th1, gam, thgam,
                      cost_fn=cost_fn,
                      precision=precision, 
                      print_output=FALSE)

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
```