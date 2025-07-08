# Two-step usage example

source("./normal/plan.R")
source("./normal/characteristics.R")
print("An example of optimal plan computation")

# input parameters:
# hypothesized values:
th0 <- 0             # theta0
th1 <- 0.21          # theta1
thgam <- c(th0, th1) # points at which ASC is weighted
gam <- c(0.5, 0.5)   # weights for weighting ASC

# Lagrangian multipliers
l0 <- 1526  # lambda0
l1 <- 326   # lambda1

M <- c(100, 100)           # group sizes
precision <- 0.05          # grid step
cost_fn <- function (x) x  # cost of a group group of x observations

# main function call for optimal plan cumputation:
test <- calculate_plan(M, l0, l1, th0, th1, gam, thgam,
                     cost_fn=cost_fn,
                     precision=precision,
                     print_output=FALSE)

alpha <- 1 - operating_characteristic(test, th0)
beta <- operating_characteristic(test, th1)
print(paste("alpha:", alpha, ", beta:", beta))
# [1] "alpha: 0.0250017679625842 , beta: 0.200056694872385"

ASN0 <- average_sampling_cost(test, th0, cost_fn=cost_fn)
ASN1 <- average_sampling_cost(test, th1, cost_fn=cost_fn)
print(paste("ASN0:", ASN0, ", ASN1:", ASN1))
# [1] "ASN0: 118.489308233436 , ASN1: 148.471027781829"


print(paste("Continuation interval for the first step:",
            head(test$data[[1]]$grid, 1),
            tail(test$data[[1]]$grid, 1)))
# [1] "Continuation interval for the first step: 8.60693365315347 23.3352164728895"

print(paste("Critical constant for the second step:",
            test$data[[2]]$const))
# [1] "Critical constant for the second step: 28.350037287989"





