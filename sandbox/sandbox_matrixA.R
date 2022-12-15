# test the function matrixAInnerADMM

library(CVN)

m <- 20

lambda1 = .4
lambda2 = .3
global_rho = 1
n1 <- global_rho * lambda1
n2 <- global_rho * lambda2

W <- matrix(1, m, m)

#W <- matrix(runif(m*m), ncol = m)
#W <- W %*% t(W) 
#W <- W / max(W) 
#W <- W 
#diag(W) <- 0

D <- CVN::create_matrix_D(W, lambda1, lambda2, global_rho)

matrix_A_inner_ADMM(m, D)
test(m,D)

CVN::matrix_A_inner_ADMM(m, D)


test <- function(m = 5, D) { 
  # create the variable
  a <- CVXR::Variable(1) 
  
  # objective function
  objective <- CVXR::Minimize(a)
  
  # constraints:
  R <- CVXR::diag(a, m, m) - t(D) %*% D
  
  constraint <- {R %>>% 0} # diag(a,m) %>>% R      # A - D'D must be positive semidefinite 
  
  # define the problem using CVXR:
  problem <- CVXR::Problem(objective, constraints = list(constraint))
  #problem <- CVXR::Problem(objective, constraints = constraints)
  
  # solve
  solution <- CVXR::solve(problem, solver = "SCS")
  
  # return the value of a
  return(solution$getValue(a))
}

test(m, D)
