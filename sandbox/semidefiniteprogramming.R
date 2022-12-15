library(CVXR)
library(CVN)
library(matrixcalc)

m <- 9

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

D <- create_matrix_D(W, lambda1, lambda2, global_rho)
c <- nrow(D)

R <- -t(D) %*% D

# Final? 

A <- CVXR::Variable(m, m)
a <- CVXR::Variable(1) 

objective <- CVXR::Minimize(abs(sum(A)) / 5)

constraint1 <- A == diag(a,m)
constraint2 <- A %>>% R

problem <- CVXR::Problem(objective, constraints = list(constraint1, constraint2))

solution <- CVXR::solve(problem, solver = "SCS")

solution$getValue(a)
Q = diag(solution$getValue(a), m) - R
is.positive.semi.definite(R)
is.positive.semi.definite(Q)






# As a vector
a <- CVXR::Variable(m, 1) 
objective <- CVXR::Minimize(sum(a) / m)
constraint1 <- -t(D) %*% D + diag(a) >= 0
constraint2 <- a >= 0

(Q = t(D) %*% D + diag(solution$getValue(a)))
is.positive.semi.definite(Q)

# As a matrix ---------------------------

A <- CVXR::Variable(m, m)#, PSD = TRUE)

a <- CVXR::Variable(1)
objective <- CVXR::Minimize(abs(sum(A)) / 5)

A %>>% R

# set constraint that it should be a diagonal matrix, a*I
R = t(D) %*% D
constraint1 <- A == diag(a,m)
constraint2 <- A %>>% R
#constraint3 <- a >= 0 

PSDConstraint(expr, id = NA_integer_)


CVXR::PSDConstraint(A - R %>>% 0, id = NA_integer_)

problem <- CVXR::Problem(objective, constraints = list(constraint1, constraint2))




is_dcp(problem)

solution <- CVXR::solve(problem, solver = "SCS")

a_tilde = solution$getValue(A)[1,1]
solution$getValue(A)

(Q = -t(D) %*% D + solution$getValue(A))

eigen(Q)

is.positive.semi.definite(Q)
is.positive.semi.definite(solution$getValue(A) + t(D) %*% D)





a <- CVXR::Variable(1)

objective <- CVXR::Minimize(a)

constraint1 <- a >= 0
constraint2 <- diag(a, m) - t(D) %*% D >= 0

problem <- CVXR::Problem(objective, constraints = list(constraint1, constraint2))

solution <- CVXR::psolve(problem)

s <- solution
s$getValue(a)

Q = -t(D) %*% D + diag(s$getValue(a), m) 

eigen(Q)

solution$value
is.positive.semi.definite(Q)
