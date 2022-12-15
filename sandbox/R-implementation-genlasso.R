altZ <- function(y, D, W, lambda1, lambda2, global_rho, diagA = 2, rho = 1, max_iter = 1000, eps = 10^-10, old = TRUE, truncate = 10^-5) { 
  
  m <- length(y)
  mm <- nrow(D)
  
  eta1 = lambda1 / global_rho
  eta2 = lambda2 / global_rho
  
  A <- diag(diagA, m) 
  
  if (!is.positive.semi.definite(A - t(D) %*% D)) { 
    stop("The matrix A - D'D must be positive semidefinite") 
  }
  
  beta_new <- rep(0, m)
  beta_old <- rep(0, m)
  
  alpha_new  <- as.matrix(rep(0, max(m, mm)), ncol = 1)
  alpha_old1 <- as.matrix(rep(0, max(m, mm)), ncol = 1)
  alpha_old2 <- as.matrix(rep(0, max(m, mm)), ncol = 1)
  
  C  <- 1 / (rho * diagA + 1) 
  Cb <- C * rho * diagA
  Cy <- C*y
  CD <- C*D 
  
  iter <- 0 
  
  repeat{ 
    
    #x <- crossprod(CD, (2*alpha_old1 - alpha_old2))
    
    alpha <- (2*alpha_old1 - alpha_old2)
    delta <- rep(0, m)
    
    
    # for (i in 1:m) {
    #   cat(sprintf("i = %d\t m = %d\tmm = %d  -------------\n", i, m, mm))
    # 
    #   delta[i] <- eta1 * alpha[i]
    #   if (i != m) {
    #     for (j in ((i+1):m)) {
    #       cat(sprintf('+ (%d, %d)\t%d\n', i,j, i+j+m-1))
    #       #print(W[i,j])
    #       #print(delta[i])
    #       #print(alpha[i + j + m - 1])
    #       delta[i] <- delta[i] + eta2*W[i,j]*alpha[i+j+m-1]
    #     }
    #   }
    # 
    #   if (i != 1) {
    #     for (j in (1:(i-1))) {
    #       cat(sprintf('- (%d, %d)\t%d\n', j, i, i + j + m - 1))
    #       delta[i] <- delta[i] - eta2*W[j,i]*alpha[i+j+m-1]
    #     }
    #   }
    # }
    
    k <- c(0, seq(m-1,2))
    for (i in 2:(m-1)) { 
      k[i] <- k[i] + k[i-1]
    }
    #print(k)
    
    # k 
    
    #k <- rev(c(seq(m-1:1))-1)
    #k <- cumsum(k)
    #k <- c(0,3,5)
    
    for (i in 1:m) {
      #cat(sprintf("i = %d\t m = %d\tmm = %d  -------------\n", i, m, mm))
      
      #cat(sprintf("k: %d\n",k[i]))
      delta[i] <- eta1 * alpha[i]
      if (i != m) {
        for (j in ((i+1):m)) {
          #   cat(sprintf('+ (%d, %d)\t%d\n', i,j, m + k[i] + (j - i)))
          #print(W[i,j])
          #print(delta[i])
          #print(alpha[i + j + m - 1])
          delta[i] <- delta[i] + eta2*W[i,j]*alpha[m + k[i] + (j - i)]
          #  cat(sprintf("%g\t%g\t%g\n", delta[i], W[i,j], alpha[m + k[i] + (j - i)]))
        }
      }
      
      if (i != 1) {
        for (j in (1:(i-1))) {
          # cat(sprintf('- (%d, %d)\t%d\n', i, j, m + k[j] - (j - i)))
          delta[i] <- delta[i] - eta2*W[j,i]*alpha[m + k[j] - (j - i)]
          #cat(sprintf("%g\t%g\t%g\n", delta[i], W[i,j], alpha[m + k[i] + (j - i)]))
        }
      }
      
    }
    
    #cat(sprintf("orig:\n"))
    #print(x)
    
    #cat(sprintf("new:\n"))
    
    if (!old) { 
      x <- C*delta
    }
    
    beta_new <- Cb*beta_old + Cy - x 
    #beta_new <- beta_old + (y - delta) / (rho*diagA)
    #  c*(rho*a*beta_old + y - x)
    
    #alpha_new <- alpha_old1 + rho * D %*% beta_new
    
    
    
    for (i in 1:m) { 
      alpha_new[i] <- alpha_old1[i] + rho * eta1 * beta_new[i]  
      if (alpha_new[i] > 1) { 
        alpha_new[i] = 1   
      } 
      if (alpha_new[i] < -1) { 
        alpha_new[i] = -1   
      } 
    }
    
    k = m+1
    #print(beta_new)
    #print(beta_old)
    #cat("\n")
    for (i in 1:(m-1)) {
      for(j in (i+1):m) { 
        #cat(sprintf("k: %d\t(i,j) = (%d, %d)\n", k, i, j))
        alpha_new[k] = alpha_old1[k] + rho * eta2*W[i,j]*(beta_new[i] - beta_new[j]) ; 
        #cat(sprintf("old: %g\t new: %g\t W: %g\teta2: %g\n", alpha_old1[k], alpha_new[k], W[i,j], eta2))
        if (alpha_new[k] > 1) { 
          alpha_new[k] = 1   
        } 
        if (alpha_new[k] < -1) { 
          alpha_new[k] = -1   
        } 
        k <- k+1  
      }
    }
    
    #alpha_new <- pmax(pmin(alpha_new, 1), -1)
    
    #print(alpha_new) 
    #print(alpha_old1) 
    #print(alpha_old2) 
    
    if (sum(abs(beta_new - beta_old)) <= eps || iter >= max_iter) { 
      break 
    }
    
    beta_old <- beta_new
    alpha_old2 <- alpha_old1
    alpha_old1 <- alpha_new
    
    #print(alpha_new)
    
    #print(beta_new)
    #if (iter == 2) { 
    #  return(beta_old) ; 
    #}
    
    iter <- iter + 1
  }
  
  
  #print(iter)
  beta_new = sapply(beta_new, function(x) 
    if (abs(x) < truncate) {
      return(0)
    } else {
      return(x)
    })
  
  beta_new
}