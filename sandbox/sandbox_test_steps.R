
for (m in 3:5) { 
  determine_steps(m)

  k <- c(0, seq(m-1,2))
  for (i in 2:(m-1)) {
    k[i] <- k[i] + k[i-1]
  }
  print(k)
}

m = 5

k <- c(0, seq(m-1,2))
for (i in 2:(m-1)) {
  k[i] <- k[i] + k[i-1]
}


for (i in 1:m) {
  if (i != m) {
    for (j in ((i+1):m)) {
      cat(sprintf("g_%d\t+ (%d,%d): a_%d \n", i, i, j, m + k[i] + (j - i)))
    }
  }
  
  if (i != 1) {
    for (j in (1:(i-1))) {
      cat(sprintf("g_%d\t- (%d,%d): a_%d \n", i, j, i, m + k[j] - (j - i)))
    }
  }
}

CVN::show_indices(m)

k = m+1
for (i in 1:(m-1)) {
  for(j in (i+1):m) { 
    cat(sprintf("%d: (%d, %d)\n", k, i, j))
    k <- k + 1
  }
 
}
