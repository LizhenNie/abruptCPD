our_method_single_run = function(d, n0 = NULL, n1 = NULL, correct = FALSE){
  
  n = nrow(d)
  if(is.null(n0)){n0 = ceiling(0.05 * n)}
  if(is.null(n1)){n1 = floor(0.95 * n)}
  
  # the scaling parameter quantities used in T_scale
  colsum = apply(d, 1, sum)
  totalsum = sum(colsum)
  s = sum((colsum * 2 / n - totalsum / n^2)^2) / 4 / n - 1 / 4 * totalsum^2 / n^4
  s = sqrt(s)
  
  # the (i,j)-th element in matrix dd is sum of all elements in the top left submatrix 
  # to the  (i,j)-th element in distance matrix d
  dd = matrix(0, nrow = n+1, ncol = n+1) 
  for (row in 2:(n+1)){
    for (col in 2:(n+1)){
      dd[row,col] = dd[row, col-1] + dd[row-1, col] - dd[row-1, col-1] + d[row-1,col-1]
    }
  }
  dd = dd[2:(n+1), 2:(n+1)]
  
  # functions to calculate average within / between group similarities (divided by index t)
  b1 = function(t){return(dd[t,t] / t^2)}
  b2 = function(t){ return( ( dd[n,n]-dd[t,n]-dd[n,t]+dd[t,t] ) / (n-t)^2 ) }
  a = function(t){ return( ( dd[n,t]-dd[t,t] ) / t / (n-t) ) }
  
  # calculating T_loc (T1), T_scale (T2) and T_joint (T3)
  bb1 = rep(0, n); bb2 = rep(0, n); aa = rep(0, n)
  S1_f = rep(0, n); S2_f = rep(0, n); S3_f = rep(0, n); correct_S2_f = rep(0, n)
  e2 = dd[n, n] / 2 / n / (n-1)
  
  for (t in n0 : n1){
    
    bb1[t] = b1(t)
    bb2[t] = b2(t)
    aa[t] = a(t)
    
    if(correct){ # if using higher-order corrections
      S1_f[t] = t * (n-t) / n * ( aa[t] - 0.5 * bb1[t] * t / (t-1) - 0.5 * bb2[t] * (n-t) / (n-t-1) )
    }else{
      S1_f[t] = t * (n-t) / n * ( aa[t] - 0.5 * bb1[t] - 0.5 * bb2[t] ) 
    }
    
    S2_f[t] = 1 / s / 2 * sqrt(t * (n-t) / n) * (bb1[t] - bb2[t])
    S3_f[t] = 1 / s^2 * n / t / (n-t) * S1_f[t]^2  + (S2_f[t] - 1 / sqrt(n) / s / sqrt( t / n * (1 - t / n) ) * e2 * (t / n * (n-t) / (n-t-1) - t / (t-1) * (n-t) / n) )^2
    
    if (correct){ # if using higher-order corrections
      correct_S2_f[t] = abs(S2_f[t] - 1 / sqrt(n) / s / sqrt( t/n * (1-t/n) ) * e2 * (t/n*(n-t)/(n-t-1) - t/(t-1)*(n-t)/n) )
    }else{
      S2_f[t] = abs(S2_f[t])
    }
    
  }
  
  
  # calculating S_loc (S1), S_sale (S2), S_joint (S3) 
  Z1 = max(S1_f[n0:n1])
  predicted1 = which.max(S1_f[n0:n1]) + n0-1
  
  if(correct){ # if using higher-order corrections
    Z2 = max(correct_S2_f[n0:n1])
    predicted2 = which.max(correct_S2_f[n0:n1]) + n0-1
  }else{
    Z2 = max(S2_f[n0:n1])
    predicted2 = which.max(S2_f[n0:n1]) + n0-1
  }
  
  Z3 = max(S3_f[n0:n1])
  predicted3 = which.max(S3_f[n0:n1]) + n0-1
  
  return (list(
    Z1 = Z1, predicted1 = predicted1, #D-loc
    Z2 = Z2, predicted2 = predicted2, #D-scale
    Z3 = Z3, predicted3 = predicted3, #D-joint
    colsum = colsum, totalsum = totalsum
  ))
  
}


our_method = function(d, n0 = NULL, n1 = NULL, correct = FALSE, permutations = 500){
  
  n = nrow(d)
  if(is.null(n0)){n0 = ceiling(0.05 * n)}
  if(is.null(n1)){n1 = floor(0.95 * n)}
  
  # --------- calculate statistic and location --------
  res = our_method_single_run(d, n0 = n0, n1 = n1, correct = correct)
  Z1 = res$Z1; predicted1 = res$predicted1
  Z2 = res$Z2; predicted2 = res$predicted2
  Z3 = res$Z3; predicted3 = res$predicted3
  colsum = res$colsum; totalsum = res$totalsum
  
  # --------- calculate p-value ------------
  permutation_pval = function(Z1, Z2, Z3, permutations){
    
    nulldist1 = rep(NA, permutations)
    nulldist2 = rep(NA, permutations)
    nulldist3 = rep(NA, permutations)
    nulldistDubey = rep(NA, permutations)
    
    for (i in 1:permutations){
      
      permu = sample(1:n, n)
      newd = d[permu, permu]
      
      newres = our_method_single_run(newd, n0 = n0, n1 = n1, correct = correct)
      nulldist1[i] = newres$Z1
      nulldist2[i] = newres$Z2
      nulldist3[i] = newres$Z3
      
    }
    
    return(
      list( 
        Z1_pval = mean(nulldist1 >= Z1), 
        Z2_pval = mean(nulldist2 >= Z2), 
        Z3_pval = mean(nulldist3 >= Z3))
    )
    
  }
  
  tmp = permutation_pval(Z1 = Z1, Z2 = Z2, Z3 = Z3, permutations = permutations)
  Z1_pval = tmp$Z1_pval
  Z2_pval = tmp$Z2_pval
  Z3_pval = tmp$Z3_pval
  
  return(
    list(
      Z1 = Z1, location1 = predicted1, Z1_pval = Z1_pval, 
      Z2 = Z2, location2 = predicted2, Z2_pval = Z2_pval, 
      Z3 = Z3, location3 = predicted3, Z3_pval = Z3_pval
    )
  )
  
}
