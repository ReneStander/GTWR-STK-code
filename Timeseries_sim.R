# Linear trend

lin_trend <- function(data, times, var_at, beta1){
  ts_lin1 <- matrix(NA, 
                    nrow = length(data), ncol = times)
  
  for(i in 1:times){
    at <- sqrt(var_at)*rnorm(length(data))
    zt <- data + beta1*i + at
    
    ts_lin1[,i] <- zt
  }
  
  return(ts_lin1)
}


# MA(1)

ma1 <- function(data, times, var_at, theta){
  ts_ma1 <- matrix(NA, 
                   nrow = length(data), ncol = times)
  
  at_1 <- rnorm(length(data))
  
  for(i in -10:times){
    at <- sqrt(var_at)*rnorm(length(data))
    zt <- data + at - theta*at_1
    
    if(i > 0){
      ts_ma1[,i] <- zt
    }
    at_1 <- at
  }
  
  return(ts_ma1)
}

# AR(1)

ar1 <- function(data, times, var_at, theta0, phi){
  ts_ar1 <- matrix(NA, 
                   nrow = length(data), ncol = times)
  
  Zt_1 <- data
  
  for(i in 1:times){
    at <- sqrt(var_at)*rnorm(length(data))
    zt <- theta0 + at + phi*Zt_1
    
    ts_ar1[,i] <- zt
    Zt_1 <- zt
  }
  
  return(ts_ar1)
}