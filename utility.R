library(esreg)
library(conquer)
library(FarmTest)
library(pracma)
library(adaHuber)
library(BB)
library(numDeriv)

hub_reg <- function(X, Y, c, intercept=TRUE)
############################################
############# Huber regression #############
############################################
## Input 
# X : n by p design matrix
# Y : n by 1 response vector
# c : positive tuning parameter in the Huber loss  
# intercept : logical flag for adding an intercept to the model; default is TRUE.

## Output
# coefficients : estimated regression coefficients 
# residuals : fitted residuals
{
  grad.huber <- function(u,c) {
    w = as.numeric((u*(u<=c & u>=0) + c*(u>c)) + (u*(u>=-c & u<0) - c*(u < -c)))
    return(w/length(u))
  }
  
  tol = 10^(-5)
  n = nrow(X)
  p = ncol(X)
  mX = colMeans(X)
  sdX = apply(X, 2, sd)
  
  if (intercept){
    X = cbind(matrix(1, n, 1), (X - matrix(mX, nr=n, nc=p, byrow=TRUE)) / sdX)
  } else {
    X = X / sdX
    }
  ols = lm(Y~X-1)
  beta = ols$coef
  res   = ols$residuals
  grad0 <- -t(X) %*% grad.huber(res, c)
  diff_beta = -grad0
  beta = beta + diff_beta
  res = Y - X %*% beta
  
  iter = 0
  repeat{
    if( max(abs(diff_beta)) < tol | iter > 500) break
    grad1 <- -t(X) %*% grad.huber(res, c)
    diff_grad = grad1 - grad0
    
    r0 = sum(diff_beta*diff_beta)
    r1 = sum(diff_grad*diff_grad)
    if (r1 == 0){
      lr = 1
    } else {
      r01 = sum(diff_beta*diff_grad)
      lr1 = r01/r1
      lr2 = r0/r01
      lr  = min(c(lr1,lr2,10))
    }
    
    grad0 = grad1
    diff_beta = -lr*grad1
    beta = beta + diff_beta
    res = Y - X %*% beta
    iter = iter + 1
  }
  
  if (intercept){
    beta[-1] = beta[-1] / sdX
    beta[1]  = beta[1] - mX %*% beta[-1]
  } else {
    beta = beta / sdX
  }
  
  outlist = list(coefficients=matrix(beta), residuals=matrix(res))
  return(outlist)
}

es_reg <- function(X, Y, alpha=0.5, robust=TRUE, method='adaptive')
####################################################################
########## Joint quantile & expected shortfall regression ########## 
####################################################################
## Input 
# X : n by p design matrix
# Y : n by 1 response vector
# alpha : quantile level between 0 and 1
# robust : logical flag for returning a robust ES estimator
# method : tuning methods ('adaptive' & 'naive') for choosing the robustification parameter

## Output
# coefficients : estimated regression coefficients 
# residuals : fitted residuals
{
  qr_fit  = conquer(X, Y, tau=alpha)
  qr_nres = qr_fit$res * (qr_fit$res <= 0)
  Ynew    = qr_nres/alpha + (Y - qr_fit$res)
  
  if (robust == FALSE){
    es_coef = lm(Ynew~X)$coef
    } else {
      if (method == 'adaptive'){
        fitting = adaHuber.reg(X, Ynew, method='adaptive')
        es_coef = fitting$coef
        true_tau <- fitting$tau
        return(true_tau)
        # es_coef = huber.reg(X, Ynew, method='adaptive')
        } else if (method == 'naive') {
          n <- nrow(X)
          p <- ncol(X)
          tau = sd(qr_nres) * sqrt(n / (p + log(n)))
          es_coef = hub_reg(alpha*X, alpha*Ynew, c=tau)$coef
          es_coef[1] = es_coef[1]/alpha
        }
    }
  outlist<-list(coef_q=matrix(qr_fit$coef), coef_e=matrix(es_coef))
  return(outlist)
}

es_norm <- function(alpha=0.5, loc=0, sd=1){
  return(loc - sd * dnorm(qnorm(alpha)) / alpha)  
}

es_t <- function(alpha=0.5, df=2){
  quan_t = qt(alpha, df) 
  return( (df + quan_t^2) * dt(quan_t, df) / ((1-df) * pt(quan_t, df)) )
}

rademacher <- function(n=1){
  return( 2*rbinom(n, 1, 0.5) - 1 )
}

unif_sphere <- function(p=1){
  tmp = rnorm(p)
  return( tmp / sum(tmp^2)^0.5 )
}


#new
# Esreg Simulation
esreg_simulation <- function(Z, y) {
  fit <- esreg(y ~ Z, alpha = 0.1)
  coefs = fit$coefficients
  return(coefs)
}


huber_loss <- function(S, tau) {
  binary_1 <- 1*(abs(S) <= tau)
  binary_2 <- 1*(abs(S) > tau)
  output_1 <- S^2 / 2
  output_2 <- tau * abs(S) - ((tau^2) / 2)
  return(binary_1*output_1 + binary_2*output_2)
}

huber_loss_gd <- function(S, tau) {
  binary_1 <- 1*(S < -tau)
  binary_2 <- 1*(abs(S) <= tau)
  binary_3 <- 1*(S > tau)
  output_1 <- -tau
  output_2 <- S
  output_3 <- tau
  return(binary_1*output_1 +  binary_2*output_2 + binary_3*output_3)
}

huber_loss_second_gd <- function(S, tau) {
  tf_arr <- abs(S) <= tau
  return(sum(tf_arr * 1))
}

# output dimension: n*1
orthogonal_score <- function(theta, beta_hat, X, y, alpha) {
  X_mul_beta <- X %*% beta_hat
  X_mul_theta <- X %*% theta
  multiplier_vec <- y <= X_mul_beta
  multiplier_vec <- multiplier_vec * 1
  
  factor_1 <- alpha * X_mul_theta
  factor_2 <- multiplier_vec * y
  factor_3 <- (multiplier_vec - alpha) * X_mul_beta
  return(factor_1 - factor_2 + factor_3)
}

solve_tau <- function(tau, w, n, p) {
  delta = 1/n
  RHS = (p + log(delta^(-1)))
  with_tau = sum(1*(abs(w) >= tau))
  with_w = sum((w*(1*(tau > abs(w))))^2)
  sol = (with_w / (RHS - with_tau))^(0.5)
  return(sol)
}

single_huber <- function(x) {
  tau=70
  if (abs(x) <= tau) {
    return(0.5*x^(0.5))
  } else {
    return(tau*(abs(x) - 0.5*tau))
  }
}

calc_residual <- function(Z, y, alpha, theta) {
  return(y/alpha - (Z %*% theta))#(y - alpha*(Z %*% theta))
}

min_max_normalize <- function(Z) {
  p <- size(Z)[3]
  for(i in 1:p) {
    col <- Z[,i]
    minimum <- min(col)
    maximum <- max(col)
    diff <- maximum - minimum
    Z[,i] <- (col - minimum)/diff
  }
  return(Z)
}

f <- function(tau){
  # print(w)
  return(abs((sum((w*(abs(w) <= tau))^2) + sum((tau*(abs(w) > tau))^2))/(n*(tau^2)) - (p+log(n))/n))
}

f_gd <- function(tau, w) {
  n = length(w)
  contant = (-2)*(w*(abs(w) < tau))^2
  return(sum(contant*(tau)^(-3)) / n)
}

tau_optim <- function(w) {
  # tau = 5
  gd_prev <- 1
  tau_prev <- 20
  alpha <- 1
  loss <- rep(0, 100)
  for (i in 1:100) {
    gd <- f_gd(tau, w)
    loss[i] = f(tau, w) - 0
  }
  plot(loss)
  return(tau)
}

# core
es_reg_new <- function(Z, y, epochs=100, alpha=0.05) {
  
  qr_fit = conquer(Z, y, tau=alpha)
  qr_nres = qr_fit$res * (qr_fit$res <= 0)
  y = qr_nres/alpha + (y - qr_fit$res)
  #y = qr_nres + alpha*(y - qr_fit$res)
  
  bq_beta <- qr_fit$coeff
  theta <- as.vector(rep(0, length(bq_beta)))
  beta_hat <- bq_beta
  losses <- c()
  early_stopping_count <- 0
  cur_min <- 10000
  x_axes <- 1
  
  # tau = 60
  
  ones = rep(1, times=length(Z[,1]))
  old_Z = Z
  Z = cbind(ones, Z)
  
  tau = es_reg(old_Z, y, alpha=alpha, method='adaptive')
  print(tau)
  
  for (i in 1:epochs) {
    # get loss --> Ortho;Huber loss
    loss_vec <- rep(0, length(Z[1,]))
    loss_gd_vec <- rep(0, length(Z[1,]))
    n <- length(Z[,1])
    gd <- rep(0, length(Z[1,]))
    
    w = calc_residual(Z, y, alpha, theta)
    
    #if (i == 1){
    #f <- function(tau){
      # print(w)
      #return((sum((w*(abs(w) <= tau))^2) + sum((tau*(abs(w) > tau))^2))/(n*(tau^2)) - (p+log(n))/n)
    #}
    #tau = spg(par = 70, fn = f, quiet = TRUE)$par}
    # print(tau)
    
    S_vec <- orthogonal_score(theta, beta_hat, Z, y, alpha)
    
    loss_vec <- huber_loss(S_vec, tau=tau)
    loss_gd_vec <- huber_loss_gd(S_vec, tau=tau)
    
    XT_X_inv_XT <- inv(t(as.matrix(Z)) %*% as.matrix(Z)) %*% t(as.matrix(Z))
    psi_z_t <- loss_gd_vec
    
    RHS <- XT_X_inv_XT %*% psi_z_t
    
    summation_gd <- huber_loss_second_gd(S_vec, tau=tau)
    avg_gd <- summation_gd / n
    
    theta <- theta - (1/avg_gd) * RHS
    
    # normalize the vector to record
    loss_per_epoch <- sum(loss_vec)
    losses <- c(losses, loss_per_epoch)
    # print(paste(c("loss by function psi:", loss_per_epoch), collapse = " "))
    x_axes <- x_axes + 1
    
    # early stopping condition
    if (loss_per_epoch+0 > cur_min) {
      early_stopping_count <- early_stopping_count + 1
    } else {
      cur_min <- loss_per_epoch
      early_stopping_count <- 0
    }
    
    if (early_stopping_count > 5) {
      print(paste(c("The iteration early stops at epoch:", x_axes-1), collapse = " "))
      break
    }
  }
  
  #plot(losses, type='l', main = 'Loss Plot', xlab = 'Epoch', ylab = 'Loss')
  # print(tau)
  return(theta)
}
