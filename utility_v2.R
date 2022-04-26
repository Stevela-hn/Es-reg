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
        print(true_tau)
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

give_tau <- function(X, Y, alpha=0.5, robust=TRUE, method='adaptive')
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
# Es-reg functions
huber_loss <- function(S, tau) {
  mins = pmin(abs(S), tau)
  return(mins*(abs(S) - mins/2))
}

huber_loss_gd <- function(S, tau) {
  return(pmin(pmax(S,-tau),tau))
}

huber_loss_second_gd <- function(S, tau) {
  return(sum(abs(S) <= tau))
}

# Neyman Orthogonal Score
orthogonal_score <- function(theta, beta_hat, X, y, alpha) {
  X_mul_beta <- X %*% beta_hat
  multiplier_vec <- y <= X_mul_beta
  factor_1 <- alpha * (X %*% theta)
  factor_2 <- multiplier_vec * y
  factor_3 <- (alpha - multiplier_vec) * X_mul_beta
  return(-factor_1 + factor_2 + factor_3)
}

calc_residual <- function(Z, y, alpha, theta) {
  return(y/alpha - (Z %*% theta))#(y - alpha*(Z %*% theta))
}

# search methodology
bisection <- function(a, b, f, n, p, res, tol, max_iter) {
  left = a
  right = b
  for (i in 1:max_iter) {
    c = (left + right) / 2
    cur = f(c, res, n, p)
    if (abs(cur) < tol) {
      break
    }
    if (f(left, res, n, p)*cur < 0) {
      right = c
    }
    if (cur * f(right, res, n, p) < 0) {
      left = c
    }
  }
  return(c)
}

# equation to solve for tau
func <- function(tau, res, n, p) {
  delta = 1/n
  RHS = p + log(delta^(-1))
  with_tau = sum(abs(res) >= tau)
  with_w = sum((res*(tau > abs(res)))^2/tau^2)
  return((with_w + with_tau)/n - RHS/n)
}

# core
es_reg_new <- function(Z, y, epochs=150, alpha=0.05) {
  n <- length(Z[,1])
  qr_fit = conquer(Z, y, tau=alpha)
  qr_nres = qr_fit$res * (qr_fit$res <= 0)
  # y = qr_nres/alpha + (y - qr_fit$res)
  y = qr_nres + alpha*(y - qr_fit$res)
  
  bq_beta <- qr_fit$coeff
  theta <- bq_beta
  beta_hat <- bq_beta
  losses <- c()
  early_stopping_count <- 0
  cur_min <- 10000
  
  ones = rep(1, times=n)
  # get ready for regression
  Z = cbind(ones, Z)
  
  for (i in 1:epochs) {
    res = calc_residual(Z, y, alpha, theta)
    
    # optimization for tau
    if (i == 1) {
      tau = bisection(3, 200, func, n, p+1, res, 0.000000001, 100)
    }
    
    # update process
    S <- orthogonal_score(theta, beta_hat, Z, y, alpha)
    loss <- huber_loss(S, tau=tau)
    loss_gd <- huber_loss_gd(S, tau=tau)
    
    update <- inv(t(Z) %*% Z) %*% t(Z) %*% loss_gd
    
    second_gd <- huber_loss_second_gd(S, tau=tau)
    
    theta <- theta + (1/(second_gd / n)) * update
    
    # sum of error -> 1-norm
    loss_sum <- sum(abs(loss))
    losses <- c(losses, loss_sum)
    
    # early stopping condition
    if (loss_sum > cur_min) {
      early_stopping_count <- early_stopping_count + 1
    } else {
      cur_min <- loss_sum
      early_stopping_count <- 0
    }
    
    if (early_stopping_count > 5) {
      print(paste(c("The iteration early stops at epoch:", i), collapse = " "))
      break
    }
  }
  
  plot(losses, type='l', main = 'Loss Plot', xlab = 'Epoch', ylab = 'Loss')
  print(tau)
  return(theta)
}