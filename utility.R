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
  factor_3 <- (alpha - multiplier_vec) * X_mul_beta# (multiplier_vec - alpha) * X_mul_beta
  return(-factor_1 + factor_2 + factor_3)
}

solve_tau <- function(tau, w, n, p) {
  delta = 1/n
  RHS = (p + log10(delta^(-1)))
  with_tau = sum(1*(abs(w) >= tau))
  with_w = sum((abs(w)*(1*(tau > abs(w))))^2)
  sol = (with_w / (RHS - with_tau))^(0.5)
  return(sol)
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

bisection <- function(a, b, f, n, p, res, tol, max_iter) {
  left = a
  right = b
  for (i in 1:max_iter) {
    c = (left + right) / 2
    cur = f(c, res, n, p)
    if (abs(cur) < tol) {
      # print('root found!')
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

func <- function(tau, res, n, p) {
  delta = 1/n
  RHS = p + log(delta^(-1))
  with_tau = sum(1*(abs(res) >= tau))
  with_w = sum((res*(1*(tau > abs(res))))^2/tau^2)
  return((with_w + with_tau)/n - RHS/n)
}

# core
es_reg_new <- function(Z, y, epochs=150, alpha=0.05) {
  # tau = give_tau(Z, y, alpha = alpha, method = 'adaptive')
  
  qr_fit = conquer(Z, y, tau=alpha)
  qr_nres = qr_fit$res * (qr_fit$res <= 0)
  # y = qr_nres/alpha + (y - qr_fit$res)
  y = qr_nres + alpha*(y - qr_fit$res)
  
  bq_beta <- qr_fit$coeff
  theta <- bq_beta #as.vector(rep(0, length(bq_beta)))
  beta_hat <- bq_beta
  losses <- c()
  early_stopping_count <- 0
  cur_min <- 10000
  x_axes <- 1
  # tau = 50 # not robust
  
  ones = rep(1, times=length(Z[,1]))

  Z = cbind(ones, Z)
  
  # return(es_reg_cf(Z, y, alpha, qr_fit$coeff))
  
  for (i in 1:epochs) {
    # get loss --> Ortho;Huber loss
    loss_vec <- rep(0, length(Z[1,]))
    loss_gd_vec <- rep(0, length(Z[1,]))
    n <- length(Z[,1])
    gd <- rep(0, length(Z[1,]))
    
    # -------------------------------------------------------------
    res = calc_residual(Z, y, alpha, theta)
    # tau = solve_tau(tau, res, n, p)
    
    # func <- function(tau) {
      # delta = 1/n
      # RHS = p + log(delta^(-1))
      # with_tau = sum(1*(abs(res) >= tau))
      # with_w = sum((res*(1*(tau > abs(res))))^2/tau^2)
      # return(abs(with_w + with_tau - RHS))
    # }
    
    if (i == 1) {
      tau = bisection(3, 200, func, n, p+1, res, 0.000000001, 100)
      # tau = BBoptim(par=100, fn=func, quiet=TRUE, lower=50)$par
      # print(tau)
    }
    
    
    # -------------------------------------------------------------
    
    S_vec <- orthogonal_score(theta, beta_hat, Z, y, alpha)
    
    loss_vec <- huber_loss(S_vec, tau=tau)
    loss_gd_vec <- huber_loss_gd(S_vec, tau=tau)
    
    XT_X_inv_XT <- inv(t(as.matrix(Z)) %*% as.matrix(Z)) %*% t(as.matrix(Z))
    psi_z_t <- loss_gd_vec
    
    RHS <- XT_X_inv_XT %*% psi_z_t
    
    summation_gd <- huber_loss_second_gd(S_vec, tau=tau)
    avg_gd <- summation_gd / n
    
    theta <- theta + (1/avg_gd) * RHS
    
    # print((1/avg_gd) * RHS)
    
    # sqrt(sum(loss_vec^2))
    
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
  
  plot(losses, type='l', main = 'Loss Plot', xlab = 'Epoch', ylab = 'Loss')
  print(tau)
  return(theta)
}

es_reg_cf <- function(X, y, alpha, beta) {
  LHS = inv(t(X) %*% X) * (1/alpha)
  cond = 1*(y <= X %*% beta)
  vec = y - X %*% beta
  # print(diag(vec))
  overall = rep(0, length(y))
  for (i in 1:length(overall)) {
    overall[i] = vec[i] * cond[i]
  }
  # overall = diag(vec) %*% cond
  RHS = t(X) %*% overall
  return(beta + LHS %*% RHS)
}