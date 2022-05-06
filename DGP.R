# Data tests on paper

# random
# not running as gamma, eta not specified
random.draw <- function() {
  z.0 <- matrix(runif(n * p, 0, 1), n, p)
  err <- rt(n, df)
  y.0 <- as.numeric(cbind(1, z.0) %*% gamma + (z.0 %*% eta[-1]) * err)
  outlist <- list(Z = z.0, y = as.vector(y.0))
  return(outlist)
}

#1
eGARCH <- function() {
  z.1 <- rt(1000,7.39)
  mdl <- ugarchspec(variance.model = list(model = 'eGARCH', garchOrder = c(1,1), mean.model = list(armaOrder = c(0,0))))
  x_fit <- ugarchfit(spec = mdl, data = z.1)
  forc <- ugarchforecast(x_fit)
  y.1 <- fitted(x_fit)
  outlist <- list(Z = z.1, y = as.vector(y.1))
  return(outlist)
}

#2
arGARCH <- function() {
  fixed.p <- list(mu = 0, # our mu (intercept)
                  ar1 = 0.1, # our phi_1 (AR(1) parameter of mu_t)
                  ma1 = 0.5, # our theta_1 (MA(1) parameter of mu_t)
                  omega = 0.01, # our alpha_0 (intercept)
                  alpha1 = 0.1, # our alpha_1 (GARCH(1) parameter of sigma_t^2)
                  beta1 = 0.85) 
  armaOrder <- c(1,1) # ARMA order
  garchOrder <- c(1,1) # GARCH order
  varModel <- list(model = "sGARCH", garchOrder = garchOrder)
  spec <- ugarchspec(varModel, mean.model = list(armaOrder = armaOrder),
                     fixed.pars = fixed.p, distribution.model = "norm")
  n <- 1000 # sample size (= length of simulated paths)
  X <- ugarchpath(spec, n.sim = n, m.sim = 1, rseed = 271)
  z.2 <- fitted(X)
  fit <- ugarchfit(spec, data = z.2)
  y.2 <- fitted(fit)
  outlist <- list(Z = z.2, y = as.vector(y.2))
  return(outlist)
}



GAS.STD <- function() {
  library(GAS)
  set.seed(123)
  n = 2500
  p = 1
  z.3 <- rt(n,7.39)
  GASSpec <- UniGASSpec(Dist = "std", ScalingType = "Identity", GASPar = list(location = TRUE, scale = TRUE, shape = TRUE))
  Fit <- UniGASFit(GASSpec, z.3)
  Sim <- UniGASSim(Fit, T.sim = n)
  y.3 = Sim@Data$vY
  outlist <- list(Z = z.3, y = as.vector(y.3))
  return(outlist)
}


GAS.SSTD <- function() {
  library(GAS)
  set.seed(123)
  n = 2500
  p = 1
  z.4 <- rt(1000,7.39)
  GASSpec <- UniGASSpec(Dist = "sstd", ScalingType = "Identity", GASPar = list(location = TRUE, scale = TRUE, skewness = TRUE, shape = TRUE))
  Fit <- UniGASFit(GASSpec, z.4)
  Sim <- UniGASSim(Fit, T.sim = n)
  y.4 = Sim@Data$vY
  outlist <- list(Z = z.4, y = as.vector(y.4))
  return(outlist)
}