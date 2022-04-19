library(MultiRNG)
library(esreg)

n <- 2500
p <- 3
alpha <- 0.05

# T-error Test
df <- 2.5
gamma <- c(-2, unif_sphere(p)) * sqrt(df/(df-2))
# eta <- 0.5*c(0, rbinom(p, 1, 0.2))
eta = 0.5*c(0, 1, 1, rep(1, p-2)) # prof's ver

beta_q <- gamma + eta * qt(alpha, df)
beta_es <- gamma + eta * es_t(alpha, df)

# calculate running time for each simulation case
num_simlations <- 100

package_diffs <- rep(0, num_simlations)
es_reg_diffs <- rep(0, num_simlations)
given_diffs <- rep(0, num_simlations)
package_time <- rep(0, num_simlations)
es_reg_time <- rep(0, num_simlations)
given_time <- rep(0, num_simlations)

for (i in 1:num_simlations) {
  # DGP
  Z <- matrix(runif(n * p, 0, 1), n, p)
  
  err <- rt(n, df)
  y <- as.numeric(cbind(1, Z) %*% gamma + (Z %*% eta[-1]) * err)
  
  # package framework
  start = Sys.time()
  coefs <- esreg(y~Z, alpha = alpha)
  end = Sys.time()
  package_time[i] <- as.numeric(difftime(end, start, units = "secs"))
  bq_beta_pkg <- coefs$coefficients_q
  be_theta_pkg <- coefs$coefficients_e
  
  
  # customized framework
  start = Sys.time()
  gd_theta <- es_reg_new(Z, y, alpha=alpha)
  end = Sys.time()
  es_reg_time[i] <- as.numeric(difftime(end, start, units = "secs"))
  
  # given framework
  start = Sys.time()
  given_theta <- es_reg(Z, y, alpha=alpha, method='adaptive')
  end = Sys.time()
  given_time[i] <- as.numeric(difftime(end, start, units = "secs"))
  
  # 2-norm difference
  pkg_benchmark_diff <- sum((beta_es - be_theta_pkg) ^ 2)^(0.5)
  simu_benchmark_diff <- sum((beta_es - gd_theta) ^ 2)^(0.5)
  given_benchmark_diff <- sum((beta_es - given_theta$coef_e) ^ 2)^(0.5)
  
  package_diffs[i] <- pkg_benchmark_diff
  es_reg_diffs[i] <- simu_benchmark_diff
  given_diffs[i] <- given_benchmark_diff
  print(paste(c("T-Error Simulation ", as.character(i)), collapse = " "))
}

out <- data.frame(package_diffs, es_reg_diffs, given_diffs)
colnames(out) <- c('package', 'new-es', 'given-es')
boxplot(out, ylab = 'l2-error', main='Error Comparison on T-Error')
time_out <- data.frame(package_time, es_reg_time, given_time)
colnames(time_out) <- c('package', 'new-es', 'given-es')
boxplot(time_out, ylab = 'runtime', main='Runtime Comparison on T-Error')