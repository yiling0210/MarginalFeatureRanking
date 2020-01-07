library(mvtnorm)
gen_data <- function(n, p, mu_0, mu_1, Sigma_0, Sigma_1) {

  d <- length(mu_0) # the dimension of X

  y <- sample(0:1, n, replace=TRUE, prob=c(1-p, p))
  x <- matrix(NA, nrow=n, ncol=d)

  x[y==0,] <- rmvnorm(n=sum(y==0), mean=mu_0, sigma=Sigma_0)
  x[y==1,] <- rmvnorm(n=sum(y==1), mean=mu_1, sigma=Sigma_1)

  return(list(x=x, y=y))
}
