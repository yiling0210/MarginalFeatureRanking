#' classical criterion for marginal feature ranking
#'
#' This function allows you to implement marginal feature ranking based on empirical risk under the classical criterion.
#' @name cc
#' @param x n by p feature matrix
#' @param y response vector of length n
#' @param B number of random splitting of data
#' @param multi.core logical if use parallel computing. Default to FALSE
#' @param ranseed random seed
#' @return a p by 2 dataframe whose first column contains classical criterion values and the second column contains rankings based on the classical criterion values. Rank = 1 has the lowest classical criterion value.
#' @author Jingyi Jessica Li, Yiling Chen (\email{yiling0210@@ucla.edu}), Xin Tong
#' @references \bold{FILL HERE}
#' @examples
#' set.seed(1)
#' gen_data <- function(n, p, mu_0, mu_1, Sigma_0, Sigma_1) {
#' d <- length(mu_0) # the dimension of X
#' y <- sample(0:1, n, replace=TRUE, prob=c(1-p, p))
#' x <- matrix(NA, nrow=n, ncol=d)
#' x[y==0,] <- rmvnorm(n=sum(y==0), mean=mu_0, sigma=Sigma_0)
#' x[y==1,] <- rmvnorm(n=sum(y==1), mean=mu_1, sigma=Sigma_1)
#' return(list(x=x, y=y))
#' }
#' n <- 400 ## sample size
#' B <- 11 ## number of random splits
#' d <- 30 ## total number of features
#' d_select <- 10 ## number of important features
#' mu_0 <- c(rep(-1.5, d_select), rnorm(d-d_select))
#' mu_1 <- c((d_select:1)*.1, mu_0[(d_select+1):d])
#' Sigma_0 <- diag(rep(2^2, d))
#' Sigma_1 <- Sigma_0
#' data = gen_data(n, p, mu_0, mu_1, Sigma_0, Sigma_1)
#' re = cc(data$x, data$y, B, multi.core=F)
#' re


cc <- function(x, y, B, multi.core=F,ranseed = 1001){
  set.seed(ranseed)
  p <- ncol(x)
  n <- nrow(x)
  idx0 <- which(y == 0)
  idx1 <- which(y == 1)
  n0 <- length(idx0)
  n1 <- length(idx1)
  x0 <- x[idx0,]
  x1 <- x[idx1,]

  if (multi.core) {
    res <- mclapply(1:B, FUN=function(b) {
      n0_ts <- ceiling(n0/2)
      n0_lo <- n0 - n0_ts
      idx0_ts <- sample(1:n0, n0_ts)
      x0_ts <- x0[idx0_ts,]
      x0_lo <- x0[-idx0_ts,]

      n1_ts <- ceiling(n1/2)
      n1_lo <- n1 - n1_ts
      idx1_ts <- sample(1:n1, n1_ts)
      x1_ts <- x1[idx1_ts,]
      x1_lo <- x1[-idx1_ts,]

      x_ts <- rbind(x0_ts, x1_ts)
      y_ts <- c(rep(0, n0_ts), rep(1, n1_ts))

      x_lo <- rbind(x0_lo, x1_lo)
      y_lo <- c(rep(0, n0_lo), rep(1, n1_lo))
      n_lo <- n0_lo + n1_lo

      sapply(1:p, FUN=function(j) {
        ################### KDE ratio approach ######################
        # for the j-th feature, construct kde functions on the train-scoring data and evaluate the function at the left-out data
        d0 <- kde(x=x0_ts[,j], eval.points=x_lo[,j])$estimate
        d0[is.na(d0) | d0<=0] <- 10^(-88)
        d1 <- kde(x=x1_ts[,j], eval.points=x_lo[,j])$estimate
        d1[is.na(d1) | d1<=0] <- 10^(-88)
        # classification scores on the left-out data
        s <- d1 / d0
        s0 <- s[1:n0_lo]
        s1 <- s[(n0_lo+1):n_lo]
        ## classical criterion
        t_cl <- n0_ts/n1_ts
        y_lo_hat <- as.numeric( s > t_cl )
        R_kde <- sum((y_lo - y_lo_hat) != 0) / n_lo
        ## NP criterion
        # classification scores on the left-out class 0 sample
        # t <- s0
        # t <- sort(t)
        # # type II errors corresponding to varying alphas
        # R1_kde_s <- sapply(alpha_s, FUN=function(alpha) {
        #   # thresholds
        #   t_np <- t[np_order_stat(n0_lo, alpha, delta)]
        #   # type II error
        #   sum(s1 <= t_np) / n1_lo
        # })

        # return(c(R_kde, R1_kde_s))
        return(R_kde)
      })
    }, mc.cores=detectCores()-1)
  } else {
    res <- lapply(1:B, FUN=function(b) {
      n0_ts <- ceiling(n0/2)
      n0_lo <- n0 - n0_ts
      idx0_ts <- sample(1:n0, n0_ts)
      x0_ts <- x0[idx0_ts,]
      x0_lo <- x0[-idx0_ts,]

      n1_ts <- ceiling(n1/2)
      n1_lo <- n1 - n1_ts
      idx1_ts <- sample(1:n1, n1_ts)
      x1_ts <- x1[idx1_ts,]
      x1_lo <- x1[-idx1_ts,]

      x_ts <- rbind(x0_ts, x1_ts)
      y_ts <- c(rep(0, n0_ts), rep(1, n1_ts))

      x_lo <- rbind(x0_lo, x1_lo)
      y_lo <- c(rep(0, n0_lo), rep(1, n1_lo))
      n_lo <- n0_lo + n1_lo

      sapply(1:p, FUN=function(j) {
        ################### KDE ratio approach ######################
        # for the j-th feature, construct kde functions on the train-scoring data and evaluate the function at the left-out data
        d0 <- kde(x=x0_ts[,j], eval.points=x_lo[,j])$estimate
        d0[is.na(d0) | d0<=0] <- 10^(-88)
        d1 <- kde(x=x1_ts[,j], eval.points=x_lo[,j])$estimate
        d1[is.na(d1) | d1<=0] <- 10^(-88)
        # classification scores on the left-out data
        s <- d1 / d0
        s0 <- s[1:n0_lo]
        s1 <- s[(n0_lo+1):n_lo]
        ## classical criterion
        t_cl <- n0_ts/n1_ts
        y_lo_hat <- as.numeric( s > t_cl )
        R_kde <- sum((y_lo - y_lo_hat) != 0) / n_lo
        ## NP criterion
        # classification scores on the left-out class 0 sample
        # t <- s0
        # t <- sort(t)
        # # type II errors corresponding to varying alphas
        # R1_kde_s <- sapply(alpha_s, FUN=function(alpha) {
        #   # thresholds
        #   t_np <- t[np_order_stat(n0_lo, alpha, delta)]
        #   # type II error
        #   sum(s1 <= t_np) / n1_lo
        # })
        #
        #
        # return(c(R_kde, R1_kde_s))
        return(R_kde)
      })
    })
  }
  for (b0 in 1:B) {
    tmp <- res[[b0]]
    if (is.numeric(tmp)) {
      break
    }
  }
  count_tmp <- 1
  if (b0 < B) {
    for (b in (b0+1):B) {
      if (is.numeric(res[[b]])) {
        tmp <- tmp + res[[b]]
        count_tmp <- count_tmp + 1
      }
    }
  }


  CC = tmp/count_tmp
  rk = val_to_rank(CC)
  out = cbind.data.frame(CC= CC, rank = rk)
  return(out)


}


