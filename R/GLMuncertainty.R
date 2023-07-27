#' translate GLM coef from the logit / log-odds scale to the probability scale
#'
#' @param x coefficient in GLM object
inv.logit <- function(x) {
  exp(x) / (1 + exp(x))
}
#' calculate difference in probabilities between two groups from GLM output
#'
#' @param m is the GLM object
#' @export
glm.prob.diff <- function(m)  {
  inv.logit(sum(stats::coef(m))) - inv.logit(stats::coef(m)[1])
}




#' Calculate the standard error of a difference between two properties
#' from a GLM analysis using the delta method
#'
#' @references Ver Hoef, J. M. (2012). Who invented the delta method? American Statistician, 66(2), 124–127. https://doi.org/10.1080/00031305.2012.687494
#' @param m is the GLM object
#' @export
glm.se.data <- function(m) {
  x <- stats::coef(m)[1]
  y <- stats::coef(m)[2]
  # intercept and contrast
  b <-
    matrix(c(
      exp(2 * x) / (1 + exp(x)) ^ 2 - exp(x) / (1 + exp(x)) - exp(2 * x + 2 * y) / (1 + exp(x + y)) ^ 2 + exp(x + y) / (1 + exp(x + y)), exp(x + y) / (1 + exp(x + y)) - exp(2 * x + 2 * y) / (1 + exp(x + y)) ^ 2), 2, 1)
  return(sqrt(as.numeric(t(b) %*% stats::vcov(m) %*% b)))
}




#' Calculate the profile likelihood-based confidence intervals for a difference 
#' in probabilities between two groups.To generate the profile likelihood CI, we 
#' fixed δp to a range of values between -0.99 and 0.99, and found the value of 
#' p0  that maximised the likelihood of observing the data, given the fixed 
#' value of δp. We then recorded the upper and lower values of δ_p for which the 
#' twice the log likelihood difference from the unconstrained model was less than 
#' the 95% quantile of a chi-squared distribution with one degree of freedom. 
#'
#' @param x is predictor variable
#' @param y is response variable
#' @references Cole, S. R., Chu, H., & Greenland, S. (2014). Maximum likelihood, profile likelihood, and penalized likelihood: A primer. American Journal of Epidemiology, 179(2), 252–260. https://doi.org/10.1093/aje/kwt245
#' @export
profCI <- function(x, y) {
  # log-likelihood at the estimated difference in probabilities
  # between the groups
  est <- tapply(y, x, mean)
  mlest <- sum(log(stats::dbinom(y, 1, est[x + 1])))
  # for any given difference in probabilities, find value of p
  # (the probability for the first group) that maximises the
  # log likelihood and return that log likelihood
  MLstar <- function(delta) {
    f <- function(p) {
      sum(log(stats::dbinom(y, 1, p + delta * x)))
    }
    low <- 0.001
    if (delta < 0)
      low = -delta + 0.001

    upp <- 0.999 - delta
    if (delta < 0)
      upp = 0.999

    ml <- stats::optim(
      0.01,
      fn = f,
      method = "L-BFGS-B",
      lower = low,
      upper = upp,
      control = list(fnscale = -1)
    )$value
  }
  # for a range of differences work out the highest likelihood
  deltaRange <- seq(-0.99, 0.99, length.out = 500)
  profLogLik <- data.frame(delta = deltaRange, ploglik =
                             unlist(lapply(deltaRange, MLstar)))
  # find and return range of values of delta values that differ
  profLogLik$logRatio <- mlest - profLogLik$ploglik
  profLogLik <- subset(profLogLik, 2 * profLogLik$logRatio < 3.9415)
  return(c(min(profLogLik$delta), max(profLogLik$delta)))
}








#' Calculate the score confidence interval. In contrast to Wald-type CI, which
#' uses the estimated standard error to calculate CIs, the score method uses
#' the null standard error. A wrapper.
#' @param x1 defined in scoreCI function below
#' @param n1 defined in scoreCI function below
#' @param x2 defined in scoreCI function below
#' @param n2  defined in scoreCI function below
#' @param conf.level defined in scoreCI function below
diffscoreciALT <- function (x1, n1, x2, n2, conf.level) {
  px = x1 / n1
  py = x2 / n2
  z = stats::qchisq(conf.level, 1)
  proot = px - py
  dp = 1 - proot
  niter = 1
  while (niter <= 50) {
    dp = 0.5 * dp
    up2 = proot + dp
    score = PropCIs::z2stat(px, n1, py, n2, up2)
    if (score < z) {
      proot = up2
    }
    niter = niter + 1
    if ((dp < 1e-07) || (abs(z - score) < 1e-06)) {
      niter = 51
      ul = up2
    }
  }
  proot = px - py
  dp = 1 + proot
  niter = 1
  while (niter <= 50) {
    dp = 0.5 * dp
    low2 = proot - dp
    score = PropCIs::z2stat(px, n1, py, n2, low2)
    if (score < z) {
      proot = low2
    }
    niter = niter + 1
    if ((dp < 1e-07) || (abs(z - score) < 1e-06)) {
      ll = low2
      niter = 51
    }
  }
  cint <- c(ll, ul)
  attr(cint, "conf.level") <- conf.level
  rval <- list(conf.int = cint)
  class(rval) <- "htest"
  return(rval)
}




#' Calculate the score confidence interval. In contrast to Wald-type CI, which
#' uses the estimated standard error to calculate CIs, the score method uses
#' the null standard error
#'
#' @param x is predictor variable
#' @param y is response variable
#' 
#' @references Wilson, E. B. . (1927). Probable Inference, the Law of Succession, and Statistical Inference. Journal of the American Statistical Association, 22(158), 209–212. https://doi.org/10.2307/2276774
#' @references Brown, L. D., Cai, T. T., & DasGupta, A. (2001). Interval estimation of binomial proportion. Statistical Science, 16(2), 101–133. https://doi.org/10.1214/ss/1009213286
#' @references Scherer, R. (2018). PropCIs: Various Confidence Interval Methods for Proportions. https://cran.r-project.org/package=PropCIs
#' @export
scoreCI <- function(x, y) {
  n1 <- table(x)["0"]
  n2 <- table(x)["1"]
  x1 <- table(x == 0 &
              y == 1)["TRUE"]
  if (is.na(x1))
    x1 <- 0

  x2 <- table(x == 1 &
              y == 1)["TRUE"]
  if (is.na(x2))
    x2 <- 0

  # group 1 and 2 may be reversed, but that is because the function
  # defines delta the other way around
  ##  as.numeric(diffscoreci(x1=x2,n1=n2,x2=x1,n2=n1,conf.level=0.95)$conf.int)
  as.numeric(diffscoreciALT(
    x1 = x2,
    n1 = n2,
    x2 = x1,
    n2 = n1,
    conf.level = 0.95
  )$conf.int)
}