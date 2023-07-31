#' Inverse of the logit link function. 
#'
#' @param x A numeric object.
#' @return A numeric object containing the inverse logit of the input value.
#' @details
#' The inverse logit is defined by \eqn{exp(x)/(1+exp(x))}. Values in x are transformed
#' from a real number (the logarithm of the odds; value between -Inf and Inf) to a probability 
#' with a value between 0 and 1.
#' @usage inv.logit(x)
#' 
#' @export
inv.logit <- function(x) {
  exp(x) / (1 + exp(x))
}




#' Calculate the difference in probability from logistic regression result.
#'
#' @param m A GLM object with a logit link function and a binomial response.
#' @return A numeric object containing the difference in probability  
#' from logistic regression result.
#' @details
#' The difference in probability between two groups is defined by 
#' \eqn{g^{-1}(GLMintercept + GLMslope) - g^{-1}(GLMintercept)}, 
#' where \eqn{g^{-1}()} is the inverse of the logit link function.
#' @usage glm.prob.diff(m)
#' @examples
#' data("startrek")
#' startrek$uniform <- ifelse(startrek$uniform == "red", 1, 0)
#' startrek$status <- ifelse(startrek$status == "dead", 1, 0)
#' glm_redshirt <- glm(status ~ uniform, data = startrek, family = "binomial")
#' glm.prob.diff(glm_redshirt)
#' 
#' @export
glm.prob.diff <- function(m)  {
  if(length(table(m$model[[2]])) != 2) 
    warning("m may not represent a 2x2 contingency table analysis")
  if(m$family$link!="logit")
    warning("m may not be a GLM object with a logit link function")
  
  inv.logit(sum(stats::coef(m))) - inv.logit(stats::coef(m)[1])
}




#' Delta Method Standard Errors
#' 
#' Calculate the standard error of the difference in probability from logistic 
#' regression result using the delta method.
#'
#' @references Ver Hoef, J. M. (2012). Who invented the delta method? American Statistician, 66(2), 124–127. \doi{https://doi.org/10.1080/00031305.2012.687494}
#' @references Lynch, M., & Walsh, B. (1998). Appendix 1: Expectations, Variances, and Covariances of Compound Variables. In Genetics and Analysis of Quantitative Traits (pp. 807–822). Sinauer Associates, Inc.
#' 
#' @param m A GLM object with a logit link function and a binomial response.
#' @return A numeric object containing the standard error for the difference
#' in probability from logistic regression result.
#' @details
#' The delta method is a common technique to approximate the standard errors for 
#' nonlinear transformations of model parameters. It is based on computing the 
#' variance for a Taylor series linearization of a function. For more detailed 
#' information, please see the references. 
#' @usage glm.se.data(m)
#' @examples
#' data("startrek")
#' startrek$uniform <- ifelse(startrek$uniform == "red", 1, 0)
#' startrek$status <- ifelse(startrek$status == "dead", 1, 0)
#' glm_redshirt <- glm(status ~ uniform, data = startrek, family = "binomial")
#' glm.se.data(glm_redshirt)
#' 
#' @export
glm.se.data <- function(m) {
  if(length(table(m$model[[2]])) != 2) 
    warning("m may not represent a 2x2 contingency table analysis")
  if(m$family$link!="logit")
    warning("m may not be a GLM object with a logit link function")
  
  x <- stats::coef(m)[1]
  y <- stats::coef(m)[2]
  b <-
    matrix(c(
      exp(2 * x) / (1 + exp(x)) ^ 2 - exp(x) / (1 + exp(x)) - exp(2 * x + 2 * y) / (1 + exp(x + y)) ^ 2 + exp(x + y) / (1 + exp(x + y)), exp(x + y) / (1 + exp(x + y)) - exp(2 * x + 2 * y) / (1 + exp(x + y)) ^ 2), 2, 1)
  return(sqrt(as.numeric(t(b) %*% stats::vcov(m) %*% b)))
}




#' Profile Likelihood-Based Confidence Intervals
#' 
#' Calculate profile likelihood-based confidence intervals of the difference in probability. 
#' 
#' @references Cole, S. R., Chu, H., & Greenland, S. (2014). Maximum likelihood, profile likelihood, and penalized likelihood: A primer. American Journal of Epidemiology, 179(2), 252–260. \doi{https://doi.org/10.1093/aje/kwt245}
#' 
#' @param x A predictor variable.
#' @param y A (binomial) response variable.
#' @return A numeric object containing the profile likelihood-based confidence 
#' intervals for the difference in probability.
#' @details
#' To generate the profile likelihood CI of the difference in probability, the difference 
#' in probabilities (δp) is 
#' fixed to a range of values between -0.99 and 0.99. The value of p0 
#' (the probability for the first group) that 
#' maximizes the likelihood of observing the data gives the fixed value of δp. 
#' Reported are the upper and lower values of δp for which twice the log 
#' likelihood difference from the unconstrained model was less than 
#' the 95% quantile of a chi-squared distribution with one degree of freedom. 
#' @usage profCI(x, y)
#' @examples
#' data("startrek")
#' startrek$uniform <- ifelse(startrek$uniform == "red", 1, 0)
#' startrek$status <- ifelse(startrek$status == "dead", 1, 0)
#' profCI(x = startrek$uniform, y = startrek$status)
#'
#' @export
profCI <- function(x, y) {
  # log-likelihood at the estimated difference in probabilities between the groups
  est <- tapply(y, x, mean)
  mlest <- sum(log(stats::dbinom(y, 1, est[x + 1])))
  # for any given difference in probabilities, find value of p0
  # that maximizes the log likelihood and return that log likelihood
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




#' Calculate score confidence intervals of the difference in probability. 
#' Internal function of scoreCI, 1.
#' 
#' @references Scherer, R. (2018). PropCIs: Various Confidence Interval Methods for Proportions. \url{https://cran.r-project.org/package=PropCIs}
#' 
#' @param p1x Ratio of success counts in sample 1 / sample size in sample 1
#' @param nx Sample size in sample 1
#' @param p1y Ratio of success counts in sample 2 / sample size in sample 2
#' @param ny  Sample size in sample 2
#' @param dif (p1x - p1y) - (0.5 * (1 + (p1x - p1y)))
#' @details
#' Adapted from the internal 'z2stat' function of 'diffscoreci' from the R package
#' {PropCIs}. For extreme differences, the 'diffscoreci' can throw up errors, which is 
#' ultimately a result of a machine error in floating point math. More specifically, 
#' the ratio \eqn{v/u^3} should never be above one, but it sometimes can be with 
#' floating point error. This internal function should correct the issue.
z2statALT<-function(p1x, nx, p1y, ny, dif){
  diff = p1x - p1y - dif
  if (abs(diff) == 0) {
    fmdiff = 0
  }
  else {
    t = ny/nx
    a = 1 + t
    b = -(1 + t + p1x + t * p1y + dif * (t + 2))
    c = dif * dif + dif * (2 * p1x + t + 1) + p1x + t * p1y
    d = -p1x * dif * (1 + dif)
    v = (b/a/3)^3 - b * c/(6 * a * a) + d/a/2
    s = sqrt((b/a/3)^2 - c/a/3)
    if (v > 0) {
      u = s
    }else {
      u = -s
    }
    if(v/u^3>1){
      w = (3.141592654 + acos(1))/3
    }else{
      w = (3.141592654 + acos(v/u^3))/3
    }
    
    p1d = 2 * u * cos(w) - b/a/3
    p2d = p1d - dif
    nxy = nx + ny
    var = (p1d * (1 - p1d)/nx + p2d * (1 - p2d)/ny) * nxy/(nxy - 
                                                             1)
    fmdiff = diff^2/var
  }
  return(fmdiff)
}




#' Calculate score confidence intervals of the difference in probability.
#' Internal function of scoreCI, 2.
#' 
#' @references Scherer, R. (2018). PropCIs: Various Confidence Interval Methods for Proportions. \url{https://cran.r-project.org/package=PropCIs}
#' 
#' @param x1 Success counts in sample 1
#' @param n1 Sample size in sample 1
#' @param x2 Success counts in sample 2
#' @param n2  Sample size in sample 2
#' @param conf.level Confidence coefficient
#' @details
#' Adapted from the internal 'z2stat' function of 'diffscoreci' from the R package
#' {PropCIs}. For extreme differences, the diffscoreci can throw up errors, which is 
#' ultimately a result of a machine error in floating point math. More specifically, 
#' the ratio \eqn{v/u^3} should never be above one, but it sometimes can be with 
#' floating point error. This internal function should correct it.
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
    score = z2statALT(px, n1, py, n2, up2)
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
    score = z2statALT(px, n1, py, n2, low2)
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



#' Score Confidence Intervals
#' 
#' Calculate score confidence intervals of the difference in probability.
#' 
#' @references Wilson, E. B. . (1927). Probable Inference, the Law of Succession, and Statistical Inference. Journal of the American Statistical Association, 22(158), 209–212. \doi{https://doi.org/10.2307/2276774}
#' @references Brown, L. D., Cai, T. T., & DasGupta, A. (2001). Interval estimation of binomial proportion. Statistical Science, 16(2), 101–133. \doi{https://doi.org/10.1214/ss/1009213286}
#' @references Scherer, R. (2018). PropCIs: Various Confidence Interval Methods for Proportions. \url{https://cran.r-project.org/package=PropCIs}
#' 
#' @param x A predictor variable.
#' @param y A (binomial) response variable.
#' @return A numeric object containing the score confidence 
#' intervals for the difference in probability.
#' @details
#' In contrast to the Wald Interval, which uses the estimated standard error to 
#' calculate CIs, the score method uses the null standard error. The method used here
#' is adapted in part from the internal 'z2stat' function of 'diffscoreci' from the R package
#' {PropCIs}. 
#' @usage scoreCI(x, y)
#' @examples
#' data("startrek")
#' startrek$uniform <- ifelse(startrek$uniform == "red", 1, 0)
#' startrek$status <- ifelse(startrek$status == "dead", 1, 0)
#' scoreCI(x = startrek$uniform, y = startrek$status)
#' 
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