% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ProbUncertainty.R
\name{scoreCI}
\alias{scoreCI}
\title{Score Confidence Intervals}
\usage{
scoreCI(x, y)
}
\arguments{
\item{x}{A predictor variable.}

\item{y}{A (binomial) response variable.}
}
\value{
A numeric object containing the score confidence
intervals for the difference in probability.
}
\description{
Calculate score confidence intervals of the difference in probability.
}
\details{
In contrast to the Wald Interval, which uses the estimated standard error to
calculate CIs, the score method uses the null standard error. The method used here
is adapted in part from the internal 'z2stat' function of 'diffscoreci' from the R package
{PropCIs}.
}
\examples{
data("startrek")
startrek$uniform <- ifelse(startrek$uniform == "red", 1, 0)
startrek$status <- ifelse(startrek$status == "dead", 1, 0)
scoreCI(x = startrek$uniform, y = startrek$status)

}
\references{
Wilson, E. B. . (1927). Probable Inference, the Law of Succession, and Statistical Inference. Journal of the American Statistical Association, 22(158), 209–212. \doi{https://doi.org/10.2307/2276774}

Brown, L. D., Cai, T. T., & DasGupta, A. (2001). Interval estimation of binomial proportion. Statistical Science, 16(2), 101–133. \doi{https://doi.org/10.1214/ss/1009213286}

Scherer, R. (2018). PropCIs: Various Confidence Interval Methods for Proportions. \url{https://cran.r-project.org/package=PropCIs}
}
