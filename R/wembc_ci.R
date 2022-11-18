#' Wilson-score interval with exact method 
#' boundary correction for extreme cases
#'
#' @description This is an efficient function to compute a 
#' binomial proportion with confidence intervals using the 
#' wilson-score interval method, with a boundary correction for when the 
#' number of successes are either 0, or equal the number of trials.
#'
#' @param x Numeric vector containing the number of successes.
#' @param n Numeric vector containing the number of trials.
#' @param alpha Numeric vector containing the level of significance.
#'
#' @return
#' An n x 3 `matrix` with the columns:
#'  \item{estimate}{point-estimate of the data.}
#'  \item{conf.low}{lower confidence bound of the estimate.}
#'  \item{conf.high}{upper confidence bound of the estimate.}
#'
#'
#' @references
#'  Sakthivel Sivam, S. M. 
#'  "Everything or Nothing-A Better Confidence Intervals for Binomial
#'  Proportion in Clinical Trial Data Analysis." 
#'  Sakthivel Sivam, Quartesian LLC,
#'  Princeton, New Jersey Subbiah Meenakshisundaram, L. N Government College, 
#'  Ponneri, India 2016 (2014).
#'  \url{https://www.lexjansen.com/pharmasug/2016/SP/PharmaSUG-2016-SP08.pdf}
#'
#' @examples
#' library(epimetrics)
#' wembc_ci(0, 0:10)
#' wembc_ci(0:10, 0:10)
#' combs <- expand.grid(0:999, 0:999)
#' wembc_ci(combs[,2], combs[,1])
#' @export
wembc_ci <- function(x, n, alpha = 0.05){
  max_l <- max(length(x), length(n), length(alpha))
  x <- numeric(max_l) + x
  n <- numeric(max_l) + n
  alpha <- numeric(max_l) + alpha
  z <- stats::qnorm(1 - (alpha/2))
  p <- x/n
  q <- 1 - p
  lb_x <- which(x == 0 & n != 0) # Lower bound extreme case
  ub_x <- which(x == n & n != 0) # Upper bound extreme case
  a <- sqrt(z^2 + 4 * n * p * q)
  b <- 2 * (n + z^2)
  lcl <- (2 * n * p + z^2 - z * a)/b
  ucl <- (2 * n * p + z^2 + z * a)/b
  lcl[lb_x] <- 0
  ucl[lb_x] <- 1 - (alpha[lb_x]/2)^(1/n[lb_x])
  lcl[ub_x] <- (alpha[ub_x]/2)^(1/n[ub_x])
  ucl[ub_x] <- 1
  ci_m <- matrix(c(p, lcl, ucl), ncol = 3, byrow = FALSE)
  colnames(ci_m) <- c("estimate", "conf.low", "conf.high")
  ci_m
}
