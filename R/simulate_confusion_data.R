#' Simulate correct and incorrect observations of diagnostic test result
#'
#' @description This function uses Monte Carlo simulated 
#' diagnostic observations based on the assumption that the 
#' number of true positives and true negatives are independent
#' binomial random variables, using the sample sensitivity and specificity
#' as estimates of the true proportions respectively.
#' @param TP True positives. Must be a numeric vector of length 1.
#' @param FP False positives. Must be a numeric vector of length 1.
#' @param FN False negatives. Must be a numeric vector of length 1.
#' @param TN True negatives. Must be a numeric vector of length 1.
#' @param R Number of Monte Carlo simulations. Numeric vector of length 1.
#'
#' @return
#' An n x 4 matrix where n is the number of simulations.
#'
#' @examples
#' # Using the results from Anu Kantele et al. (2021) as an example
#' library(epimetrics)
#' TP <- 392
#' FP <- 75
#' FN <- 35
#' TN <- 807
#' sim_cm <- simulate_confusion_data(TP, FP, FN, TN, R = 1000)
#' sim_cm
#' ppv(x = sim_cm)
#' TP/(TP + FP)
#' hist(ppv(x = sim_cm)[,1])
#' @export
simulate_confusion_data <- function(TP, FP, FN, TN, R = 10^5){
  stopifnot(length(TP) == 1)
  stopifnot(length(FP) == 1)
  stopifnot(length(FN) == 1)
  stopifnot(length(TN) == 1)
  stopifnot(length(R) == 1)
  n_ones <- TP + FN
  n_zeros <- TN + FP
  TPR <- TP/(TP + FN)
  TNR <- TN/(TN + FP)
  TPs <- stats::rbinom(R, n_ones, TPR)
  TNs <- stats::rbinom(R, n_zeros, TNR)
  FPs <- n_zeros - TNs
  FNs <- n_ones - TPs
  CM <- matrix(c(TPs, FPs, FNs, TNs), byrow = FALSE, ncol = 4, nrow = R)
  # rownames(CM) <- seq_len(nrow(CM))
  colnames(CM) <- c("TP", "FP", "FN", "TN")
  CM
}
