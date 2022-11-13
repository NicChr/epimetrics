#' Diagnostic metrics to measure the performance of a clinical diagnostic test
#'
#' @description
#'  Standard 2x2 confusion matrix:
#'  \tabular{lcc}{
#'           \tab  Disease \tab  No Disease  \cr
#'  Positive \tab  a      \tab  b    \cr
#'  Negative \tab  c      \tab  d   \cr
#'  }
#'
#'  \bold{Important note} \cr
#' You can either supply `sensitivity`, `specificity` and 
#' `prevalence` \bold{OR} `x`.
#' You cannot supply both as this will throw an error.
#'
#' @param sensitivity Sensitivity or true positive rate of 
#' diagnostic test.
#' @param specificity Specificity or true negative rate of 
#' diagnostic test.
#' @param prevalence Prevalence or proportion of population 
#' who have the disease,
#' @param x An n x 4 matrix containing the numbers of
#' (TP) true positives, (FP) false positives,
#' (FN) false negatives and (TN) true negatives
#' in order from left to right. \cr
#' This can easily be used with `simulate_confusion_data()`,
#' see examples for more detail.\cr
#' You can supply anything that can be coerced to a numeric matrix
#' on a by-row basis.
#' @param cm A 2x2 confusion matrix containing data of 
#' aggregate results and outcomes.
#' Accepted inputs are:
#' \itemize{
#'  \item a vector or factor of the format \code{c(a,b,c,d)} 
#'  based on the table above.
#'  \item a 2x2 table, matrix, or data.frame like the one above.
#'  }
#' @param conf.int Logical. Should confidence intervals be added?
#' @param conf.level (1-alpha) Significance level for 
#' confidence interval calculations.
#' @param R Number of simulation replicates. \cr
#'
#'
#' @return The function `diagnostic_metrics()`
#' returns a `data.frame` containing the variables:
#'  \item{characteristic}{The diagnostic test metric.}
#'  \item{abbreviation}{A 3-letter abbreviation of the characteristic.}
#'  \item{estimate}{Point-estimate of the test characteristic.}
#'  \item{conf.low}{Lower bound of the confidence interval.}
#'  \item{conf.high}{Upper bound of the confidence interval.}
#'
#' The other standalone metrics return either a 
#' numeric vector estimate or matrix
#' containing an estimate and (optionally) confidence bounds. \cr
#' They can be called either by their lower-case 3-letter abbreviation,
#' or by the full name of the metric. \cr
#'  If `x` is supplied then a `matrix` containing the
#'  below 3 variables above is returned:
#'
#'  \item{estimate}{point-estimate of the data.}
#'  \item{conf.low}{lower confidence bound of the estimate.}
#'  \item{conf.high}{upper confidence bound of the estimate.}
#'
#'  If sensitivity, specificity and prevalence are supplied
#'  then a numeric vector containing only the estimate is returned.
#'
#' The following metrics which are all returned by `diagnostic_metrics()`
#' or produced individually using the standalone functions are described:
#'  \item{TPR}{Sensitivity, or true positive rate.
#'  This is the proportion of diseased individuals that test positive.
#'  Pr(+Result|+Outcome)}
#'  \item{TNR}{Specificity, or true negative rate.
#'  This is the proportion of not diseased individuals that test negative.
#'  Pr(-Result|-Outcome)}
#'  \item{ACC}{Accuracy. Proportion of correctly classified individuals.}
#'  \item{BAC}{Balanced accuracy. Median of Sensitivity and Specificity.}
#'  \item{PPV}{Positive predictive value.
#'  Proportion of individuals that have the disease
#'  after receiving a positive result.
#'  Pr(+Outcome|+Result)}
#'  \item{NPV}{Negative predictive value.
#'  Proportion of individuals that don't have the disease
#'  after receiving a negative result.
#'  Pr(-Outcome|-Result)}
#'  \item{PLR}{Positive likelihood ratio.
#'  Ratio of probability of true positive over false positive or
#'  the odds of a positive result being correct.}
#'  \item{NLR}{Negative likelihood ratio.
#'  Ratio of probability of false negative over true negative or the
#'  odds of negative result being incorrect.}
#'  \item{DOR}{Diagnostic odds ratio.
#'  Odds of testing positive if individual has the disease relative to the
#'  odds of testing positive if the individual does not have the disease.}
#'  \item{TPV}{True prevalence. 
#'  Proportion of individuals with the disease.}
#'  \item{APV}{Apparent prevalence. 
#'  Proportion of individuals testing positive.}
#'  \item{PRT}{Prevalence threshold.
#'  Threshold defining the inflection point such that 
#'  the PPV of a diagnostic test drops
#'  precipitously as prevalence decreases below this point.}
#'  \item{YIX}{Youden index.
#'  A diagnostic metric of performance which is the difference
#'  between the rate of correctly classifying one group
#'  against incorrectly classifying the other.
#'  TPR - FPR = TNR - FNR = TPR + TNR - 1.}
#'  \item{NND}{Number needed to diagnose. The number of patients 
#'  that need to be examined
#'  in order to correctly detect one person with the disease of interest.}
#'
#'  All metrics/characteristics are accompanied with confidence intervals.
#'  Wilson score intervals are used for most metrics,
#'  except for the likelihood ratios and diagnostic odds ratio which use
#'  an asymptotic normal approximation,
#'  and the prevalence threshold which uses
#'  Monte Carlo simulated percentile confidence intervals.
#'  All wilson-score intervals are produced based on the
#'  wilson-score method with exact boundary correction for extreme cases
#'  except for the Youden index and number needed to diagnose metrics
#'  which use a standard wilson-score (without continuity correction).
#'
#'  \bold{Notes on the Youden index}
#'
#'  Generally, any diagnostic test with a Youden index less than 0.5
#'  is considered inadequate for diagnostic screening.
#'  The Youden index can be used as a way of choosing a decision threshold in
#'  binary classifier analyses,
#'  where maximising the index across all 
#'  (sensitivity, 1-specificity) pairs on the ROC
#'  curve results in the "optimal" cutoff.
#'  This is mathematically equivalent to maximising
#'  the vertical distance between the diagonal chance line and the curve.
#'  Furthermore, cutoffs based on the Youden index implicitly assume that the
#'  cost of a false positive is equal to the cost of a 
#'  false negative and hence may not be suitable in 
#'  real-world applications where type I errors may incur a
#'  greater clinical, economic, or otherwise 
#'  cost relative to type II errors, or vice-versa.
#'
#' @references
#'  Shan G. Improved Confidence Intervals for the Youden Index.
#'  PLoS One. 2015 Jul 1;10(7):e0127272.
#'  doi: 10.1371/journal.pone.0127272
#'
#'  Bender R. Calculating confidence intervals for the number needed to treat.
#'  Control Clin Trials. 2001 Apr;22(2):102-10.
#'  doi: s0197-2456(00)00134-3
#'
#'  Sakthivel Sivam, S. M. 
#'  "Everything or Nothing-A Better Confidence Intervals for Binomial
#'  Proportion in Clinical Trial Data Analysis." 
#'  Sakthivel Sivam, Quartesian LLC,
#'  Princeton, New Jersey Subbiah Meenakshisundaram, L. N Government College, 
#'  Ponneri, India 2016 (2014).
#'  \url{https://www.lexjansen.com/pharmasug/2016/SP/PharmaSUG-2016-SP08.pdf}
#'  
#'  Balayla J. Prevalence threshold and the geometry of screening curves.
#'  PLoS One. 2020 Oct 7;15(10):e0240215.
#'  doi: 10.1371/journal.pone.0240215
#'
#'  Smits, N. A note on Youden's J and its cost ratio.
#'  BMC Med Res Methodol 10, 89 (2010).
#'  doi: 10.1186/1471-2288-10-89
#'
#' Kantele A, Paajanen J, Turunen S, Pakkanen SH, Patjas A, Itkonen L,
#' Heiskanen E, Lappalainen M, Desquilbet L, Vapalahti O,
#' Hielm-Björkman A.
#' Scent dogs in detection of COVID-19: triple-blinded randomised trial
#' and operational real-life screening in airport setting.
#' BMJ Glob Health. 2022 May;7(5):e008024.
#' doi: 10.1136/bmjgh-2021-008024
#'
#' @examples
#' library(epimetrics)
#' # Using the results from Anu Kantele et al. (2021),
#' # we can calculate the diagnostic accuracy of scent dogs 
#' # in detection of COVID-19
#' # in the randomised trial experiment compared 
#' # to a gold standard rt-PCR test.
#'
#' # Using the results from the sniffed samples of all dogs, 
#' # we add the samples not sniffed and
#' # regard them as test-negatives
#' TP <- 392
#' FP <- 75
#' FN <- 35
#' TN <- 807
#' FN2 <- 456 - 427
#' TN2 <- 1224 - 882
#' confusion <- matrix(c(TP, FP, (FN + FN2),
#'                       (TN + TN2)), nrow = 2, ncol = 2, byrow = TRUE)
#' rownames(confusion) <- c("+", "-")
#' colnames(confusion) <- c("rt-PCR+", "rt-PCR-")
#' confusion
#' # All diagnostic metrics with CIs
#' diagnostic_metrics(confusion)
#' 
#' # You can use standalone metrics too
#' TPR <- sensitivity(x = confusion)
#' TNR <- specificity(x = confusion)
#' PPV <- ppv(x = confusion)
#' TPR
#' TNR
#' PPV
#' 
#' # power sample size calculation (based on normal-approximation)
#' # to achieve 90% sensitivity and specificity
#' # Based on an expected 80% probability that the lower bound of
#' # the 95% CI will be greater than 80%
#' p0 <- 0.8
#' p1 <- 0.9
#' alpha <- 0.05
#' power <- 0.8
#' beta <- 1 - power
#' z <- qnorm(1-(alpha/2))
#' 
#' n <- (p0*(1-p0)*(z + qnorm(1 - beta)*sqrt((p1*(1-p1)) / 
#'                                             ((p0*(1-p0)))))^2 ) / (p1 - p0)^2
#' cat("The minimal number of PCR+ and PCR- is:", ceiling(n))
#' 
#' # We can plot the changing predictive values as prevalence changes
#' 
#' prevalence <- seq(0, 1, 0.01)
#' 
#' plot(prevalence*100, ppv(TPR[,"estimate"],
#'                          TNR[,"estimate"],
#'                          prevalence) * 100,
#'      col = "blue",
#'      xlab = "Prevalence %", ylab = "PPV %",
#'      main = "Predictive value as prevalence changes")
#' points(prevalence * 100, npv(TPR[,"estimate"],
#'                              TNR[,"estimate"],
#'                              prevalence) * 100,
#'        col = "orange")
#' legend(60, 40,
#'        legend = c("PPV", "NPV", "Prevalence threshold"),
#'        col = c("blue", "orange", "black"),
#'        lty = c(2,2,3))
#' abline(v = prevalence_threshold(TPR[,"estimate"],
#'                                 TNR[,"estimate"]) * 100,
#'        lty = 3)
#' 
#' n1 <- TP + FN + FN2 # rt-PCR+
#' n2 <- TN + TN2 + FP # rt-PCR-
#' n <- n1 + n2 # Number of samples
#' 
#' # Use our estimate of TPR and TNR to get expected numbers of outcomes
#' n1 <- round(n * prevalence)
#' n2 <- n - n1
#' 
#' TPR <- TPR[, "estimate"]
#' TNR <- TNR[, "estimate"]
#' 
#' TPs <- round(n1 * TPR)
#' FNs <- n1 - TPs
#' TNs <- round(n2 * TNR)
#' FPs <- n2 - TNs
#' 
#' x <- matrix(c(TPs, FPs, FNs, TNs), byrow = FALSE, ncol = 4)
#' sensitivities <- sensitivity(x = x)
#' sensitivities <- cbind(sensitivities, prevalence * 100)
#' # True Sensitivity doesn't change as it's not a function of prevalence
#' # Uncertainty in the data that produced our
#' # estimate increases as prevalences decreases
#' sensitivities
#' 
#' ppvs <- ppv(x = x)
#' ppvs <- cbind(ppvs, prevalence * 100)
#' ppvs
#' # We can see that in a low prevalence scenario (1%), the ppv
#' # drops to as low as 13% (95%CI: 8%-20%)
#' 
#' # Use simulation to estimate uncertainty in any metric.
#' # Here we choose sensitivity and ppv
#' sim_cm <- simulate_confusion_data(TP, FP, (FN + FN2), (TN + TN2),
#'                                   R = 10^4)
#' sim_sensitivities <- sensitivity(x = sim_cm, conf.int = FALSE)[,"estimate"]
#' hist(sim_sensitivities)
#' sim_sens_lcl <- quantile(sim_sensitivities, 0.05/2)
#' sim_sens_ucl <- quantile(sim_sensitivities, 1 - (0.05/2))
#' abline(v = sim_sens_lcl)
#' abline(v = sim_sens_ucl)
#' # Set conf.int = F for faster estimate
#' sim_ppvs <- as.numeric(ppv(x = sim_cm, conf.int = FALSE))
#' hist(sim_ppvs)
#' sim_ppv_lcl <- quantile(sim_ppvs, 0.05/2)
#' sim_ppv_ucl <- quantile(sim_ppvs, 1 - (0.05/2))
#' abline(v = sim_ppv_lcl)
#' abline(v = sim_ppv_ucl)
#' 
#' # Compare to the usual wilson-score method
#' c(sim_sens_lcl, sim_sens_ucl) # Percentile method
#' sensitivity(x = c(TP, FP, FN + FN2, TN + TN2)) # Wilson-score method
#' c(sim_ppv_lcl, sim_ppv_ucl) # Percentile method
#' ppv(x = c(TP, FP, FN + FN2, TN + TN2)) # Wilson-score method
#' @rdname diagnostic_metrics
#' @export
diagnostic_metrics <- function(cm, conf.int = TRUE, 
                               conf.level = 0.95, R = 10^5){
  if (length(as.vector(t(cm))) != 4) {
    stop(print(knitr::kable(
      matrix(c("a", "b", "c", "d"),
             nrow = 2,
             byrow = TRUE,
             dimnames = list(c("+Test", "-Test"),
                             c("+Outcome", "-Outcome"))))),
      "cm must represent a 2x2 confusion matrix of above form.\n
      Try cm = c(a, b, c, d)")
    }
  x <- flatten_confusion(cm)
  TP <- x[,"TP"][[1]]
  FP <- x[,"FP"][[1]]
  FN <- x[,"FN"][[1]]
  TN <- x[,"TN"][[1]]
  get_characteristics(TP, FP, FN, TN, conf.int = conf.int, 
                      alpha = 1 - conf.level, R = R)
}
#' @rdname diagnostic_metrics
#' @export
sensitivity <- function(sensitivity, specificity, prevalence,
                        x, conf.int = TRUE, conf.level = 0.95){
  x_exists <- contains_empirical(sensitivity, specificity, prevalence, x)
  if (x_exists){
    x <- flatten_confusion(x)
    TP <- x[, "TP", drop = TRUE]
    FN <- x[, "FN", drop = TRUE]
    a <- TP
    b <- TP + FN
    if (conf.int){
      est <- wembc_ci(a, b, alpha = 1 - conf.level)
    } else {
      est <- est_matrix(a / b)
    }
    rownames(est) <- sprintf("(%.0f/%.0f)", a, b)
    est
  } else {
   sensitivity
  }
}
#' @rdname diagnostic_metrics
#' @export
#' @usage NULL
tpr <- sensitivity
#' @rdname diagnostic_metrics
#' @export
specificity <- function(sensitivity, specificity, prevalence,
                        x, conf.int = TRUE, conf.level = 0.95){
  x_exists <- contains_empirical(sensitivity, specificity, prevalence, x)
  if (x_exists){
    x <- flatten_confusion(x)
    TN <- x[, "TN", drop = TRUE]
    FP <- x[, "FP", drop = TRUE]
    a <- TN
    b <- TN + FP
    if (conf.int){
      est <- wembc_ci(a, b, alpha = 1 - conf.level)
    } else {
      est <- est_matrix(a / b)
    }
    rownames(est) <- sprintf("(%.0f/%.0f)", a, b)
    est
  } else {
    specificity
  }
}
#' @rdname diagnostic_metrics
#' @export
#' @usage NULL
tnr <- specificity
#' @rdname diagnostic_metrics
#' @export
accuracy <- function(sensitivity, specificity, prevalence,
                     x, conf.int = TRUE, conf.level = 0.95){
  x_exists <- contains_empirical(sensitivity, specificity, prevalence, x)
  if (x_exists){
    x <- flatten_confusion(x)
    n <- rowSums(x)
    TP <- x[, "TP", drop = TRUE]
    TN <- x[, "TN", drop = TRUE]
    a <- TP + TN
    b <- n
    if (conf.int){
      est <- wembc_ci(a, b, alpha = 1 - conf.level)
    } else {
      est <- est_matrix(a / b)
    }
    rownames(est) <- sprintf("(%.0f/%.0f)", a, b)
    est
  } else {
    (sensitivity * prevalence) + ( (specificity) * (1 - prevalence) )
  }
}
#' @rdname diagnostic_metrics
#' @export
#' @usage NULL
acc <- accuracy
#' @rdname diagnostic_metrics
#' @export
balanced_accuracy <- function(sensitivity, specificity, prevalence,
                              x, conf.int = TRUE, conf.level = 0.95){
  x_exists <- contains_empirical(sensitivity, specificity, prevalence, x)
  if (x_exists){
    x <- flatten_confusion(x)
    TP <- x[, "TP", drop = TRUE]
    FP <- x[, "FP", drop = TRUE]
    FN <- x[, "FN", drop = TRUE]
    TN <- x[, "TN", drop = TRUE]
    if (conf.int){
      est <- (youden_ci(TP, FP, FN, TN, alpha = 1 - conf.level) + 1)/2
    } else {
      TPR <- TP / (TP + FN)
      TNR <- TN / (TN + FP)
      est <- est_matrix((TPR + TNR) / 2)
    }
    est
  } else {
    (sensitivity + specificity)/2
  }
}
#' @rdname diagnostic_metrics
#' @export
#' @usage NULL
bac <- balanced_accuracy
#' @rdname diagnostic_metrics
#' @export
positive_predictive_value <- function(sensitivity, specificity, prevalence,
                x, conf.int = TRUE, conf.level = 0.95){
  x_exists <- contains_empirical(sensitivity, specificity, prevalence, x)
  if (x_exists){
    x <- flatten_confusion(x)
    TP <- x[, "TP", drop = TRUE]
    FP <- x[, "FP", drop = TRUE]
    a <- TP
    b <- TP + FP
    if (conf.int){
      est <- wembc_ci(a, b, alpha = 1 - conf.level)
    } else {
      est <- est_matrix(a / b)
    }
    rownames(est) <- sprintf("(%.0f/%.0f)", a, b)
    est
  } else {
    (sensitivity * prevalence) / 
      ((sensitivity * prevalence) + 
         (1 - specificity) * (1 - prevalence) )
  }
}
#' @rdname diagnostic_metrics
#' @export
#' @usage NULL
ppv <- positive_predictive_value
#' @rdname diagnostic_metrics
#' @export
negative_predictive_value <- function(sensitivity, specificity, prevalence,
                                      x, conf.int = TRUE, conf.level = 0.95){
  x_exists <- contains_empirical(sensitivity, specificity, prevalence, x)
  if (x_exists){
    x <- flatten_confusion(x)
    FN <- x[, "FN", drop = TRUE]
    TN <- x[, "TN", drop = TRUE]
    a <- TN
    b <- TN + FN
    if (conf.int){
      est <- wembc_ci(a, b, alpha = 1 - conf.level)
    } else {
      est <- est_matrix(a / b)
    }
    rownames(est) <- sprintf("(%.0f/%.0f)", a, b)
    est
  } else {
    (specificity * (1 - prevalence) ) / 
      (specificity * (1 - prevalence) + 
         ( (1 - sensitivity) * prevalence ) )
  }
}
#' @rdname diagnostic_metrics
#' @export
#' @usage NULL
npv <- negative_predictive_value
#' @rdname diagnostic_metrics
#' @export
positive_likelihood_ratio <- function(sensitivity, specificity, prevalence,
                x, conf.int = TRUE, conf.level = 0.95){
  x_exists <- contains_empirical(sensitivity, specificity, prevalence, x)
  if (x_exists){
    x <- flatten_confusion(x)
    TP <- x[, "TP", drop = TRUE]
    FP <- x[, "FP", drop = TRUE]
    FN <- x[, "FN", drop = TRUE]
    TN <- x[, "TN", drop = TRUE]
    if (conf.int){
      est <- plr_ci(TP, FP, FN, TN, alpha = 1 - conf.level)
    } else {
      est <- est_matrix((TP / (TP + FN) ) / (FP / (TN + FP)))
    }
    est
  } else {
    # if (missing(prevalence)){
    #   sensitivity / (1 - specificity)
    # } else {
    #   (prevalence*sensitivity) / ( (1 - prevalence)*(1 - specificity) )
    # }
    sensitivity / (1 - specificity)
  }
}
#' @rdname diagnostic_metrics
#' @export
#' @usage NULL
plr <- positive_likelihood_ratio
#' @rdname diagnostic_metrics
#' @export
negative_likelihood_ratio <- function(sensitivity, specificity, prevalence,
                x, conf.int = TRUE, conf.level = 0.95){
  x_exists <- contains_empirical(sensitivity, specificity, prevalence, x)
  if (x_exists){
    x <- flatten_confusion(x)
    TP <- x[, "TP", drop = TRUE]
    FP <- x[, "FP", drop = TRUE]
    FN <- x[, "FN", drop = TRUE]
    TN <- x[, "TN", drop = TRUE]
    if (conf.int){
      est <- nlr_ci(TP, FP, FN, TN, alpha = 1 - conf.level)
    } else {
      est <- est_matrix((FN / (TP + FN) ) / (TN / (TN + FP)))
    }
    est
  } else {
    # if (missing(prevalence)){
    #   (1 - sensitivity) / specificity
    # } else {
    #   ( (prevalence)*(1 - sensitivity) ) / ( (1 - prevalence)*(specificity) )
    # }
    (1 - sensitivity) / specificity
  }
}
#' @rdname diagnostic_metrics
#' @export
#' @usage NULL
nlr <- negative_likelihood_ratio
#' @rdname diagnostic_metrics
#' @export
diagnostic_odds_ratio <- function(sensitivity, specificity, prevalence,
                x, conf.int = TRUE, conf.level = 0.95){
  x_exists <- contains_empirical(sensitivity, specificity, prevalence, x)
  if (x_exists){
    x <- flatten_confusion(x)
    TP <- x[, "TP", drop = TRUE]
    FP <- x[, "FP", drop = TRUE]
    FN <- x[, "FN", drop = TRUE]
    TN <- x[, "TN", drop = TRUE]
    if (conf.int){
      est <- dor_ci(TP, FP, FN, TN, alpha = 1 - conf.level)
    } else {
      PLR <- (TP / (TP + FN) ) / (FP / (TN + FP))
      NLR <- (FN / (TP + FN) ) / (TN / (TN + FP))
      DOR <- PLR / NLR
      est <- est_matrix(DOR)
    }
    est
  } else {
    (sensitivity * specificity) / ((1 - sensitivity) * (1 - specificity))
  }
}
#' @rdname diagnostic_metrics
#' @export
#' @usage NULL
dor <- diagnostic_odds_ratio
#' @rdname diagnostic_metrics
#' @export
true_prevalence <- function(sensitivity, specificity, prevalence,
                            x, conf.int = TRUE, conf.level = 0.95){
  x_exists <- contains_empirical(sensitivity, specificity, prevalence, x)
  if (x_exists){
    x <- flatten_confusion(x)
    TP <- x[, "TP", drop = TRUE]
    FN <- x[, "FN", drop = TRUE]
    n <- rowSums(x)
    a <- TP + FN
    b <- n
    if (conf.int){
      est <- wembc_ci(a, b, alpha = 1 - conf.level)
    } else {
      est <- est_matrix(a / b)
    }
    rownames(est) <- sprintf("(%.0f/%.0f)", a, b)
    est
  } else {
    prevalence
  }
}
#' @rdname diagnostic_metrics
#' @export
#' @usage NULL
tpv <- true_prevalence
#' @rdname diagnostic_metrics
#' @export
apparent_prevalence <- function(sensitivity, specificity, prevalence,
                                x, conf.int = TRUE, conf.level = 0.95){
  x_exists <- contains_empirical(sensitivity, specificity, prevalence, x)
  if (x_exists){
    x <- flatten_confusion(x)
    TP <- x[, "TP", drop = TRUE]
    FP <- x[, "FP", drop = TRUE]
    FN <- x[, "FN", drop = TRUE]
    TN <- x[, "TN", drop = TRUE]
    n <- rowSums(x)
    a <- TP + FP
    b <- n
    if (conf.int){
      est <- wembc_ci(a, b, alpha = 1 - conf.level)
    } else {
      est <- est_matrix(a / b)
    }
    rownames(est) <- sprintf("(%.0f/%.0f)", a, b)
    est
  } else {
    (sensitivity * prevalence) + ( (1 - specificity) * (1 - prevalence) )
  }
}
#' @rdname diagnostic_metrics
#' @export
#' @usage NULL
apv <- apparent_prevalence
#' @rdname diagnostic_metrics
#' @export
prevalence_threshold <- function(sensitivity, specificity, prevalence,
                                 x, conf.int = TRUE, 
                                 conf.level = 0.95, R = 10^5){
  x_exists <- contains_empirical(sensitivity, specificity, prevalence, x)
  if (x_exists){
    x <- flatten_confusion(x)
    TP <- x[, "TP", drop = TRUE]
    FP <- x[, "FP", drop = TRUE]
    FN <- x[, "FN", drop = TRUE]
    TN <- x[, "TN", drop = TRUE]
    if (conf.int){
      est <- prevalence_threshold_ci(TP, FP, FN, TN, 
                                     alpha = 1 - conf.level, R = R)
    } else {
      TPR <- TP / (TP + FN)
      TNR <- TN / (TN + FP)
      a <- (sqrt(TPR * (1 - TNR)) + TNR - 1)
      b <- (TPR + TNR - 1)
      prev_thresh <- a/b
      # Special case
      prev_thresh[a == 0 & b == 0] <- 0.5
      prev_thresh
      est <- est_matrix(prev_thresh)
    }
    est
  } else {
    a <- (sqrt(sensitivity * (1 - specificity)) + specificity - 1)
    b <- (sensitivity + specificity - 1)
    prev_thresh <- a/b
    # Special case
    prev_thresh[a == 0 & b == 0] <- 0.5
    prev_thresh
  }
}
#' @rdname diagnostic_metrics
#' @export
#' @usage NULL
prt <- prevalence_threshold
#' @rdname diagnostic_metrics
#' @export
youden_index <- function(sensitivity, specificity, prevalence,
                         x, conf.int = TRUE, conf.level = 0.95){
  x_exists <- contains_empirical(sensitivity, specificity, prevalence, x)
  if (x_exists){
    x <- flatten_confusion(x)
    TP <- x[, "TP", drop = TRUE]
    FP <- x[, "FP", drop = TRUE]
    FN <- x[, "FN", drop = TRUE]
    TN <- x[, "TN", drop = TRUE]
    if (conf.int){
      est <- youden_ci(TP, FP, FN, TN, alpha = 1 - conf.level)
    } else {
      TPR <- TP / (TP + FN)
      TNR <- TN / (TN + FP)
      est <- est_matrix(TPR + TNR - 1)
    }
    est
  } else {
    sensitivity + specificity - 1
  }
}
#' @rdname diagnostic_metrics
#' @export
#' @usage NULL
yix <- youden_index
#' @rdname diagnostic_metrics
#' @export
number_needed_to_diagnose <- function(sensitivity, specificity, prevalence,
                                      x, conf.int = TRUE, 
                                      conf.level = 0.95, R = 10^5){
  x_exists <- contains_empirical(sensitivity, specificity, prevalence, x)
  if (x_exists){
    x <- flatten_confusion(x)
    TP <- x[, "TP", drop = TRUE]
    FP <- x[, "FP", drop = TRUE]
    FN <- x[, "FN", drop = TRUE]
    TN <- x[, "TN", drop = TRUE]
    if (conf.int){
      est <- nnd_ci(TP, FP, FN, TN, alpha = 1 - conf.level, R = R)
    } else {
      TPR <- TP / (TP + FN)
      TNR <- TN / (TN + FP)
      est <- est_matrix(1/(TPR + TNR - 1))
    }
    est
  } else {
    1/(sensitivity + specificity - 1)
  }
}
#' @rdname diagnostic_metrics
#' @export
#' @usage NULL
nnd <- number_needed_to_diagnose

contains_empirical <- function(sensitivity, specificity, prevalence, x){
  x_exists <- !missing(x)
  rates_exist <- !missing(
    sensitivity) || 
    !missing(specificity) || 
    !missing(prevalence
    )
  if (x_exists && rates_exist){
    stop("Please supply either sensitivity, specificity 
         and prevalence OR a 2x2 confusion matrix x")
  }
  return(x_exists)
}
est_matrix <- function(x){
  M <- matrix(x, ncol = 1, byrow = FALSE)
  colnames(M) <- "estimate"
  M
}
get_characteristics <- function(TP, FP, FN, TN, 
                                conf.int = TRUE, alpha = 0.05, R = 10^5){
  binom_ci <- function(...) wembc_ci(..., 
                                     alpha = alpha)[1,, drop = TRUE][c(2,3)]
  stopifnot(length(TP) == 1)
  stopifnot(length(FP) == 1)
  stopifnot(length(FN) == 1)
  stopifnot(length(TN) == 1)
  stopifnot(length(alpha) == 1)
  stopifnot(length(R) == 1)
  x <- c(TP, FP, FN, TN)
  n <- sum(x)
  if (n == 0) stop("Cannot calculate characteristics with 0 observations")
  TPR <- TP / (TP + FN)
  FNR <- FN / (TP + FN)
  FPR <- FP / (TN + FP)
  TNR <- TN / (TN + FP)
  PPV <- TP / (TP + FP)
  NPV <- TN / (TN + FN)
  PLR <- TPR / (1 - TNR)
  NLR <- (1 - TPR) / TNR
  DOR <- PLR / NLR
  ACC <- (TP + TN) / n
  BAC <- (TPR + TNR)/2
  PRT <- prevalence_threshold(TPR, TNR)
  TPV <- (TP + FN) / n
  APV <- (TP + FP) / n
  YIX <- TPR + TNR - 1
  # YIX <- (TN*TP - FN*FP)/((TN + FP) * (TP + FN))
  NND <- 1 / YIX
  characteristics <- c("sensitivity",
                       "specificity",
                       "accuracy",
                       "balanced_accuracy",
                       "positive_predictive_value",
                       "negative_predictive_value",
                       "positive_likelihood_ratio",
                       "negative_likelihood_ratio",
                       "diagnostic_odds_ratio",
                       "true_prevalence",
                       "apparent_prevalence",
                       "prevalence_threshold",
                       "youden_index",
                       "number_needed_to_diagnose")
  characteristic_abbrv <- c("TPR",
                            "TNR",
                            "ACC",
                            "BAC",
                            "PPV",
                            "NPV",
                            "PLR",
                            "NLR",
                            "DOR",
                            "TPV",
                            "APV",
                            "PRT",
                            "YIX",
                            "NND")
  characteristic_est <- c(TPR, TNR, ACC, BAC, PPV, NPV,
                          PLR, NLR, DOR, TPV, APV, PRT, YIX, NND)
  metrics_tbl <- data.frame(characteristic = characteristics,
                            abbreviation = characteristic_abbrv,
                            estimate = characteristic_est)
  if (conf.int){
    outcomes <- c(rep_len(1, TP + FN), rep_len(0, TN + FP))
    predictions <- c(rep_len(1, TP), rep_len(0, FN),
                     rep_len(1, FP), rep_len(0, TN))
    # Binomial CI using wilson-score
    critical_val <- stats::qnorm(1 - (alpha/2))
    TPR_ci <- binom_ci(TP, (TP + FN))
    TNR_ci <- binom_ci(TN, (TN + FP))
    ACC_ci <- binom_ci(TP + TN, n)
    PPV_ci <- binom_ci(TP, (TP + FP))
    NPV_ci <- binom_ci(TN, (TN + FN))
    TPV_ci <- binom_ci(TP + FN, n)
    APV_ci <- binom_ci(TP + FP, n)
    YIX_ci <- youden_ci(TP, FP, FN, TN, 
                        alpha = alpha)[1,,drop = TRUE][c(2,3)]
    NND_ci <- nnd_ci(TP, FP, FN, TN, 
                     alpha = alpha, R = R)[1,,drop = TRUE][c(2,3)]
    PLR_ci <- plr_ci(TP, FP, FN, TN, 
                     alpha = alpha)[1,,drop = TRUE][c(2,3)]
    NLR_ci <- nlr_ci(TP, FP, FN, TN, 
                     alpha = alpha)[1,,drop = TRUE][c(2,3)]
    # Odds ratio CI
    DOR_ci <- dor_ci(TP, FP, FN, TN, 
                     alpha = alpha)[1,,drop = TRUE][c(2,3)]
    # Simulated CI
    PRT_ci <- prevalence_threshold_ci(TP, FP, FN, TN, 
                                      alpha = alpha, 
                                      R = R)[1,,drop = TRUE][c(2,3)]
    BAC_ci <- (YIX_ci + 1)/2
    characteristic_lcl <- vapply(list(TPR_ci, TNR_ci, ACC_ci, BAC_ci, 
                                      PPV_ci, NPV_ci, PLR_ci, NLR_ci, 
                                      DOR_ci, TPV_ci, APV_ci, 
                                      PRT_ci, YIX_ci, NND_ci),
                                         function(x) x[[1]], 
                                 FUN.VALUE = numeric(1))
    characteristic_ucl <- vapply(list(TPR_ci, TNR_ci, ACC_ci, BAC_ci, 
                                      PPV_ci, NPV_ci, PLR_ci, NLR_ci, 
                                      DOR_ci, TPV_ci, APV_ci, 
                                      PRT_ci, YIX_ci, NND_ci),
                                 function(x) x[[2]], 
                                 FUN.VALUE = numeric(1))
    metrics_tbl[["conf.low"]] <- characteristic_lcl
    metrics_tbl[["conf.high"]] <- characteristic_ucl
  }
  metrics_tbl
}
# Youden index confidence interval using wilson-score interval
youden_ci <- function(TP, FP, FN, TN, alpha = 0.05){
  # As described in Guogen Shan (2015) using wilson-score interval
  n <- TP + FN
  m <- TN + FP
  p1 <- FP/m
  p2 <- TP/n

  # J <- p2 - p1 # Jouden index
  J <- (TP / (TP + FN)) + (TN / (TN + FP)) - 1
  z <- stats::qnorm(1-(alpha/2))
  se_wilson1 <- sqrt( ((p1 * (1-p1))/m) + (z^2)/(4 * m^2))
  se_wilson2 <- sqrt( ((p2 * (1-p2))/n) + (z^2)/(4 * n^2))
  l1 <- (1/(1 + ((z^2)/m))) * (p1 + (z^2/(2*m)) - z*se_wilson1)
  u1 <- (1/(1 + ((z^2)/m))) * (p1 + (z^2/(2*m)) + z*se_wilson1)
  l2 <- (1/(1 + ((z^2)/n))) * (p2 + (z^2/(2*n)) - z*se_wilson2)
  u2 <- (1/(1 + ((z^2)/n))) * (p2 + (z^2/(2*n)) + z*se_wilson2)
  # J_lcl <- J - z * sqrt( ((l1 * (1-l1))/m) + ((u2 * (1-u2))/n))
  # J_ucl <- J + z * sqrt( ((u1 * (1-u1))/m) + ((l2 * (1-l2))/n))
  delta <- z * sqrt( ((u1 * (1-u1))/m) + ((l2 * (1-l2))/n))
  epsilon <- z * sqrt( ((l1 * (1-l1))/m) + ((u2 * (1-u2))/n))
  J_lcl <- J - delta
  J_ucl <- J + epsilon
  J_m <- matrix(c(J, J_lcl, J_ucl), ncol = 3, byrow = FALSE)
  colnames(J_m) <- c("estimate", "conf.low", "conf.high")
  J_m
}
nnd_ci <- function(TP, FP, FN, TN, alpha = 0.05, R = 10^5){
  max_l <- max(length(TP), length(FP), length(FN), length(TN), length(alpha))
  TP <- numeric(max_l) + TP
  FP <- numeric(max_l) + FP
  FN <- numeric(max_l) + FN
  TN <- numeric(max_l) + TN
  alpha <- numeric(max_l) + alpha
  n <- TP + FN
  m <- TN + FP
  p1 <- FP/m
  p2 <- TP/n
  # J <- p2 - p1 # Jouden index
  J <- (TP / (TP + FN)) + (TN / (TN + FP)) - 1
  z <- stats::qnorm(1-(alpha/2))
  se_wilson1 <- sqrt( ((p1 * (1-p1))/m) + (z^2)/(4 * m^2))
  se_wilson2 <- sqrt( ((p2 * (1-p2))/n) + (z^2)/(4 * n^2))
  l1 <- (1/(1 + ((z^2)/m))) * (p1 + (z^2/(2*m)) - z*se_wilson1)
  u1 <- (1/(1 + ((z^2)/m))) * (p1 + (z^2/(2*m)) + z*se_wilson1)
  l2 <- (1/(1 + ((z^2)/n))) * (p2 + (z^2/(2*n)) - z*se_wilson2)
  u2 <- (1/(1 + ((z^2)/n))) * (p2 + (z^2/(2*n)) + z*se_wilson2)
  delta <- z * sqrt( ((u1 * (1-u1))/m) + ((l2 * (1-l2))/n))
  epsilon <- z * sqrt( ((l1 * (1-l1))/m) + ((u2 * (1-u2))/n))
  J_lcl <- J - delta
  J_ucl <- J + epsilon
  nnd_lcl <- 1/J_ucl
  nnd_ucl <- 1/J_lcl
  nnd_m <- matrix(c(1/J, nnd_lcl, nnd_ucl), ncol = 3, byrow = FALSE)
  nnd_m <- cbind(nnd_m, rowSums(sign(nnd_m[,c(2, 3), drop = FALSE])))
  which_span0 <- which(nnd_m[,4, drop = TRUE] == 0)
  nnd_m <- nnd_m[, c(1, 2, 3), drop = FALSE]
  nnd_m_span0 <- nnd_m[which_span0,, drop = FALSE]
  if (length(which_span0 > 0)){
    stopifnot(length(R) == 1)
    message("NND CI spans 0. The largest bounds produced by the 
            collective wilson-score, 
            substitution, 
            and simulated percentile methods will be used.")
    confusion_m <- matrix(c(TP[which_span0], FP[which_span0], 
                            FN[which_span0], TN[which_span0],
                            alpha[which_span0]),
                          byrow = FALSE, ncol = 5)
    w_nnd_lcl <- nnd_m_span0[,3, drop = TRUE]
    w_nnd_ucl <- nnd_m_span0[,2, drop = TRUE]
    # Substitution method
    nnd_se <- sqrt( ((p2*(1 - p2))/n) + ((p1 * (1 - p1))/m))
    exact_nnd_lcl <- 1/( (p2 - p1) - z * nnd_se )
    exact_nnd_ucl <- 1/( (p2 - p1) + z * nnd_se )
    exact_nnd_lcl <- exact_nnd_lcl[which_span0]
    exact_nnd_ucl <- exact_nnd_ucl[which_span0]
    # Percentile Method
    if ((nrow(confusion_m) * R) >= 1e07) {
      warning("N simulations >= 1e07, consider reducing R")
      }
    sim_nnd_ci <- t(apply(confusion_m, 1, function(x){
      sim_nnd <- 1/sim_youden_index(x[1], x[2], x[3], x[4], R = R)
      stats::quantile(sim_nnd, c(x[5]/2, 1 - (x[5]/2)), na.rm = TRUE)
    }))
    sim_nnd_lcl <- sim_nnd_ci[,1, drop = TRUE]
    sim_nnd_ucl <- sim_nnd_ci[,2, drop = TRUE]
    nnd_lcl_span0 <- pmin(w_nnd_lcl, exact_nnd_lcl, sim_nnd_lcl, na.rm = TRUE)
    nnd_ucl_span0 <- pmax(w_nnd_ucl, exact_nnd_ucl, sim_nnd_ucl, na.rm = TRUE)
    nnd_m[which_span0, 2] <- nnd_lcl_span0
    nnd_m[which_span0, 3] <- nnd_ucl_span0
  }
  colnames(nnd_m) <- c("estimate", "conf.low", "conf.high")
  nnd_m
}
sim_youden_index <- function(TP, FP, FN, TN, R = 10^5){
  stopifnot(length(TP) == 1)
  stopifnot(length(FP) == 1)
  stopifnot(length(FN) == 1)
  stopifnot(length(TN) == 1)
  stopifnot(length(R) == 1)
  n_ones <- TP + FN
  n_zeros <- TN + FP
  TPR <- TP / (TP + FN)
  TNR <- TN / (TN + FP)
  TPs <- stats::rbinom(R, n_ones, TPR)
  TNs <- stats::rbinom(R, n_zeros, TNR)
  FPs <- n_zeros - TNs
  FNs <- n_ones - TPs
  TPRs <- TPs / (TPs + FNs)
  TNRs <- TNs / (TNs + FPs)
  TPRs + TNRs - 1
}
plr_ci <- function(TP, FP, FN, TN, alpha = 0.05){
  z <- stats::qnorm(1 - (alpha / 2))
  PLR <- (TP / (TP + FN) ) / (FP / (TN + FP))
  PLR_se <- sqrt((1 / TP) - (1 / (TP + FN)) +
                   (1 / FP) - (1 / (FP + TN)))
  PLR_lcl <- exp(log(PLR) - z * PLR_se)
  PLR_ucl <- exp(log(PLR) + z * PLR_se)
  PLR_m <- matrix(c(PLR, PLR_lcl, PLR_ucl),
                  byrow = FALSE, ncol = 3)
  colnames(PLR_m) <- c("estimate", "conf.low", "conf.high")
  PLR_m
}
nlr_ci <- function(TP, FP, FN, TN, alpha = 0.05){
  z <- stats::qnorm(1 - (alpha/2))
  NLR <- (FN / (TP + FN) ) / (TN / (TN + FP))
  NLR_se <- sqrt((1/FN) - (1/(TP + FN)) +
                   (1/TN) - (1/(FP + TN)))
  NLR_lcl <- exp(log(NLR) - z * NLR_se)
  NLR_ucl <- exp(log(NLR) + z * NLR_se)
  NLR_m <- matrix(c(NLR, NLR_lcl, NLR_ucl),
                  byrow = FALSE, ncol = 3)
  colnames(NLR_m) <- c("estimate", "conf.low", "conf.high")
  NLR_m
}
dor_ci <- function(TP, FP, FN, TN, alpha = 0.05){
  z <- stats::qnorm(1 - (alpha/2))
  PLR <- (TP / (TP + FN) ) / (FP / (TN + FP))
  NLR <- (FN / (TP + FN) ) / (TN / (TN + FP))
  DOR <- PLR/NLR
  DOR_se <- sqrt((1/TP) + (1/TN) + (1/FP) + (1/FN))
  DOR_lcl <- exp(log(DOR) - z * DOR_se)
  DOR_ucl <- exp(log(DOR) + z * DOR_se)
  DOR_m <- matrix(c(DOR, DOR_lcl, DOR_ucl),
                  byrow = FALSE, ncol = 3)
  colnames(DOR_m) <- c("estimate", "conf.low", "conf.high")
  DOR_m
}
sim_prevalence_threshold <- function(TP, FP, FN, TN, R = 10^5){
  stopifnot(length(TP) == 1)
  stopifnot(length(FP) == 1)
  stopifnot(length(FN) == 1)
  stopifnot(length(TN) == 1)
  stopifnot(length(R) == 1)
  n_ones <- TP + FN
  n_zeros <- TN + FP
  TPR <- TP / (TP + FN)
  TNR <- TN / (TN + FP)
  TPs <- stats::rbinom(R, n_ones, TPR)
  TNs <- stats::rbinom(R, n_zeros, TNR)
  FPs <- n_zeros - TNs
  FNs <- n_ones - TPs
  TPRs <- TPs / (TPs + FNs)
  TNRs <- TNs / (TNs + FPs)
  a <- (sqrt(TPRs * (1 - TNRs)) + TNRs - 1)
  b <- (TPRs + TNRs - 1)
  prev_thresh <- a/b
  # Special case
  prev_thresh[a == 0 & b == 0] <- 0.5
  prev_thresh
}
prevalence_threshold_ci <- function(TP, FP, FN, TN, alpha = 0.05, R = 10^5){
  stopifnot(length(R) == 1)
  TPR <- TP / (TP + FN)
  TNR <- TN / (TN + FP)
  a <- (sqrt(TPR * (1 - TNR)) + TNR - 1)
  b <- (TPR + TNR - 1)
  prev_thresh <- a/b
  # Special case
  prev_thresh[a == 0 & b == 0] <- 0.5
  # Simulation
  alpha <- numeric(length(prev_thresh)) + alpha
  n_ones <- TP + FN
  n_zeros <- TN + FP
  confusion_m <- matrix(c(TP, FP, FN, TN, alpha),
                        ncol = 5, byrow = FALSE)
  if ((nrow(confusion_m) * R) >= 1e07) {
    warning("N simulations >= 1e07, consider reducing R")
  }
  sim_prev_thresh_ci <- t(apply(confusion_m, 1, function(x){
    sim_prev_thresh <- sim_prevalence_threshold(x[1], x[2], x[3], x[4], R = R)
    stats::quantile(sim_prev_thresh, c(x[5]/2, 1 - (x[5]/2)), na.rm = TRUE)
  }))
  sim_prev_thresh_lcl <- sim_prev_thresh_ci[,1, drop = TRUE]
  sim_prev_thresh_ucl <- sim_prev_thresh_ci[,2, drop = TRUE]
  prev_thresh_m <- matrix(c(prev_thresh,
                            sim_prev_thresh_lcl,
                            sim_prev_thresh_ucl),
                          byrow = FALSE, ncol = 3)
  colnames(prev_thresh_m) <- c("estimate", "conf.low", "conf.high")
  prev_thresh_m
}

