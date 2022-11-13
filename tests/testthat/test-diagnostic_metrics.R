
testthat::test_that("Confusion matrix tests", {
  set.seed(1)
  x1 <- sample(c(0, 1), size = 50, replace = TRUE)
  x2 <- sample(c(0, 1), size = 50, replace = TRUE)
  TP <- sum(x1 == 1 & x2 == 1)
  FP <- sum(x1 == 1 & x2 == 0)
  FN <- sum(x1 == 0 & x2 == 1)
  TN <- sum(x1 == 0 & x2 == 0)
  CM <- matrix(c(TP, FP, FN, TN), byrow = TRUE, nrow = 2, ncol = 2,
               dimnames = list("prediction" = c(1, 0),
                               "outcome" = c(1, 0)))
  testthat::expect_identical(CM, confusion_matrix(x2, x1))
  testthat::expect_identical(flatten_confusion(1:4),
                             matrix(as.numeric(1:4), byrow = TRUE, nrow = 1,
                                    dimnames = list("", c("TP", "FP", "FN", "TN"))))
  testthat::expect_error(flatten_confusion(1:5))
})

testthat::test_that("Simulation tests", {
  set.seed(1)
  x1 <- simulate_confusion_data(0, 1, 2, 3, R = 10^5)
  testthat::expect_identical(dim(x1), c(100000L, 4L))
  testthat::expect_true(all(rowSums(x1) == 6))
  testthat::expect_error(simulate_confusion_data(NULL, R = 100))
  testthat::expect_error(simulate_confusion_data(1))
  testthat::expect_error(simulate_confusion_data("1", "2", "3", "4"))
})


testthat::test_that("Compare standalone metrics to diagnostic_metrics", {
  v1 <- c(20, 180, 10, 1820)
  v2 <- factor(c(20, 180, 10, 1820))
  m1 <- matrix(v1, byrow = TRUE, ncol = 2)
  df1 <- as.data.frame(m1)
  tab1 <- as.table(m1)

  data1 <- data.frame(disease = c(rep(1, 30),rep(0, 2000)),
                      result = c(rep(1, 20),rep(0, 10), rep(0, 1820), rep(1, 180)))
  data2 <- data.frame(disease = factor(c(rep(1, 30),rep(0, 2000))),
                      result = factor(c(rep(1, 20),rep(0, 10), rep(0, 1820), rep(1, 180))))
  set.seed(1)
  t1 <- diagnostic_metrics(v1)
  set.seed(1)
  t2 <- diagnostic_metrics(v2)
  set.seed(1)
  t3 <- diagnostic_metrics(m1)
  set.seed(1)
  t4 <- diagnostic_metrics(df1)
  set.seed(1)
  t5 <- diagnostic_metrics(confusion_matrix(data1$disease, data1$result))
  set.seed(1)
  t6 <- diagnostic_metrics(confusion_matrix(data2$disease, data2$result))
  set.seed(1)
  t7 <- diagnostic_metrics(tab1)
  objs <- mget(c("t1", "t2", "t3", "t4", "t5", "t6", "t7"))
  all_equal_pairs <- outer(objs, objs, Vectorize(all.equal))
  all_the_same <- all(all_equal_pairs)
  testthat::expect_true(all_the_same)
})

testthat::test_that("diagnostic_metrics errors", {
  set.seed(42)
  x <- matrix(sample(0:10, size = 4 * 20, replace = TRUE), ncol = 4)
  test <- c(44, 20, 11, 23)
  expected_tbl <- data.frame(characteristic = c(
    'sensitivity',
    'specificity',
    'accuracy',
    'balanced_accuracy',
    'positive_predictive_value',
    'negative_predictive_value',
    'positive_likelihood_ratio',
    'negative_likelihood_ratio',
    'diagnostic_odds_ratio',
    'true_prevalence',
    'apparent_prevalence',
    'prevalence_threshold',
    'youden_index',
    'number_needed_to_diagnose'),
    abbreviation = c('TPR',
      'TNR',
      'ACC',
      'BAC',
      'PPV',
      'NPV',
      'PLR',
      'NLR',
      'DOR',
      'TPV',
      'APV',
      'PRT',
      'YIX',
      'NND'),
    estimate = c(
      0.8,
      0.534883720930233,
      0.683673469387755,
      0.667441860465116,
      0.6875,
      0.676470588235294,
      1.72,
      0.373913043478261,
      4.6,
      0.561224489795918,
      0.653061224489796,
      0.432621812306111,
      0.334883720930232,
      2.98611111111111),
    conf.low = c(
      0.676351208427361,
      0.38915640472314,
      0.586161217162896,
      0.571883665701014,
      0.566076642972519,
      0.508428945939082,
      1.21607523460162,
      0.205711804084279,
      1.88528627155179,
      0.462510021091726,
      0.554661429907658,
      0.384207922611085,
      0.143767331402029,
      2.00641401698805),
    conf.high = c(
      0.884477850902388,
      0.674889422958782,
      0.767329398710772,
      0.749200810882781,
      0.787689334605179,
      0.808683564359745,
      2.43274422159346,
      0.67964482984115,
      11.2237596588358,
      0.655320184199497,
      0.739914083320398,
      0.471490343538977,
      0.498401621765563,
      6.95568311832689
    )
  )
  set.seed(1)
  tbl_res1 <- diagnostic_metrics(test, R = 10^4)
  tbl_res2 <- diagnostic_metrics(test, conf.int = FALSE)
  testthat::expect_equal(tbl_res1, expected_tbl)
  testthat::expect_equal(tbl_res2, expected_tbl[, c("characteristic", 
                                                    "abbreviation", 
                                                    "estimate"),
                                                    drop = FALSE])
  
  testthat::expect_error(diagnostic_metrics(x))
  testthat::expect_error(diagnostic_metrics(NULL))
  testthat::expect_error(diagnostic_metrics(NA))
  testthat::expect_error(diagnostic_metrics(rep_len(0, 4)))
  testthat::expect_error(diagnostic_metrics(1:4, R = NULL))
  testthat::expect_error(diagnostic_metrics(1:4, R = c(1, 2)))
  testthat::expect_error(diagnostic_metrics(x))
})
testthat::test_that("Compare standalone metrics to each other", {
  set.seed(42)
  x <- matrix(sample(10:20, size = 4 * 20, replace = TRUE), ncol = 4)
  sens <- cbind(sensitivity(x = x),
                tpr(x = x, conf.int = FALSE))
  spec <- cbind(specificity(x = x),
                tnr(x = x, conf.int = FALSE))
  acc <- cbind(accuracy(x = x),
               acc(x = x, conf.int = FALSE))
  bac <- cbind(balanced_accuracy(x = x),
               bac(x = x, conf.int = FALSE))
  ppv <- cbind(positive_predictive_value(x = x),
               ppv(x = x, conf.int = FALSE))
  npv <- cbind(negative_predictive_value(x = x),
                npv(x = x, conf.int = FALSE))
  plr <- cbind(positive_likelihood_ratio(x = x),
                plr(x = x, conf.int = FALSE))
  nlr <- cbind(negative_likelihood_ratio(x = x),
               nlr(x = x, conf.int = FALSE))
  dor <- cbind(diagnostic_odds_ratio(x = x),
               dor(x = x, conf.int = FALSE))
  tpv <- cbind(true_prevalence(x = x),
               tpv(x = x, conf.int = FALSE))
  apv <- cbind(apparent_prevalence(x = x),
               apv(x = x, conf.int = FALSE))
  prt <- cbind(prevalence_threshold(x = x),
               prt(x = x, conf.int = FALSE))
  yix <- cbind(youden_index(x = x),
               yix(x = x, conf.int = FALSE))
  testthat::expect_message(nnd(x = x))
  nnd <- suppressMessages(cbind(number_needed_to_diagnose(x = x),
               nnd(x = x, conf.int = FALSE)))
  tol <- 1e-14
  testthat::expect_identical(vapply(list(sens, spec, acc, bac,
                                         ppv, npv, plr, nlr,
                                         dor, tpv, apv, prt,
                                         yix, nnd),
                                    dim, FUN.VALUE = integer(2)),
                             matrix(c(rep_len(20L, 14),
                                      rep_len(4L, 14)),
                                    nrow = 2,
                                    byrow = TRUE))
  testthat::expect_identical(vapply(list(sens, spec, acc, bac,
                                         ppv, npv, plr, nlr,
                                         dor, tpv, apv, prt,
                                         yix, nnd),
                                    function(x) sum((x[,1] - x[,4]) <= tol),
                                    FUN.VALUE = integer(1)),
                             rep_len(20L, 14))
  testthat::expect_error(sensitivity(0.8, 0.9, 0.05, x = 1:4))
  testthat::expect_error(specificity(0.8, 0.9, 0.05, x = 1:4))
  testthat::expect_error(ppv(0.8, 0.9, 0.05, x = 1:4))
  testthat::expect_error(npv(0.8, 0.9, 0.05, x = 1:4))
  testthat::expect_identical(sensitivity(0.8, 0.9, 0.05), 0.8)
  testthat::expect_identical(specificity(0.8, 0.9, 0.05), 0.9)
  testthat::expect_identical(true_prevalence(0.8, 0.9, 0.05), 0.05)
  TPR <- tpr(x = x[1,])[1,1][[1]]
  TNR <- tnr(x = x[1,])[1,1][[1]]
  PREV <- tpv(x = x[1,])[1,1][[1]]
  testthat::expect_equal(accuracy(TPR, TNR, PREV),
                         acc(x = x[1,])[1,1][[1]])
  testthat::expect_equal(balanced_accuracy(TPR, TNR, PREV),
                         bac(x = x[1,])[1,1][[1]])
  testthat::expect_equal(positive_predictive_value(TPR, TNR, PREV),
                             ppv(x = x[1,])[1,1][[1]])
  testthat::expect_equal(negative_predictive_value(TPR, TNR, PREV),
                         npv(x = x[1,])[1,1][[1]])
  testthat::expect_equal(positive_likelihood_ratio(TPR, TNR, PREV),
                         plr(x = x[1,])[1,1][[1]])
  testthat::expect_equal(negative_likelihood_ratio(TPR, TNR, PREV),
                         nlr(x = x[1,])[1,1][[1]])
  testthat::expect_equal(diagnostic_odds_ratio(TPR, TNR, PREV),
                         dor(x = x[1,])[1,1][[1]])
  testthat::expect_equal(true_prevalence(TPR, TNR, PREV),
                         tpv(x = x[1,])[1,1][[1]])
  testthat::expect_equal(apparent_prevalence(TPR, TNR, PREV),
                         apv(x = x[1,])[1,1][[1]])
  testthat::expect_equal(prevalence_threshold(TPR, TNR, PREV),
                         prt(x = x[1,])[1,1][[1]])
  testthat::expect_equal(youden_index(TPR, TNR, PREV),
                         yix(x = x[1,])[1,1][[1]])
  testthat::expect_equal(suppressMessages(number_needed_to_diagnose(TPR, TNR, PREV)),
                         suppressMessages(nnd(x = x[1,])[1,1][[1]]))
})
testthat::test_that("Expect all outputs to be equal", {
  set.seed(42)
  x <- sample(1:100, size = 4)
  result1 <- diagnostic_metrics(x)
  set.seed(42)
  x <- sample(1:100, size = 4)
  sens <- sensitivity(x = x)
  spec <- specificity(x = x)
  acc <- accuracy(x = x)
  bac <- balanced_accuracy(x = x)
  ppv <- ppv(x = x)
  npv <- npv(x = x)
  plr <- plr(x = x)
  nlr <- nlr(x = x)
  dor <- dor(x = x)
  tpv <- true_prevalence(x = x)
  apv <- apparent_prevalence(x = x)
  prt <- prevalence_threshold(x = x)
  yi <- youden_index(x = x)
  nnd <- number_needed_to_diagnose(x = x)

  result2 <- rbind(sens, spec, acc, bac, ppv, npv, plr, nlr, dor, tpv, apv, prt, yi, nnd)

  tol <- 1e-14
  testthat::expect_true(all(abs(result1[,"estimate"] - result2[, "estimate"]) <= tol))
  testthat::expect_true(all(abs(result1[,"conf.low"] - result2[, "conf.low"]) <= tol))
  testthat::expect_true(all(abs(result1[,"conf.high"] - result2[, "conf.high"]) <= tol))
  testthat::expect_warning(prevalence_threshold(x = simulate_confusion_data(10,20,30,40, R = 10^4),
                                                R = 10^3))
  testthat::expect_warning(
    suppressMessages(
      number_needed_to_diagnose(x = simulate_confusion_data(10,20,30,40, 
                                                            R = 10^4),
                                R = 2 * 10^3)
      )
    )
})
