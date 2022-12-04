
# epimetrics

<!-- badges: start -->

[![R-CMD-check](https://github.com/NicChr/epimetrics/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/NicChr/epimetrics/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/NicChr/epimetrics/branch/main/graph/badge.svg)](https://app.codecov.io/gh/NicChr/epimetrics?branch=main)
<!-- badges: end -->

## A package for easy calculation of diagnostic test metrics with confidence intervals

## Installation

``` r
remotes::install_github("NicChr/epimetrics")
```

``` r
library(epimetrics)
```

You can easily calculate metrics from data, which is typically compared
to a gold standard

``` r
# Using a 2x2 confusion matrix example from wikipedia 
# https://en.wikipedia.org/wiki/Sensitivity_and_specificity
TP <- 20
FP <- 180
FN <- 10
TN <- 1820
CM <- matrix(c(TP, FP, FN, TN), byrow = TRUE, ncol = 2)
rownames(CM) <- c("+Test", "-Test")
colnames(CM) <- c("+Disease", "-Disease")
CM
#>       +Disease -Disease
#> +Test       20      180
#> -Test       10     1820
set.seed(42)
epimetrics(CM)
#>               characteristic abbreviation    estimate   conf.low   conf.high
#> 1                sensitivity          TPR  0.66666667 0.48780052  0.80769502
#> 2                specificity          TNR  0.91000000 0.89665917  0.92176885
#> 3                   accuracy          ACC  0.90640394 0.89295467  0.91831800
#> 4          balanced_accuracy          BAC  0.78833333 0.69865185  0.85909261
#> 5  positive_predictive_value          PPV  0.10000000 0.06567045  0.14940581
#> 6  negative_predictive_value          NPV  0.99453552 0.98997008  0.99702909
#> 7  positive_likelihood_ratio          PLR  7.40740741 5.54896809  9.88826817
#> 8  negative_likelihood_ratio          NLR  0.36630037 0.22078856  0.60771246
#> 9      diagnostic_odds_ratio          DOR 20.22222222 9.32219270 43.86717642
#> 10           true_prevalence          TPV  0.01477833 0.01037124  0.02101836
#> 11       apparent_prevalence          APV  0.09852217 0.08630637  0.11225456
#> 12      prevalence_threshold          PRT  0.26869764 0.24313957  0.30216948
#> 13              youden_index          YIX  0.57666667 0.39730369  0.71818522
#> 14 number_needed_to_diagnose          NND  1.73410405 1.39239846  2.51696630
```

0r from known or theoretical estimates of sensitivity, specificity and
prevalence

``` r
sens <- 0.75 # Sensitivity
spec <- 0.95 # Specificity
prev <- seq(0, 1, 0.01) # Prevalence
ppv <- ppv(sens, spec, prev)
npv <- npv(sens, spec, prev)
plot(prev, ppv, col = "blue")
points(prev, npv, col = "orange")
legend(0.6, 0.4, legend = c("PPV", "NPV"),
       col = c("blue", "orange"),
       pch = c(1, 1))
```

![](man/figures/README-predictive_value-1.png)<!-- -->

You can also use standalone metrics

``` r
sensitivity(x = CM)
#>          estimate  conf.low conf.high
#> (20/30) 0.6666667 0.4878005  0.807695
specificity(x = CM)
#>             estimate  conf.low conf.high
#> (1820/2000)     0.91 0.8966592 0.9217688
ppv(x = CM, conf.level = 0.95)
#>          estimate   conf.low conf.high
#> (20/200)      0.1 0.06567045 0.1494058
```

Easily simulate confusion matrix data

``` r
simulated_data <- simulate_confusion_data(TP, FP, FN, TN, R = 10)
simulated_data
#>       TP  FP FN   TN
#>  [1,] 22 174  8 1826
#>  [2,] 22 165  8 1835
#>  [3,] 23 183  7 1817
#>  [4,] 22 159  8 1841
#>  [5,] 22 169  8 1831
#>  [6,] 23 191  7 1809
#>  [7,] 17 182 13 1818
#>  [8,] 22 187  8 1813
#>  [9,] 16 163 14 1837
#> [10,] 21 198  9 1802
ppv(x = simulated_data, conf.int = FALSE)
#>            estimate
#> (22/196) 0.11224490
#> (22/187) 0.11764706
#> (23/206) 0.11165049
#> (22/181) 0.12154696
#> (22/191) 0.11518325
#> (23/214) 0.10747664
#> (17/199) 0.08542714
#> (22/209) 0.10526316
#> (16/179) 0.08938547
#> (21/219) 0.09589041
ppv(x = simulated_data)
#>            estimate   conf.low conf.high
#> (22/196) 0.11224490 0.07530253 0.1640945
#> (22/187) 0.11764706 0.07898878 0.1716981
#> (23/206) 0.11165049 0.07555952 0.1619601
#> (22/181) 0.12154696 0.08165367 0.1771706
#> (22/191) 0.11518325 0.07730682 0.1682336
#> (23/214) 0.10747664 0.07269141 0.1561055
#> (17/199) 0.08542714 0.05401935 0.1325375
#> (22/209) 0.10526316 0.07054724 0.1542279
#> (16/179) 0.08938547 0.05576870 0.1402561
#> (21/219) 0.09589041 0.06357524 0.1421381
```

If your data are in non-aggregate form, you can use `confusion_matrix`

``` r
outcomes <- c(rep(1, TP + FN), # With disease
              rep(0, FP + TN)) # Without disease
predictions <- c(rep(1, TP), rep(0, FN), # Disease group predictions
                 rep(1, FP), rep(0, TN)) # Non-disease group predictions
(cm <- confusion_matrix(outcomes, predictions))
#>           outcome
#> prediction  1    0
#>          1 20  180
#>          0 10 1820
set.seed(42)
epimetrics(cm)
#>               characteristic abbreviation    estimate   conf.low   conf.high
#> 1                sensitivity          TPR  0.66666667 0.48780052  0.80769502
#> 2                specificity          TNR  0.91000000 0.89665917  0.92176885
#> 3                   accuracy          ACC  0.90640394 0.89295467  0.91831800
#> 4          balanced_accuracy          BAC  0.78833333 0.69865185  0.85909261
#> 5  positive_predictive_value          PPV  0.10000000 0.06567045  0.14940581
#> 6  negative_predictive_value          NPV  0.99453552 0.98997008  0.99702909
#> 7  positive_likelihood_ratio          PLR  7.40740741 5.54896809  9.88826817
#> 8  negative_likelihood_ratio          NLR  0.36630037 0.22078856  0.60771246
#> 9      diagnostic_odds_ratio          DOR 20.22222222 9.32219270 43.86717642
#> 10           true_prevalence          TPV  0.01477833 0.01037124  0.02101836
#> 11       apparent_prevalence          APV  0.09852217 0.08630637  0.11225456
#> 12      prevalence_threshold          PRT  0.26869764 0.24313957  0.30216948
#> 13              youden_index          YIX  0.57666667 0.39730369  0.71818522
#> 14 number_needed_to_diagnose          NND  1.73410405 1.39239846  2.51696630
```
