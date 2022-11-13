
## A package for easy calculation of diagnostic test metrics with confidence intervals

## Installation

``` r
remotes::install_github("Nic-Chr/epimetrics")
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
diagnostic_metrics(CM)
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
#> 12      prevalence_threshold          PRT  0.26869764 0.24320128  0.30216948
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
#>  [1,] 24 192  6 1808
#>  [2,] 17 184 13 1816
#>  [3,] 22 180  8 1820
#>  [4,] 21 147  9 1853
#>  [5,] 20 172 10 1828
#>  [6,] 20 170 10 1830
#>  [7,] 23 172  7 1828
#>  [8,] 20 180 10 1820
#>  [9,] 17 180 13 1820
#> [10,] 25 196  5 1804
ppv(x = simulated_data, conf.int = FALSE)
#>            estimate
#> (24/216) 0.11111111
#> (17/201) 0.08457711
#> (22/202) 0.10891089
#> (21/168) 0.12500000
#> (20/192) 0.10416667
#> (20/190) 0.10526316
#> (23/195) 0.11794872
#> (20/200) 0.10000000
#> (17/197) 0.08629442
#> (25/221) 0.11312217
ppv(x = simulated_data)
#>            estimate   conf.low conf.high
#> (24/216) 0.11111111 0.07581156 0.1600014
#> (17/201) 0.08457711 0.05347491 0.1312604
#> (22/202) 0.10891089 0.07303049 0.1593885
#> (21/168) 0.12500000 0.08323007 0.1835359
#> (20/192) 0.10416667 0.06844892 0.1554131
#> (20/190) 0.10526316 0.06918068 0.1569911
#> (23/195) 0.11794872 0.07989412 0.1707652
#> (20/200) 0.10000000 0.06567045 0.1494058
#> (17/197) 0.08629442 0.05457499 0.1338396
#> (25/221) 0.11312217 0.07780627 0.1616578
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
diagnostic_metrics(cm)
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
#> 12      prevalence_threshold          PRT  0.26869764 0.24358070  0.30188803
#> 13              youden_index          YIX  0.57666667 0.39730369  0.71818522
#> 14 number_needed_to_diagnose          NND  1.73410405 1.39239846  2.51696630
```
