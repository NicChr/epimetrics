#' Confusion matrix
#'
#' @description Wrapper around `table()` that
#' creates an n x m confusion matrix of the number of
#' correct and incorrect classes of a multi-class variable where
#' n is the number of distinct predicted classes and m is the
#' number of distinct outcome classes.
#' The only difference between `confusion_matrix()` and `table()`
#' is that it sorts the rows and columns
#' of the predicted and actual classes respectively in descending order,
#' similar to the way that confusion matrices are usually presented.
#'
#' @param outcome Vector of outcome/response/actual classes. e.g. 0, 1.
#' @param prediction Vector of predicted classes. e.g. "A", "B".
#' 
#' @return
#' An n x m `matrix` where `n` is the number of 
#' distinct predicted classes and `m` is the number of 
#' distinct outcome classes.
#'
#' @examples
#' library(epimetrics)
#' library(MASS)
#' 
#' pima_train <- Pima.tr # Training data set
#' pima_test <- Pima.te # Test data set
#' 
#' rf_train <- glm(type ~ ., data = pima_train, family = "binomial")
#' rf_prob_scores <- predict(rf_train, newdata = pima_test,
#'                           type = "response")
#' # Typically you will use an ROC curve to determine best cutoff
#' rf_predictions <- ifelse(rf_prob_scores >= 0.5, 1, 0)
#' outcomes <- pima_test$type
#' cm <- confusion_matrix(outcomes, rf_predictions)
#' cm
#' epimetrics(cm)
#' @export
confusion_matrix <- function(outcome, prediction){
  cm <- table(outcome, prediction)
  cm <- apply(cm, 2, rev)
  cm <- apply(cm, 1, rev)
  cm
}
flatten_confusion <- function(x){
  y <- as.vector(t(x))
  if ( (length(y) %% 4) != 0) {
    stop("number of rows and columns of x must be a multiple of 4")
  }
  flat_cm <- matrix(as.numeric(y), ncol = 4, byrow = TRUE)
  rownames(flat_cm) <- rep_len("", nrow(flat_cm))
  colnames(flat_cm) <- c("TP", "FP", "FN", "TN")
  flat_cm
}

