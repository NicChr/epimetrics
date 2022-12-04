recycle_args <- function(..., length, set_names = FALSE){
  dots <- list(...)
  missing_length <- missing(length)
  if (missing_length){
    recycle_length <- max(lengths(dots))
  } else {
    recycle_length <- length
  }
  if (missing_length && base::length(unique(lengths(dots))) == 1){
    out <- dots
  } else {
    out <- lapply(dots, function(x) rep_len(x, recycle_length))
  }
  if (set_names){
    arg_names <- vapply(match.call(expand.dots = FALSE)[["..."]],
                        FUN = deparse,
                        FUN.VALUE = character(1))
    names(out) <- arg_names
  }
  out
}
max_lengths <- function(...){
  dots <- list(...)
  max(lengths(dots))
}
assign2 <- function(x, values){
  stopifnot(is.list(values))
  for (i in seq_len(length(values))){
    assign(x[i], values[[i]], envir = parent.frame(n = 1))
  }
}