#' Binarize matrix
#'
#' @importFrom Matrix summary
#'
#' @keywords internal
#'
binarize <- function(x) {

  x.summ <- Matrix::summary(x)
  y <- x.summ[, 3]
  y[which(y > 0)] <- 1
  xx <- Matrix::sparseMatrix(
    i = x.summ[, 1],
    j = x.summ[, 2],
    x = y,
    dims = dim(x)
  )
  rownames(xx) <- rownames(x)
  colnames(xx) <- colnames(x)


  xx
}



#' Suppress the messages
#'
#' @keywords internal
#'
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}



#' Reverse the \code{is.null} function
#'
#' @keywords internal
#'
is.not.null <- function(x) !is.null(x) # define a function
