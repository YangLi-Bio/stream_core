# Script directory
code.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/R/"
data.dir <- "/fs/ess/PCON0022/liyang/STREAM/Codes/stream/data/"


# Binarize matrix
#' @import Matrix
binarize <- function(x) {

  # Libraries
  library(Matrix)


  x.summ <- summary(x)
  y <- x.summ[, 3]
  y[which(y > 0)] <- 1
  # x.summ[, 3] <- y

  xx <- sparseMatrix(
    i = x.summ[, 1],
    j = x.summ[, 2],
    x = y,
    dims = dim(x)
  )

  rownames(xx) <- rownames(x)
  colnames(xx) <- colnames(x)

  xx
}


# Suppress the messages
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}


# Reverse the boolean variables
`%!in%` <- Negate(`%in%`) # define the negation of %in%
is.not.null <- function(x) !is.null(x) # define a function
