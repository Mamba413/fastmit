#'@title kNN Mutual Information Estimators
#'
#'@description Estimate mutual information based on the distribution of
#'nearest neighborhood distances. The kNN method is described by Kraskov, et. al (2004).
#'
#'@param x A numeric vector, matrix, data.frame or \code{\link{dist}} object.
#'@param y A numeric vector, matrix, data.frame or \code{\link{dist}} object.
#'@param k Order of neighborhood to be used in the kNN method.
#'@param distance Bool flag for considering \code{x} and \code{y} as distance matrices or not.
#'If \code{distance = TRUE}, \code{x} and \code{y} would be considered as distance matrices,
#'otherwise, these arguments are treated as data and
#'Euclidean distance would be implemented for the samples in \code{x} and \code{y}.
#'Default: \code{distance = FALSE}.
#'
#'@return
#'\item{\code{mi}}{The estimated mutual information.}
#'
#'@details
#'If two samples are passed to arguments \code{x} and \code{y}, the sample sizes
#'(i.e. number of rows of the matrix or length of the vector) must agree.
#'Moreover, data being passed to \code{x} and \code{y} must not contain missing or infinite values.
#'
#'@references Kraskov, A., Stögbauer, H., & Grassberger, P. (2004). Estimating mutual information. Physical review E 69(6): 066138.
#'
#'@useDynLib fastmit, .registration = TRUE
#'@export
#'
#'@import Rcpp
#'@importFrom stats dist
#'
#'@examples
#'library(fastmit)
#'set.seed(1)
#'x <- rnorm(100)
#'y <- x + rnorm(100)
#'mi(x, y, k = 5, distance = FALSE)
#'
#'set.seed(1)
#'x <- rnorm(100)
#'y <- 100 * x + rnorm(100)
#'distx <- dist(x)
#'disty <- dist(y)
#'mi(distx, disty, k = 5, distance = TRUE)
mi <- function(x, y, k = 5, distance = FALSE) {
  stopifnot(k == round(k, 0) && k > 0)
  stopifnot(all(is.finite(x)))
  stopifnot(all(is.finite(y)))

  if (class(x)[1] == "dist" && class(y)[1] == "dist") {
    distance <- TRUE
  }
  if (distance) {
    x <- as.matrix(x)
    y <- as.matrix(y)
    stopifnot(nrow(x) == nrow(x) && nrow(y) == ncol(y) && nrow(x) == nrow(y))
    mi <- knn_mi(x, y, k)
  } else {
    x <- as.matrix(x)
    y <- as.matrix(y)
    stopifnot(nrow(x) == nrow(y))
    distx <- as.matrix(dist(x))
    disty <- as.matrix(dist(y))
    mi <- knn_mi(distx, disty, k)
  }
  return(mi)
}
