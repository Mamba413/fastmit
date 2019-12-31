#'@title Mutual Information Test
#'
#'@description Mutual Information test of independence.
#'Mutual Information are generic dependence measures in Banach spaces.
#'
#'@inheritParams mi
#'@param num.permutations The number of permutation replications.
#'If \code{num.permutations = 0}, the function just returns the Mutual Information statistic.
#'Default: \code{num.permutations = 99}.
#'@param seed The random seed. Default: \code{seed = 1}.
#'
#'@return If \code{num.permutations > 0}, \code{mi.test} returns a \code{htest}
#'class object containing the following components:
#'\item{\code{statistic}}{Mutual Information statistic.}
#'\item{\code{p.value}}{The p-value for the test.}
#'\item{\code{replicates}}{Permutation replications of the test statistic.}
#'\item{\code{size}}{Sample size.}
#'\item{\code{alternative}}{A character string describes the alternative hypothesis.}
#'\item{\code{method}}{A character string indicates what type of test was performed.}
#'\item{\code{data.name}}{Description of data.}
#'If \code{num.permutations = 0}, \code{mi.test} returns a statistic value.
#'
#'@details
#'If two samples are passed to arguments \code{x} and \code{y}, the sample sizes
#'(i.e. number of rows of the matrix or length of the vector) must agree.
#'Moreover, data being passed to \code{x} and \code{y} must not contain missing or infinite values.
#'
#'\code{mi.test} utilizes the Mutual Information statistics (see \code{\link{mi}})
#'to measure dependence and derive a \eqn{p}-value via replicating the random permutation \code{num.permutations} times.
#'
#'@export
#'
#'@examples
#' library(fastmit)
#' set.seed(1)
#' error <- runif(50, min = -0.3, max = 0.3)
#' x <- runif(50, 0, 4*pi)
#' y <- cos(x) + error
#' # plot(x, y)
#' res <- mi.test(x, y)
mi.test <- function(x, y, k = 5, distance = FALSE,
                    num.permutations = 99, seed = 1) {
  stopifnot(k == round(k, 0) && k > 0)
  stopifnot(all(is.finite(x)))
  stopifnot(all(is.finite(y)))
  stopifnot(num.permutations == round(num.permutations, 0) && num.permutations >= 0)

  data_name <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  if (length(data_name) > 1) {
    data_name <- ""
  }

  if (class(x)[1] == "dist" && class(y)[1] == "dist") {
    distance <- TRUE
  }
  if (distance) {
    distx <- as.matrix(x)
    disty <- as.matrix(y)
    stopifnot(nrow(x) == nrow(x) && nrow(y) == ncol(y) && nrow(x) == nrow(y))
  } else {
    x <- as.matrix(x)
    y <- as.matrix(y)
    stopifnot(nrow(x) == nrow(y))
    distx <- as.matrix(dist(x))
    disty <- as.matrix(dist(y))
  }

  ini_mi <- knn_mi(distx, disty, k)

  if(num.permutations == 0) return(ini_mi)

  # permuted_mi_value <- c()
  # for(r in 1:num.permutations)
  # {
  #   set.seed(r + seed - 1)
  #   targetindex <- sample(1:nrow(disty), nrow(disty), replace = FALSE)
  #   distypert <- disty[targetindex,targetindex]
  #   permuted_mi_value[r] <- knn_mi(distx, distypert, k)
  # }
  # pvalue <- (sum(abs(permuted_mi_value) > abs(ini_mi)) + 1) / (1 + num.permutations)
  
  pvalue <- mi_test(distx, disty, k, num.permutations, ini_mi)

  alternative_message <- "random variables are dependent"
  test_method <- "Mutual Information test of independence"

  data_name <- paste0(data_name,"\nnumber of observations = ", nrow(distx))
  data_name <- paste0(data_name, "\nreplicates = ", num.permutations)
  names(ini_mi) <- "MI"
  res <- list(estimate = ini_mi,
              p.value = pvalue,
              replicates = num.permutations,
              size = nrow(distx),
              alternative = alternative_message,
              method = test_method,
              data.name = data_name)
  class(res) <- "htest"
  return(res)
}
