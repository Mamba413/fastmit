library(testthat)
library(fastmit)

context("mi and mi.test function")
skip_on_cran()

test_that("Error if computation result for mutual information is wrong!", {
  target_value <- 0.3856349
  names(target_value) <- "MI"
  expect_equal(mi(1:10, 1:10), target_value)
  expect_equal(mi.test(1:10, 1:10, num.permutations = 0), target_value)
  dx <- dist(1:10)
  dy <- dist(1:10)
  expect_equal(mi(dx, dy, distance = TRUE), target_value)
  expect_equal(mi.test(dx, dy, distance = TRUE, num.permutations = 0), target_value)
})
