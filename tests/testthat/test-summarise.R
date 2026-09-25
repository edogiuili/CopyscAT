test_that("cutAverage averages only values inside the configured quantiles", {
  session <- local_copyscat_session()
  # lowerTrim 0.5, upperTrim 0.8 -> keep values between the 50th and 80th pct.
  x <- 1:100
  qs <- stats::quantile(x, c(0.5, 0.8), na.rm = TRUE, type = 8)
  expect_equal(cutAverage(x), mean(x[x >= qs[1] & x <= qs[2]]))
})

test_that("cutAverage returns NA for an all-NA vector rather than erroring", {
  session <- local_copyscat_session()
  expect_true(is.na(cutAverage(as.numeric(c(NA, NA, NA)))))
})

test_that("cutAverage ignores NAs among real values", {
  session <- local_copyscat_session()
  expect_false(is.na(cutAverage(c(1, 2, 3, NA, 5))))
})

test_that("makeFakeCells draws the requested number of cells around the mean", {
  set.seed(1)
  cells <- makeFakeCells(10, num = 500, sd_fact = 0.1)
  expect_length(cells, 500)
  expect_equal(mean(cells), 10, tolerance = 0.05)
  expect_equal(stats::sd(cells), 1, tolerance = 0.15)
})

test_that("makeFakeCells with zero SD factor returns a constant vector", {
  expect_equal(unique(makeFakeCells(5, num = 10, sd_fact = 0)), 5)
})
