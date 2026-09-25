test_that("testOverlap detects overlapping intervals and merges their span", {
  # chr1:1-5 and chr1:3-8 overlap; merged span is 1-8.
  res <- testOverlap("chr1_1.chr1_5", "chr1_3.chr1_8")
  expect_equal(res[3], "TRUE")
  expect_equal(res[1], "1")
  expect_equal(res[2], "8")
})

test_that("testOverlap reports FALSE for disjoint intervals", {
  res <- testOverlap("chr1_1.chr1_2", "chr1_50.chr1_60")
  expect_equal(res[3], "FALSE")
})

test_that("testOverlap is symmetric in its arguments", {
  a <- testOverlap("chr1_10.chr1_20", "chr1_15.chr1_30")
  b <- testOverlap("chr1_15.chr1_30", "chr1_10.chr1_20")
  expect_equal(a[3], b[3])
})

test_that("testOverlap treats a contained interval as overlapping", {
  res <- testOverlap("chr1_1.chr1_100", "chr1_10.chr1_20")
  expect_equal(res[3], "TRUE")
})
