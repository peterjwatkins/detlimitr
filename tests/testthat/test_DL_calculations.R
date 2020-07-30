context("Testing DL calculation functions")
library("testthat")

data("paraquat")
tst <- round(dl_miller(paraquat$x, paraquat$y)[1],1)
test_that("dl_miller test", {
  expect_equal(tst, expected = 0.19 , tolerance = 0.1)
})

data("vhtab6")
d <- adjustcolnames(vhtab6)
tst <- round(dl_vogelhad(d$x, d$y)[1], 1)
test_that("dl_vogelhad test", {
  expect_equal(tst, expected = 8.7 , tolerance = 0.5)
})

data("mtbe")
tst <- round(dl_hubertvos(mtbe$x, mtbe$y)[1], 1)
test_that("dl_hubertvos test", {
  expect_equal(tst, expected = 1.2 , tolerance = 0.5)
})
