context("Testing DL calculation functions")
library("testthat")

#-----------------------------------------------------
data("paraquat")
tst <- round(dl_miller(paraquat$x, paraquat$y),1)
test_that("dl_miller test", {
  expect_equal(tst, expected = 0.19 , tolerance = 0.1)
})

# Check object length
test_that("dl_miller object length", {
  expect_length(dl_miller(paraquat$x, paraquat$y), 1)
})

#----------------------------------------------------
data("vhtab6")
d <- adjustdf(vhtab6)
tst <- round(dl_vogelhad(d$x, d$y), 1)
test_that("dl_vogelhad test", {
  expect_equal(tst, expected = 8.7 , tolerance = 0.2)
})

# Check object length
test_that("dl_vogelhad object length", {
  expect_length(dl_vogelhad(d$x, d$y), 1)
})

#----------------------------------------------------
data("mtbe")
tst <- round(dl_hubertvos(mtbe$x, mtbe$y), 2)
test_that("dl_hubertvos test", {
  expect_equal(tst, expected = 1.2 , tolerance = 0.5)
})

# Check object length
test_that("dl_hubertvos object length", {
  expect_length(dl_hubertvos(d$x, d$y), 1)
})
#----------------------------------------------------

