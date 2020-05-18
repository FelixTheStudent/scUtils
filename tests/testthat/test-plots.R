context("plot functions")



test_that("zero expression is allowed, negative is not", {
  expect_error(feat(data.frame(u1=1:5, u2=5:1), -1:3))
  p <- feat(data.frame(u1=1:5, u2=5:1), 0:4)
  expect_true(min(p$data$expression) > 0)
})


test_that("closed_breaks_log output is as expected",{
  expect_equal( closed_breaks_log(c(.9, 200.1), base= 2), c(.9, 2, 8, 32, 200.1) )
  expect_equal(closed_breaks_log(c(.1, 15), base=2), c(.1, .25, 1, 4, 15))
  expect_equal(closed_breaks_log(c(.1, 500), base=2), c(.1, .5, 4, 32, 500))
  expect_equal(closed_breaks_log(c(.12444, 977.223), base=2),
               c(.12444, .5, 4, 32, 256, 977.223))
  expect_equal( closed_breaks_log(c(.038, 14.8), base= 2),
                c(.038, .125, .5, 2, 14.8) )
  # zeros throw error; excluding them is not job of the breaks function
  expect_error( closed_breaks_log(c(0, 4.9), base= 2) )
})

test_that("closed_labels output is as expected",{
  expect_equal( closed_labels(c(.9, 2, 8, 32, 200.1)), c(".9", "2", "8", "32", "200.1") )
  expect_equal(closed_labels(c(.1, .25, 1, 4, 15)), c(".1", ".25", "1", "4", "15"))
  expect_equal(closed_labels(c(.1, .5, 4, 32, 500)), c(".1", ".5", "4", "32", "500") )
  expect_equal(closed_labels(c(.12444, .5, 4, 32, 256, 977.223)),
               c(".12", ".5", "4", "32", "256", "977.2"))
  # note that .125 is rounded to .12 on Linux. On OS, this might be different.
  expect_equal( closed_labels(c(.038, .125, .5, 2, 14.8)), c(".04", ".12", ".5", "2", "14.8") )
  expect_equal( closed_labels(c(426, 10000, 161333)), c("426", "1e+04", "1.6e+05") )
})

test_that("in some cases closed_breaks_log returns less than 5 breaks",{
  # not sure yet if this is pleasant behavior, more breaks would be nicer I think.
  # It's just when I implemented the closed_breaks_log function at 2 am, this
  # is where I drew the line and told myself to accept it.
  expect_equal( closed_breaks_log(c(.9, 200.1), base=10), c(.9, 10, 200.1) )
  expect_equal(closed_breaks_log(c(.244, 101.223), base=2),
               c(.244, 1, 8, 101.223))

})


test_that("falling back to log_sub_breaks works", {
  expect_equal(closed_breaks_log(c(.1, .244), base=2), c(0.0625, 0.125, 0.25))
  expect_equal(closed_breaks_log(c(100, 150), base=2), c(64, 128, 256))
})


test_that("all kinds of colnames are allowed", {
  # u1/2 and something that is not u1/2
  p <- feat(data.frame(u1=rnorm(150, c(1, 8)), u2=rnorm(150, c(-1, 1))),
       expression = rlnorm(150, meanlog=c(0, 3)))
  expect_silent(print(p))
  p <- feat(data.frame(d1=rnorm(150, c(1, 8)), d2=rnorm(150, c(-1, 1))),
       expression = rlnorm(150, meanlog=c(0, 3)))
  expect_silent(print(p))
  # special characters in dim names
  p <- feat(data.frame(p_1_2=rnorm(150, c(1, 8)), b..2=rnorm(150, c(-1, 1))),
       expression = rlnorm(150, meanlog=c(0, 3)))
  expect_silent(print(p))
  # colnames are NULL
  df <- data.frame(p_1_2=rnorm(150, c(1, 8)), b..2=rnorm(150, c(-1, 1)))
  colnames(df) <- NULL
  p <- feat(df, expression = rlnorm(150, meanlog=c(0, 3)))
  expect_silent(print(p))
  expect_equal(p$labels, list(y="Dim2", x="Dim1", colour="expression"))
})


test_that("All kinds of 2D embeddings are allowed", {
  # matrix is allowed
  p<-feat(matrix(c(rnorm(150, c(1,8)), rnorm(150, c(-1,1))), ncol=2,
                dimnames = list(as.character(1:150), c("m1", "m2"))),
       expression =rlnorm(150, meanlog=c(0, 3)))
  expect_silent(print(p)); rm(p)
  # matrix with empty dimnames:
  p <- feat(matrix(c(rnorm(150, c(1,8)), rnorm(150, c(-1,1))), ncol=2,
                dimnames = list(NULL)),
       expression =rlnorm(150, meanlog=c(0, 3)))
  expect_silent(print(p)); rm(p)
  # tibble is allowed
  if(requireNamespace("tibble", quietly = TRUE)) {
   p <- feat(tibble::tibble(d1=rnorm(150, c(1, 8)), d2=rnorm(150, c(-1, 1))),
       expression = rlnorm(150, meanlog=c(0, 3)))
   expect_silent(print(p)); rm(p)
  }

})
