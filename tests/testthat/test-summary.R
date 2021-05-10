test_that("To see whether the conf.int argument does not occur an error", {
  ex <- EMSS(response = cigs_intervals ~ educ,
             selection = smoker ~ educ + age,
             data = Smoke)
  expect_equal(summary(ex, conf.int = FALSE), summary(ex, conf.int = TRUE))
})
