testthat::test_that(
  desc = "get_exactmass works", 
  code = {
    inp <- c("C6H12O6", "Na", "H1", "x")
    out <- InterpretMSSpectrum::get_exactmass(x = inp)
    testthat::expect_length(out, length(inp))
    testthat::expect_equal(names(out), inp)
    testthat::expect_equal(unname(out), c(180.063388, 22.989768, 1.007825, NA))
  }
)
