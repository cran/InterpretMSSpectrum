testthat::test_that(
  desc = "GetIsotopeDistribution works", 
  code = {
    inp <- "C12H40O2S2Si3"
    out <- InterpretMSSpectrum:::GetIsotopeDistribution(fml = inp)
    testthat::expect_true(is.matrix(out))
    testthat::expect_length(out, 6)
    testthat::expect_equal(unname(out[2,]), c(0.6528, 0.1945, 0.1527))
    
    out2 <- InterpretMSSpectrum:::GetIsotopeDistribution(fml = inp, res=100000)
    testthat::expect_equal(unname(out2[2,]), c(0.7206000, 0.1316000, 0.1478000))
  }
)
