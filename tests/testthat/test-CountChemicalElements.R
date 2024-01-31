testthat::test_that(
  desc = "CountChemicalElements returns expected result",
  code = {
    testthat::expect_equal(unname(InterpretMSSpectrum::CountChemicalElements("CHC")), c(2,1))
  }
)
