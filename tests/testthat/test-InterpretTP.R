testthat::test_that(
  desc = "InterpretTP returns expected result",
  code = {
    # simple calculations
    inp <- InterpretMSSpectrum::apci_spectrum
    out <- InterpretMSSpectrum::InterpretTP(fml = "C14H33NO4Si3", spec=inp, param="APCIpos", silent = TRUE)
    testthat::expect_true(is.list(out))
    inp <- InterpretMSSpectrum::esi_spectrum
    out <- InterpretMSSpectrum::InterpretTP(fml = "C14H33NO4Si3", spec=inp, param="ESIneg", silent = TRUE)
    testthat::expect_true(is.list(out))
    out <- InterpretMSSpectrum::InterpretTP(fml = "C14H33NO4Si3", spec=inp, param="ESIpos", silent = TRUE)
    testthat::expect_true(is.list(out))
  }
)
