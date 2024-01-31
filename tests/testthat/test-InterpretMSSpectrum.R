testthat::test_that(
  desc = "InterpretMSSpectrum works", 
  code = {
    inp <- InterpretMSSpectrum::apci_spectrum
    cp <- "Glutamic acid (3TMS), C14H33NO4Si3, 364.1790"
    mdb <- data.frame(
      "Name" = c("Glutamic acid (3TMS)", "other peak with same sum formula"),
      "Formula" = c("C14H33NO4Si3", "C14H33NO4Si3"),
      "M+H" = c(364.179, 364.179), stringsAsFactors = FALSE, check.names = FALSE
    )
    fdb <- system.file("extdata", "APCI_min.db", package = "InterpretMSSpectrum")

    out <- InterpretMSSpectrum::InterpretMSSpectrum(spec=inp, correct_peak=cp, met_db=mdb, formula_db=fdb, silent=TRUE)
    
    testthat::expect_length(out, 10)
    testthat::expect_equal(attr(out, "stats")[,"mz"], c(202.0893, 246.1339, 274.1288, 348.1477, 364.1789))

    # test the main output as snapshot test
    #testthat::expect_snapshot(InterpretMSSpectrum::InterpretMSSpectrum(spec=inp, correct_peak=cp, met_db=mdb, formula_db=fdb, silent=FALSE))
    
    inp <- InterpretMSSpectrum::esi_spectrum    
    out <- InterpretMSSpectrum::InterpretMSSpectrum(spec=inp, precursor = 263.0534, param = "ESIneg", dppm = 1, silent=TRUE)
    testthat::expect_length(out, 6)
    
  }
)
