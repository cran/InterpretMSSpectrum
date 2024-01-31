testthat::test_that(
  desc = "check_correct_peak works", 
  code = {
    inp <- "Name, C6H12O6, 0"
    testthat::expect_equal(as.character(InterpretMSSpectrum:::check_correct_peak(x = inp)), inp)
    # return NULL in case that input formula contains other elements than allowed by isotopes parameter
    testthat::expect_null(InterpretMSSpectrum:::check_correct_peak(x = inp, matrix("C", ncol=1)))
    testthat::expect_warning(InterpretMSSpectrum:::check_correct_peak(x = inp, matrix("C", ncol=1), silent=FALSE))
    testthat::expect_warning(InterpretMSSpectrum:::check_correct_peak(x = "wrong, string", silent=FALSE))
    testthat::expect_warning(InterpretMSSpectrum:::check_correct_peak(x = "wrong, string, really", silent=FALSE))
  }
)

testthat::test_that(
  desc = "check_met_db works", 
  code = {
    # minimal valid input
    inp <- data.frame("Name" = "", "Formula" = "C", "mz" = 0)
    # convert 'mz' to 'M+H' in colnames
    testthat::expect_equal(InterpretMSSpectrum:::check_met_db(x = inp), setNames(inp, c("Name", "Formula", "M+H")))
    # return NULL in case that Formula entries are invalid
    testthat::expect_null(InterpretMSSpectrum:::check_met_db(x = data.frame("Name" = "", "Formula" = "", "mz" = 0)))
    testthat::expect_null(InterpretMSSpectrum:::check_met_db(x = data.frame("Name" = "", "Formula" = "X", "mz" = 0)))
    
    # process a more realistic version
    inp <- data.frame(
      "Name" = c("Glutamic acid (3TMS)", "other peak with same sum formula"),
      "Formula" = c("C14H33NO4Si3", "C14H33NO4Si3"),
      "M+H" = c(364.179, 364.179), stringsAsFactors = FALSE, check.names = FALSE
    )
    testthat::expect_true(is.data.frame(InterpretMSSpectrum:::check_met_db(x = inp)))
    # colnames are somewhat case insensitive
    colnames(inp)[1:2] <- tolower(colnames(inp)[1:2])
    testthat::expect_equal(colnames(InterpretMSSpectrum:::check_met_db(x = inp))[1:2], c("Name", "Formula"))
    # omit additional columns
    inp <- cbind(inp, inp)
    testthat::expect_length(InterpretMSSpectrum:::check_met_db(x = inp), 3)
    # return NULL in case that important column cant be found
    inp <- inp[,!colnames(inp)=="name"]
    testthat::expect_null(InterpretMSSpectrum:::check_met_db(x = inp))
  }
)
