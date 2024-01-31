testthat::test_that(
	desc = "ReadSpecClipboard works as expected", 
	code = {
	  fn <- tempfile(fileext = ".txt")
    x <- as.matrix(InterpretMSSpectrum::apci_spectrum)
    rownames(x) <- 1:nrow(x)
    write.table(x = x, file = fn, row.names = FALSE)
    testthat::expect_equal(InterpretMSSpectrum::ReadSpecClipboard(con = fn), x)
	}
)
