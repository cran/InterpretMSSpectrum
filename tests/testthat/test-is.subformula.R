testthat::test_that(
  desc = "is.subformula works", 
  code = {
    out <- InterpretMSSpectrum:::is.subformula(f_sub = "C4H8O5", f_main = "C6H12O6")
    testthat::expect_true(out)
    testthat::expect_equal(names(out), "C2H4O1")
    out2 <- InterpretMSSpectrum:::is.subformula(f_sub = "CH", f_main = "C2Na")
    testthat::expect_false(out2)
    testthat::expect_equal(names(out2), "C1Na1")
    out3 <- InterpretMSSpectrum:::is.subformula(f_sub = "CNa", f_main = "CH", substitutions = data.frame("s1"=c("Na1"),"s2"=c("H1")))
    testthat::expect_true(out3)
    testthat::expect_equal(names(out3), "")
  }
)
