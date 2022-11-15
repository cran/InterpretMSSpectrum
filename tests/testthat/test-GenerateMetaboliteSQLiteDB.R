testthat::test_that(
  desc = "GenerateMetaboliteSQLiteDB will create a data frame for APCI", 
  code = {
    mr <- c(100, 105)
    db <- GenerateMetaboliteSQLiteDB(dbfile = NULL, ionization = "APCI", mass_range = mr, ncores = 1)
    testthat::expect_equal(nrow(db), 140)
  }
)

testthat::test_that(
  desc = "GenerateMetaboliteSQLiteDB will create a data frame for ESI", 
  code = {
    mr <- c(100.6, 101.6)
    db <- GenerateMetaboliteSQLiteDB(dbfile = NULL, ionization = "ESI", mass_range = mr, ncores = 1)
    testthat::expect_equal(nrow(db), 170)
  }
)
