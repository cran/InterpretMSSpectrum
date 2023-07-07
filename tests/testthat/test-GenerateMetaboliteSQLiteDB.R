testthat::test_that(
  desc = "GenerateMetaboliteSQLiteDB will create a data frame for APCI", 
  code = {
    mr <- c(100, 105)
    db <- GenerateMetaboliteSQLiteDB(dbfile = NULL, ionization = "APCI", mass_range = mr, ncores = 1)
    if (requireNamespace("Rdisop", quietly = TRUE)) {
      testthat::expect_equal(nrow(db), 140)
    } else {
      testthat::expect_null(db)
    }
  }
)

testthat::test_that(
  desc = "GenerateMetaboliteSQLiteDB will create a data frame for ESI", 
  code = {
    mr <- c(100.6, 101.6)
    db <- GenerateMetaboliteSQLiteDB(dbfile = NULL, ionization = "ESI", mass_range = mr, ncores = 1)
    if (requireNamespace("Rdisop", quietly = TRUE)) {
      testthat::expect_equal(nrow(db), 170)
    } else {
      testthat::expect_null(db)
    }
  }
)
