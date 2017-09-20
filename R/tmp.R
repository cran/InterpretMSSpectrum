# formula_db <- "D:/Bruker/R/Rpackage_InterpretMSSpectrum/ESI.db"
# test_fml <- "C45H74O17"
# test_mz <- Rdisop::getMolecule(test_fml)$exactmass-0.0005
# dmz <- 0.002
# db_con <- DBI::dbConnect(RSQLite::SQLite(), formula_db)
# dbq <- DBI::dbSendQuery(db_con, "SELECT * FROM sfdb WHERE Mass > (:x) AND Mass < (:y);", data.frame(x=test_mz-dmz, y=test_mz+dmz))
# out <- DBI::dbFetch(dbq, -1)
# DBI::dbClearResult(dbq)
# out[out$Formula==test_fml,]
# 
# 
# test_fml <- "C12H2O2N15"
# test_mz <- Rdisop::getMolecule(test_fml)$exactmass-0.0005
# dmz <- 0.002
# # will be FALSE
# PlausibleFormula(x=test_fml, ruleset="ESI")
# # will generate a small test.db
# GenerateMetaboliteSQLiteDB(dbfile="test.db", ionization="ESI", mass_range=test_mz+c(-1,1)*2*dmz, ncores=6)
# formula_db <- "D:/Bruker/R/Rpackage_InterpretMSSpectrum/test.db"
# # repeat above statements to find it not contained in 'test.db'
# 
# 
