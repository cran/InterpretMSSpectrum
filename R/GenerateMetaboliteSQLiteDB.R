#'@title GenerateMetaboliteSQLiteDB.
#'
#'@description
#'\code{GenerateMetaboliteSQLiteDB} will set up a SQLite data base containing potential metabolite formulas, their masses and isotopic distribution for use with \link{InterpretMSSpectrum}.
#'
#'@details
#'The process takes a long time for larger masses (>400 Da). Parallel processing with 8 cores is highly recommended. Alternatively pre-processed versions can be downloaded on request to \email{jan.lisec@charite.de}. To process a 1 Da range (from 900 to 901) for ESI does take ~5minutes on 8 cores.
#'
#'@param dbfile path and filename of the final SQLiteDB.
#'@param ionization Has to be specified to account for different plausibility rules and elemental composition.
#'@param mass_range For testing use default range, otherwise use your measurement range.
#'@param ncores Number of cores. Use as many as possible.
#'
#'@return
#'Returns NULL invisible. Will write an SQL_DB as specified by 'dbfile'.
#'
#'@examples
#'#this would be relatively fast, but for higher masses it is getting much slower
#'\donttest{
#'GenerateMetaboliteSQLiteDB(dbfile="APCI.db", ionization="APCI", mass_range=c(80,140), ncores=1)
#'GenerateMetaboliteSQLiteDB(dbfile="ESI.db", ionization="ESI", mass_range=c(900,901), ncores=8)
#'}
#'
#'@export
#'
GenerateMetaboliteSQLiteDB <- function(dbfile="SQLite_APCI.db", ionization=c("APCI","ESI")[1], mass_range=c(80,160), ncores=2) {
  if (ionization=="ESI") {
    allowed_elements <- c("C","H","N","O","P","S","Cl")
    maxElements="P4S4Cl2"
  }
  if (ionization=="APCI") {
    allowed_elements <- c("C","H","N","O","P","S","Si")
    maxElements="P2S2"
  }
  em <- 0.00055
  elements <- Rdisop::initializeElements(allowed_elements)
  mmin <- mass_range[1]
  mmax <- mass_range[2]
  dmz <- 0.001*2^ifelse(mmin<100,14,ifelse(mmin<300,10,ifelse(mmin<1000,1,ifelse(mmin<1400,0))))

  # process small fractions of mass range, either
  # PARALLEL
  if (ncores>1) {
    doParallel::registerDoParallel(parallel::makeCluster(ncores))
    dmz <- 0.001
    isotopes <- matrix(allowed_elements,ncol=1)
    empty_result <- data.frame("Formula"=I(character(0)), "Valid"=character(0), "Mass"=numeric(0),
                               "m0"=numeric(0), "m1"=numeric(0), "m2"=numeric(0), "m3"=numeric(0), 
                               "a0"=numeric(0), "a1"=numeric(0), "a2"=numeric(0), "a3"=numeric(0))
    suppressWarnings(
    res <- plyr::ldply(seq(mmin+dmz, mmax-dmz, by=2*0.001), function(mz) {
      molecules <- Rdisop::decomposeMass(mass=mz, mzabs=dmz, ppm=0, z=1, maxisotopes=4, elements=elements, minElements="C1", maxElements=maxElements)
      if (!is.null(molecules)) {
        out <- data.frame("Formula"=I(molecules$formula), "Valid"=molecules$valid, "Mass"=molecules$exactmass)
        out <- cbind(out, plyr::ldply(molecules$isotopes,function(x){matrix(c(x[1,],x[2,]),nrow=1,dimnames=list(NULL,paste0(rep(c("m","a"),each=4),0:3)))}))
        out <- out[sapply(out[,"Formula"], InterpretMSSpectrum::PlausibleFormula, ruleset=ionization),,drop=FALSE]
        out[,"Formula"] <- sapply(out[,"Formula"], function(fml) { enviPat::check_chemform(isotopes=isotopes, fml)[,"new_formula"] })
      } else {
        out <- empty_result
      }
      return(out)
    }, .parallel = TRUE)
    )
  # SERIAL
  } else {
    res <- NULL
    isotopes <- matrix(allowed_elements,ncol=1)
    os <- utils::object.size(res)
    mtmp <- mmin+dmz
    while (mmin<mmax && os<4*1024^3) {
      molecules <- Rdisop::decomposeMass(mass=mtmp, mzabs=dmz, ppm=0, z=1, maxisotopes=4, elements=elements, minElements="C1", maxElements=maxElements)
      cat(paste("\nmz =", mtmp, "; mzabs =", dmz, "; nmolecules =", length(molecules[[1]]), "; ntotal =", nrow(res), "; ", format(utils::object.size(res), units="Gb")))
      mmin <- mtmp+dmz
      if (!is.null(molecules)) {
        out <- data.frame("Formula"=I(molecules$formula), "Valid"=molecules$valid, "Mass"=molecules$exactmass)
        out <- cbind(out, plyr::ldply(molecules$isotopes,function(x){matrix(c(x[1,],x[2,]),nrow=1,dimnames=list(NULL,paste0(rep(c("m","a"),each=4),0:3)))}))
        out <- out[sapply(out[,"Formula"], InterpretMSSpectrum::PlausibleFormula, ruleset=ionization),,drop=FALSE]
        out[,"Formula"] <- sapply(out[,"Formula"], function(fml) { enviPat::check_chemform(isotopes=isotopes, fml)[,"new_formula"] })
        if (nrow(out)>=1)  res <- rbind(res,out)
        if (length(molecules[[1]])>100000) dmz <- 0.5*dmz
      }
      mtmp <- mmin+dmz
      os <- utils::object.size(res)
    }
  }
  #save(res, file="res.RData")

# remove redundancies
  if (any(duplicated(res[,1]))) {
    cat("\nRemoving redundancies...")
    res <- res[!duplicated(res[,1]),]
  }

# stop cluster
  if (ncores>1) {
    cat("\nStopping cluster and quit...")
    doParallel::stopImplicitCluster()
  }
  
# make SQLite-DB out of res
  cat("\nWriting DB File...\n\n")
  con <- DBI::dbConnect(RSQLite::SQLite(), dbfile)
  DBI::dbWriteTable(conn=con, name="sfdb", res, append=F, overwrite=T)
  db <- DBI::dbSendQuery(con, "CREATE INDEX idx ON sfdb (Mass)")
  DBI::dbClearResult(db)
  DBI::dbDisconnect(con)
  
  invisible(NULL)

}
