#'@title PlausibleFormula.
#'
#'@description
#'\code{PlausibleFormula} will count the chemical elements within a provided formula and apply a rule set to evaluate if this formula is plausible.
#'It is currently exported to allow the user to modify/substitute the result but may become an internal function in the future.
#'
#'@details
#'The rules are developed based on comparison with ~1000 compounds found in the Golm Metabolome Database (GMD) for silylated metabolites detectable in GC-MS.
#'The following elements will be counted and tested: c("C","H","N","O","P","S","Si").
#'
#'@param x A chemical formula. No validity checks are performed.
#'@param ruleset Currently either "GCMS" or "LCMS" as a keyword. In future potentially a link to a texfile containing the ruleset.
#'
#'@return
#'Either TRUE or FALSE. Use sapply(vector, PlausibleFormula) if necessary.
#'
#'@examples
#'PlausibleFormula(x="C6H12O6")
#'PlausibleFormula(x="C5H13O14P3")
#'
#'@export
#'
#'@keywords internal

PlausibleFormula <- function(x, ruleset=c("GCMS","LCMS")[1]) {
  if (ruleset=="GCMS") {
    if (length(x)==1 && is.character(x)) {
      out_name <- x
      # check if a certain SumFormula is likely to appear in GCMS-APCI data
      x <- CountChemicalElements(x, ele=c("C","H","N","O","P","S","Si"))
      out <- TRUE
      # test maximal number of elemental occurences (empirical cutoff evaluation based on GMD entries)
      ulim <- list("C"=60, "H"=100, "N"=10, "O"=27, "P"=4, "S"=3, "Si"=12)
      out <- out & all(sapply(names(ulim), function(y) { x[y]<=ulim[[y]] }), na.rm=T)
      # test minimal number of elemental occurences (empirical cutoff evaluation based on GMD entries)
      llim <- list("C"=1, "H"=4)
      out <- out & all(sapply(names(llim), function(y) { x[y]>=llim[[y]] }), na.rm=T)
      # test individual element ratios
      out <- out & x["N"]<=x["C"]
      out <- out & x["O"]<=x["C"]
      out <- out & 6*x["P"]<=x["C"]
      out <- out & 4*x["S"]<=x["C"]
      out <- out & sum(x[c("N","O","P","S")])<=x["C"]
      # test Si to O+N+S ratio (assumption: O, N and S Groups will most likely be Silylated, which may not be always True)
      # allow aditionally non Silylated  Compounds as some of those exist in GMD
      #out <- out & ifelse(x["Si"]==0, TRUE, !(6*x["Si"] < x["O"]+x["N"]+x["S"]+x["P"]))
      out <- out & ifelse(x["Si"]==0, TRUE, !(8*x["Si"] < x["O"]+x["N"]+x["S"]+x["P"]))
      out <- out & !(x["Si"] > 2*(x["O"]+x["N"]+x["S"]+x["P"]))
      #out <- out & !(x["Si"] >= x["O"]+x["N"]+x["S"])
      # test Si to C ratio (assumption: Si is still present as SiC3H8 Group, which may not be always True)
      # lowest value allowed for C is 3*Si-2 to account for e.g. CH4 Loss of Phosphate-3TMS or TMS-Dimer (mz=147, F=C5H15OSi2)
      out <- out & (3*x["Si"]-2 <= x["C"])
      # test Si to H ratio (assumption: Si is still present as SiC3H8 Group, which may not be always True)
      out <- out & (8*x["Si"]-3 <= x["H"])
      # S and P are rare to occur together (though they sometimes may)
      out <- out & x["S"]*x["P"]==0
      # if S=2 than Si>2 (...)
      out <- out & ifelse(x["S"]==2, x["S"]<x["Si"], TRUE)
      # if P>1 than Si>2P (...)
      out <- out & ifelse(x["P"]>0, 2*x["P"]<=x["Si"], TRUE)
      # if P>1 than O>=4 (...)
      out <- out & ifelse(x["P"]>0, x["O"]>=4, TRUE)
      names(out) <- out_name
    } else {
      out <- NA
      names(out) <- paste(x, collapse="_", sep="_")
    }
  }
  if (ruleset=="LCMS") {
    if (length(x)==1 && is.character(x)) {
      out_name <- x
      # check if a certain SumFormula is likely to appear in LCMS data
      x <- CountChemicalElements(x, ele=c("C","H","N","O","P","S"))
      m <- sum(x*c(12.0, 1.007825, 14.003074, 15.994915, 30.973762, 31.972071))
      out <- TRUE
      # test maximal number of elemental occurences (empirical cutoff evaluation based on GMD entries)
      ulim <- list("C"=2+0.07*m, "H"=6+0.14*m, "N"=5+0.01*m, "O"=3+0.03*m, "P"=min(c(8, 1+0.01*m)), "S"=4)
      out <- out & all(sapply(names(ulim), function(y) { x[y]<=ulim[[y]] }), na.rm=T)
      llim <- list("C"=-15+0.04*m, "H"=-10+0.03*m, "N"=0, "O"=-25+0.03*m, "P"=0, "S"=0)
      #out <- out & all(sapply(names(llim), function(y) { x[y]>=llim[[y]] }), na.rm=T)
      # test individual element ratios
      out <- out & x["N"]<=x["C"]
      out <- out & (x["O"]-4*x["P"])<=(1+x["C"])
      out <- out & x["P"]<=x["C"]
      out <- out & x["S"]<=x["C"]
      out <- out & (x["N"]+x["P"]+x["O"]-4*x["P"]+x["S"])<=(2+x["C"])
      # S and P are rare to occur together (though they sometimes may)
      out <- out & x["S"]*x["P"]==0
      # if P>1 than O>=4 (...)
      out <- out & ifelse(x["P"]>0, x["O"]>=4, TRUE)
      names(out) <- out_name
    } else {
      out <- NA
      names(out) <- paste(x, collapse="_", sep="_")
    }
  }
  return(out)
}
