#' @title RestrictByTypicalLosses.
#'
#' @description \code{RestrictByTypicalLosses} will remove all formulas from a group which are inconsistent with an observed typical loss.
#'
#' @details Not exported.
#'
#' @param rdisop_res Internal result structure of InterpretMSSpectrum.
#' @param tl Typical loss, a mass with the according sum formula as name.
#' @param neutral_loss_cutoff Cutoff in mDa for accepting an internal mass difference as a given neutral loss.
#' @param punish If NULL remove fragments, if numeric 0..1 don't remove but punish the score.
#' @param substitutions Potential substitutions like H/Na in ESI mode, provided as two column matrix.
#'
#' @return Modified 'rdisop_res' object.
#'
#' @keywords internal
#' @noRd
#'
RestrictByTypicalLosses <- function(rdisop_res=NULL, tl=NULL, neutral_loss_cutoff=0.5, punish=0.1, substitutions=NULL) {
  for (i in 1:(length(rdisop_res)-1)) {
    for (j in (i+1):length(rdisop_res)) {
      # is tl found between these fragment groups?
      im <- rdisop_res[[i]][,"Mass"]
      ifm <- rdisop_res[[i]][,"Formula"]
      jm <- rdisop_res[[j]][,"Mass"]
      jfm <- rdisop_res[[j]][,"Formula"]
      ind <- data.frame("i"=rep(1:nrow(rdisop_res[[i]]),each=nrow(rdisop_res[[j]])), 
                        "j"=rep(1:nrow(rdisop_res[[j]]),times=nrow(rdisop_res[[i]])), 
                        "k"=abs(rep(jm, times=length(im))-rep(im, each=length(jm))-tl)<(neutral_loss_cutoff/1000))
      # test for correct sumformula combinations
      if (any(ind[,"k"])) {
        ind[,"f"] <- F
        for (k in which(ind[,"k"])) {
          # Rdisop version causes memory overflow and was substituted by enviPat version
          #if (Rdisop::subMolecules(jfm[ind[k,"j"]], names(tl))$formula == ifm[ind[k,"i"]]) ind[k,"f"] <- TRUE
          #if (enviPat::subform(jfm[ind[k,"j"]], names(tl)) == ifm[ind[k,"i"]]) ind[k,"f"] <- TRUE
          ind[k,"f"] <- is.subformula(f_sub=ifm[ind[k,"i"]], f_main=jfm[ind[k,"j"]], substitutions=substitutions)
        }
        # filter for correct sumformula combinations if present
        if (any(ind[,"f"])) {
          #ind[ind[,"k"]&ind[,"f"],]
          #browser()
          if (is.null(punish)) {
            rdisop_res[[i]] <- rdisop_res[[i]][ind[ind[,"f"],"i"],]
            rdisop_res[[j]] <- rdisop_res[[j]][ind[ind[,"f"],"j"],]
          } else {
            rdisop_res[[i]][-ind[ind[,"f"],"i"],"Score"] <- punish*rdisop_res[[i]][-ind[ind[,"f"],"i"],"Score"]
            rdisop_res[[j]][-ind[ind[,"f"],"j"],"Score"] <- punish*rdisop_res[[j]][-ind[ind[,"f"],"j"],"Score"]
          }
        }
      }
    }
  }
  return(rdisop_res)
}