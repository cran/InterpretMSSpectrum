#'@title ScoreFormulaCombination.
#'
#'@description
#'\code{ScoreFormulaCombination} takes a rdisop formula seach result as a list. This list will contain formula suggestions for certain masses (including Score information) and masses should be in ascending order.
#'
#'@details
#'Not exported.
#'
#'@param rdisop_res Internal result structure of InterpretMSSpectrum.
#'@param nl_vec Named vector of neutral losses or NULL.
#'@param punish_invalid Lower formulas tagged invalid by Rdisop by factor 1-punish_invalid, e.g. punish_invalid=0.5 would lower a score of 80 to 40.
#'@param punish_S Lower formulas containing Sufor by factor 1-punish_S*n_S, e.g. punish_S=0.2 would lower a score of 100 to 60 if 2 Sulfor are contained.
#'@param punish_nonplausible Check for all potential fragments not only if they are a subformula but if this neutral loss is a plausible formula itself and lower according to
#'@param return_rank Integer, will return the n-th best combination; if NA will return a ranked list of all found combinations.
#'@param neutral_loss_cutoff Cutoff in mDa for accepting an internal mass difference as a given neutral loss.
#'@param silent Print some stats or not.

#'@return
#'Output will be the most likely combination of fragments by evaluating the largest mean score combination of all fragments.
#'
#'@import plyr
#'@import enviPat
#'
#'@keywords internal
#'
ScoreFormulaCombination <- function(rdisop_res, nl_vec=c("C1H4"=16.0313), punish_invalid=0.5, punish_S=0.2, punish_nonplausible=0.2, return_rank=1, neutral_loss_cutoff=0.5, silent=TRUE) {

  # weight formulas tagged "invalid" by 0.5
  if (is.numeric(punish_invalid) && abs(punish_invalid-0.5)<=0.5) {
    rdisop_res <- lapply(rdisop_res, function(x) {
      x[,"Score"] <- x[,"Score"]*c(1-punish_invalid,1)[as.numeric(x[,"Valid"])]
      return(x)
    })
  }
  
  # weight formulas containing Sulfor
  if (is.numeric(punish_S) && abs(punish_S-0.5)<=0.5) {
    rdisop_res <- lapply(rdisop_res, function(x) {
      nS <- sapply(x[,"Formula"], CountChemicalElements, ele="S")
      x[,"Score"] <- x[,"Score"]*(1-punish_S*nS)
      return(x)
    })
  }
  
  # start checking at top-mass levels
  top_formulas <- rdisop_res[[length(rdisop_res)]]

  # evaluate mean score of potential fragment group
  top_formulas_list <- lapply(1:nrow(top_formulas), function(k) {
    plyr::ldply(1:length(rdisop_res), function(m) {
      #print(paste(k,m))
      test_subformula <- sapply(rdisop_res[[m]][,"Formula"], function(y) { is.subformula(f_sub=y, f_main=top_formulas[k,"Formula"]) })
      # enviPat alternative is slower :(
      #test_subformula <- sapply(rdisop_res[[m]][,"Formula"], function(y) { enviPat::check_ded(top_formulas[k,"Formula"], y) })=="FALSE"
      if (sum(test_subformula)>=1 && m!=length(rdisop_res)) {
        #if (sum(test_subformula)>=2) print(paste(top_formulas[k,1], rdisop_res[[m]][test_subformula,"Formula"]))
        # exclude subformulas which do not fit a neutral loss if detectable
        test_nl <- matrix(abs(top_formulas[k,"Mass"]-rep(rdisop_res[[m]][test_subformula,"Mass"],each=length(nl_vec))-nl_vec) <= (neutral_loss_cutoff/1000),nrow=length(nl_vec), ncol=sum(test_subformula), dimnames=list(names(nl_vec),names(which(test_subformula))))
        if (!is.null(nl_vec) && any(test_nl)) {
          for (n in which(test_subformula)) {
            this_loss <- which(test_nl[,names(test_subformula)[n]])
            if (length(this_loss)!=1) {
              test_subformula[n] <- FALSE
            } else {
              if (enviPat::subform(top_formulas[k,"Formula"], names(nl_vec)[this_loss]) != rdisop_res[[m]][n,"Formula"]) test_subformula[n] <- FALSE
            }
          }
        }
      }
      if (any(test_subformula)) {
        sf <- which(test_subformula)
        if (is.numeric(punish_nonplausible) && abs(punish_nonplausible-0.5)<=0.5 && m!=length(rdisop_res)) {
          for (n in sf) {
            rdisop_res[[m]][n,"Score"] <- round(ifelse(PlausibleFormula(x=enviPat::subform(top_formulas[k,"Formula"], rdisop_res[[m]][n,"Formula"])), 1, 1-punish_nonplausible)*rdisop_res[[m]][n,"Score"],1)
          }
        }
        return(rdisop_res[[m]][sf[which.max(rdisop_res[[m]][sf,"Score"])],])
      } else {
        dummy <- rdisop_res[[m]][1,,drop=F]
        dummy[1,] <- NA
        return(dummy)
      }
    })
  })

  # get mean score for top formulas
  top_formulas_mean_score <- sapply(top_formulas_list, function(x) { mean(x[,"Score"],na.rm=T) })

  # print results overview
  if (!silent) {
    ncol_print <- min(c(8, nrow(rdisop_res[[length(rdisop_res)]])))
    tmp <- matrix(NA, nrow=length(rdisop_res), ncol=ncol_print, dimnames=list(round(sapply(rdisop_res,attr,"M0"),4), round(sort(top_formulas_mean_score,decreasing=TRUE),2)[1:ncol_print]))
    tf_rnk <- order(top_formulas_mean_score,decreasing=TRUE)
    for (i in 1:ncol(tmp)) tmp[,i] <- lapply(top_formulas_list, function(x) { x[,"Formula"] })[[tf_rnk[i]]]
    cat(paste("\n\nTop",ifelse(ncol_print<nrow(rdisop_res[[length(rdisop_res)]]), paste(ncol_print,"of",nrow(rdisop_res[[length(rdisop_res)]])), ncol_print),"remaining Formula combinations (ordered after mean score of combination)...\n"))
    print(tmp)
  }

  # return ordered result list
  if (length(return_rank)==1 & (is.na(return_rank) || return_rank<=length(top_formulas_mean_score))) {
    if (is.na(return_rank)) {
      return(top_formulas_list[order(top_formulas_mean_score, decreasing=TRUE)])
    } else {
      k <- order(top_formulas_mean_score, decreasing=TRUE)[return_rank]
      out <- top_formulas_list[[k]]
      return(out[!is.na(out[,1]),])
    }
  }

}