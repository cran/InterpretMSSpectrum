#'@export
print.findMAIN <- function (x, ...)
  {
    nres <- length(x)
    i <- 1
    scores <- summary(x)
    scores_i <- attr(x[[i]], "scores")
    nprec <- length(unique(sapply(x, function(x) attr(x, "scores")[,"adductmz"])))
    nadduct <- length(unique(sapply(x, function(x) attr(x, "scores")[,"adducthyp"])))
    message(
      sprintf(
        "Analyzed %d neutral mass hypotheses (%d peaks * %d adducts), kept %d",
        nprec * nadduct,
        nprec,
        nadduct,
        nres
      )
    )
    message(
      sprintf(
        "Selected m/z=%.4f as %s adduct of neutral mass %.4f with score %.2f.",
        scores_i[, 1],
        rownames(scores_i)[1],
        scores_i[, 3],
        scores_i[, 10]
      )
    )
    print(x[[i]])
  }