#'@export
print.findMAIN <- function (x, ...)
  {
    nres <- length(x)
    i = 1
    scores <- summary(x)
    scores1 <- attr(x[[i]], "scores")
    nprec <- length(unique(sapply(x, function(x)
      attr(x, "scores")[,
                        "adductmz"])))
    nadduct <-
      length(unique(sapply(x, function(x)
        attr(x, "scores")[,
                          "adducthyp"])))
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
        "Selected %.4f as [M%+.2f], neutral %.4f, score %.2f",
        scores1[, 1],
        scores1[, 2],
        scores1[, 3],
        scores1[, 10]
      )
    )
    print(x[[i]])
  }