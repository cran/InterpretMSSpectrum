#'@export
summary.findMAIN <-
  function (object, ...)
  {
    attrname <- "scores"
    out <- data.frame(matrix(ncol = ncol(attr(object[[1]], attrname)), nrow = length(object)))
    colnames(out) <- colnames(attr(object[[1]], attrname))
    for (i in 1:length(object)) {
      out[i,] <- attr(object[[i]], attrname)
    }
    for (i in c(1, 2, 3)) out[, i] <- round(out[, i], 4)
    for (i in c(5, 6, 7, 8, 10))  out[, i] <- round(out[, i], 2)
    out[,"adducthyp"] <- sapply(object, function(x) {rownames(attr(x, "scores"))[1]})
    return(out)
  }