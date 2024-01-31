#' @title ReadSpecClipboard.
#' @description Read a mass spectrum from the windows clipboard.
#' @param con A connection other than 'clipboard' can be provided.
#' @return A spectrum as two-column matrix.
#' @export
#' @examples
#' \dontrun{
#'   if (length(grep("Windows", utils::sessionInfo()$running))==1) {
#'     x <- InterpretMSSpectrum::apci_spectrum
#'     write.table(x, "clipboard", sep="\t", row.names=FALSE)
#'     InterpretMSSpectrum::ReadSpecClipboard()
#'   }
#' }
ReadSpecClipboard <- function(con = "clipboard") {
  if (length(grep("Windows", utils::sessionInfo()$running))!=1 && con == "clipboard") return(NULL)
  # source could be Excel (German/English) or DA directly
  spec <- readLines(con = con)
  spec <- gsub("\t", " ", spec) # replace Tabs
  if (length(grep("[^[:digit:],/. ]", spec[1])) == 1) spec <- spec[-1] # strip header if present
  spec <- gsub(",", ".", spec) # replace Colons
  spec <- gsub(" +$", "", spec) # trim white space end
  spec <- gsub("^ +", "", spec) # trim white space start
  # convert to numeric matrix
  spec <- as.matrix(
    t(sapply(spec, function(x) {
      as.numeric(strsplit(x, " ")[[1]])
    }))
  )
  if (ncol(spec) >= 3) spec <- spec[, 2:3]
  rownames(spec) <- 1:nrow(spec)
  colnames(spec) <- c("mz", "int")
  return(spec)
}
