ExtractChemicalElements <-
function(x) {
    # returns all elements within a given sum-formula
    # all elements start with a large letter...
    p <- gregexpr("[[:upper:]]", x)[[1]]
    # split initial string at the large letter positions
    out <- sapply(1:length(p), function(i) {
        substr(x, p[i], ifelse(i == length(p), nchar(x), p[i + 1] - 1))
    })
    # remove all non letters (digits, brackets, charges...)
    out <- gsub("[^[:alpha:]]", "", out)
    # give a warning in case that elements were found repeatedly
    if (any(duplicated(out))) warning(paste("Element", paste(out[duplicated(out)], collapse="; "), "was found several times."))
    # and return
    return(unique(out))
}
