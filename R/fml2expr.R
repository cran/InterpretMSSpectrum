#' @title Turn text with chemical formulas into expression vector.
#'
#' @description It is often helpful to annotate a Figure with chemical formulas.
#'     However, to increase readability of chemical formulas, certain conventions
#'     have to be met. These concern, among others, superscript and subscript 
#'     text which, in R, can only be provided via expressions. `fml2expr` will
#'     convert any character vector into an expression vector, aiming to identify
#'     and format potential contained chemical formulas.
#' 
#' @param x Character vector containing chemical formulas.
#' @param expr Will return vector of characters when FALSE and expressions when TRUE.
#'
#' @returns A vector of expressions.
#' @export
#'
#' @examples
#' texts <- c(
#'   "C6H12Cl6",         # simple formula
#'   "H3Na",             # simple formula
#'   "Text with blank",  # simple text
#'   "A1B2",             # No valid elements but looks like formula
#'   "Fe(NO3)3",         # group with index
#'   "SO4^2-",           # charged molecule
#'   "Na+", "Cl-",       # simple ions
#'   paste0("H", intToUtf8(0x2082), "O"), # unicode-subscript
#'   "Ca(OH)2-"          # group with index and charge
#' )
#' 
#' # Check that all examples can be converted to expression
#' all(sapply(texts , function(x) { !inherits(try(fml2expr(x), silent=TRUE), "try-error") }))
#' exprs <- fml2expr(texts)
#' str(exprs)
#' 
#' # Plot with legend
#' plot(1:10, 1:10, pch = 16, col = 1:10, xlim = c(0, 11), ylim = c(0, 11))
#' legend("topright", legend = exprs, col = 1:10, pch = 16, title = "Formeln")
#' 
#' # Careful! text() does not accept an expression vector
#' text(3:4,3:4,labels=exprs[3:4])
#' for (i in 5:6) text(i,i,labels=exprs[[i]])
#' 
#' # you can also return a named character vector (names are input, values are modified)
#' fml2expr(texts, expr = FALSE)

fml2expr <- function(x, expr = TRUE) {
  
  replace_unicode <- function(text) {
    # Tiefgestellte Ziffern 0-9 und Zeichen +, -
    subscripts <- intToUtf8(0x2080:0x2089, multiple = TRUE)
    subscripts_extra <- intToUtf8(c(0x208A, 0x208B), multiple = TRUE)
    subscripts_replacements <- c(as.character(0:9), "+", "-")
    
    # Hochgestellte Ziffern 0-9 und Zeichen +, -
    superscripts <- c(intToUtf8(0x2070), intToUtf8(0x00B9), intToUtf8(0x00B2), intToUtf8(0x00B3), intToUtf8(0x2074:0x2079, multiple = TRUE))
    superscripts_extra <- intToUtf8(c(0x207A, 0x207B), multiple = TRUE) 
    superscripts_replacements <- c(as.character(0:9), "+", "-")
    
    # Kombiniere alle zu ersetzenden Zeichen
    all_unicode <- c(subscripts, subscripts_extra, superscripts, superscripts_extra)
    all_replacements <- c(subscripts_replacements, superscripts_replacements)
    
    # Ersetze alle Unicode-Zeichen durch normale Zeichen
    for (i in seq_along(all_unicode)) {
      text <- gsub(all_unicode[i], all_replacements[i], text, fixed = TRUE)
    }
    
    return(text)
  }
  
  process_element <- function(text) {
    # 1. replace unicode characters
    text <- replace_unicode(text)
    
    # 2. replace blanks
    text <- gsub(" ", "~", text)
    
    # 3. replace groups with brackets and index: (OH)2 -> (OH)[2]
    text <- gsub("\\(([A-Za-z0-9]+)\\)(\\d{1,3})", "(\\1)[\\2]", text, perl = TRUE)
    
    # 4. replace standard elements: LETTERS{1,1}letters{0,1}:digit:{1,3}
    text <- gsub("([A-Z][a-z]?)(\\d{1,3})", "\\1[\\2]", text, perl = TRUE)
    
    # 5. replace charges:
    # a) oxidation ^2-, ^3+: ^'2-', ^'3+'
    text <- gsub("\\^(\\d+[+-])", "^'\\1'", text, perl = TRUE)
    
    # b) simple + or - at end:  ^'+', ^'-'
    text <- gsub("([A-Za-z0-9\\]\\)])([+-])$", "\\1^'\\2'", text, perl = TRUE)
    
    # c) simple + or - inside: ^'+', ^'-'
    #text <- gsub("([A-Za-z0-9\\]\\)])([+-])(?![0-9])", "\\1^'\\2'", text, perl = TRUE)
    
    # 6. add ~ between elements:
    text <- gsub("\\](?=\\()", "]~", text, perl = TRUE)
    text <- gsub("\\](?=[A-Z])", "]~", text, perl = TRUE)
    
    # 7. remove ~ from end if present
    text <- sub("~$", "", text)
    
    # return value
    if (expr) {
      parse(text = text)[[1]]
    } else {
      text
    }
  }
  
  # return as list when length of input >=2
  n <- length(x)
  if (expr) {
    if (n==1) process_element(x) else lapply(x, process_element)
  } else {
    if (n==1) process_element(x) else vapply(x, process_element, FUN.VALUE = character(1))
  }
}
