#' @title Render SMILES into 2D image for plotting via rcdk.
#' @description This function uses the rcdk to parse the smiles into a mol, with
#'     options to switch kekulise (aromaticity detection) on or off and to define
#'     desired coordinates. Output requires that plot.new has been called, i.e.
#'     this is designed to be used directly during plotting.
#'     This function uses default depiction options.
#' @param smiles A valid SMILES code for rendering (e.g. \code{"c1ccccc1"}).
#' @param kekulise If \code{TRUE}, performs CDK aromaticiy detection, which is
#'     recommended. Setting \code{FALSE} can force rendering of invalid SMILES
#'     with undefined bond orders. Older rcdk versions may behave differently.
#' @param coords This is used to control the size of the image within the plot.
#'     Values \code{c(xmin,xmax,ymin,ymax)} in user coordinates.
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#' @details More information about aromaticity: \url{https://github.com/CDK-R/cdkr/issues/49}
#' @examples
#' smiles <- "OS(=O)(=O)c1ccc(cc1)C(CC(=O)O)CC(=O)O"
#' plot.new()
#' plot.window(xlim=c(0,200), ylim=c(0,100))
#' renderSMILES(smiles,kekulise=FALSE)
#' renderSMILES(smiles,kekulise=TRUE)
#' renderSMILES(smiles, coords = c(100,150,0,50))
#' rect(100,0,150,50)
#' @return Returns an image for use during plotting
#' @noRd
#' @keywords internal
renderSMILES <- function(smiles, kekulise=TRUE, coords=c(0,100,0,100), gp_usr = NULL) {
  if (is.null(gp_usr)) gp_usr <- graphics::par("usr")
  if (nchar(smiles)>1) {
    mol <- rcdk::parse.smiles(smiles,kekulise=kekulise)[[1]]
    mol <- rcdk::generate.2d.coordinates(mol)
    #depictor <- rcdk::get.depictor(width = 200, height = 200, zoom = 1.0, fillToFit = FALSE)
    depictor <- rcdk::get.depictor(width = 300, height = 300, zoom = 1.0, fillToFit = FALSE)
    img <- tryCatch({
      (rcdk::view.image.2d(mol, depictor = depictor))
    }, error = function(e) {
      img <- ""
      message(paste("Invalid SMILES not rendered: ", smiles, sep=""))
    })
    if (length(img)<=2 && kekulise) {
      message("Replotting with kekulise=FALSE")
      mol <- rcdk::parse.smiles(smiles,kekulise=FALSE)[[1]]
      mol <- rcdk::parse.smiles(smiles,kekulise=kekulise)[[1]]
      mol <- rcdk::generate.2d.coordinates(mol)
      x_len <- diff(coords[1:2])
      y_len <- diff(coords[3:4])
      depictor <- rcdk::get.depictor(width = x_len, height = y_len, zoom = 1.0, fillToFit = FALSE)
      img <- tryCatch({
        (rcdk::view.image.2d(mol, depictor = depictor))
      }, error = function(e) {
        img <- ""
        message(paste("Invalid SMILES not plotted without kekulise either: ", smiles, sep=""))
      })
    }
    if (length(img)>2) {
      
      # make points of white color transparent
      img[,,4] <- matrix(abs(as.numeric(img[,,1]==1 & img[,,2]==1 & img[,,3]==1)-1), nrow = nrow(img[,,1]))
      
      xx <- mean(coords[c(1,2)])
      x_len <- diff(coords[1:2])
      yy <- mean(coords[c(3,4)])
      y_len <- diff(coords[3:4])
      
      # remove empty rows and apply y-shift
      empty_rows <- apply(img[,,4], 1, function (x) { all(x==0) })
      strip_height_percent <- 1-(table(diff(which(empty_rows)))[1])/nrow(img[,,4])
      y_len <- y_len * strip_height_percent
      # that would be strict with respect to upper/lower half
      # but IMSS is determining sum formula annotation according to 0.9*max(y)
      # hence this is adjusted here to 0.75*max
      peak_in_upper_half <- yy > mean(gp_usr[c(3:4)])
      
      ## target center adjustment based on location:
      ## y-location of subplot
      lh <- diff(gp_usr[3:4])/30
      #browser()
      dont_account_for_lh <- peak_in_upper_half && yy < (gp_usr[4] - 0.1 * diff(gp_usr[3:4]))
      # get adjusted bottom border to account for stripping
      yy <- yy + ifelse(peak_in_upper_half, -1, 1) * ifelse(dont_account_for_lh, 0.5, 1.5) * lh + ifelse(peak_in_upper_half, -1, 0) * y_len
      coords[3] <- yy
      coords[4] <- yy + y_len
      stripped_img <- img[!empty_rows,,]
      
      # remove empty cols and apply x-shift
      empty_cols <- apply(stripped_img[, , 4], 2, function(x) { all(x == 0) })
      if (!all(empty_cols)) {
        strip_width_percent <- sum(!empty_cols) / length(empty_cols)
        x_len <- diff(coords[1:2]) * strip_width_percent
        coords[1] <- xx - x_len/2
        coords[2] <- xx + x_len/2
        if (coords[1] < gp_usr[1]) {
          shft <- gp_usr[1] - coords[1]
          coords[1] <- coords[1] + shft
          coords[2] <- coords[2] + shft
        }
        if (coords[2] > gp_usr[2]) {
          shft <- gp_usr[2] - coords[2]
          coords[1] <- coords[1] + shft
          coords[2] <- coords[2] + shft
        }
        stripped_img <- stripped_img[, !empty_cols, , drop = FALSE]
      }
      
      # graphics::rect(coords[1], coords[3], coords[2], coords[4])
      graphics::rasterImage(stripped_img, coords[1], coords[3], coords[2], coords[4])
    }
  }
  invisible(img)
}

#' @title Determine square_subplot_coord.
#' @description Compute square subplot coordinates in user space
#' @param x x center coordinate of square.
#' @param y y center coordinate of square.
#' @param gp_usr numeric(4); c(xlim, ylim) of spec.
#' @param w Relative proportion of figure region to cover.
#' @return Returns coordinates for a square shape sub-plot.
#' @examples
#' square_subplot_coord(x=50, y=50, gp_usr=c(0,100,0,100))
#'
#' @noRd
#' @keywords internal
#'
square_subplot_coord <- function(x, y, gp_usr = NULL, w = 0.2) {
  
  if (is.null(gp_usr)) gp_usr <- graphics::par("usr")
  xlim <- gp_usr[1:2]
  ylim <- gp_usr[3:4]
  
  ## device size in inch
  dev_s <- grDevices::dev.size("in")
  if (!all(is.finite(dev_s))) dev_s <- c(1,1)
  
  ## subplot size in inch
  subplot_size_in <- w * min(dev_s)
  
  ## subplot size in user coords
  x_len <- subplot_size_in * diff(xlim) / dev_s[1]
  y_len <- subplot_size_in * diff(ylim) / dev_s[2]
  
  ## final coords
  x_start <- x - x_len/2
  x_end   <- x + x_len/2
  y_start <- y - y_len/2
  y_end   <- y + y_len/2
  
  c(x_start, x_end, y_start, y_end)
}
