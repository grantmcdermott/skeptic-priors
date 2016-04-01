#######################################
#######################################

### Decimal function (to make sure, e.g. three decimals places are 
### always printed in tables) ###

decimals <- function(x, k) {
  as.double(format(round(x, k), nsmall = k))
}

#######################################
#######################################

### Match short prior names to long prior names
match_priors <- function(x) {
  x <- gsub("ni", "Noninformative", x)
  x <- gsub("lukemod", "Moderate Lukewarmer", x)
  x <- gsub("lukestrong", "Strong Lukewarmer", x)
  x <- gsub("denmod", "Moderate Denier", x)
  x <- gsub("denstrong", "Strong Denier", x)
  return(x)
}

#######################################
#######################################

### Match short coeficient names to long coeficient names
match_coefs <- function(x) {
  x <- gsub("alpha", "Constant", x)
  x <- gsub("beta", "Total radiative forcing", x)
  x <- gsub("gamma", "Stratospheric aerosols", x)
  x <- gsub("delta", "SOI", x)
  x <- gsub("eta", "AMO", x)
  return(x)
}


#######################################
#######################################

### Match short RCP names to long RCP names
match_rcps <- function(x) {
  x <- gsub("rcp26", rcp_names[1], x)
  x <- gsub("rcp45", rcp_names[2], x)
  x <- gsub("rcp60", rcp_names[3], x)
  x <- gsub("rcp85", rcp_names[4], x)
  return(x)
}

#######################################
#######################################

### GGPLOT2 predictions theme
theme_coefs <-
  cowplot::theme_cowplot() +
  theme(
    # axis.line.x = element_line(linetype = 1), ## Temporary bug(?) in cowplot theme: missing axis line
    # axis.line.y = element_line(linetype = 1), ## Ditto
    text = element_text(family = font_type),
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size=14),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size=14),
    legend.position = "none",
    strip.text = element_text(size = 18, colour = "black"),
    strip.background = element_rect(fill = "white"), ## Facet strip
    panel.margin = unit(2, "lines") ## Increase gap between facet panels
    ) 

theme_pred <-
  cowplot::theme_cowplot() +
  theme(
    # axis.line.x = element_line(linetype = 1), ## Temporary bug(?) in cowplot theme: missing axis line
    # axis.line.y = element_line(linetype = 1), ## Ditto
    text = element_text(family = font_type),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20, angle = 0),
    axis.text  = element_text(size=18),
    # panel.grid.major.x = element_blank(),
    # panel.grid.major.y = element_line(colour = "grey90", size = 1),
    # legend.position = c(.16, .81),
    legend.position = c(.2, .75),
    legend.title = element_blank(), # switch off the legend title
    legend.text = element_text(size=18),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.key.width = unit(3.75, "line"),
    legend.key.height = unit(2.25, "line"),
    legend.key.size = unit(2, "line")
  ) 


#######################################
#######################################

### Decimal/character align with xtable and LaTeX
## Based on https://gist.github.com/jbryer/4458674
library(xtable) 

#' Prints a LaTeX table with numeric columns aligned on their decimal points.
#' 
#' This function wraps the \code{\link{xtable}} and \code{\link{print.xtable}}
#' functions in the \code{xtable} package so that numeric columns are aligned
#' on their decimal place.
#' 
#' See \url{http://jason.bryer.org/posts/2013-01-04/xtable_with_aligned_decimals.html}
#' for more information.
#' 
#' @author Jason Bryer <jason@@bryer.org>
#' @param x a data frame to create a LaTeX table from.
#' @param cols a numeric vector indicating which columns should be aligned on
#'        decimal points. It defaults to all columns of type numeric.
#' @param colAlignment named character vector where each element name corresponds to a
#         column name and the value is the LaTeX alignment (i.e. l, r, or c).
#' @param tocharFun the function used to convert the numeric vecotr to a character
#'        vector. This defaults to \code{\link{prettyNum}}, but other possible
#'        options are \code{\link{as.character}}, \code{\link{format}}, 
#'        \code{\link{formatC}}, or some other custom function.
#' @param ... other parameters passed to \code{tocharFun}, \code{\link{xtable}},
#'        and \code{\link{print.xtable}}.
#' @seealso xtable
#' @export
xtable.decimal <- function(x, 
                           cols=which(lapply(x, class) == 'numeric'), 
                           colAlignment, 
                           tocharFun=prettyNum,
                           ...) {
  splitCol <- function(x, ...) {
    s <- strsplit(tocharFun(unlist(x), ...), split='.', fixed=TRUE)
    right <- sapply(s, FUN=function(x) { ifelse(length(x) == 2, x[2], '0') })
    left <- sapply(s, FUN=function(x) { x[1] })
    data.frame(left=left, right=right, stringsAsFactors=FALSE)
  }
  
  cols <- cols[order(cols, decreasing=TRUE)]
  colnames <- names(x)
  for(i in cols) {
    if(i == 1) {
      tmp <- cbind(splitCol(x[,1], ...), x[,2:ncol(x)])
      names(tmp)[1:2] <- paste(names(tmp)[1], c('left','right'), sep='.')
      names(tmp)[3:ncol(x)] <- names(x)[2:ncol(x)]
      x <- tmp
    } else if(i == ncol(x)) {
      tmp <- cbind(x[,1:(ncol(x)-1)], splitCol(x[,ncol(x)], ...))
      names(tmp)[1:(ncol(tmp)-2)] <- names(x)[1:(ncol(x)-1)]
      names(tmp)[(ncol(tmp)-1):ncol(tmp)] <- paste(names(x)[ncol(x)], 
                                                   c('left','right'), sep='.')
      x <- tmp
    } else {
      tmp <- cbind(x[,1:(i-1)], splitCol(x[,i], ...), x[,(i+1):ncol(x)])
      names(tmp)[1:(i-1)] <- names(x)[1:(i-1)]
      names(tmp)[i:(i+1)] <- paste(names(x)[i], c('left','right'), sep='.')
      names(tmp)[(i+2):ncol(tmp)] <- names(x)[(i+1):ncol(x)]
      x <- tmp
    }
  }
  
  colnames[cols] <- paste('\\multicolumn{2}{c}{', colnames[cols], '}', sep='')
  colnames <- paste(colnames, collapse=' & ')
  
  addtorow <- list()
  addtorow$pos <- list()
  addtorow$pos[[1]] <- c(0)
  addtorow$command <- paste( colnames, ' \\\\ ', sep='')
  
  align <- rep('l', ncol(x))
  if(!missing(colAlignment)) {
    for(i in seq_along(colAlignment)) {
      align[names(x) == names(colAlignment)[i]] <- colAlignment[i]
    }
  }
  align[grep('.left$', names(x), perl=TRUE)] <- 'r@{.}'
  align <- c('l', align) #Add an alignment for row names
  
  xtab <- xtable(x, align=align, ...)
  args <- list(...)
  args <- args[names(args) != 'caption']
  do.call(print, 
          c(list(xtab, add.to.row=addtorow, include.rownames=FALSE, include.colnames=FALSE), 
            args)
          )
}
