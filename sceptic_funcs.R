## Load packages ##
library(readr) ## For reading in data files
library(LearnBayes) ## Mostly for simulating noninformative prior (using random multivarite normal command)
# library(arm)
library(rjags) ## For running the MCMC (Gibbs) sampler
# library(coda) ## For converting MCMC objects and diagnostics. Loads as rjags dependency.
library(dclone) ## Allows parallel updating of JAGS models
library(snow) ## Allows clusters: i.e. subsidiary R programmes running separately on your computer's different CPUs
library(jagstools) ## devtools::install_github("johnbaums/jagstools") For extracting summary statistics from MCMC chain
library(ggplot2)
library(cowplot) ## For cowplot ggplot theme
library(ggthemes) ## For additional ggplot2 themes (e.g. "few") 
library(RColorBrewer)
library(grid) ## To adjust legend key width and size in ggplot2 themes that don't naturally support a grid
library(gridExtra) ## Facilitates easier labelling in ggplot2 (e.g. annote with extrafont fonts)
library(extrafont) ## For additional fonts in ggplot2
library(stargazer) ## For nice LaTeX tables
library(dplyr) ## For manipulating and munging data frames
library(tidyr) ## For tidying data frames
library(purrr) ## For manipulating vectors and functions (complements dplyr)


#######################################
#######################################
## Choose font type for graphs (note extrafont package installation instructions)
font_type <- choose_font("Palatino Linotype") ## Will use ggplot2 default font if not available

## Assign colours and names for later graphs ##
rcp_names <- c(expression("RCP 2.6 (420 ppmv CO"[2]*")"),
               expression("RCP 4.5 (540 ppmv CO"[2]*")"),
               expression("RCP 6.0 (670 ppmv CO"[2]*")"),
               expression("RCP 8.5 (940 ppmv CO"[2]*")"))  
# rcp_cols <- c("limegreen", "orchid", "orange", "red2")
rcp_cols <- c("darkgreen", "darkorchid", "darkorange2", "darkred")
rcp_fills <- c("lightgreen", "orchid", "orange", "red")

# prior_cols <- c("dodgerblue2", "limegreen", "orange", "red2", "gray20")
prior_cols <- c(brewer.pal(12, "Paired")[c(2, 4, 8, 6)], "#000000")
prior_names <- c("Strong Denier", "Moderate Denier", 
                 "Strong Lukewarmer", "Moderate Lukewarmer", 
                 "Noninformative")

#######################################
#######################################

## Negate version of %in% function
"%nin%" <- Negate("%in%")

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

### GGPLOT2 themes

theme_coefs <-
  # cowplot::theme_cowplot() +
  theme(
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
  # cowplot::theme_cowplot() +
  theme(
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

theme_2100 <-
  # cowplot::theme_cowplot() +
  theme(
    text = element_text(family = font_type),
    legend.position = "bottom",
    legend.title = element_blank(),
    # strip.text = element_text(size = 18, colour = "black"),
    strip.background = element_rect(fill = "white"), ## Facet strip
    panel.margin = unit(2, "lines") ## Increase gap between facet panels
  ) 

theme_tcr <-
  theme(
    text = element_text(family = font_type),
    legend.position = "bottom",
    legend.title = element_blank()
  )


#######################################
#######################################

#' ### Decimal/character align with xtable and LaTeX
#' ## Based on https://gist.github.com/jbryer/4458674
#' library(xtable) 
#' 
#' #' Prints a LaTeX table with numeric columns aligned on their decimal points.
#' #' 
#' #' This function wraps the \code{\link{xtable}} and \code{\link{print.xtable}}
#' #' functions in the \code{xtable} package so that numeric columns are aligned
#' #' on their decimal place.
#' #' 
#' #' See \url{http://jason.bryer.org/posts/2013-01-04/xtable_with_aligned_decimals.html}
#' #' for more information.
#' #' 
#' #' @author Jason Bryer <jason@@bryer.org>
#' #' @param x a data frame to create a LaTeX table from.
#' #' @param cols a numeric vector indicating which columns should be aligned on
#' #'        decimal points. It defaults to all columns of type numeric.
#' #' @param colAlignment named character vector where each element name corresponds to a
#' #         column name and the value is the LaTeX alignment (i.e. l, r, or c).
#' #' @param tocharFun the function used to convert the numeric vecotr to a character
#' #'        vector. This defaults to \code{\link{prettyNum}}, but other possible
#' #'        options are \code{\link{as.character}}, \code{\link{format}}, 
#' #'        \code{\link{formatC}}, or some other custom function.
#' #' @param ... other parameters passed to \code{tocharFun}, \code{\link{xtable}},
#' #'        and \code{\link{print.xtable}}.
#' #' @seealso xtable
#' #' @export
#' xtable.decimal <- function(x, 
#'                            cols=which(lapply(x, class) == 'numeric'), 
#'                            colAlignment, 
#'                            tocharFun=prettyNum,
#'                            ...) {
#'   splitCol <- function(x, ...) {
#'     s <- strsplit(tocharFun(unlist(x), ...), split='.', fixed=TRUE)
#'     right <- sapply(s, FUN=function(x) { ifelse(length(x) == 2, x[2], '0') })
#'     left <- sapply(s, FUN=function(x) { x[1] })
#'     data.frame(left=left, right=right, stringsAsFactors=FALSE)
#'   }
#'   
#'   cols <- cols[order(cols, decreasing=TRUE)]
#'   colnames <- names(x)
#'   for(i in cols) {
#'     if(i == 1) {
#'       tmp <- cbind(splitCol(x[,1], ...), x[,2:ncol(x)])
#'       names(tmp)[1:2] <- paste(names(tmp)[1], c('left','right'), sep='.')
#'       names(tmp)[3:ncol(x)] <- names(x)[2:ncol(x)]
#'       x <- tmp
#'     } else if(i == ncol(x)) {
#'       tmp <- cbind(x[,1:(ncol(x)-1)], splitCol(x[,ncol(x)], ...))
#'       names(tmp)[1:(ncol(tmp)-2)] <- names(x)[1:(ncol(x)-1)]
#'       names(tmp)[(ncol(tmp)-1):ncol(tmp)] <- paste(names(x)[ncol(x)], 
#'                                                    c('left','right'), sep='.')
#'       x <- tmp
#'     } else {
#'       tmp <- cbind(x[,1:(i-1)], splitCol(x[,i], ...), x[,(i+1):ncol(x)])
#'       names(tmp)[1:(i-1)] <- names(x)[1:(i-1)]
#'       names(tmp)[i:(i+1)] <- paste(names(x)[i], c('left','right'), sep='.')
#'       names(tmp)[(i+2):ncol(tmp)] <- names(x)[(i+1):ncol(x)]
#'       x <- tmp
#'     }
#'   }
#'   
#'   colnames[cols] <- paste('\\multicolumn{2}{c}{', colnames[cols], '}', sep='')
#'   colnames <- paste(colnames, collapse=' & ')
#'   
#'   addtorow <- list()
#'   addtorow$pos <- list()
#'   addtorow$pos[[1]] <- c(0)
#'   addtorow$command <- paste( colnames, ' \\\\ ', sep='')
#'   
#'   align <- rep('l', ncol(x))
#'   if(!missing(colAlignment)) {
#'     for(i in seq_along(colAlignment)) {
#'       align[names(x) == names(colAlignment)[i]] <- colAlignment[i]
#'     }
#'   }
#'   align[grep('.left$', names(x), perl=TRUE)] <- 'r@{.}'
#'   align <- c('l', align) #Add an alignment for row names
#'   
#'   xtab <- xtable(x, align=align, ...)
#'   args <- list(...)
#'   args <- args[names(args) != 'caption']
#'   do.call(print, 
#'           c(list(xtab, add.to.row=addtorow, include.rownames=FALSE, include.colnames=FALSE), 
#'             args)
#'           )
#' }
