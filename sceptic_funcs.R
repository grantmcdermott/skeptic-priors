### Facet_wrap labeller for parsing strip titles and text in ggplot2 wrapped facets ###
## Especially useful for Geek letters and math
## See: http://stackoverflow.com/a/33488476/4115816

facet_wrap_labeller <- function(gg.plot, labels = NULL, labeller = label_value) {
  #works with R 3.1.2 and ggplot2 1.0.1
  require(gridExtra)
  
  # old labels
  g <- ggplotGrob(gg.plot)
  gg <- g$grobs      
  strips <- grep("strip_t", names(gg))
  modgrobs <- lapply(strips, function(i) {
    getGrob(gg[[i]], "strip.text", grep=TRUE, global=TRUE)
  })
  old_labels <- sapply(modgrobs, function(i) i$label)
  
  # find new labels
  if (is.null(labels)) # no labels given, use labeller function
    new_labels <- labeller(names(gg.plot$facet$facets), old_labels)
  else if (is.null(names(labels))) # unnamed list of labels, take them in order
    new_labels <- as.list(labels)
  else { # named list of labels, go by name where provided, otherwise keep old
    new_labels <- sapply(as.list(old_labels), function(i) {
      if (!is.null(labels[[i]])) labels[[i]] else i
    })
  }
  
  # replace labels
  for(i in 1:length(strips))  {
    gg[[strips[i]]]$children[[modgrobs[[i]]$name]] <- 
      editGrob(modgrobs[[i]], label=new_labels[[i]])
  }
  
  g$grobs <- gg
  class(g) = c("arrange", "ggplot",class(g))
  return(g) 
}

## Additional code for getting facet_wrap_labeller function to work with R 3.2.2
## Needed because of changes to gridExtra
print.arrange <- function(x){
  grid::grid.draw(x)
}


#######################################
#######################################

### Decimal function (to make sure, e.g. three decimals places are always printed in tables) ###

decimals <- function(x, k) {
  format(round(x, k), nsmall = k)
  }
