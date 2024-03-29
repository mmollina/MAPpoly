#' Plot a genetic map 
#' 
#' This function plots a genetic linkage map(s) generated by \code{MAPpoly}.
#' The map(s) should be passed as a single object or a list of objects of class \code{mappoly.map}.
#'
#' @param  map.list A list of objects or a single object of class \code{mappoly.map}
#'
#' @param horiz logical. If FALSE, the maps are plotted vertically with the first map to the left. 
#'              If TRUE  (default), the maps are plotted horizontally with the first at the bottom
#'
#' @param col a vector of colors for each linkage group.  (default = 'lightgray')
#'            \code{ggstyle} produces maps using the default \code{ggplot} color palette. 
#'            
#' @param title a title (string) for the maps (default = 'Linkage group')
#'
#' @return A \code{data.frame} object containing the name of the markers and their genetic position
#' 
#' @examples
#'  ## hexafake map
#'  plot_map_list(maps.hexafake, horiz = FALSE)
#'  plot_map_list(maps.hexafake, col = c("#999999", "#E69F00", "#56B4E9"))
#'  
#'  ## solcap map
#'  plot_map_list(solcap.dose.map, col = "ggstyle")
#'  plot_map_list(solcap.dose.map, col = "mp_pallet3", horiz = FALSE)
#'  
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \doi{10.1534/g3.119.400378}
#'
#' @export plot_map_list
#'
plot_map_list <- function(map.list, horiz = TRUE, 
                          col = "lightgray", 
                          title = "Linkage group"){
  if(inherits(map.list, "mappoly.map"))
    map.list <- list(map.list)
  if (any(!sapply(map.list, inherits, "mappoly.map"))) 
    stop("All elemnts in 'map.list' should be of class 'mappoly.map'")
  if(all(col  ==  "ggstyle"))
    col  <- gg_color_hue(length(map.list))
  if(all(col  ==  "mp_pallet1"))
    col  <- mp_pallet1(length(map.list))
  if(all(col  ==  "mp_pallet2"))
    col  <- mp_pallet2(length(map.list))
  if(all(col  ==  "mp_pallet3"))
    col  <- mp_pallet3(length(map.list))
  if(length(col) == 1)
    col <- rep(col, length(map.list))
  z <- NULL
  if(is.null(names(map.list)))
    names(map.list) <- 1:length(map.list)
  max.dist <- max(sapply(map.list, function(x) sum(imf_h(x$maps[[1]]$seq.rf))))
  if(horiz){
    plot(0, 
         xlim = c(0, max.dist), 
         ylim = c(0,length(map.list)+1), 
         type = "n", axes = FALSE, 
         xlab = "Map position (cM)", 
         ylab = title)
    axis(1)
    for(i in 1:length(map.list)){
      d <- extract_map(map.list[[i]])
      z <- rbind(z, data.frame(mrk = map.list[[i]]$info$mrk.names, 
                             LG = names(map.list)[i], pos = d))
      plot_one_map(d, i = i, horiz = TRUE, col = col[i])   
    }
    axis(2, at = 1:length(map.list), labels = names(map.list), lwd = 0, las = 2)
  } else{
    plot(0, 
         ylim = c(-max.dist, 0), 
         xlim = c(0,length(map.list)+1), 
         type = "n", axes = FALSE, 
         ylab = "Map position (cM)", 
         xlab = title)
    x <- axis(2, labels = FALSE, lwd = 0)
    axis(2, at = x, labels = abs(x))
    for(i in 1:length(map.list)){
      d <- extract_map(map.list[[i]])
      z <- rbind(z, data.frame(mrk = map.list[[i]]$info$mrk.names, 
                             LG = names(map.list)[i],pos = d))
      plot_one_map(d, i = i, horiz = FALSE, col = col[i])  
    }
    axis(3, at = 1:length(map.list), labels = names(map.list), lwd = 0, las = 2)
  }
  invisible(z)
}

#' Extract the maker position from an object of class 'mappoly.map'
#'
#' @param input.map An object of class \code{mappoly.map}
#' @param phase.config which phase configuration should be used. "best" (default) 
#'                     will choose the maximum likelihood configuration
#' @examples
#'  x <- maps.hexafake[[1]]$info$genome.pos/1e6
#'  y <- extract_map(maps.hexafake[[1]])
#'  plot(y~x, ylab = "Map position (cM)", xlab = "Genome Position (Mbp)")
#' @export
extract_map <- function(input.map, phase.config = "best")
{
  if (!inherits(input.map, "mappoly.map")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  }
  ## choosing the linkage phase configuration
  LOD.conf <- get_LOD(input.map, sorted = FALSE)
  if(phase.config  ==  "best") {
    i.lpc <- which.min(LOD.conf)
  } else if (phase.config > length(LOD.conf)) {
    stop("invalid linkage phase configuration")
  } else i.lpc <- phase.config
  x <- cumsum(c(0, imf_h(input.map$maps[[i.lpc]]$seq.rf)))
  names(x) <- input.map$info$mrk.names
  x
}

#' plot a single linkage group with no phase
#' @keywords internal
plot_one_map <- function(x, i = 0, horiz = FALSE, col = "lightgray")
{
  if(horiz)
  {
    rect(xleft = x[1], ybottom = i-0.25, 
         xright = tail(x,1), ytop = i+0.25,
         col = col)
    for(j in 1:length(x))
      lines(x = c(x[j], x[j]), y = c(i-0.25, i+0.25), lwd = .5)
  } else {
    x <- -rev(x)
    rect(xleft = i-0.25, ybottom = x[1], 
         xright = i+0.25, ytop = tail(x,1),
         col = col)
    for(j in 1:length(x))
      lines(y = c(x[j], x[j]), x = c(i-0.25, i+0.25), lwd = .5)
  }
}
