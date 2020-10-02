#' Homolog probabilities 
#' 
#' Compute homolog probabilities for all individuals in the full-sib
#' population given a map and conditional genotype probabilities. 
#'
#' @param input.genoprobs an object of class \code{mappoly.genoprob}
#'
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
#' 
#'@examples
#'    \donttest{
#'      ## tetraploid example
#'      w1 <- calc_genoprob(solcap.dose.map[[1]])
#'      h.prob <- calc_homoprob(w1)
#'      print(h.prob)
#'      plot(h.prob, ind = 5, use.plotly = FALSE)
#'      ## using error modeling (removing noise)
#'      w2 <- calc_genoprob_error(solcap.err.map[[1]])
#'      h.prob2 <- calc_homoprob(w2)
#'      print(h.prob2)
#'      plot(h.prob2, ind = 5, use.plotly = FALSE)
#'   }
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari M., Olukolu B. A.,  Pereira G. da S., 
#'     Khan A., Gemenet D., Yencho G. C., Zeng Z-B. (2020), 
#'     Unraveling the Hexaploid Sweetpotato Inheritance 
#'     Using Ultra-Dense Multilocus Mapping, 
#'     _G3: Genes, Genomes, Genetics_. 
#'     \url{https://doi.org/10.1534/g3.119.400620} 
#'     
#' @importFrom ggplot2 ggplot geom_density ggtitle facet_grid theme_minimal ylab xlab aes vars
#' @importFrom plotly ggplotly
#' @export calc_homoprob
#' 
calc_homoprob<-function(input.genoprobs, verbose = TRUE){
  if(inherits(input.genoprobs, "mappoly.genoprob")) 
    input.genoprobs <- list(input.genoprobs)
  if(!inherits(input.genoprobs, "list"))
    stop(deparse(substitute(input.genoprobs)), 
         " is not an object of class 'mappoly.genoprob' neither a list containing 'mappoly.genoprob' objects.")
  if (any(!sapply(input.genoprobs, inherits, "mappoly.genoprob"))) 
    stop(deparse(substitute(input.genoprobs)), 
         " is not an object of class 'mappoly.genoprob' neither a list containing 'mappoly.genoprob' objects.")
  df.res <- NULL
  for(j in 1:length(input.genoprobs)){
    if (verbose) cat("\nLinkage group ", j, "...")
    stt.names<-dimnames(input.genoprobs[[j]]$probs)[[1]] ## state names
    mrk.names<-dimnames(input.genoprobs[[j]]$probs)[[2]] ## mrk names
    ind.names<-dimnames(input.genoprobs[[j]]$probs)[[3]] ## individual names
    v<-c(2,4,6,8,10,12)
    names(v)<-choose(v,v/2)^2
    m<-v[as.character(length(stt.names))]
    hom.prob<-array(NA, dim = c(m*2, length(mrk.names), length(ind.names)))
    dimnames(hom.prob)<-list(letters[1:(2*m)], mrk.names, ind.names)
    for(i in letters[1:(2*m)])
      hom.prob[i,,] <- apply(input.genoprobs[[j]]$probs[grep(stt.names, pattern = i),,], c(2,3), function(x) round(sum(x, na.rm = TRUE),4))
    df.hom<-reshape2::melt(hom.prob)
    map<-data.frame(map.position = input.genoprobs[[j]]$map, marker = names(input.genoprobs[[j]]$map))
    colnames(df.hom)<-c("homolog", "marker", "individual", "probability")
    df.hom<-merge(df.hom, map, sort = FALSE)
    df.hom$LG <- j
    df.res<-rbind(df.res, df.hom)
  }
  if (verbose) cat("\n")
  structure(list(info = list(m = m, nind = length(ind.names)) , homoprob = df.res), class = "mappoly.homoprob")
}

#' @export
print.mappoly.homoprob<-function(x, ...){
  head(x$homoprob, 20)
}

#' Plots mappoly.homoprob
#' 
#' @param x an object of class \code{mappoly.homoprob}
#' 
#' @param stack logical. If \code{TRUE}, probability profiles of all homologues
#'              are stacked in the plot (default = FALSE)
#'              
#' @param lg indicates which linkage group should be plotted. If \code{NULL} 
#'           (default), it plots the first linkage group. If 
#'           \code{"all"}, it plots all linkage groups
#'           
#' @param ind indicates which individuals should be plotted. It can be the 
#'            position of the individuals in the dataset or it's name. 
#'            If \code{NULL} (default), the function plots the first 
#'            individual
#'            
#' @param use.plotly if \code{TRUE} (default), it uses plotly interactive graphics
#'
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
#' 
#' @param ... unused arguments
#' @export
plot.mappoly.homoprob<-function(x, stack = FALSE, lg = NULL, 
                                ind = NULL, use.plotly = TRUE, verbose = TRUE,  ...){
  all.ind<-as.character(unique(x$homoprob$individual))
  #### Individual handling ####
  if(length(ind) > 1){
      if (verbose) message("More than one individual provided: using the first one")
    ind<-ind[1]
  }
  if(is.null(ind)){
    ind<-as.character(all.ind[1])
    df.pr1<-subset(x$homoprob, individual == ind)  
  } else if(is.numeric(ind)) {
    if(ind > length(all.ind))
      stop("Please chose an individual number between 1 and ", length(all.ind))
    ind<-as.character(all.ind[ind])
    df.pr1<-subset(x$homoprob, individual == ind)  
  } else if (is.character(ind)){
    if(!ind%in%all.ind)
      stop("Invalid individual name")
  } else stop("Invalid individual name")
  #### LG handling ####
  if(is.null(lg))
    lg <- 1
  if(all(lg=="all"))
    lg <- unique(x$homoprob$LG)
  LG<-individual<-map.position<-probability<-homolog<-NULL
  if(length(lg) > 1 & !stack)
  {
    if (verbose) message("Using 'stack = TRUE' to plot multiple linkage groups")
    stack <- TRUE
  }
  if(stack){
    ##subset linkage group
    if(!is.null(lg)){
      df.pr1<-subset(x$homoprob, LG%in%lg)
      df.pr1<-subset(df.pr1, individual == ind)
    } else 
      df.pr1<-subset(x$homoprob, individual == ind)
    p <- ggplot2::ggplot(df.pr1, ggplot2::aes(x = map.position, y = probability, fill = homolog, color  = homolog)) +
      ggplot2::geom_density(stat = "identity", alpha = 0.7, position = "stack") + 
      ggplot2::ggtitle(ind) + 
      ggplot2::facet_grid(rows = ggplot2::vars(LG)) + 
      ggplot2::ylab(label = "Homologs probabilty") +
      ggplot2::xlab(label = "Map position")
  } else {
    ##subset linkage group
    if(is.null(lg)){
      lg<-1
      df.pr1<-subset(x$homoprob, LG %in% lg)
    } else df.pr1<-subset(x$homoprob, LG %in% lg)
    df.pr1<-subset(df.pr1, individual == ind)  
    p <- ggplot2::ggplot(df.pr1, ggplot2::aes(x = map.position, y = probability, fill = homolog, color  = homolog)) +
      ggplot2::geom_density(stat = "identity", alpha = 0.7) + 
      ggplot2::ggtitle(paste(ind, "   LG", lg)) + 
      ggplot2::facet_grid(rows = ggplot2::vars(homolog)) + 
      ggplot2::theme_minimal() + 
      ggplot2::ylab(label = "Homologs probabilty") +
      ggplot2::xlab(label = "Map position")
  }
  if(use.plotly)
    p <- plotly::ggplotly(p)
  return(p)
}
