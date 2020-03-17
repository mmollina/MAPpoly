#' Physical versus genetic distance 
#' 
#' This function plots scatterplot(s) of physical distance (in Mb) versus the genetic 
#' distance (in cM). Map(s) should be passed as a single object or a list of objects 
#' of class \code{mappoly.map}.
#'
#' @param  map.list A list or a single object of class \code{mappoly.map}
#' 
#' @param config A vector containing the position of the
#'     configuration(s) to be plotted. If \code{'best'} (default), plots the configuration
#'     with the highest likelihood for all elements in \code{'map.list'} 
#'     
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \url{https://doi.org/10.1534/g3.119.400378}
#'
#' @export plot_map_list
#'
plot_genome_vs_map<-function(map.list, config = "best"){
  if(class(map.list) == "mappoly.map")  
    map.list<-list(map.list)
  if(any(sapply(map.list, class)!="mappoly.map"))
    stop("All elemnts in 'map.list' should be of class 'mappoly.map'")
  if(length(config)!=length(map.list))
    config<-rep(config[1], length(map.list))
  geno.vs.map <- NULL
  for(i in 1:length(map.list)){
    ## Get linkage phase configuration
    LOD.conf <- get_LOD(map.list[[i]], sorted = FALSE)
    if(config[i] == "best") {
      i.lpc <- which.min(LOD.conf)
    } else if(is.numeric(config[i])){
      i.lpc <- config[i]
      if(i.lpc > length(LOD.conf))
        stop("invalid linkage phase configuration")
      } else stop("invalid linkage phase configuration")
    LG <- genomic.pos <- map.pos <- NULL
    geno.vs.map<-rbind(geno.vs.map,
                       data.frame(mrk.names = map.list[[i]]$info$mrk.names,
                                  map.pos = cumsum(imf_h(c(0, map.list[[i]]$maps[[i.lpc]]$seq.rf))),
                                  genomic.pos = map.list[[i]]$info$sequence.pos/1e6, 
                                  LG = paste0("LG_", i),
                                  chr = map.list[[i]]$info$sequence))
  }
  p<-ggplot2::ggplot(geno.vs.map, ggplot2::aes(genomic.pos, map.pos)) +
    ggplot2::geom_point(alpha = 1/5, ggplot2::aes(colour = LG)) +
    ggplot2::facet_grid(chr~LG) +  ggplot2::xlab("Genome Position (Mb)") + ggplot2::ylab("Map Distance (cM)")
  p
}
