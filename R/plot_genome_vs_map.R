#' Physical versus genetic distance 
#' 
#' This function plots scatterplot(s) of physical distance (in Mbp) versus the genetic 
#' distance (in cM). Map(s) should be passed as a single object or a list of objects 
#' of class \code{mappoly.map}.
#'
#' @param  map.list A list or a single object of class \code{mappoly.map}
#' 
#' @param phase.config A vector containing which phase configuration should be
#'  plotted. If \code{'best'} (default), plots the configuration
#'  with the highest likelihood for all elements in \code{'map.list'} 
#'     
#' @param same.ch.lg Logical. If \code{TRUE} displays only the scatterplots between the 
#'   chromosomes and linkage groups with the same number. Default is \code{FALSE}.   
#'    
#' @examples
#'   ## tetraploid example
#'   plot_genome_vs_map(solcap.mds.map)
#'   plot_genome_vs_map(solcap.mds.map, same.ch.lg = TRUE)
#'   
#'   ## hexaploid example
#'   plot_genome_vs_map(maps.hexafake)
#'   plot_genome_vs_map(maps.hexafake, same.ch.lg = TRUE)
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
#' @export plot_genome_vs_map
plot_genome_vs_map<-function(map.list, phase.config = "best", same.ch.lg = FALSE){
  if(inherits(map.list, "mappoly.map")) 
    map.list <- list(map.list)
  if(!inherits(map.list, "list"))
    stop("All elemnts in 'map.list' should be of class 'mappoly.map'")
  if (any(!sapply(map.list, inherits, "mappoly.map"))) 
    stop("All elemnts in 'map.list' should be of class 'mappoly.map'")
  if(length(phase.config)!=length(map.list))
    phase.config<-rep(phase.config[1], length(map.list))
  if(any(sapply(map.list, function(x) is.null(x$info$sequence.pos))))
    stop("All elemnts in 'map.list' should have the genomic position of the markers")
  geno.vs.map <- NULL
  for(i in 1:length(map.list)){
    ## Get linkage phase configuration
    LOD.conf <- get_LOD(map.list[[i]], sorted = FALSE)
    if(phase.config[i] == "best") {
      i.lpc <- which.min(LOD.conf)
    } else if(is.numeric(phase.config[i])){
      i.lpc <- phase.config[i]
      if(i.lpc > length(LOD.conf))
        stop("invalid linkage phase configuration")
    } else stop("invalid linkage phase configuration")
    LG <- genomic.pos <- map.pos <- NULL
    geno.vs.map<-rbind(geno.vs.map,
                       data.frame(mrk.names = map.list[[i]]$info$mrk.names,
                                  map.pos = cumsum(imf_h(c(0, map.list[[i]]$maps[[i.lpc]]$seq.rf))),
                                  genomic.pos = map.list[[i]]$info$sequence.pos/1e6, 
                                  LG = as.factor(i),
                                  chr = as.factor(map.list[[i]]$info$sequence)))
  }
  geno.vs.map$chr <- factor(geno.vs.map$chr, levels = sort(levels(geno.vs.map$chr))) 
  if(same.ch.lg){
    p<-ggplot2::ggplot(geno.vs.map, ggplot2::aes(genomic.pos, map.pos)) +
      ggplot2::geom_point(alpha = 1/5, ggplot2::aes(colour = LG)) +
      ggplot2::facet_wrap(~LG, nrow = floor(sqrt(length(map.list)))) +  
      ggplot2::labs(subtitle = "Linkage group", x = "Genome position (Mbp)", y = "Map position (cM)") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position="none", plot.subtitle = ggplot2::element_text(hjust = 0.5)) 
    
  } else {
    p<-ggplot2::ggplot(geno.vs.map, ggplot2::aes(genomic.pos, map.pos)) +
      ggplot2::geom_point(alpha = 1/5, ggplot2::aes(colour = LG)) +
      ggplot2::facet_grid(LG~chr) + 
      ggplot2::labs(x = "Genome position (Mbp)", y = "Map position (cM)") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position="none") 
  }
  p
}
