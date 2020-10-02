#'  Subset recombination fraction matrices
#'
#'  Get a subset of an object of class \code{mappoly.rf.matrix}, i.e.
#'  recombination fraction and LOD score matrices based in a
#'  sequence of markers.
#'
#' @param input.mat an object of class \code{mappoly.rf.matrix}
#'
#' @param input.seq an object of class \code{mappoly.sequence}, with 
#' a sequence of markers contained in \code{input.mat}
#'
#' @return an object of class \code{mappoly.rf.matrix},
#' which is a subset of \code{'input.mat'}. 
#' See \code{\link[mappoly]{rf_list_to_matrix}} for details
#'     
#' @examples
#'     # sequence with 20 markers
#'     mrk.seq<-make_seq_mappoly(hexafake, 1:20)
#'     mrk.pairs<-est_pairwise_rf(input.seq = mrk.seq,
#'                                verbose=TRUE)
#'     ## Full recombination fraction matrix
#'     mat<-rf_list_to_matrix(input.twopt=mrk.pairs)
#'     plot(mat)
#'     ## Matrix subset
#'     id <- make_seq_mappoly(hexafake, 1:10)
#'     mat.sub<-make_mat_mappoly(mat, id)
#'     plot(mat.sub)
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
#' @export
#' 
make_mat_mappoly<-function(input.mat, input.seq){
  ## checking for correct object
  input_classes1 <- c("mappoly.sequence")
  if (!inherits(input.seq, input_classes1)) {
    stop(deparse(substitute(input.seq)), " is not an object of class 'mappoly.sequence'")
  }
  input_classes2 <-c("mappoly.rf.matrix")
  if (!inherits(input.mat, input_classes2)) {
    stop(deparse(substitute(input.mat)), " is not an object of class 'mappoly.rf.matrix'")
  }
  if(input.mat$cl == "poly.haplo.est.two.pts.pairwise")
    stop("This function does not work for recombination fractions matrices originated from blocks of markers")
  input.mat$thresh.LOD.ph <- NULL
  input.mat$thresh.LOD.rf <- NULL
  input.mat$thresh.rf <- NULL
  id <- get(input.seq$data.name, pos =1)$mrk.names[input.seq$seq.num]
  input.mat$rec.mat <- input.mat$rec.mat[id, id]
  input.mat$lod.mat <- input.mat$lod.mat[id, id]
  input.mat
}
