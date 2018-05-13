#'  Subset recombination fraction matrices
#'
#'  Get a subset of an object of class \code{'mappoly.rf.matrix'}, i.e.,
#'  recombination fraction and LOD score matrices based in a
#'  sequence of markers.
#'
#' @param input.mat an object of class \code{mappoly.rf.matrix}.
#'
#' @param input.seq an object of class \code{mappoly.sequence}.
#'     It must be contained in 'input.mat'
#'
#' @return a subset of \code{'input.mat'}. This object is also
#'     of class \code{mappoly.rf.matrix}.
#' @examples
#'  \dontrun{
#'     data(hexafake)
#'     all.mrk<-make_seq_mappoly(hexafake, 'all')
#'     red.mrk<-elim_redundant(all.mrk)
#'     unique.mrks<-make_seq_mappoly(red.mrk)
#'     counts.web<-cache_counts_twopt(unique.mrks, get.from.web = TRUE)
#'     all.pairs<-est_pairwise_rf(input.seq = unique.mrks,
#'                                count.cache = counts.web,
#'                                n.clusters = 16,
#'                                verbose=TRUE)
#'
#'     ## Full recombination fraction matrix
#'     mat.full<-rf_list_to_matrix(input.twopt=all.pairs)
#'     plot(mat.full)
#'
#'     lgs <- group_mappoly(input.mat = mat.full,
#'                          input.seq = unique.mrks,
#'                          expected.groups = 3,
#'                          verbose=TRUE)
#'     lgs
#'     plot(lgs)
#'     lg1 <- make_seq_mappoly(lgs, 1)
#'     lg2 <- make_seq_mappoly(lgs, 2)
#'     lg3 <- make_seq_mappoly(lgs, 3)
#'
#'     ##Plot matrices
#'     m1<-make_mat_mappoly(input.seq = lg1, input.mat = mat.full)
#'     m2<-make_mat_mappoly(input.seq = lg2, input.mat = mat.full)
#'     m3<-make_mat_mappoly(input.seq = lg3, input.mat = mat.full)
#'     op<-par(mfrow = c(1,3), pty = "s")
#'     plot(m1, main.text = "LG1")
#'     plot(m2, main.text = "LG2")
#'     plot(m3, main.text = "LG3")
#'     par(op)
#'    }
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2017) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _submited_
#'
#' @export

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
  input.mat$rec.mat <- input.mat$rec.mat[as.character(input.seq$seq.num), as.character(input.seq$seq.num)]
  input.mat$lod.mat <- input.mat$lod.mat[as.character(input.seq$seq.num), as.character(input.seq$seq.num)]
  input.mat
}
