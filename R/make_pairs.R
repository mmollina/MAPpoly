#'  Subset pairwise recombination fractions
#'
#'  Get a subset of an object of class \code{poly.est.two.pts.pairwise} (i.e.
#'  recombination fraction) and LOD score statistics for all possible linkage
#'  phase combinations based on a sequence of markers.
#'
#' @param input.twopt an object of class \code{poly.est.two.pts.pairwise}
#'
#' @param input.seq an object of class \code{mappoly.sequence}, with 
#' a sequence of markers contained in \code{input.twopt}
#'
#' @return an object of class \code{poly.est.two.pts.pairwise} which is a
#'  subset of \code{input.twopt}. 
#'  See \code{\link[mappoly]{est_pairwise_rf}} for details
#'
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
#'     mat.full<-rf_list_to_matrix(input.twopt = all.pairs)
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
#'     p1<-make_pairs_mappoly(input.seq = lg1, input.twopt = all.pairs)
#'     p2<-make_pairs_mappoly(input.seq = lg2, input.twopt = all.pairs)
#'     p3<-make_pairs_mappoly(input.seq = lg3, input.twopt = all.pairs)
#'
#'     m1<-rf_list_to_matrix(input.twopt = p1)
#'     m2<-rf_list_to_matrix(input.twopt = p2)
#'     m3<-rf_list_to_matrix(input.twopt = p3)
#'
#'     op<-par(mfrow = c(1,3), pty = "s")
#'     plot(m1, main.text = "LG1")
#'     plot(m2, main.text = "LG2")
#'     plot(m3, main.text = "LG3")
#'     par(op)
#'    }
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

make_pairs_mappoly<-function(input.twopt, input.seq){
  ## checking for correct object
  input_classes1 <- c("mappoly.sequence")
  if (!inherits(input.seq, input_classes1)) {
    stop(deparse(substitute(input.seq)), " is not an object of class 'mappoly.sequence'")
  }
  input_classes2 <-c("poly.est.two.pts.pairwise")
  if (!inherits(input.twopt, input_classes2)) {
    stop(deparse(substitute(input.twopt)), " is not an object of class 'poly.est.two.pts.pairwise'")
  }
  if(!all(input.seq$seq.num%in%input.twopt$seq.num))
    stop(deparse(substitute(input.twopt)), " does not contain (some) markers present in ", deparse(substitute(input.seq)))

  w<-combn(sort(input.seq$seq.num), 2)
  idx<-apply(w, 2, function(x) paste(sort(x), collapse="-"))
  input.twopt$n.mrk <- length(input.seq$seq.num)
  input.twopt$seq.num <- sort(input.seq$seq.num)
  input.twopt$pairwise<-input.twopt$pairwise[idx]
  input.twopt
}
