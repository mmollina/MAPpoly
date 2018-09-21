#'  Remove markers that do not meet a LOD criteria
#'
#'  Remove markers that do not meet a LOD and recombination fraction
#'  criteria for at least a percentage of the pairwise markers
#'  combinations.
#'
#' \code{thresh_LOD_ph} should be set in order to only selects
#'     recombination fractions which have LOD scores associated to the
#'     linkage phase configuration bigger than \code{thresh_LOD_ph}
#'     for the second most likely linkage phase configuration.
#'     Notice that eliminated markers are usually unlinked to the
#'     set of markers analyzed.
#'
#' @param input.twopt an object of class \code{poly.est.two.pts.pairwise}
#'
#' @param thresh.LOD.ph LOD score threshold for linkage phase
#'     configuration.
#'
#' @param thresh.LOD.rf LOD score threshold for recombination phase
#'
#' @param thresh.rf recombination fraction threshold
#'
#' @param n.clusters number of parallel processes to spawn
#'
#' @param thresh.perc threshold for the percentage of the pairwise markers
#'  combinations that should be considered in order to
#'  kept the marker. A \code{thresh.perc = 0.05} means that, at least
#'  5% of the pairwise combinations should be present in order to
#'  kept the marker.
#'
#' @return A filtered object of class \code{mappoly.sequence}
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
#'                          inter = TRUE,
#'                          comp.mat = TRUE, #this data has physical information
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
#'
#'     ## Removing disruptive SNPs
#'
#'     lg1.filt<-rf_snp_filter(p1, 5, 5, 0.15, thresh.perc = 0.05)
#'     lg2.filt<-rf_snp_filter(p2, 5, 5, 0.15, thresh.perc = 0.05)
#'     lg3.filt<-rf_snp_filter(p3, 5, 5, 0.15, thresh.perc = 0.05)
#'
#'     p1.filt<-make_pairs_mappoly(input.seq = lg1.filt, input.twopt = all.pairs)
#'     p2.filt<-make_pairs_mappoly(input.seq = lg2.filt, input.twopt = all.pairs)
#'     p3.filt<-make_pairs_mappoly(input.seq = lg3.filt, input.twopt = all.pairs)
#'
#'     m1.filt<-rf_list_to_matrix(input.twopt = p1.filt)
#'     m2.filt<-rf_list_to_matrix(input.twopt = p2.filt)
#'     m3.filt<-rf_list_to_matrix(input.twopt = p3.filt)
#'
#'     op<-par(mfrow = c(2,3), pty = "s")
#'     plot(m1, main.text = "LG1")
#'     plot(m2, main.text = "LG2")
#'     plot(m3, main.text = "LG3")
#'     plot(m1.filt, main.text = "LG1.filt")
#'     plot(m2.filt, main.text = "LG2.filt")
#'     plot(m3.filt, main.text = "LG3.filt")
#'     par(op)
#'    }
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2018) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _submited_. \url{https://doi.org/10.1101/415232}
#'
#' @export rf_snp_filter
#'
rf_snp_filter<-function(input.twopt,
                        thresh.LOD.ph = 5,
                        thresh.LOD.rf = 5,
                        thresh.rf = 0.15,
                        thresh.perc = 0.05,
                        n.clusters = 1)
{
  ## checking for correct object
  input_classes <-c("poly.est.two.pts.pairwise")
  if (!inherits(input.twopt, input_classes)) {
    stop(deparse(substitute(input.twopt)),
         " is not an object of class 'poly.est.two.pts.pairwise'")
  }
  thresh.missing<-1-thresh.perc
  rf_mat<- rf_list_to_matrix(input.twopt = input.twopt, thresh.LOD.ph = thresh.LOD.ph,
                             thresh.LOD.rf = thresh.LOD.rf, thresh.rf = thresh.rf,
                             n.clusters = n.clusters, verbose = FALSE)
  ## Removing markers that have too many large recombination fractions
  x <- apply(rf_mat$rec.mat, 1, function(x) sum(is.na(x)))
  o <- names(which(x < quantile(x, probs = thresh.missing)))
  ch_filt<-make_seq_mappoly(input.obj = get(input.twopt$data.name), arg = o, data.name = input.twopt$data.name)
  ch_filt
}
