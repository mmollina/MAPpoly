#'  Remove markers that do not meet a LOD criteria
#'
#'  Remove markers that do not meet a LOD and recombination fraction
#'  criteria for at least a percentage of the pairwise marker
#'  combinations. It also removes markers with strong evidence of
#'  linkage across the whole linkage group (false positive).
#'
#' \code{thresh.LOD.ph} should be set in order to only select
#'     recombination fractions that have LOD scores associated to the
#'     linkage phase configuration higher than \code{thresh_LOD_ph}
#'     when compared to the second most likely linkage phase configuration.
#'     That action usually eliminates markers that are unlinked to the
#'     set of analyzed markers.
#'
#' @param input.twopt an object of class \code{poly.est.two.pts.pairwise}
#'
#' @param thresh.LOD.ph LOD score threshold for linkage phase configuration (default = 5)
#'
#' @param thresh.LOD.rf LOD score threshold for recombination fraction (default = 5)
#'
#' @param thresh.rf threshold for recombination fractions (default = 0.15) 
#'
#' @param probs indicates the probability corresponding to the filtering quantiles. (defaul = c(0.05, 1))
#' 
#' @param ncpus number of parallel processes (i.e. cores) to spawn (default = 1)
#' 
#' @param diagnostic.plot if \code{TRUE} produces a diagnostic plot
#' 
#' @return A filtered object of class \code{mappoly.sequence}. 
#' See \code{\link[mappoly]{make_seq_mappoly}} for details
#' 
#' @examples
#'     all.mrk<-make_seq_mappoly(hexafake, 1:20)
#'     red.mrk<-elim_redundant(all.mrk)
#'     unique.mrks<-make_seq_mappoly(red.mrk)
#'     all.pairs<-est_pairwise_rf(input.seq = unique.mrks,
#'                                ncpus = 1,
#'                                verbose=TRUE)
#'
#'     ## Full recombination fraction matrix
#'     mat.full<-rf_list_to_matrix(input.twopt=all.pairs)
#'     plot(mat.full)
#'
#'     ## Removing disruptive SNPs
#'     tpt.filt<-rf_snp_filter(all.pairs, 2, 2, 0.07, probs = c(0.15, 1))
#'     p1.filt<-make_pairs_mappoly(input.seq = tpt.filt, input.twopt = all.pairs)
#'     m1.filt<-rf_list_to_matrix(input.twopt = p1.filt)
#'     plot(mat.full, main.text = "LG1")
#'     plot(m1.filt, main.text = "LG1.filt")
#'    
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu} with updates by Gabriel Gesteira, \email{gabrielgesteira@usp.br}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \url{https://doi.org/10.1534/g3.119.400378}
#'
#' @export rf_snp_filter
#' @importFrom ggplot2 ggplot geom_histogram aes scale_fill_manual xlab ggtitle
#' 
rf_snp_filter<-function(input.twopt,
                        thresh.LOD.ph = 5,
                        thresh.LOD.rf = 5,
                        thresh.rf = 0.15,
                        probs = c(0.05, 1),
                        ncpus = 1L,
                        diagnostic.plot = TRUE)
{
    ## checking for correct object
    input_classes <-c("poly.est.two.pts.pairwise")
    if (!inherits(input.twopt, input_classes)) {
        stop(deparse(substitute(input.twopt)),
             " is not an object of class 'poly.est.two.pts.pairwise'")
    }
    probs <- range(probs)
    ## Getting filtered rf matrix
    rf_mat<- rf_list_to_matrix(input.twopt = input.twopt, thresh.LOD.ph = thresh.LOD.ph,
                               thresh.LOD.rf = thresh.LOD.rf, thresh.rf = thresh.rf,
                               ncpus = ncpus, verbose = FALSE)
    x <- apply(rf_mat$rec.mat, 1, function(x) sum(!is.na(x)))
    th <- quantile(x, probs = probs)
    rem <- c(which(x < th[1]), which(x > th[2]))
    ids <- names(which(x >= th[1] & x <= th[2]))
    value <- type <- NULL
    if(diagnostic.plot){
        d<-rbind(data.frame(type = "original", value = x),
                 data.frame(type = "filtered", value = x[ids]))
        p<-ggplot2::ggplot(d, ggplot2::aes(value)) +
            ggplot2::geom_histogram(ggplot2::aes(fill = type),
                                    alpha = 0.4, position = "identity", binwidth = 30) +
            ggplot2::scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
            ggplot2::ggtitle( paste0("Filtering probs: [", probs[1], " : ", probs[2], "]")) +
            ggplot2::xlab(paste0("Non 'NA' values at LOD.ph = ", thresh.LOD.ph, ", LOD.rf = ", thresh.LOD.rf, ", and thresh.rf = ", thresh.rf))
        print(p)
    }
    ## Returning sequence object
    ch_filt<-make_seq_mappoly(input.obj = get(input.twopt$data.name, pos = 1), arg = ids, data.name = input.twopt$data.name)
    ch_filt
}
