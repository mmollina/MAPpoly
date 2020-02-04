#'  Remove markers that do not meet a LOD criteria
#'
#'  Remove markers that do not meet a LOD and recombination fraction
#'  criteria for at least a percentage of the pairwise marker
#'  combinations. It also removes markers with strong evidence of
#' linkage across the whole linkage group (false positive).
#'
#' \code{thresh.LOD.ph} should be set in order to only select
#'     recombination fractions that have LOD scores associated to the
#'     linkage phase configuration higher than \code{thresh_LOD_ph}
#'     when compared to the second most likely linkage phase configuration.
#'     Please notice that eliminated markers are usually unlinked to the
#'     set of markers analyzed.
#'
#' @param input.twopt an object of class \code{poly.est.two.pts.pairwise}
#'
#' @param thresh.LOD.ph LOD score threshold for linkage phase
#'     configuration (default = 5)
#'
#' @param thresh.LOD.rf LOD score threshold for recombination fraction (default = 5)
#'
#' @param thresh.rf threshold for recombination fractions (default = 0.15) 
#'#'
#' @param thresh.perc threshold for the percentage of the pairwise marker
#'  combinations that should be considered in order to
#'  keep the marker. A \code{thresh.perc = 0.05} means that, at least
#'  5\% of the pairwise combinations should be present in order to
#'  keep the marker (default = 0.05)
#'
#' @param remove.fp numeric value from 0.0 to 0.5 (default = NULL). When defined,
#' this parameter identifies and removes markers that presented more than 90\%
#' of its pairwise recombination fractions below \code{remove.fp} value throughout
#' the linkage group
#' 
#' @param n.clusters number of parallel processes (i.e. cores) to spawn (default = 1)
#'
#' @return A filtered object of class \code{mappoly.sequence}. 
#' See \code{\link[mappoly]{make_seq_mappoly}} for details
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
#'
rf_snp_filter<-function(input.twopt,
                        thresh.LOD.ph = 5,
                        thresh.LOD.rf = 5,
                        thresh.rf = 0.15,
                        thresh.perc = 0.05,
                        remove.fp = NULL,
                        n.clusters = 1)
{
    ## checking for correct object
    input_classes <-c("poly.est.two.pts.pairwise")
    if (!inherits(input.twopt, input_classes)) {
        stop(deparse(substitute(input.twopt)),
             " is not an object of class 'poly.est.two.pts.pairwise'")
    }

    ## Identifying false positives
    if (!is.null(remove.fp)){
        if (remove.fp < 0.0 || remove.fp > 0.5)
            stop("'remove.fp' value should be defined between 0.0 and 0.5.")
        rf_mat_fp = rf_list_to_matrix(input.twopt = input.twopt, n.clusters = n.clusters, verbose = FALSE)
        xo = apply(rf_mat_fp$rec.mat, 2, function(x) sum(x < remove.fp, na.rm = T))
        o2 = names(which(xo > 0.90*(length(xo))))
        o1 = !(names(xo) %in% o2)
    }

    ## Getting filtered rf matrix
    thresh.missing<-1-thresh.perc
    rf_mat<- rf_list_to_matrix(input.twopt = input.twopt, thresh.LOD.ph = thresh.LOD.ph,
                               thresh.LOD.rf = thresh.LOD.rf, thresh.rf = thresh.rf,
                               n.clusters = n.clusters, verbose = FALSE)

    ## Removing false positives
    if (!is.null(remove.fp))
        rf_mat$rec.mat = rf_mat$rec.mat[o1,o1]

    ## Removing markers that presented too many values above threshold 
    x <- apply(rf_mat$rec.mat, 1, function(x) sum(is.na(x)))
    o <- names(which(x < quantile(x, probs = thresh.missing)))

    ## Checking for remaining markers with rfs below threshold
    rem = names(x)[!(names(x) %in% o)]
    remaining = rf_mat$rec.mat[rem,rem]
    ev.mrks = rownames(remaining[which(remaining[,ncol(remaining)] < thresh.rf),])
    if (!is.null(remove.fp) && length(o2) > 0)
        cat('The following markers presented more than 90% of recombination fractions below', remove.fp, 'and were removed:', paste0(o2), '\n')
    if (length(ev.mrks) > 0)
        cat('The following markers were also removed, but presented one or more recombination fractions below', thresh.rf,':', paste0(ev.mrks), '\n')

    ## Returning sequence object
    ch_filt<-make_seq_mappoly(input.obj = get(input.twopt$data.name), arg = o, data.name = input.twopt$data.name)
    ch_filt
}
