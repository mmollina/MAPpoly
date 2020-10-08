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
#' @param thresh.LOD.ph LOD score threshold for linkage phase
#'     configuration (default = 5)
#'
#' @param thresh.LOD.rf LOD score threshold for recombination fraction (default = 5)
#'
#' @param thresh.rf threshold for recombination fractions (default = 0.15) 
#'#'
#' @param thresh.perc threshold for the percentage of the pairwise marker
#'  combinations that should be considered in order to
#'  keep the marker. For example, \code{thresh.perc = 0.05} means that at least
#'  5\% of the pairwise combinations should be present in order to
#'  keep the marker (default = 0.05)
#'
#' @param remove.fp numeric value from 0.0 to 0.5 (default = NULL). When defined,
#' this parameter identifies and removes markers that presented more than 90\%
#' of its pairwise recombination fractions below \code{remove.fp} value throughout
#' the linkage group
#' 
#' @param ncpus number of parallel processes (i.e. cores) to spawn (default = 1)
#'
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
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
#'     tpt.filt<-rf_snp_filter(all.pairs, 2, 2, 0.07, thresh.perc = 0.1)
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
#'
rf_snp_filter<-function(input.twopt,
                        thresh.LOD.ph = 5,
                        thresh.LOD.rf = 5,
                        thresh.rf = 0.15,
                        thresh.perc = 0.05,
                        remove.fp = NULL,
                        ncpus = 1L,
                        verbose = TRUE)
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
        rf_mat_fp = rf_list_to_matrix(input.twopt = input.twopt, ncpus = ncpus, verbose = FALSE)
        xo = apply(rf_mat_fp$rec.mat, 2, function(x) sum(x < remove.fp, na.rm = T))
        o2 = names(which(xo > 0.90*(length(xo))))
        o1 = !(names(xo) %in% o2)
    }

    ## Getting filtered rf matrix
    thresh.missing<-1-thresh.perc
    rf_mat<- rf_list_to_matrix(input.twopt = input.twopt, thresh.LOD.ph = thresh.LOD.ph,
                               thresh.LOD.rf = thresh.LOD.rf, thresh.rf = thresh.rf,
                               ncpus = ncpus, verbose = FALSE)

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
    if (verbose && !is.null(remove.fp) && length(o2) > 0)
        cat('The following markers presented more than 90% of recombination fractions below', remove.fp, 'and were removed:', paste0(o2), '\n')
    if (verbose && length(ev.mrks) > 0)
        cat('The following markers were also removed, but presented one or more recombination fractions below', thresh.rf,':', paste0(ev.mrks), '\n')

    ## Returning sequence object
    ch_filt<-make_seq_mappoly(input.obj = get(input.twopt$data.name), arg = o, data.name = input.twopt$data.name)
    ch_filt
}
