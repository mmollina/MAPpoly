#'  Subset pairwise recombination fractions
#'
#'  Get a subset of an object of class \code{mappoly.twopt} or \code{mappoly.twopt2} (i.e.
#'  recombination fraction) and LOD score statistics for all possible linkage
#'  phase combinations based on a sequence of markers.
#'
#' @param input.twopt an object of class \code{mappoly.twopt}
#'
#' @param input.seq an object of class \code{mappoly.sequence}, with 
#' a sequence of markers contained in \code{input.twopt}
#'
#' @return an object of class \code{mappoly.twopt} which is a
#'  subset of \code{input.twopt}. 
#'  See \code{\link[mappoly]{est_pairwise_rf}} for details
#'
#' @examples
#'     ## selecting some markers along the genome
#'     some.mrk <- make_seq_mappoly(hexafake, seq(1, 1500, 30))
#'     all.pairs <- est_pairwise_rf(input.seq = some.mrk)
#'     mat.full <- rf_list_to_matrix(input.twopt = all.pairs)
#'     plot(mat.full)
#'     
#'     ## selecting two-point information for chromosome 1
#'     mrks.1 <- make_seq_mappoly(hexafake, names(which(some.mrk$chrom == 1)))
#'     p1 <- make_pairs_mappoly(input.seq = mrks.1, input.twopt = all.pairs)
#'     m1 <- rf_list_to_matrix(input.twopt = p1)
#'     plot(m1, main.text = "LG1")
#'    
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \doi{10.1534/g3.119.400378} 
#'
#' @export

make_pairs_mappoly <- function(input.twopt, input.seq){
  ## checking for correct object
  input_classes1 <- c("mappoly.sequence")
  if (!inherits(input.seq, input_classes1)) {
    stop(deparse(substitute(input.seq)), " is not an object of class 'mappoly.sequence'")
  }
  if (inherits(input.twopt, "mappoly.twopt")){
    w <- combn(sort(input.seq$seq.num), 2)
    idx <- apply(w, 2, function(x) paste(sort(x), collapse = "-"))
    input.twopt$n.mrk <- length(input.seq$seq.num)
    input.twopt$seq.num <- sort(input.seq$seq.num)
    input.twopt$pairwise <- input.twopt$pairwise[idx]
    return(input.twopt)
  } 
  else if (inherits(input.twopt, "mappoly.twopt2")){
    id <- match(input.seq$seq.mrk.names, input.twopt$seq.mrk.names)
    input.twopt$n.mrk <- length(id)
    input.twopt$seq.mrk.names <- input.twopt$seq.mrk.names[id]
    input.twopt$chisq.pval.thres <- input.twopt$chisq.pval.thres
    input.twopt$chisq.pval <- input.twopt$chisq.pval[id]
    input.twopt$pairwise$rf <- input.twopt$pairwise$rf[id,id] 
    input.twopt$pairwise$LOD.rf <- input.twopt$pairwise$LOD.rf[id,id] 
    input.twopt$pairwise$LOD.ph <- input.twopt$pairwise$LOD.ph[id,id]
    input.twopt$pairwise$Sh.P1 <- input.twopt$pairwise$Sh.P1[id,id]
    input.twopt$pairwise$Sh.P2 <- input.twopt$pairwise$Sh.P2[id,id]
    return(input.twopt)
  } 
  else{
    stop(deparse(substitute(input.twopt)),
         " is not an object of class 'mappoly.twopt' or 'mappoly.twopt2'")
  }
}
