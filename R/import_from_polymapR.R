#' Import data from polymapR
#'
#' Function to import datasets from polymapR
#'
#' @param input.data  a \code{polymapR} dataset
#' @param ploidy the ploidy level     
#' @param parent1 name of parent 1
#' @param parent2 name of parent 2
#' @param filter.non.conforming if \code{TRUE} (default) exclude samples with non 
#'     expected genotypes under no double reduction
#'     
#' @examples
#' require(polymapR)
#' data("screened_data3")
#' mappoly.data <- import_data_from_polymapR(screened_data3, 4)
#' plot(mappoly.data)
#'
#' @author Marcelo Mollinari \email{mmollin@ncsu.edu}
#'
#' @references
#'     Bourke PM et al: (2019) PolymapR — linkage analysis and genetic map 
#'     construction from F1 populations of outcrossing polyploids. 
#'     _Bioinformatics_ 34:3496–3502.
#'     \url{https://doi.org/10.1093/bioinformatics/bty1002}
#' 
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \url{https://doi.org/10.1534/g3.119.400378}
#'     
#' @export import_data_from_polymapR
import_data_from_polymapR <- function(input.data, 
                                      ploidy, 
                                      parent1 = "P1", 
                                      parent2 = "P2", 
                                      filter.non.conforming = TRUE){
  geno.dose <- input.data[,-match(c(parent1, parent2), colnames(input.data)), drop = FALSE]
  mappoly.data <- structure(list(m = ploidy,
                                 n.ind = ncol(geno.dose),
                                 n.mrk = nrow(geno.dose),
                                 ind.names = colnames(geno.dose),
                                 mrk.names = rownames(geno.dose),
                                 dosage.p = input.data[,parent1],
                                 dosage.q = input.data[,parent2],
                                 sequence = NA,
                                 sequence.pos = NA,
                                 seq.ref = NULL,
                                 seq.alt = NULL,
                                 all.mrk.depth = NULL,
                                 prob.thres = NULL,
                                 geno.dose = geno.dose,
                                 nphen = 0,
                                 phen = NULL,
                                 kept = NULL,
                                 elim.correspondence = NULL),
                            class = "mappoly.data")
  if(filter.non.conforming){
    mappoly.data<-filter_non_conforming_classes(mappoly.data)
    Ds <- array(NA, dim = c(ploidy+1, ploidy+1, ploidy+1))
    for(i in 0:ploidy)
      for(j in 0:ploidy)
        Ds[i+1,j+1,] <- segreg_poly(m = ploidy, dP = i, dQ = j)
    Dpop<-cbind(mappoly.data$dosage.p, mappoly.data$dosage.q)
    M<-t(apply(Dpop, 1, function(x) Ds[x[1]+1, x[2]+1,]))
    dimnames(M)<-list(mappoly.data$mrk.names, c(0:ploidy))
    M<-cbind(M, mappoly.data$geno.dose)
    mappoly.data$chisq.pval<-apply(M, 1, mrk_chisq_test, m = ploidy)
  }
  mappoly.data
}


#' Import phased map list from polymapR
#'
#' Function to import phased map lists from polymapR
#' 
#' @param maplist a list of phased maps obtained using function 
#' \code{create_phased_maplist} from package \code{polymapR} 
#' @param mappoly.data a dataset used to obtain \code{maplist}, 
#' converted into class \code{mappoly.data}
#' @param ploidy the ploidy level     
#'     
#' @examples
#' require(polymapR)
#' ## Loading polymapR example
#' data("integrated.maplist", "screened_data3", "marker_assignments_P1","marker_assignments_P2")
#' maplist <- create_phased_maplist(maplist = integrated.maplist,
#'                                  dosage_matrix.conv = screened_data3,
#'                                  marker_assignment.1=marker_assignments_P1,
#'                                  marker_assignment.2=marker_assignments_P2,
#'                                  ploidy = 4)
#'  ## Importing polymapR dataset                                
#'  mappoly.data <- import_data_from_polymapR(screened_data3, 4)
#'  plot(mappoly.data) 
#'  
#'  ## Importing polymapR phased maplist
#'  mappoly.maplist <- import_phased_maplist_from_polymapR(maplist, mappoly.data)
#'  plot_map_list(mappoly.maplist)
#'  ## plot phased map
#'  plot(mappoly.maplist[[1]])
#'  ## plot a segment of phased map (from 0 to 20 cM)
#'  plot(mappoly.maplist[[1]], mrk.names = TRUE, left.lim = 0, right.lim = 20, cex = .7)
#'
#' @author Marcelo Mollinari \email{mmollin@ncsu.edu}
#'
#' @references
#'     Bourke PM et al: (2019) PolymapR — linkage analysis and genetic map 
#'     construction from F1 populations of outcrossing polyploids. 
#'     _Bioinformatics_ 34:3496–3502.
#'     \url{https://doi.org/10.1093/bioinformatics/bty1002}
#' 
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \url{https://doi.org/10.1534/g3.119.400378}
#'     
#' @export import_phased_maplist_from_polymapR
import_phased_maplist_from_polymapR <- function(maplist, 
                                                mappoly.data, 
                                                ploidy = NULL){
  input_classes <- c("list")
  if (!inherits(maplist, input_classes)) {
    stop(deparse(substitute(maplist)), " is not a list of phased maps.")
  }
  X <- maplist[[1]]
  if(is.null(ploidy))
    m <- (ncol(X)-2)/2
  MAPs <- vector("list", length(maplist))
  for(i in 1:length(MAPs)){
    X <- maplist[[i]]
    seq.num <- match(X$marker, mappoly.data$mrk.names)
    seq.rf <- mf_h(diff(X$position)) ## Using haldane
    seq.rf[seq.rf <= 1e-05] <- 1e-4
    P = ph_matrix_to_list(X[,3:(m+2)])
    Q = ph_matrix_to_list(X[,3:(m+2) + m])
    names(P) <- names(Q) <- seq.num
    seq.ph <- list(P = P, Q = Q)
    maps <- vector("list", 1)
    maps[[1]] <- list(seq.num = seq.num, seq.rf = seq.rf, seq.ph = seq.ph, loglike = 0)
    MAPs[[i]] <- structure(list(info = list(m = (ncol(X)-2)/2,
                                            n.mrk = nrow(X),
                                            seq.num = seq.num,
                                            mrk.names = as.character(X$marker),
                                            seq.dose.p = mappoly.data$dosage.p[as.character(X$marker)],
                                            seq.dose.q = mappoly.data$dosage.q[as.character(X$marker)],
                                            sequence = rep(i, nrow(X)),
                                            sequence.pos = NULL,
                                            seq.ref = NULL,
                                            seq.alt = NULL,
                                            chisq.pval = mappoly.data$chisq.pval[as.character(X$marker)],
                                            data.name = as.character(sys.call())[3], 
                                            ph.thresh = NULL),
                                maps = maps),
                           class = "mappoly.map")
    MAPs[[i]] <- loglike_hmm(MAPs[[i]], mappoly.data)
  }
   MAPs
}
