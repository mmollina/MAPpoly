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
#' \dontrun{
#' require(polymapR)
#' data("screened_data3")
#' mappoly.data <- import_data_from_polymapR(screened_data3, 4)
#' plot(mappoly.data)
#'}
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
#' \dontrun{
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
#'  plot(mappoly.maplist[[1]], mrk.names = T, left.lim = 0, right.lim = 20, cex = .7)
#'  plot(mappoly.maplist[[2]])
#'  plot(mappoly.maplist[[3]])
#'  plot(mappoly.maplist[[4]])
#'  plot(mappoly.maplist[[5]])
#'  
#'  ## Computing conditional genotype probabilities
#'  genoprob0 <- lapply(mappoly.maplist, calc_genoprob, step = 1)
#'  
#'  ## Computing preferential pairing profiles
#'  pref.pair0 <- calc_prefpair_profiles(genoprob0)
#'  plot(pref.pair0, min.y.prof = .25, max.y.prof = 0.4, P = "P1", Q = "P2")
#'  
#'  ## Computing homolog probabilities
#'  h.prob0<-calc_homoprob(genoprob0)
#'  plot(h.prob0, ind = "F1_030") ## plot haplotype of individual "F1_030"
#'  
#'  #### Computing conditional genotype probabilities including error
#'  genoprob1 <- lapply(mappoly.maplist, calc_genoprob_error, step = 1, error = 0.05)
#'  
#'  ## Computing preferential pairing profiles
#'  pref.pair1 <- calc_prefpair_profiles(genoprob1)
#'  plot(pref.pair1, min.y.prof = .25, max.y.prof = 0.4, P = "P1", Q = "P2")
#'  
#'  ## Computing homolog probabilities
#'  h.prob1<-calc_homoprob(genoprob1)
#'  plot(h.prob1, ind = "F1_030") ## plot haplotype of individual "F1_030" 
#'  
#'  #### Reestimating recombination fractions using HMM
#'  cl <- parallel::makeCluster(5)
#'  parallel::clusterEvalQ(cl, require(mappoly))
#'  parallel::clusterExport(cl,  "mappoly.data")
#'  reest.maps <- parallel::parLapply(cl, mappoly.maplist, 
#'                                    est_full_hmm_with_global_error, 
#'                                    error = 0.05)
#'  parallel::stopCluster(cl)
#'  
#'  ## Computing conditional genotype probabilities
#'  genoprob2 <- lapply(reest.maps, calc_genoprob_error, step = 1, error = 0.05)
#'  
#'  ## Computing preferential pairing profiles
#'  pref.pair2 <- calc_prefpair_profiles(genoprob2)
#'  plot(pref.pair2, min.y.prof = .25, max.y.prof = 0.4, P = "P1", Q = "P2")
#'  
#'  ## Computing homolog probabilities
#'  h.prob2<-calc_homoprob(genoprob2)
#'  plot(h.prob2, ind = "F1_030") 
#'  
#'  #### Reconstructing the map using MAPpoly
#'  s <- make_seq_mappoly(mappoly.data, "all")
#'  tpt <- est_pairwise_rf(input.seq = s, n.clusters = 7)
#'  mat <- rf_list_to_matrix(make_pairs_mappoly(tpt, s))
#'  grs <- group_mappoly(input.mat = mat,
#'                       expected.groups = 5,
#'                       inter = TRUE)
#'  grs
#'  LG <- vector("list", 5)
#'  op <- par(mfrow = c(2,3))
#'  for(i in 1:5){
#'    s.temp <-  make_seq_mappoly(grs, arg = i)
#'    tpt.temp <- make_pairs_mappoly(tpt, s.temp)
#'    sf<-rf_snp_filter(input.twopt = tpt.temp, 
#'                      thresh.LOD.ph = 1, 
#'                      thresh.LOD.rf = 1, 
#'                      thresh.perc = 0.02)
#'    M <- make_mat_mappoly(input.mat = mat, sf)
#'    o <- mds_mappoly(M)
#'    so<-make_seq_mappoly(o)
#'    plot(M, ord = so$seq.mrk.names, main.text = paste("LG", i), index = FALSE)
#'    LG[[i]] <- list(s = so, tpt = tpt.temp)
#'    cat("\n")
#'  }
#'  par(op) 
#'  MAPs <- vector("list", 5)
#'  for(i in 1:5){
#'  MAPs[[i]] <- est_rf_hmm_sequential(input.seq = LG[[i]]$s,
#'                                     start.set = 6,
#'                                     thres.twopt = 10,
#'                                     thres.hmm = 50,
#'                                     extend.tail = 30,
#'                                     twopt = LG[[i]]$tpt,
#'                                     verbose = TRUE,
#'                                     tol = 10e-2,
#'                                     tol.final = 10e-4,
#'                                     phase.number.limit = 20,
#'                                     sub.map.size.diff.limit =  5,
#'                                     info.tail = TRUE,
#'                                     reestimate.single.ph.configuration = TRUE)
#'  }
#'  cl <- parallel::makeCluster(5)
#'  parallel::clusterEvalQ(cl, require(mappoly))
#'  parallel::clusterExport(cl,  "mappoly.data")
#'  recons.maps <- parallel::parLapply(cl, MAPs, 
#'                                     est_full_hmm_with_global_error, 
#'                                     error = 0.05)
#'  parallel::stopCluster(cl)
#'  
#'  ## Comparing resulting maps
#'  ## polymapR
#'  summary_maps(mappoly.maplist) 
#'  
#'  ## MAPpoly
#'  summary_maps(recons.maps) 
#'
#'  ## Computing conditional genotype probabilities
#'  genoprob3 <- lapply(recons.maps, 
#'                      calc_genoprob_error, 
#'                      step = 1, 
#'                      error = 0.05)
#'  
#'  ## Computing preferential pairing profiles
#'  pref.pair3 <- calc_prefpair_profiles(genoprob3)
#'  plot(pref.pair3, min.y.prof = .25, max.y.prof = 0.4, P = "P1", Q = "P2")
#'  
#'  ## Comparing homolog probabilities with different mapping approaches
#'  h.prob3<-calc_homoprob(genoprob3)
#'  ## plot haplotype of individual 10 (polymapR)
#'  plot(h.prob0, ind = "F1_030", use.plotly = FALSE) 
#'  ## plot haplotype of individual 10 (polymapR + HMM error modeling) 
#'  plot(h.prob1, ind = "F1_030", use.plotly = FALSE)  
#'  ## plot haplotype of individual 10 (reestimated: MAPpoly)
#'  plot(h.prob2, ind = "F1_030", use.plotly = FALSE) 
#'  ## plot haplotype of individual 10 (reconstructed: MAPpoly)
#'  plot(h.prob3, ind = "F1_030", use.plotly = FALSE) 
#'}
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