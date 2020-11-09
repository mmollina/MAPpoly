#' Multipoint analysis using Hidden Markov Models in autopolyploids
#'
#' Performs the multipoint analysis proposed by \cite{Mollinari and
#' Garcia (2019)} in a sequence of markers
#'
#'  This function first enumerates a set of linkage phase configurations
#'  based on two-point recombination fraction information using a threshold
#'  provided by the user (argument \code{thresh}). After that, for each
#'  configuration, it reconstructs the genetic map using the
#'  HMM approach described in Mollinari and Garcia (2019). As result, it returns 
#'  the multipoint likelihood for each configuration in form of LOD Score comparing 
#'  each configuration to the most likely one. It is recommended to use a small number 
#'  of markers (e.g. 50 markers for hexaploids) since the possible linkage 
#'  phase combinations bounded only by the two-point information can be huge. 
#'  Also, it can be quite sensible to small changes in \code{'thresh'}. 
#'  For a large number of markers, please see \code{\link[mappoly]{est_rf_hmm_sequential}}.
#
#' @param input.seq an object of class \code{mappoly.sequence}
#'
#' @param input.ph an object of class \code{two.pts.linkage.phases}. 
#' If not available (default = NULL), it will be computed
#'
#' @param thres LOD Score threshold used to determine if the linkage phases
#'     compared via two-point analysis should be considered. Smaller 
#'     values will result in smaller number of linkage phase 
#'     configurations to be evaluated by the multipoint algorithm. 
#'
#' @param twopt an object of class \code{poly.est.two.pts.pairwise}
#'     containing two-point information
#'
#' @param verbose if \code{TRUE}, current progress is shown; if
#'     \code{FALSE} (default), no output is produced
#'
#' @param tol the desired accuracy (default = 1e-04)
#'
#' @param est.given.0.rf logical. If TRUE returns a map forcing all
#' recombination fractions equals to 0 (1e-5, for internal use only. 
#' Default = FALSE)
#'
#' @param reestimate.single.ph.configuration logical. If \code{TRUE}
#' returns a map without re-estimating the map parameters for cases
#' where there is only one possible linkage phase configuration. 
#' This argument is intended to be used in a sequential map construction
#' 
#' @param high.prec logical. If \code{TRUE} (default) uses high precision 
#' long double numbers in the HMM procedure
#' 
#' @param x an object of the class \code{mappoly.map}
#' 
#' @param detailed logical. if TRUE, prints the linkage phase configuration and the marker 
#' position for all maps. If FALSE (default), prints a map summary 
#' 
#' @param left.lim the left limit of the plot (in cM, default = 0). 
#' 
#' @param right.lim the right limit of the plot (in cM, default = Inf, i.e., 
#'                  will print the entire map)
#' 
#' @param phase logical. If \code{TRUE} (default) plots the phase configuration
#'  for both parents 
#' 
#' @param mrk.names if TRUE, marker names are displayed (default = FALSE)
#' 
#' @param cex The magnification to be used for marker names
#' 
#' @param config should be \code{'best'} or the position of the
#'     configuration to be plotted. If \code{'best'}, plot the configuration
#'     with the highest likelihood
#'
#' @param P a string containing the name of parent P
#' 
#' @param Q a string containing the name of parent Q
#' 
#' @param ... currently ignored
#'
#' @return A list of class \code{mappoly.map} with two elements: 
#' 
#' i) info:  a list containing information about the map, regardless of the linkage phase configuration:
#' \item{m}{the ploidy level}
#' \item{n.mrk}{number of markers}
#' \item{seq.num}{a vector containing the (ordered) indices of markers in the map, 
#'                according to the input file}
#' \item{mrk.names}{the names of markers in the map}
#' \item{seq.dose.p}{a vector containing the dosage in parent 1 for all markers in the map}
#' \item{seq.dose.q}{a vector containing the dosage in parent 2 for all markers in the map}
#' \item{sequence}{a vector indicating the sequence (usually chromosome) each marker belongs 
#'                 as informed in the input file. If not available, 
#'                 \code{sequence = NULL}}
#' \item{sequence.pos}{physical position (usually in megabase) of the markers into the sequence}
#' \item{seq.ref}{reference base used for each marker (i.e. A, T, C, G). If not available, 
#'                 \code{seq.ref = NULL}}                 
#' \item{seq.alt}{alternative base used for each marker (i.e. A, T, C, G). If not available, 
#'                 \code{seq.ref = NULL}}
#' \item{chisq.pval}{a vector containing p-values of the chi-squared test of Mendelian 
#'                   segregation for all markers in the map}                 
#' \item{data.name}{name of the dataset of class \code{mappoly.data}}
#' \item{ph.thres}{the LOD threshold used to define the linkage phase configurations to test}
#' 
#' ii) a list of maps with possible linkage phase configuration. Each map in the list is also a 
#'    list containing
#' \item{seq.num}{a vector containing the (ordered) indices of markers in the map, 
#'                according to the input file}
#' \item{seq.rf}{a vector of size (\code{n.mrk - 1}) containing a sequence of recombination 
#'               fraction between the adjacent markers in the map}
#' \item{seq.ph}{linkage phase configuration for all markers in both parents}
#' \item{loglike}{the hmm-based multipoint likelihood}
#'
#' @examples
#'     mrk.subset<-make_seq_mappoly(hexafake, 1:10)
#'     red.mrk<-elim_redundant(mrk.subset)
#'     unique.mrks<-make_seq_mappoly(red.mrk)
#'     subset.pairs<-est_pairwise_rf(input.seq = unique.mrks,
#'                                   ncpus = 1,
#'                                   verbose=TRUE)
#'
#'     ## Estimating subset map with a low tolerance for the E.M. procedure
#'     ## for CRAN testing purposes
#'     subset.map <- est_rf_hmm(input.seq = unique.mrks,
#'                              thres = 2,
#'                              twopt = subset.pairs,
#'                              verbose = TRUE,
#'                              tol = 0.1,
#'                              est.given.0.rf = FALSE)
#'     subset.map
#'     ## linkage phase configuration with highest likelihood
#'     plot(subset.map, mrk.names = TRUE, config = "best")
#'     ## the second one
#'     plot(subset.map, mrk.names = TRUE, config = 2)
#' 
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     https://doi.org/10.1534/g3.119.400378 
#' @rdname est_rf_hmm
#' @export est_rf_hmm
est_rf_hmm <- function(input.seq, input.ph = NULL,
                       thres = 0.5, twopt = NULL,
                       verbose = FALSE, 
                       tol = 1e-04,
                       est.given.0.rf=FALSE,
                       reestimate.single.ph.configuration = TRUE,
                       high.prec = TRUE) {
  ## checking for correct object
  input_classes <- c("mappoly.sequence", "two.pts.linkage.phases")
  if (!inherits(input.seq, input_classes[1])) {
    stop(deparse(substitute(input.seq)), " is not an object of class 'mappoly.sequence'")
  }
  if(length(input.seq$seq.num) == 1)
    stop("Input sequence contains only one marker.", call. = FALSE)
  if (is.null(input.ph)) {
    if (is.null(twopt))
      stop("Two-point analysis not found!")
    if (verbose)
      cat("\nListing all configurations under threshold", thres, "using two-point information...\n")
    input.ph <- ls_linkage_phases(input.seq = input.seq, thres = thres, twopt = twopt)
  }
  if (!inherits(input.ph, input_classes[2])) {
    stop(deparse(substitute(input.ph)), " is not an object of class 'two.pts.linkage.phases'")
  }
  if(length(input.seq$seq.num) == 2){
    maps<-vector("list", 1)
    res.temp <- twopt$pairwise[[paste(sort(input.seq$seq.num), collapse = "-")]]
    dp <- get(input.seq$data.name, pos =1)$dosage.p[input.seq$seq.num]
    dq <- get(input.seq$data.name, pos =1)$dosage.q[input.seq$seq.num]
    sh <- as.numeric(unlist(strsplit(rownames(res.temp)[1], split = "-")))
    res.ph <- list(P = c(list(0),list(0)), Q = c(list(0),list(0)))
    if(dp[1] != 0)
      res.ph$P[[1]] <- 1:dp[1]
    if(dq[1] != 0)
      res.ph$Q[[1]] <- 1:dq[1]
    if(dp[2] != 0)
    {
      v<-(rev(res.ph$P[[1]])[1]+1):(dp[2]+rev(res.ph$P[[1]])[1])
      res.ph$P[[2]]<-v-sh[1]
    }
    if(dq[2] != 0)
    {
      v<-(rev(res.ph$Q[[1]])[1]+1):(dq[2]+rev(res.ph$Q[[1]])[1])
      res.ph$Q[[2]]<-v-sh[2]
    }
    names(res.ph$P)<-names(res.ph$Q)<-input.seq$seq.num
    maps[[1]] <- list(seq.num = input.seq$seq.num,
                      seq.rf = res.temp[1,"rf"],
                      seq.ph = res.ph,
                      loglike = 0)
    return(structure(list(info = list(m = input.seq$m,
                                      n.mrk = length(input.seq$seq.num),
                                      seq.num = input.seq$seq.num,
                                      mrk.names = input.seq$seq.mrk.names,
                                      seq.dose.p = get(input.seq$data.name, pos =1)$dosage.p[input.seq$seq.mrk.names],
                                      seq.dose.q = get(input.seq$data.name, pos =1)$dosage.q[input.seq$seq.mrk.names],
                                      sequence = get(input.seq$data.name, pos =1)$sequence[input.seq$seq.mrk.names],
                                      sequence.pos = get(input.seq$data.name, pos =1)$sequence.pos[input.seq$seq.mrk.names],
                                      seq.ref = get(input.seq$data.name, pos =1)$seq.ref[input.seq$seq.mrk.names],
                                      seq.alt = get(input.seq$data.name, pos =1)$seq.alt[input.seq$seq.mrk.names],
                                      chisq.pval = get(input.seq$data.name, pos =1)$chisq.pval[input.seq$seq.mrk.names],
                                      data.name = input.seq$data.name,
                                      ph.thresh = abs(res.temp[2,1])),
                          maps = maps),
                     class = "mappoly.map"))
  }
  n.ph <- length(input.ph$config.to.test)
  ret.map.no.rf.estimation <- FALSE
  if(n.ph == 1 && !reestimate.single.ph.configuration)
    ret.map.no.rf.estimation <- TRUE
  #if (verbose) {
  #  cat("\n    Number of linkage phase configurations: ")
    #cat("\n---------------------------------------------\n|\n|--->")
  #}
  maps <- vector("list", n.ph)
  if (verbose){
    txt<-paste0("       ",n.ph, " phase(s): ")
    cat(txt)
  }
  for (i in 1:n.ph) {
    if (verbose) {
      if((i+1)%%40==0)
        cat("\n", paste0(rep(" ", nchar(txt)), collapse = ""), sep="")
      cat(". ")
    }
    if(est.given.0.rf){
      rf.temp <- rep(1e-5, length(input.seq$seq.num) - 1)
      tol=1
    }
    else{
      rf.temp <- rep(0.001, length(input.seq$seq.num) - 1)
    }
    ph <- list(P = input.ph$config.to.test[[i]]$P,
               Q = input.ph$config.to.test[[i]]$Q)
    maps[[i]] <- est_rf_hmm_single(input.seq = input.seq,
                                   input.ph.single = ph,
                                   rf.temp = rf.temp,
                                   tol = tol,
                                   verbose = FALSE,
                                   ret.map.no.rf.estimation = ret.map.no.rf.estimation, 
                                   high.prec = high.prec)
  }
  #cat("\n")
  id<-order(sapply(maps, function(x) x$loglike), decreasing = TRUE)
  maps<-maps[id]
  if (verbose)
    cat("\n")
  structure(list(info = list(m = input.seq$m,
                             n.mrk = length(input.seq$seq.num),
                             seq.num = input.seq$seq.num,
                             mrk.names = input.seq$seq.mrk.names,
                             seq.dose.p = get(input.seq$data.name, pos =1)$dosage.p[input.seq$seq.mrk.names],
                             seq.dose.q = get(input.seq$data.name, pos =1)$dosage.q[input.seq$seq.mrk.names],
                             sequence = get(input.seq$data.name, pos =1)$sequence[input.seq$seq.mrk.names],
                             sequence.pos = get(input.seq$data.name, pos =1)$sequence.pos[input.seq$seq.mrk.names],
                             seq.ref = get(input.seq$data.name, pos =1)$seq.ref[input.seq$seq.mrk.names],
                             seq.alt = get(input.seq$data.name, pos =1)$seq.alt[input.seq$seq.mrk.names],
                             chisq.pval = get(input.seq$data.name, pos =1)$chisq.pval[input.seq$seq.mrk.names],
                             data.name = input.seq$data.name, ph.thresh = input.ph$thres),
                 maps = maps),
            class = "mappoly.map")
}


#' Multipoint analysis using Hidden Markov Models: Sequential phase elimination
#'
#' Performs the multipoint analysis proposed by \cite{Mollinari and
#'  Garcia (2019)} in a sequence of markers removing unlikely phases
#' using sequential multipoint information.
#'
#' This function sequentially includes markers into a map given an
#' ordered sequence. It uses two-point information to eliminate
#' unlikely linkage phase configurations given \code{thres.twopt}. The
#' search is made within a window of size \code{extend.tail}. For the
#' remaining configurations, the HMM-based likelihood is computed and
#' the ones that pass the HMM threshold (\code{thres.hmm}) are eliminated.
#'
#' @param input.seq an object of class \code{mappoly.sequence}
#' 
#' @param twopt an object of class \code{poly.est.two.pts.pairwise}
#'     containing the two-point information
#'  
#' @param start.set number of markers to start the phasing procedure (default = 4)      
#'     
#' @param thres.twopt the LOD threshold used to determine if the linkage
#'     phases compared via two-point analysis should be considered 
#'     for the search space reduction (A.K.A. \eqn{\eta} in 
#'     \cite{Mollinari and Garcia (2019)}, default = 5)
#'     
#' @param thres.hmm the LOD threshold used to determine if the linkage
#'     phases compared via hmm analysis should be evaluated in the 
#'     next round of marker inclusion (default = 50)
#'     
#' @param extend.tail the length of the chain's tail that should
#'     be used to calculate the likelihood of the map. If \code{NULL} (default), 
#'     the function uses all markers positioned. Even if \code{info.tail = TRUE}, 
#'     it uses at least \code{extend.tail} as the tail length
#'     
#' @param phase.number.limit the maximum number of linkage phases of the sub-maps defined 
#'     by arguments \code{info.tail} and \code{extend.tail}. Default is 20. If the
#'     size exceeds this limit, the marker will not be inserted. If
#'     \code{Inf}, then it will insert all markers.
#'     
#' @param sub.map.size.diff.limit the maximum accepted length
#'     difference between the current and the previous sub-map defined 
#'     by arguments \code{info.tail} and \code{extend.tail}. If the
#'     size exceeds this limit, the marker will not be inserted. If
#'     \code{NULL}(default), then it will insert all markers.
#'     
#' @param info.tail if \code{TRUE} (default), it uses the complete informative tail
#'     of the chain (i.e. number of markers where all homologous 
#'     (\eqn{ploidy x 2}) can be distinguished) to calculate the map likelihood 
#'     
#'@param reestimate.single.ph.configuration logical. If \code{FALSE} (default)
#'     returns a map without re-estimating the map parameters in cases
#'     where there are only one possible linkage phase configuration
#'      
#' @param tol the desired accuracy during the sequential phase (default = 10e-02)
#'     
#' @param tol.final the desired accuracy for the final map (default = 10e-04)
#'     
#' @param verbose If \code{TRUE} (default), current progress is shown; if
#'     \code{FALSE}, no output is produced
#'     
#' @param detailed.verbose If \code{TRUE}, the expansion of the current 
#'     submap is shown; 
#'     
#' @param high.prec logical. If \code{TRUE} uses high precision 
#' (long double) numbers in the HMM procedure implemented in C++,
#' which can take a long time to perform (default = FALSE)
#' 
#' @return A list of class \code{mappoly.map} with two elements: 
#' 
#' i) info:  a list containing information about the map, regardless of the linkage phase configuration:
#' \item{m}{the ploidy level}
#' \item{n.mrk}{number of markers}
#' \item{seq.num}{a vector containing the (ordered) indices of markers in the map, 
#'                according to the input file}
#' \item{mrk.names}{the names of markers in the map}
#' \item{seq.dose.p}{a vector containing the dosage in parent 1 for all markers in the map}
#' \item{seq.dose.q}{a vector containing the dosage in parent 2 for all markers in the map}
#' \item{sequence}{a vector indicating the sequence (usually chromosome) each marker belongs 
#'                 as informed in the input file. If not available, 
#'                 \code{sequence = NULL}}
#' \item{sequence.pos}{physical position (usually in megabase) of the markers into the sequence}
#' \item{seq.ref}{reference base used for each marker (i.e. A, T, C, G). If not available, 
#'                 \code{seq.ref = NULL}}
#' \item{seq.alt}{alternative base used for each marker (i.e. A, T, C, G). If not available, 
#'                 \code{seq.ref = NULL}}
#' \item{chisq.pval}{a vector containing p-values of the chi-squared test of Mendelian 
#'                   segregation for all markers in the map}
#' \item{data.name}{name of the dataset of class \code{mappoly.data}}
#' \item{ph.thres}{the LOD threshold used to define the linkage phase configurations to test}
#' 
#' ii) a list of maps with possible linkage phase configuration. Each map in the list is also a 
#'    list containing
#' \item{seq.num}{a vector containing the (ordered) indices of markers in the map, 
#'                according to the input file}
#' \item{seq.rf}{a vector of size (\code{n.mrk - 1}) containing a sequence of recombination 
#'               fraction between the adjacent markers in the map}
#' \item{seq.ph}{linkage phase configuration for all markers in both parents}
#' \item{loglike}{the hmm-based multipoint likelihood}
#'
#' @examples
#'  \donttest{
#'     mrk.subset<-make_seq_mappoly(hexafake, 1:20)
#'     red.mrk<-elim_redundant(mrk.subset)
#'     unique.mrks<-make_seq_mappoly(red.mrk)
#'     subset.pairs<-est_pairwise_rf(input.seq = unique.mrks,
#'                                   ncpus = 1,
#'                                   verbose=TRUE)
#'     subset.map <- est_rf_hmm_sequential(input.seq = unique.mrks,
#'                                         thres.twopt = 5,
#'                                         thres.hmm = 10,
#'                                         extend.tail = 10,
#'                                         tol = 0.1,
#'                                         tol.final = 10e-3,
#'                                         phase.number.limit = 5,
#'                                         twopt = subset.pairs,
#'                                         verbose = TRUE)
#'      print(subset.map, detailed = TRUE)
#'      plot(subset.map)
#'      plot(subset.map, left.lim = 0, right.lim = 1, mrk.names = TRUE)
#'      plot(subset.map, phase = FALSE)
#'      
#'      ## Retrieving simulated linkage phase
#'      ph.P <- maps.hexafake[[1]]$maps[[1]]$seq.ph$P
#'      ph.Q <- maps.hexafake[[1]]$maps[[1]]$seq.ph$Q
#'      ## Estimated linkage phase
#'      ph.P.est <- subset.map$maps[[1]]$seq.ph$P
#'      ph.Q.est <- subset.map$maps[[1]]$seq.ph$Q
#'      compare_haplotypes(m = 6, h1 = ph.P[names(ph.P.est)], h2 = ph.P.est)
#'      compare_haplotypes(m = 6, h1 = ph.Q[names(ph.Q.est)], h2 = ph.Q.est)
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
#' @importFrom utils head
#' @importFrom cli rule
#' @export est_rf_hmm_sequential
est_rf_hmm_sequential<-function(input.seq,
                                twopt,
                                start.set = 4,
                                thres.twopt = 5,
                                thres.hmm = 50,
                                extend.tail = NULL,
                                phase.number.limit = 20,
                                sub.map.size.diff.limit = Inf,
                                info.tail = TRUE,
                                reestimate.single.ph.configuration = FALSE,
                                tol = 10e-2,
                                tol.final = 10e-4,
                                verbose = TRUE,
                                detailed.verbose = FALSE,
                                high.prec = FALSE)
{
  ## checking for correct object
  input_classes <- c("mappoly.sequence", "poly.est.two.pts.pairwise")
  if (!inherits(input.seq, input_classes[1])) {
    stop(deparse(substitute(input.seq)), 
         " is not an object of class 'mappoly.sequence'")
  }
  if (!inherits(twopt, input_classes[2])) {
    stop(deparse(substitute(twopt)), 
         " is not an object of class 'poly.est.two.pts.pairwise'")
  }
  ## Information about the map
  if(verbose) {
    cli::cat_line("Number of markers: ", length(input.seq$seq.num))
    msg("Initial sequence", line = 2)
  }
  if(is.null(extend.tail))
    extend.tail<-length(input.seq$seq.num)
  ##### Map in case of two markers #####
  if(length(input.seq$seq.num) == 2)
  {
    cur.map <- est_rf_hmm(input.seq = input.seq, thres = thres.twopt, 
                          twopt = twopt, tol = tol.final, high.prec = high.prec, 
                          verbose = verbose, 
                          reestimate.single.ph.configuration = reestimate.single.ph.configuration)
    return(filter_map_at_hmm_thres(cur.map, thres.hmm))
  } 
  ##### More than 3 markers ####
  ##Starting sequential algorithm for maps with more than 3 markers
  ## First step: test all possible phase configurations under 
  ## a 'thres.twopt' threshold for a sequence of size 'start.set'
  if(start.set > length(input.seq$seq.num))
    start.set <- length(input.seq$seq.num)
  if(verbose) cat(start.set, " markers...\n", sep = "")
  cte<-0
  rf.temp<-c(Inf, Inf)
  while(any(rf.temp >= sub.map.size.diff.limit))
  {
    cte<-cte+1
    if(length(rf.temp)==1) {
      ## cat("Impossible build a map using the given thresholds\n")
      ## stop()
      stop("Impossible build a map using the given thresholds\n")
    }
    if(verbose) cat(cli::symbol$bullet, "   Trying sequence:", cte:(start.set+cte-1), ":\n")
    cur.seq <- make_seq_mappoly(input.obj = get(input.seq$data.name, pos=1), na.omit(input.seq$seq.num[cte:(start.set+cte-1)]), data.name = input.seq$data.name)
    input.ph <- ls_linkage_phases(input.seq = cur.seq,
                                  thres = thres.twopt,
                                  twopt = twopt)
    cur.map <- est_rf_hmm(input.seq = cur.seq, thres = thres.twopt, 
                          twopt = twopt, tol = tol, high.prec = high.prec, 
                          verbose = verbose, input.ph = input.ph,
                          reestimate.single.ph.configuration = reestimate.single.ph.configuration)
    cur.map <- filter_map_at_hmm_thres(cur.map, thres.hmm)
    cur.map <- reest_rf(cur.map, tol=tol.final, verbose = FALSE)
    rf.temp<-imf_h(cur.map$maps[[1]]$seq.rf)
  }
  if(verbose)
    msg("Done with initial sequence", line = 2)
  if(start.set >= length(input.seq$seq.num)){
    if(verbose) cat("\nDone phasing", cur.map$info$n.mrk, "markers\nReestimating final recombination fractions.")
    return(reest_rf(cur.map, tol=tol.final))
  }
  ##### More than 'start.set' markers #####
  ## For maps with more markers than 'start.set', include
  ## next makres in a sequential fashion
  ct <- start.set+cte
  all.ph <- update_ph_list_at_hmm_thres(cur.map, thres.hmm)
  if(sub.map.size.diff.limit!=Inf & !reestimate.single.ph.configuration){
    if (verbose) message("Making 'reestimate.single.ph.configuration = TRUE' to use map expansion")
    reestimate.single.ph.configuration <- TRUE
  }
  while(ct <= length(input.seq$seq.num))
  {
    ## Number of multipoint phases evaluated in the 
    ## previous round of marker insertion
    hmm.phase.number<-length(cur.map$maps)
    ## Get the map tail containing sufficient markers to 
    ## distinguish among the m*2 possible homologs.
    if(info.tail)
      tail.temp <- get_full_info_tail(cur.map)
    input.ph <- NULL
    ## for each linkage phase inherited from HMM 
    ## analysis from previous round, evaluate the
    ## possible linkage phases for the marker inserted 
    ## in the current round using two-point information
    ## under the threshold of thres.twopt
    for(i in 1:length(all.ph$config.to.test))
    {
      seq.num <- NULL
      if(info.tail)
        seq.num <- tail.temp$maps[[i]]$seq.num 
      if(length(seq.num) < extend.tail)
        seq.num <- tail(cur.map$maps[[i]]$seq.num, extend.tail)
      ## If the tail do not contain the marker responsable for carrying 
      ## multiple linkage phases through the rounds of insertion,
      ## extend the tail so it contains that marker
      if(max(check_ls_phase(all.ph)) >= length(seq.num)){
        seq.num <- tail(cur.map$maps[[i]]$seq.num, max(check_ls_phase(all.ph)) + start.set)        
      }
      seq.cur <- make_seq_mappoly(input.obj = get(input.seq$data.name, pos=1),
                                  arg = seq.num,
                                  data.name = input.seq$data.name)
      all.ph.temp <- get_ph_list_subset(all.ph, seq.num, i)
      input.ph.temp <- ls_linkage_phases(input.seq = seq.cur,
                                         thres = thres.twopt,
                                         twopt = twopt,
                                         mrk.to.add = input.seq$seq.num[ct],
                                         prev.info = all.ph.temp)
      input.ph <- concatenate_ph_list(input.ph, input.ph.temp)
    }
    ## This is the number of linkage phases to be tested 
    ## combining HMM phases from previous round of mrk
    ## insertion and the linkage phases from the two-point
    ## information obtained in the current round
    twopt.phase.number<-length(input.ph$config.to.test)
    ## If this number is higher than phase.number.limit,
    ## procede to the next iteration
    if(length(input.ph$config.to.test) > phase.number.limit) {
      if(verbose)
        cat(crayon::italic$yellow(paste0(ct ,": not included (linkage phases)\n", sep = "")))
      ct <- ct + 1
      next()
    }
    ## Appending the marker to the numeric 
    ## sequence and makeing a new sequence
    seq.num <- c(seq.num, input.seq$seq.num[ct])
    seq.cur <- make_seq_mappoly(input.obj = get(input.seq$data.name, pos=1),
                                arg = seq.num,
                                data.name = input.seq$data.name)
    if(verbose)
      cat_phase(input.seq, input.ph, all.ph, ct, seq.num, twopt.phase.number, hmm.phase.number)
    ## Evaluation all maps using HMM
    cur.map.temp <- est_rf_hmm(input.seq = seq.cur,
                               input.ph = input.ph,
                               twopt = twopt,
                               tol = tol,
                               verbose = FALSE,
                               reestimate.single.ph.configuration = reestimate.single.ph.configuration,
                               high.prec = high.prec)
    ## Filtering linkage phase configurations under a HMM LOD threshold
    cur.map.temp <- filter_map_at_hmm_thres(cur.map.temp, thres.hmm)
    ## Gathering linkage phases of the current map, excluding the marker inserted 
    ## in the current round
    ph.new<-lapply(cur.map.temp$maps, function(x) list(P = head(x$seq.ph$P, -1), 
                                                       Q = head(x$seq.ph$Q, -1)))
    ## Gathering linkage phases of the previous map, excluding the marker inserted 
    ## in the current round
    ph.old<-lapply(cur.map$maps, function(x, id) list(P = x$seq.ph$P[id], 
                                                      Q = x$seq.ph$Q[id]), id = names(ph.new[[1]]$P))
    ## Check in which whole phase configurations the new 
    ## HMM tail should be appended
    MQ<-MP<-matrix(NA, length(ph.old), length(ph.new))
    for(j1 in 1:nrow(MP)){
      for(j2 in 1:ncol(MP)){
        MP[j1, j2]<-identical(ph.old[[j1]]$P, ph.new[[j2]]$P)
        MQ[j1, j2]<-identical(ph.old[[j1]]$Q, ph.new[[j2]]$Q)
      }
    }
    M<-which(MP & MQ, arr.ind = TRUE)
    colnames(M)<-c("old", "new")
    ## Checking for quality (map size and LOD Score)
    submap.length.old <- sapply(cur.map$maps, function(x) sum(imf_h(x$seq.rf)))
    last.dist.old <- sapply(cur.map$maps, function(x) tail(imf_h(x$seq.rf),1))
    submap.length.new <- sapply(cur.map.temp$maps, function(x) sum(imf_h(x$seq.rf)))
    last.dist.new <- sapply(cur.map.temp$maps, function(x) tail(imf_h(x$seq.rf),1))
    last.mrk.expansion <- submap.expansion <- numeric(nrow(M))
    for(j1 in 1:nrow(M)){
      submap.expansion[j1] <- submap.length.new[M[j1,2]] - submap.length.old[M[j1,1]]
      last.mrk.expansion[j1] <- last.dist.new[M[j1,2]] - last.dist.old[M[j1,1]]
    }
    LOD <- get_LOD(cur.map.temp, sorted = FALSE)
    if(sub.map.size.diff.limit!=Inf){
      selected.map <- submap.expansion < sub.map.size.diff.limit & LOD < thres.hmm
      if(detailed.verbose){
        x <- round(cbind(submap.length.new, submap.expansion, last.mrk.expansion),2)
        for(j1 in 1:nrow(x)){
          cat("    ", x[j1,1], ": (", x[j1,2],  "/", x[j1,3],")", sep = "")
          if(j1!=nrow(x)) cat("\n")
        }
      }
      if(detailed.verbose){
        if(all(!selected.map)) cat(paste0(crayon::red(cli::symbol$cross), "\n"))
        else cat(paste0(crayon::green(cli::symbol$tick), "\n"))
      }
      if(all(!selected.map)){
        if(verbose) cat(crayon::italic$yellow(paste0(ct ,": not included (map extension)\n", sep = "")))
        ct <- ct + 1
        next()
      }
    } else { selected.map <- LOD < thres.hmm }
    all.ph.temp <- update_ph_list_at_hmm_thres(cur.map.temp, Inf)
    cur.map.temp$maps <- cur.map.temp$maps[selected.map]
    cur.map <- cur.map.temp
    all.ph <- add_mrk_at_tail_ph_list(all.ph, all.ph.temp, M[selected.map,,drop=FALSE])
    ct <- ct + 1
  }
  ##### Reestimating final map ####
  ## Reestimating final map with higher tolerance
  seq.final <- make_seq_mappoly(input.obj = get(input.seq$data.name, pos=1), 
                                arg = as.numeric(names(all.ph$config.to.test[[1]]$P)), 
                                data.name = input.seq$data.name)
  #msg(paste("Done phasing", length(seq.final$seq.num), "markers"))
  if (verbose){
      msg("Reestimating final recombination fractions", line = 2)
      cat("")
  }
  final.map <- est_rf_hmm(input.seq = seq.final,
                          input.ph = all.ph,
                          twopt = twopt,
                          tol = tol.final,
                          verbose = FALSE,
                          reestimate.single.ph.configuration = TRUE,
                          high.prec = TRUE)
  if(verbose) {
    #cat("\n------------------------------------------")
    cat("Markers in the initial sequence: ", length(input.seq$seq.num), sep = "")
    cat("\nMapped markers                  : ", final.map$info$n.mrk, " (", 
        round(100*final.map$info$n.mrk/length(input.seq$seq.num),1) ,"%)\n", sep = "")
    msg("", line = 2)
  }
  return(final.map)
}

#' @rdname est_rf_hmm
#' @export 

print.mappoly.map <- function(x, detailed = FALSE, ...) {
  cat("This is an object of class 'mappoly.map'\n")
  cat("    Ploidy level:\t", x$info$m, "\n")
  cat("    No. individuals:\t", get(x$info$data.name)$n.ind, "\n")
  cat("    No. markers:\t", x$info$n.mrk, "\n")
  cat("    No. linkage phases:\t", length(x$maps), "\n")
  l <- sapply(x$maps, function(x) x$loglike)
  LOD <- round(l - max(l), 2)
  if(detailed)
  {
    for (j in 1:length(x$maps)) {
      cat("\n    ---------------------------------------------\n")
      cat("    Linkage phase configuration: ", j)
      cat("\n       log-likelihood:\t", x$map[[j]]$loglike)
      cat("\n       LOD:\t\t", format(round(LOD[j], 2), digits = 2))
      cat("\n\n")
      M <- matrix("|", x$info$n.mrk, x$info$m * 2)
      M <- cbind(M, "    ")
      M <- cbind(M, format(round(cumsum(c(0, imf_h(x$maps[[j]]$seq.rf))), 1), digits = 1))
      for (i in 1:x$info$n.mrk) {
        if (all(x$maps[[j]]$seq.ph$Q[[i]] != 0))
          M[i, c(x$maps[[j]]$seq.ph$P[[i]], x$maps[[j]]$seq.ph$Q[[i]] + x$info$m)] <- "o" else M[i, x$maps[[j]]$seq.ph$P[[i]]] <- "o"
      }
      M <- cbind(get(x$info$data.name)$mrk.names[x$maps[[j]]$seq.num], M)
      big.name<-max(nchar(M[,1]))
      format_name<-function(y, big.name){
        paste0(y, paste0(rep(" ", big.name-nchar(y)), collapse = ""))
      }
      M<-rbind(c("", letters[1:(ncol(M)-3)], "", ""), M)
      format(apply(M, 1, function(y) cat(c("\t", format_name(y[1], big.name), "\t", y[2:(x$info$m + 1)], rep(" ", 4), y[(x$info$m + 2):(x$info$m * 2 + 3)], "\n"), collapse = "")))
    }
  } else {
    cat("\n    ---------------------------------------------\n")
    cat("    Number of linkage phase configurations: ", length(x$maps))
    cat("\n    ---------------------------------------------\n")
    for (j in 1:length(x$maps)) {
      cat("    Linkage phase configuration: ", j)
      cat("\n       map length:\t", format(round(sum(imf_h(x$map[[j]]$seq.rf)), 2)))
      cat("\n       log-likelihood:\t", format(round(x$map[[j]]$loglike, 2)))
      cat("\n       LOD:\t\t", format(round(LOD[j], 2), digits = 2))
      cat("\n    ~~~~~~~~~~~~~~~~~~\n")
    }
  }
}

#' @rdname est_rf_hmm
#' @importFrom grDevices rgb
#' @importFrom graphics rect
#' @export 
plot.mappoly.map <- function(x, left.lim = 0, right.lim = Inf,
                             phase = TRUE, mrk.names = FALSE, 
                             cex = 1, config = "best", P = "Parent 1",
                             Q = "Parent 2", ...) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  if(phase){
    map.info <- prepare_map(x, config)
    if(any(map.info$ph.p=="B")){
      var.col <- c("black", "darkgray")
      var.col <- var.col[1:2]
      names(var.col) <- c("A", "B")
    } else {
      var.col <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
      names(var.col) <- c("A", "T", "C", "G")
    }
    m <- map.info$m
    if(m == 2) {
      d.col <- c(NA, "#1B9E77", "#D95F02")
    }else if(m == 4){
      d.col <- c(NA, "#1B9E77", "#D95F02", "#7570B3", "#E7298A")
    }else if(m == 6){
      d.col <- c(NA, "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02")
    }else if(m == 8){
      d.col <- c(NA, "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666") 
    } else d.col <- c(NA, gg_color_hue(m))
    names(d.col) <- 0:m
    d.col[1]<-NA
    x <- map.info$map
    lab <- names(x)
    zy <- seq(0, 0.5, length.out = m) + 1.5
    pp <- map.info$ph.p
    pq <- map.info$ph.q
    dp <- map.info$dp
    dq <- map.info$dq
    x1<-abs(left.lim - x)
    x2<-abs(right.lim - x)
    id.left<-which(x1==min(x1))[1]
    id.right<-rev(which(x2==min(x2)))[1]
    par(mai = c(1,0.15,0,0), mar=c(4.5,2,1,2))
    curx<-x[id.left:id.right]
    layout(mat =matrix(c(4,2,3, 1), ncol = 2), heights = c(2, 10), widths = c(1, 10))
    plot(x = curx,
         y = rep(.5,length(curx)),
         type = "n" , 
         ylim = c(.25, 4.5), 
         axes = FALSE, 
         xlab = "Distance (cM)", 
         ylab = "")
    lines(c(x[id.left], x[id.right]), c(.5, .5), lwd=15, col = "gray")
    points(x = curx,
           y = rep(.5,length(curx)),
           xlab = "", ylab = "", 
           pch = "|", cex=1.5, 
           ylim = c(0,2))
    axis(side = 1)
    #Parent 2
    x1 <- seq(x[id.left], x[id.right], length.out = length(curx))
    x.control <- diff(x1[1:2])/2
    if(length(x1) < 150)
      x.control <- x.control * .8
    if(length(x1) < 100)
      x.control <- x.control * .8
    if(length(x1) < 75)
      x.control <- x.control * .8
    if(length(x1) < 50)
      x.control <- x.control * .8
    if(length(x1) < 25)
      x.control <- x.control * .8
    for(i in 1:m)
    {
      lines(range(x1), c(zy[i], zy[i]), lwd=8, col = "lightgray")
      y1 <- rep(zy[i], length(curx))
      pal <- var.col[pq[id.left:id.right,i]]
      rect(xleft = x1 - x.control, ybottom = y1 -.05, xright = x1 + x.control, ytop = y1 +.05, col = pal, border = NA)
    }
    #connecting allelic variants to markers 
    for(i in 1:length(x1))
      lines(c(curx[i], x1[i]), c(0.575, zy[1]-.05), lwd=0.2)
    points(x = x1,
           y = zy[m]+0.05+dq[id.left:id.right]/20,
           col = d.col[as.character(dq[id.left:id.right])],
           pch = 19, cex = .7)
    #Parent 1
    zy<-zy + 1.1
    for(i in 1:m)
    {
      lines(range(x1), c(zy[i], zy[i]), lwd=8, col = "gray")
      y1 <- rep(zy[i], length(curx))
      pal <- var.col[pp[id.left:id.right,i]]
      rect(xleft = x1 - x.control, ybottom = y1 -.05, xright = x1 + x.control, ytop = y1 +.05, col = pal, border = NA)
    }
    points(x = x1,
           y = zy[m]+0.05+dp[id.left:id.right]/20,
           col = d.col[as.character(dp[id.left:id.right])],
           pch = 19, cex = .7)
   
    if(mrk.names)
      text(x = x1,
           y = rep(zy[m]+0.05+.3, length(curx)),
           labels = names(curx),
           srt=90, adj = 0, cex = cex)
    par(mar = c(4.5,1,1,0), xpd = TRUE) 
    plot(x = 0,
         y = 0,
         type = "n" ,
         axes = FALSE, 
         ylab = "",
         xlab = "",
         ylim = c(.25, 4.5))
    zy <- zy - 1.1
    mtext(text = Q, side = 4, at = mean(zy), line = -1, font = 4)
    for(i in 1:m)
      mtext(letters[(m+1):(2*m)][i], line = 0, at = zy[i], side = 4)
    zy <- zy + 1.1
    mtext(text = P, side = 4, at = mean(zy), line = -1, font = 4)
    for(i in 1:m)
      mtext(letters[1:m][i],  line = 0, at = zy[i], side = 4)
    par(mar = c(0,1,2,4), xpd=FALSE)
    plot(x = curx,
         y = rep(.5,length(curx)),
         type = "n" , 
         axes = FALSE, 
         xlab = "", 
         ylab = "")
    if(any(map.info$ph.p=="B")){
      legend("bottomleft", legend=c("A", "B"), 
             fill =c(var.col), title = "Variants",
             box.lty=0, bg="transparent", ncol = 2)
    } else {    
      legend("bottomleft", legend=c("A", "T", "C", "G", "-"), 
             fill =c(var.col, "white"), title = "Nucleotides",
             box.lty=0, bg="transparent", ncol = 4)
    }
    legend("bottomright", legend=names(d.col)[-1], title = "Doses" ,
           col = d.col[-1], ncol = m/2, pch = 19,
           box.lty=0)
  } else {
    plot_map_list(x) 
  }
}

#' prepare maps for plot 
#' @param void internal function to be documented
#' @keywords internal
prepare_map<-function(input.map, config = "best"){
  if (!inherits(input.map, "mappoly.map")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  }
  ## Choosing the linkage phase configuration
  LOD.conf <- get_LOD(input.map, sorted = FALSE)
  if(config == "best") {
    i.lpc <- which.min(LOD.conf)
  } else if(config == "all"){
    i.lpc <- seq_along(LOD.conf) } else if (config > length(LOD.conf)) {
      stop("invalid linkage phase configuration")
    } else i.lpc <- config
  ## Gathering marker positions
  map <- cumsum(imf_h(c(0, input.map$maps[[i.lpc]]$seq.rf)))
  names(map) <- input.map$info$mrk.names
  ## 
  ph.p <- ph_list_to_matrix(input.map$maps[[i.lpc]]$seq.ph$P, input.map$info$m)
  ph.q <- ph_list_to_matrix(input.map$maps[[i.lpc]]$seq.ph$Q, input.map$info$m)
  dimnames(ph.p) <- list(names(map), letters[1:input.map$info$m])
  dimnames(ph.q) <- list(names(map), letters[(1+input.map$info$m):(2*input.map$info$m)])
  if(is.null(input.map$info$seq.alt))
  {
    ph.p[ph.p==1] <- ph.q[ph.q==1] <- "A"
    ph.p[ph.p==0] <- ph.q[ph.q==0] <- "B"  
  } else {
    for(i in input.map$info$mrk.names){
      ph.p[i, ph.p[i,]==1] <- input.map$info$seq.alt[i]
      ph.p[i, ph.p[i,]==0] <- input.map$info$seq.ref[i]
      ph.q[i, ph.q[i,]==1] <- input.map$info$seq.alt[i]
      ph.q[i, ph.q[i,]==0] <- input.map$info$seq.ref[i]
    }
  }
  dp <- input.map$info$seq.dose.p
  dq <- input.map$info$seq.dose.q
  list(m = input.map$info$m, map = map, ph.p = ph.p, ph.q = ph.q, dp = dp, dq = dq)
}

#' Get the tail of a marker sequence up to the point where the markers
#' provide no additional information.
#'
#' @param void internal function to be documented
#' @keywords internal
get_full_info_tail <- function(input.obj, extend = NULL) {
  ## checking for correct object
  if(!inherits(input.obj, "mappoly.map"))
    stop(deparse(substitute(input.obj)), " is not an object of class 'mappoly.map''")
  if (!is.null(extend))
    if (extend > input.obj$info$n.mrk)
      return(input.obj)
  m <- input.obj$info$m
  i <- m/2
  while (i < input.obj$info$n.mrk) {
    wp <- ph_list_to_matrix(tail(input.obj$maps[[1]]$seq.ph$P, i), m)
    xp <- apply(wp, 2, paste, collapse = "-")
    wq <- ph_list_to_matrix(tail(input.obj$maps[[1]]$seq.ph$Q, i), m)
    xq <- apply(wq, 2, paste, collapse = "-")
    if (length(unique(xp)) == m && length(unique(xq)) == m)
      (break)()
    i <- i + 1
  }
  if (!is.null(extend))
    if (i < extend)
      i <- extend
  input.obj$info$n.mrk <- i
  for (j in 1:length(input.obj$maps)) {
    input.obj$maps[[j]]$loglike <- 0
    input.obj$maps[[j]]$seq.ph$P <- tail(input.obj$maps[[j]]$seq.ph$P, n = i)
    input.obj$maps[[j]]$seq.ph$Q <- tail(input.obj$maps[[j]]$seq.ph$Q, n = i)
    input.obj$maps[[j]]$seq.num <- tail(input.obj$maps[[j]]$seq.num, n = i)
    input.obj$maps[[j]]$seq.rf <- tail(input.obj$maps[[j]]$seq.rf, n = (i - 1))
  }
  return(input.obj)
}

#' remove maps under a certain threshold
#' @param void internal function to be documented
#' @keywords internal
#' @export filter_map_at_hmm_thres
filter_map_at_hmm_thres <- function(map, thres.hmm){
  map$info$ph.thresh <- NULL
  map$maps<-map$maps[get_LOD(map, sorted = FALSE) < thres.hmm]
  map
}

#' makes a phase list from map, selecting only 
#' configurations under a certain threshold
#' @param void internal function to be documented
#' @keywords internal
update_ph_list_at_hmm_thres <- function(map, thres.hmm){
  temp.map <- filter_map_at_hmm_thres(map, thres.hmm)
  config.to.test <- lapply(temp.map$maps, function(x) x$seq.ph)
  rf.vec <- t(sapply(temp.map$maps, function(x) x$seq.rf))
  rownames(rf.vec) <- names(config.to.test) <- paste("Conf", 1:length(config.to.test), sep = "-")
  structure(list(config.to.test = config.to.test, rec.frac = rf.vec, 
                 m = map$info$m, seq.num = map$maps[[1]]$seq.num, 
                 thres = map$info$ph.thresh, data.name = map$info$data.name, 
                 thres.hmm = thres.hmm),
            class = "two.pts.linkage.phases")
}

#' subset of a linkage phase list
#' @param void internal function to be documented
#' @keywords internal
get_ph_list_subset<-function(ph.list, seq.num, conf){
  config.to.test <- list(lapply(ph.list$config.to.test[[conf]], function(x, seq.num) x[as.character(seq.num)], seq.num))
  rf.vec <- ph.list$rec.frac[conf, , drop = FALSE]
  names(config.to.test) <- rownames(rf.vec)
  structure(list(config.to.test = config.to.test, rec.frac = rf.vec, 
                 m = ph.list$m, seq.num = ph.list$seq.num, 
                 thres = ph.list$ph.thresh, data.name = ph.list$data.name, 
                 thres.hmm = ph.list$thres.hmm),
            class = "two.pts.linkage.phases")
}

#' concatenate two linkage phase lists
#' @param void internal function to be documented
#' @keywords internal
concatenate_ph_list<-function(ph.list.1, ph.list.2){
  if(length(ph.list.1)==0)
    return(ph.list.2)
  config.to.test <- c(ph.list.1$config.to.test, ph.list.2$config.to.test)
  rf.vec <- rbind(ph.list.1$rec.frac, ph.list.2$rec.frac)
  rownames(rf.vec) <- names(config.to.test) <- paste("Conf", 1:length(config.to.test), sep = "-")
  structure(list(config.to.test = config.to.test, rec.frac = rf.vec, 
                 m = ph.list.1$m, seq.num = ph.list.1$seq.num, 
                 thres = ph.list.1$ph.thresh, data.name = ph.list.1$data.name, 
                 thres.hmm = ph.list.1$thres.hmm),
            class = "two.pts.linkage.phases")
}

#' add a single marker at the tail of a linkage phase list
#' @param void internal function 
#' @keywords internal
add_mrk_at_tail_ph_list <- function(ph.list.1, ph.list.2, cor.index){
  config.to.test <- vector("list", length = nrow(cor.index))
  for(i in 1:nrow(cor.index)){
    config.to.test[[i]] <- list(P = c(ph.list.1$config.to.test[[cor.index[i,1]]]$P, tail(ph.list.2$config.to.test[[cor.index[i,2]]]$P,1)),
                                Q = c(ph.list.1$config.to.test[[cor.index[i,1]]]$Q, tail(ph.list.2$config.to.test[[cor.index[i,2]]]$Q,1)))
  }
  structure(list(config.to.test = config.to.test, rec.frac = NULL, 
                 m = ph.list.1$m, seq.num = ph.list.1$seq.num, 
                 thres = ph.list.1$ph.thresh, data.name = ph.list.1$data.name, 
                 thres.hmm = ph.list.1$thres.hmm),
            class = "two.pts.linkage.phases")
}

#' Compare a list of linkage phases and return the 
#' markers for which they are different.
#' @param void internal function to be documented
#' @keywords internal
check_ls_phase<-function(ph){
  if(length(ph$config.to.test) == 1) return(0)
  id <- rep(1, length(ph$config.to.test[[1]]$P))
  for(i in 1:length(ph$config.to.test[[1]]$P)){
    w <- ph$config.to.test[[1]]$P[i]
    for(j in 2:length(ph$config.to.test))
      id[i] <- id[i] + identical(w, ph$config.to.test[[j]]$P[i])
  }
  names(id) <- names((ph$config.to.test[[1]]$P))
  w <- length(ph$config.to.test[[1]]$P) - which(id < length(ph$config.to.test)) 
  if(length(w)==0) return(0)
  w
}

#' cat for graphical representation of the phases
#' @param void internal function to be documented
#' @keywords internal
print_ph<-function(input.ph){
  phs.P<-lapply(input.ph$config.to.test, 
                function(x, m) {
                  M <- matrix("|", nrow = 1, ncol = m) 
                  M[unlist(tail(x$P, 1))] <- crayon::red(cli::symbol$bullet)
                  paste(M, collapse = "")}, 
                m = input.ph$m) 
  phs.Q<-lapply(input.ph$config.to.test, 
                function(x, m) {
                  M <- matrix("|", nrow = 1, ncol = m) 
                  M[unlist(tail(x$Q, 1))] <- crayon::cyan(cli::symbol$bullet)
                  paste(M, collapse = "")}, 
                m = input.ph$m) 
  if(length(phs.P) == 1)
    return(paste(unlist(phs.P)[1], unlist(phs.Q)[1], "                   "))
  if(length(phs.P) == 2)
    return(paste(unlist(phs.P)[1], unlist(phs.Q)[1], "     ", unlist(phs.P)[2], unlist(phs.Q)[2]))
  if(length(phs.P) > 2)
    return(paste(unlist(phs.P)[1], unlist(phs.Q)[1], " ... ", unlist(phs.P)[2], unlist(phs.Q)[2]))
}

#' cat for phase information
#' @param void internal function to be documented
#' @keywords internal
cat_phase <- function(input.seq,
                      input.ph,
                      all.ph,
                      ct,
                      seq.num,
                      twopt.phase.number,
                      hmm.phase.number){
  pc <- round(100*ct/length(input.seq$seq.num),1)
  xmax11 <- nchar(length(input.seq$seq.num))
  x11 <- ct
  x11 <- paste0(x11, paste(rep(" ", xmax11-nchar(x11)), collapse = ""))
  x12 <- length(all.ph$config.to.test[[1]]$P) + 1
  x12 <- paste0(x12, paste(rep(" ", xmax11-nchar(x12)), collapse = ""))
  x13 <- pc
  x13 <- paste0(":(",pc,"%)", paste(rep(" ", 5-nchar(x13)), collapse = ""))
  x1 <- paste0(x12,"/",x11, x13)
  xmax21 <- nchar(max(input.seq$seq.num)) 
  x21 <- input.seq$seq.num[ct]
  x21 <- paste0(x21, paste(rep(" ", xmax21-nchar(x21)), collapse = ""))
  x22 <- length(input.ph$config.to.test)
  x22 <- paste0(x22, " ph ", paste(rep(" ", 5-nchar(x22)), collapse = ""))
  x23 <- paste0("(", hmm.phase.number, "/", twopt.phase.number,")  ")
  x23 <- paste0(x23, paste(rep(" ", 10-nchar(x23)), collapse = ""))
  x2 <- paste0( x21, ": ", x22, x23)
  x31 <- length(seq.num)-1
  x31 <- paste0(x31, paste(rep(" ", xmax11-nchar(x31)), collapse = ""))
  x3 <- paste0(" -- tail: ", x31)
  cat(x1, x2, x3, print_ph(input.ph), "\n")
}
