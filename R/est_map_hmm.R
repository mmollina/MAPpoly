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
#'  For higher number of markers, please see \code{\link[mappoly]{est_rf_hmm_sequential}}.
#
#' @param input.seq an object of class \code{mappoly.sequence}
#'
#' @param input.ph an object of class \code{two.pts.linkage.phases}. 
#' If not available (default = NULL), it will be computed
#'
#' @param thres the LOD threshold used to determine if the linkage phases
#'     compared via two-point analysis should be considered. Smaller 
#'     values will result in smaller number of linkage phase 
#'     configurations
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
#' recombination fractions equals to 0 (1e-5, for internal use only. Default = FALSE)
#'
#' @param reestimate.single.ph.configuration logical. If \code{TRUE}
#' returns a map without reestimating the map parameters for cases
#' where there is only one possible linkage phase configuration. 
#' This argument is intended to be used in a sequential map contruction
#' 
#' @param high.prec logical. If \code{TRUE} (default) uses high precision 
#' long double numbers in the HMM procedure
#' 
#' @param x an object of the class \code{mappoly.map}
#' 
#' @param detailed logical. if TRUE, prints the linkage phase configuration and the marker 
#' position for all maps. if FALSE (default), prints a map summary 
#'
#' @param phase logical. If \code{TRUE} (default) plots the phase configuration
#'  for both parents 
#'
#' @param col.cte.P a single character string or a vector of strings with size 
#' equal to the number of markers in the map indicating the color of the allelic
#'  variants for parent P (default = 'red')
#' 
#' @param col.cte.Q a single character string or a vector of strings with size 
#' equal to the number of markers in the map indicating the color of the allelic
#'  variants for parent Q (default = 'blue')
#' 
#' @param mrk.names if TRUE, marker names are displayed (default = FALSE)
#' 
#' @param cex The magnification to be used for marker names
#' 
#' @param config should be \code{'best'} or the position of the
#'     configuration to be plotted. If \code{'best'}, plot the configuration
#'     with the highest likelihood
#'
#' @param ... currently ignored
#'
#' @return An object of class \code{mappoly.map} with the following structure:
#' \item{m}{the ploidy level}
#' \item{mrk.names}{the names of markers present in the sequence}
#' \item{data.name}{name of the dataset of class \code{mappoly.data}}
#' \item{ph.thres}{the LOD threshold used to define the linkage phase configurations to test}
#' \item{maps}{a list containing the sequence of markers, their recombination fractions,
#' the linkage phase configuration for all markers in both parents P and Q and the 
#' map's joint likelihood}
#'
#' @examples
#'  \dontrun{
#'     data(hexafake)
#'     mrk.subset<-make_seq_mappoly(hexafake, 1:50)
#'     red.mrk<-elim_redundant(mrk.subset)
#'     unique.mrks<-make_seq_mappoly(red.mrk)
#'     counts.web<-cache_counts_twopt(unique.mrks, get.from.web = TRUE)
#'     subset.pairs<-est_pairwise_rf(input.seq = unique.mrks,
#'                                   count.cache = counts.web,
#'                                   n.clusters = 1,
#'                                   verbose=TRUE)
#'
#'     ## Estimating subset map with a low tolerance for the E.M. procedure
#'     subset.map <- est_rf_hmm(input.seq = unique.mrks,
#'                              thres = 2,
#'                              twopt = subset.pairs,
#'                              verbose = TRUE,
#'                              tol = 0.1,
#'                              est.given.0.rf = FALSE)
#'
#'     ## Re-estimating the map with the most likely configuration
#'     subset.map1 <- est_rf_hmm_single(input.seq = unique.mrks,
#'                                     input.ph.single = subset.map$maps[[1]]$seq.ph,
#'                                     tol = 10e-3,
#'                                     verbose = TRUE)
#'
#'     subset.map$maps[[1]]$seq.ph <- subset.map1$seq.ph
#'     
#'     plot(subset.map)
#'
#'      ## Retrieving simulated linkage phase
#'      ph.P <- maps.hexafake[[1]]$maps[[1]]$seq.ph$P
#'      ph.Q <- maps.hexafake[[1]]$maps[[1]]$seq.ph$Q
#'
#'      ## Estimated linkage phase
#'      ph.P.est <- subset.map$maps[[1]]$seq.ph$P
#'      ph.Q.est <- subset.map$maps[[1]]$seq.ph$Q
#'
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
#'     https://doi.org/10.1534/g3.119.400378 
#'
#' @export est_rf_hmm
#'
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
                                      mrk.names = get(input.seq$data.name, pos =1)$mrk.names[input.seq$seq.num],
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
                             mrk.names = get(input.seq$data.name, pos =1)$mrk.names[input.seq$seq.num],
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
#' @param thres.hmm the threshold used to determine if the linkage
#'     phases compared via hmm analysis should be evaluated in the 
#'     next round of marker inclusion (default = 50)
#'     
#' @param extend.tail the length of the chain's tail that should
#'     be used to calculate the likelihood of the map. If \code{NULL} (default), 
#'     the function uses all markers positioned. Even if \code{info.tail = TRUE}, 
#'     it uses at least \code{extend.tail} as the tail length
#'     
#' @param phase.number.limit the maximum number of linkage phases of the sub-maps defined 
#'     by arguments \code{info.tail} and \code{extend.tail}. If the
#'     size exceeds this limit, the marker will not be inserted. If
#'     \code{NULL}, then it will insert all markers (default = Inf)
#'     
#' @param sub.map.size.diff.limit the maximum accepted length
#'     difference between the current and the previous sub-map defined 
#'     by arguments \code{info.tail} and \code{extend.tail}. If the
#'     size exceeds this limit, the marker will not be inserted. If
#'     \code{NULL}, then it will insert all markers (default = Inf)
#'     
#' @param info.tail if \code{TRUE} (default), it uses the complete informative tail
#'     of the chain (i.e. number of markers where all homologous 
#'     (\eqn{ploidy x 2}) can be distinguished) to calculate the map likelihood 
#'     
#'@param reestimate.single.ph.configuration logical. If \code{FALSE} (default)
#'     returns a map without reestimating the map parameters in cases
#'     where there are only one possible linkage phase configuration
#'      
#' @param tol the desired accuracy during the sequential phase (default = 10e-02)
#'     
#' @param tol.final the desired accuracy for the final map (default = 10e-04)
#'     
#' @param verbose If \code{TRUE} (default), current progress is shown; if
#'     \code{FALSE}, no output is produced
#'     
#' @param high.prec logical. If \code{TRUE} uses high precision 
#' (long double) numbers in the HMM procedure implemented in C++,
#' which can take a long time to perform (default = FALSE)
#' 
#' @return An object of class \code{mappoly.map} with the following structure:
#' \item{m}{the ploidy level}
#' \item{mrk.names}{the names of markers present in the sequence}
#' \item{data.name}{name of the dataset of class \code{mappoly.data}}
#' \item{ph.thres}{the LOD threshold used to define the linkage phase configurations to test}
#' \item{maps}{a list containing the sequence of markers, their recombination fractions,
#' the linkage phase configuration for all markers in both parents P and Q and the 
#' map's joint likelihood}
#'
#' @examples
#'  \dontrun{
#'     data(hexafake)
#'     mrk.subset<-make_seq_mappoly(hexafake, 1:50)
#'     red.mrk<-elim_redundant(mrk.subset)
#'     unique.mrks<-make_seq_mappoly(red.mrk)
#'     counts.web<-cache_counts_twopt(unique.mrks, get.from.web = TRUE)
#'     subset.pairs<-est_pairwise_rf(input.seq = unique.mrks,
#'                                   count.cache = counts.web,
#'                                   n.clusters = 1,
#'                                   verbose=TRUE)
#'     subset.map <- est_rf_hmm_sequential(input.seq = unique.mrks,
#'                                         thres.twopt = 5,
#'                                         thres.hmm = 10,
#'                                         extend.tail = 10,
#'                                         tol = 0.1,
#'                                         tol.final = 10e-3,
#'                                         twopt = subset.pairs,
#'                                         verbose = TRUE)
#'      print(subset.map, detailed = TRUE)
#'      plot(subset.map)
#'      plot(subset.map, phase = FALSE)
#'      
#'      ## Retrieving simulated linkage phase
#'      ph.P <- maps.hexafake[[1]]$maps[[1]]$seq.ph$P
#'      ph.Q <- maps.hexafake[[1]]$maps[[1]]$seq.ph$Q
#'      ## Estimated linkage phase
#'      ph.P.est <- subset.map$maps[[1]]$seq.ph$P
#'      ph.Q.est <- subset.map$maps[[1]]$seq.ph$Q
#'      ##Notice that two estimated homologous in parent P are different
#'      ##from the simulated ones
#'      compare_haplotypes(m = 6, h1 = ph.P[names(ph.P.est)], h2 = ph.P.est)
#'      compare_haplotypes(m = 6, h1 = ph.Q[names(ph.Q.est)], h2 = ph.Q.est)
#'      
#'      
#'      ############# Autotetraploid example
#'     data(tetra.solcap)
#'     s1<-make_seq_mappoly(tetra.solcap, 'seq1')
#'     red.mrk<-elim_redundant(s1)
#'     s1.unique.mrks<-make_seq_mappoly(red.mrk)
#'     counts.web<-cache_counts_twopt(s1.unique.mrks, get.from.web = TRUE)
#'     s1.pairs<-est_pairwise_rf(input.seq = s1.unique.mrks,
#'                                   count.cache = counts.web,
#'                                   n.clusters = 10,
#'                                   verbose=TRUE)
#'     unique.gen.ord<-get_genomic_order(s1.unique.mrks)
#'     ## Selecting a subset of 100 markers at the beginning of chromosome 1 
#'     s1.gen.subset<-make_seq_mappoly(tetra.solcap, rownames(unique.gen.ord)[1:100])
#'     s1.gen.subset.map <- est_rf_hmm_sequential(input.seq = s1.gen.subset,
#'                                                start.set = 10,
#'                                                thres.twopt = 10, 
#'                                                thres.hmm = 10,
#'                                                extend.tail = 50,
#'                                                info.tail = TRUE, 
#'                                                twopt = s1.pairs,
#'                                                sub.map.size.diff.limit = 5, 
#'                                                phase.number.limit = 40,
#'                                                reestimate.single.ph.configuration = TRUE,
#'                                                tol = 10e-3,
#'                                                tol.final = 10e-5)
#'      print(s1.gen.subset.map, detailed = TRUE)
#'      plot(s1.gen.subset.map)
#'      plot(s1.gen.subset.map, phase = FALSE)
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
                                phase.number.limit = Inf,
                                sub.map.size.diff.limit = Inf,
                                info.tail = TRUE,
                                reestimate.single.ph.configuration = FALSE,
                                tol = 10e-2,
                                tol.final = 10e-4,
                                verbose = TRUE,
                                high.prec = FALSE)
{
  ## checking for correct object
  if (!any(class(input.seq) == "mappoly.sequence"))
    stop(deparse(substitute(input.seq)),
         " is not an object of class 'mappoly.sequence'")
  if (!any(class(twopt) == "poly.est.two.pts.pairwise"))
    stop(deparse(substitute(twopt)),
         " is not an object of class 'poly.est.two.pts.pairwise'")
  ## Information about the map
  if(verbose) {
    cli::cat_line("Number of markers: ", length(input.seq$seq.num))
    msg("Initial sequence", line = 2)
  }
  if(is.null(extend.tail))
    extend.tail<-length(input.seq$seq.num)
  #####
  ## Map in case of two markers
  #####
  if(length(input.seq$seq.num) == 2)
  {
    cur.map <- est_rf_hmm(input.seq = input.seq, thres = thres.twopt, 
                          twopt = twopt, tol = tol.final, high.prec = high.prec, 
                          verbose = verbose, 
                          reestimate.single.ph.configuration = reestimate.single.ph.configuration)
    return(filter_map_at_hmm_thres(cur.map, thres.hmm))
  } 
  #####
  ##Starting sequential algorithm for maps with more than 3 markers
  ## Fisrt step: test all possible phase configurations under 
  ## a 'thres.twopt' threshold for a sequence of size 'start.set'
  #####
  if(start.set > length(input.seq$seq.num))
    start.set <- length(input.seq$seq.num)
  if(verbose) cat(start.set, " markers...\n", sep = "")
  cte<-0
  rf.temp<-c(Inf, Inf)
  while(any(rf.temp >= sub.map.size.diff.limit))
  {
    cte<-cte+1
    if(length(rf.temp)==1) {
      cat("Impossible build a map using the given thresholds\n")
      stop()
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
  #####
  ## For maps with more markers than 'start.set', include
  ## next makres in a sequential fashion
  #####
  ct <- start.set+cte
  all.ph <- update_ph_list_at_hmm_thres(cur.map, thres.hmm)
  if(sub.map.size.diff.limit!=Inf & !reestimate.single.ph.configuration){
    message("Making 'reestimate.single.ph.configuration = TRUE' to use map expansion")
    reestimate.single.ph.configuration <- TRUE
  }
  while(ct <= length(input.seq$seq.num))
  {
    hmm.phase.number<-length(cur.map$maps)
    if(info.tail)
      tail.temp <- get_full_info_tail(cur.map)
    input.ph <- NULL
    for(i in 1:length(all.ph$config.to.test))
    {
      seq.num <- NULL
      if(info.tail)
        seq.num <- tail.temp$maps[[i]]$seq.num 
      if(length(seq.num) < extend.tail)
        seq.num <- tail(cur.map$maps[[i]]$seq.num, extend.tail)
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
    twopt.phase.number<-length(input.ph$config.to.test)
    if(length(input.ph$config.to.test) > phase.number.limit) {
      if(verbose)
        cat(crayon::red(paste0(ct ,": not included (too many linkage phases)\n", sep = "")))
      ct <- ct + 1
      next()
    }
    seq.num <- c(seq.num, input.seq$seq.num[ct])
    seq.cur <- make_seq_mappoly(input.obj = get(input.seq$data.name, pos=1),
                                arg = seq.num,
                                data.name = input.seq$data.name)
    if(verbose)
    {
      pc <- round(100*ct/length(input.seq$seq.num),1)
      cat(ct,": ", pc,"% (", input.seq$seq.num[ct] ,"): ", length(input.ph$config.to.test), 
          " ph(s) : (", hmm.phase.number, "/", twopt.phase.number,")", 
          "--t: ", length(seq.num),"\n",
          sep = "")
    }
    cur.map.temp <- est_rf_hmm(input.seq = seq.cur,
                               input.ph = input.ph,
                               twopt = twopt,
                               tol = tol,
                               verbose = FALSE,
                               reestimate.single.ph.configuration = reestimate.single.ph.configuration,
                               high.prec = high.prec)
    cur.map.temp <- filter_map_at_hmm_thres(cur.map.temp, thres.hmm)
    ph.new<-lapply(cur.map.temp$maps, function(x) list(P = head(x$seq.ph$P, -1), 
                                                       Q = head(x$seq.ph$Q, -1)))
    ph.old<-lapply(cur.map$maps, function(x, id) list(P = x$seq.ph$P[id], 
                                                      Q = x$seq.ph$Q[id]), id = names(ph.new[[1]]$P))
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
      if(verbose){
        x <- round(cbind(submap.length.new, submap.expansion, last.mrk.expansion),2)
        for(j1 in 1:nrow(x)){
          cat("    ", x[j1,1], ": (", x[j1,2],  "/", x[j1,3],")", sep = "")
          if(j1!=nrow(x)) cat("\n")
        }
      }
      if(verbose){
        if(all(!selected.map)) cat(paste0(crayon::red(cli::symbol$cross), "\n"))
        else cat(paste0(crayon::green(cli::symbol$tick), "\n"))
      }
      if(all(!selected.map)){
        if(verbose) cat(crayon::red(paste0(ct ,": not included (map expansion)\n", sep = "")))
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
  #####
  ## Reestimating final map with higher tolerance
  #####
  seq.final <- make_seq_mappoly(input.obj = get(input.seq$data.name, pos=1), 
                                arg = as.numeric(names(all.ph$config.to.test[[1]]$P)), 
                                data.name = input.seq$data.name)
  #msg(paste("Done phasing", length(seq.final$seq.num), "markers"))
  msg("Reestimating final recombination fractions", line = 2)
  cat("")
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
    cat("\nMaped markers                  : ", final.map$info$n.mrk, " (", 
        round(100*final.map$info$n.mrk/length(input.seq$seq.num),1) ,"%)\n", sep = "")
    msg("", line = 2)
  }
  return(final.map)
}

#' @rdname est_rf_hmm
#' @keywords internal
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
#' @importFrom grid grid.roundrect  grid.rect grid.newpage grid.lines viewport pushViewport upViewport grid.text gpar unit
#' @importFrom grDevices rgb
#' @keywords internal
#' @export
plot.mappoly.map <- function(x,
                             phase = TRUE,
                             col.P = "#e41a1c",
                             col.Q = "#377eb8",
                             mrk.names = FALSE, 
                             cex = 1, 
                             config = "best", ...) {
  
  if(phase){
    grid::grid.newpage()
    m <- x$info$m
    vp1 <- grid::viewport(x = 0, y = 0.8, width = 1, height = 0.15, just = c("left", "bottom"))
    vp2 <- grid::viewport(x = 0, y = 0.6, width = 1, height = 0.15, just = c("left", "bottom"))
    vp3 <- grid::viewport(x = 0.04, y = 0.4, width = 0.92, height = 0.2, just = c("left", "bottom"))
    grid::pushViewport(vp1)
    draw_homologous(m, y.pos = 0.1, h.names = letters[1:m], parent = "P1")
    if (config == "best")
      config <- which.min(abs(get_LOD(x)))
    if (!is.numeric(config))
      stop("config must be numeic")
    if (length(x$maps) < config) {
      warning("config should be equal or less than", length(x$maps), "\nplotting the best configuration")
      config <- which.min(abs(get_LOD(x, sorted = FALSE)))
    }
    P <- x$maps[[config]]$seq.ph$P
    if(length(col.P)==1)
      col.P<-rep(col.P, length(P))
    y <- seq(from = 0.05, to = 0.95, length.out = length(P))
    for (i in 1:length(P))
      draw_alleles(m, y[i], P[[i]], col.cte = col.P[i], y.pos = 0.1)
    grid::upViewport()
    grid::pushViewport(vp2)
    draw_homologous(m, h.names = letters[(m+1):(2*m)], parent = "P2")
    Q <- x$maps[[config]]$seq.ph$Q
    if(length(col.Q)==1)
      col.Q<-rep(col.Q, length(Q))
    y <- seq(from = 0.05, to = 0.95, length.out = length(Q))
    for (i in 1:length(Q))
      draw_alleles(m, y[i], Q[[i]], col.cte = col.Q[i])
    grid::upViewport()
    grid::pushViewport(vp3)
    d <- cumsum(c(0, imf_h(x$maps[[config]]$seq.rf)))
    x1 <- d/max(d)
    grid::grid.roundrect(x = -0.01, y = 0.5, width = 1.02, height = 0.07, r = grid::unit(2, "mm"), gp = grid::gpar(fill = "white", col = "black", lwd = 0.7), just = c("left",                                                                                                                                           "center"))
    x2 <- seq(from = 0.01125, to = 0.99, length.out = length(x1))
    id<-x$maps[[config]]$seq.num
    id<-paste0(get(x$info$data.name, pos =1)$mrk.names[id], " - (", id, ")")
    for (i in 1:length(x1)) {
      grid::grid.lines(x = x1[i], y = c(0.45, 0.55), gp = grid::gpar(lwd = 1))
      grid::grid.lines(x = c(x1[i], x2[i]), y = c(0.55, 1.19), gp = grid::gpar(lwd = 0.5, lty = 1, col = "gray"))
      if(mrk.names){
        grid::grid.text(label = id[i], x = x2[i], y = 1.22,
                        rot = 90, just = "right", gp = grid::gpar(cex = cex))      
      } 
    }
    dmax <- max(d)
    d1<-pretty(seq(0, dmax, length.out = 10), n = 10)
    if(rev(d1)[1] > dmax) d1<-d1[-length(d1)]
    x3 <- seq(from = 0, to = max(d1)/dmax, length.out = length(d1))
    grid::grid.lines(x = c(0, max(x3)), y = c(0.2, 0.2), gp = grid::gpar(lwd = 0.5))
    for (i in 1:length(x3)) {
      grid::grid.lines(x = x3[i], y = c(0.2, 0.17), gp = grid::gpar(lwd = 2))
      grid::grid.text(label = round(d1[i], 1), x = x3[i], y = 0.155, rot = 90, just = "right", gp = grid::gpar(cex = 0.8))
    }
    grid::grid.text(label = "centimorgans (cM)", x = 0.5, y = -0.5)
    upViewport()
  } else {
    plot_map_list(x) 
  }
}

#' draw homologous
#' @param void interfunction to be documented
#' @keywords internal
draw_homologous <- function(m, y.pos = 0.5, h.names = letters[1:m], parent = "P") {
  st <- y.pos - m/40
  s <- seq(from = st, by = 0.1, length.out = m)
  ct<-1
  for (i in s){
    grid::grid.roundrect(x = 0.04, y = i, width = 0.92, height = 0.09, r = grid::unit(1, "mm"), gp = grid::gpar(fill = "gray", col = "darkgray", lwd = 0.1),
                         just = c("left", "center"))
    grid::grid.text(label = h.names[ct], x = 0.03, y = i,
                    just = "right", gp = grid::gpar(cex =.7))  
    ct<-ct+1
  } 
  grid::grid.text(label = parent, x = 0.055, y = i+0.2,
                  just = "right", gp = grid::gpar(cex = 1))  
}

#' draw alleles
#' @param void interfunction to be documented
#' @keywords internal
draw_alleles <- function(m, x, v, col.cte = 2, y.pos = 0.5) {
  st <- y.pos - m/40
  s <- seq(from = st, by = 0.1, length.out = m)
  u <- is.na(match(1:m, v)) + 1
  gp <- list(grid::gpar(fill = col.cte, col = col.cte, lineend = "square", linejoin = 	"mitre"),
             grid::gpar(fill = "white", col = "white", lineend = "square", linejoin = 	"mitre"),
             just = c("center", "center"))
  for (i in 1:length(s)) grid::grid.rect(x = x, y = s[i], width = 0.006, height = 0.075, gp = gp[[u[i]]])
}

#' Get the tail of a marker sequence up to the point where the markers
#' provide no additional infomation.
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export get_full_info_tail
get_full_info_tail <- function(input.obj, extend = NULL) {
  ## checking for correct object
  if (all(is.na(match(class(input.obj), c("mappoly.map", "mappoly.map.haplo")))))
    stop(deparse(substitute(input.obj)), " is not an object of class 'mappoly.map' or 'mappoly.map.haplo'")
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
#' @param void interfunction to be documented
#' @keywords internal
#' @export filter_map_at_hmm_thres
filter_map_at_hmm_thres <- function(map, thres.hmm){
  map$info$ph.thresh <- NULL
  map$maps<-map$maps[get_LOD(map, sorted = FALSE) < thres.hmm]
  map
}
#' makes a phase list from map, selecting only 
#' configurations under a certain threshold
#' @param void interfunction to be documented
#' @keywords internal
#' @export update_ph_list_at_hmm_thres
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
#' @param void interfunction to be documented
#' @keywords internal
#' @export get_ph_list_subset
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
#' @param void interfunction to be documented
#' @keywords internal
#' @export concatenate_ph_list
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
#' @param void interfunction to be documented
#' @keywords internal
#' @export add_mrk_at_tail_ph_list
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

