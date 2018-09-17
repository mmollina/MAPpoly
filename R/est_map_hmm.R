#' Multipoint analysis using Hidden Markov Models in autopolyploids
#'
#' Performs the multipoint analysis proposed by \cite{Mollinari and
  #'  Garcia  (2018)} in a sequence of markers
#'
#'  This function first enumerates a set of linkage phase configurations
#'  based on two-point recombination fraction information using a threshold
#'  provided by the user (argument \code{'thresh'}). After that, for each
#'  one of the configurations, it reconstructs the genetic map for using the
#'  HMM approach. As result, it returns the multipoint likelihood for each
#'  configuration in form of LOD Score comparing each configuration to the most
#'  likely one. It is recommended to use a small number of markers (e.g. 50 markers
#'  for hexaploids) since the possible linkage phase combinations bounded only
#'  by the two-point information can be huge. Also, it can be quite sensible
#'  to small changes in \code{'thresh'}. Thus it is strongly recommended use a small
#'  \code{'thresh'} value in the first run (0.1 - 0.5) and raises it to find a
#'  reasonable number of linkage phase configurations to be evaluated by the
#'  HMM approach. For higher number of markers, please see  \code{'est_rf_hmm_sequential'}
#'  and \code{'est_map_haplotype'}
#
#' @param input.seq an object of class \code{mappoly.sequence}.
#'
#' @param input.ph an object of class \code{two.pts.linkage.phases}.
#'
#' @param thres the threshold used to determine if the linkage phases
#'     compared via two-point analysis should be considered
#'
#' @param twopt an object of class \code{poly.est.two.pts.pairwise}
#'     containing the two-point information
#'
#' @param verbose if \code{TRUE}, current progress is shown; if
#'     \code{FALSE}, no output is produced.
#'
#' @param rf.lim limit of recombination fraction during convergence.
#'     This is specially useful for marker blocks.
#'
#' @param tol the desired accuracy.
#'
#' @param est.given.0.rf logical. If TRUE returns a map forcing all
#'     recombination fractions equals to 0 (1e-5)
#'
#' @param reestimate.single.ph.configuration logical. If \code{TRUE}
#' returns a map without reestimating the map parameters for cases
#' where there are only one possible linkage phase configuration. This argument
#' is intended to be used in a sequential map contruction.
#'
#' @param x an object of one of the classes \code{mappoly.map}
#'
#' @param detailed if TRUE print the linkage phase configuration and the marker 
#' position for all maps. if FALSE print a map summary 
#'
#' @param col.cte a single value or a vector of with size equal to the number of 
#' markers in the map indicating the color of the allelic variants. The default is \code{red}
#' 
#' @param config should be \code{'best'} or the position of the
#'     configuration to be plotted. If \code{'best'}, plot the configuration
#'     with the highest likelihood.
#'
#' @param ... currently ignored
#'
#' @return An object of class 'mappoly.map'
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
#'                              rf.lim = 0.5,
#'                              tol = 0.1,
#'                              est.given.0.rf = FALSE)
#'
#'     ## Re-estimating the map with the most likely configuration
#'     subset.map1 <- est_rf_hmm_single(input.seq = unique.mrks,
#'                                     input.ph.single = subset.map$maps[[1]]$seq.ph,
#'                                     tol = 10e-4,
#'                                     verbose = TRUE)
#'
#'     subset.map$maps[[1]]$seq.ph <- subset.map1$seq.ph
#'
#'      ## Retrieving simulated linkage phase
#'      sim.ph.P.file <- system.file('doc', 'phase_sim_hexa_P.csv', package = 'mappoly')
#'      ph.P <- read.csv2(sim.ph.P.file)
#'      ph.P <- ph_matrix_to_list(ph.P[1:50,-1])
#'      sim.ph.Q.file <- system.file('doc', 'phase_sim_hexa_Q.csv', package = 'mappoly')
#'      ph.Q <- read.csv2(sim.ph.Q.file)
#'      ph.Q <- ph_matrix_to_list(ph.Q[1:50,-1])
#'
#'      ## Estimated linkage phase
#'      ph.P.est <- subset.map$maps[[1]]$seq.ph$P
#'      ph.Q.est <- subset.map$maps[[1]]$seq.ph$Q
#'
#'      compare_haplotypes(m = 6, h1 = ph.P[names(ph.P.est)], h2 = ph.P.est)
#'      compare_haplotypes(m = 6, h1 = ph.Q[names(ph.Q.est)], h2 = ph.Q.est)
#'
#'      try(dev.off(), silent = TRUE)
#'      plot(subset.map, 1)
#'
#'    }
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2018) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _submited_. \url{https://doi.org/10.1101/415232}
#'
#' @export est_rf_hmm
#'
est_rf_hmm <- function(input.seq, input.ph = NULL,
                       thres = 0.5, twopt = NULL,
                       verbose = FALSE, rf.lim = 0.5,
                       tol = 1e-04,
                       est.given.0.rf=FALSE,
                       reestimate.single.ph.configuration = TRUE) {
  ## checking for correct object
  input_classes <- c("mappoly.sequence", "two.pts.linkage.phases")
  if (!inherits(input.seq, input_classes[1])) {
    stop(deparse(substitute(input.seq)), " is not an object of class 'mappoly.sequence'")
  }
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
    dp <- get(twopt$data.name, pos =1)$dosage.p[input.seq$seq.num]
    dq <- get(twopt$data.name, pos =1)$dosage.q[input.seq$seq.num]
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
    return(structure(list(info = list(m = input.seq$m, n.mrk = length(input.seq$seq.num),
                                      data.name = input.seq$data.name, ph.thresh = abs(res.temp[2,1])),
                          maps = maps),
                     class = "mappoly.map"))
  }
  n.ph <- length(input.ph$config.to.test)
  ret.map.no.rf.estimation <- FALSE
  if(n.ph == 1 && !reestimate.single.ph.configuration)
    ret.map.no.rf.estimation <- TRUE
  if (verbose) {
    cat("\nNumber of linkage phase configurations: ", n.ph)
    cat("\n---------------------------------------------\n")
  }
  maps <- vector("list", n.ph)
  if (verbose)
    cat(n.ph, "phase(s): ")
  for (i in 1:n.ph) {
    if (verbose) {
      cat(i, " ")
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
                                   ret.map.no.rf.estimation = ret.map.no.rf.estimation)
  }
  id<-order(sapply(maps, function(x) x$loglike), decreasing = TRUE)
  maps<-maps[id]
  if (verbose)
    cat("\n")
  structure(list(info = list(m = input.seq$m, n.mrk = length(input.seq$seq.num),
                             data.name = input.seq$data.name, ph.thresh = input.ph$thres),
                 maps = maps),
            class = "mappoly.map")
}


#' Multipoint analysis using Hidden Markov Models - Sequential phase elimination
#'
#' Performs the multipoint analysis proposed by \cite{Mollinari and
  #'  Garcia et al.  (2017)} in a sequence of markers removing unlikely phases
#' using sequential multipoint information
#'
#' This function sequentially includes markers into a map given an
#' ordered sequence. It uses two-point information to eliminate
#' unlikely linkage phase configurations given \code{thres.twopt}. The
#' search is made within a window of size \code{extend.tail}. For the
#' remaining configurations the HMM-based likelihood is computed and
#' the ones that pass \code{thres.hmm} are eliminated.
#'
#' @param input.seq an object of class \code{mappoly.sequence}.
#'
#' @param thres.twopt the threshold used to determine if the linkage
#'     phases compared via two-point analysis should be considered 
#'     for the search space reduction. (A.K.A. \deqn{\eta} in 
#'     Mollinari and Garcia 2018). 
#'     
#' @param thres.hmm the threshold used to determine if the linkage
#'     phases compared via hmm analysis should be evaluated in the 
#'     next round of marker inclusion.   
#'
#' @param info.tail if \code{TRUE} uses the complete informative tail
#'     of the chain (i.e. ploidy x 2 homologous can be distinguished) 
#'     to calculate the likelihood of the linkage phases
#'
#' @param extend.tail the length of the tail of the chain that should
#'     be used to calculate the likelihood of the linkage phases. If
#'     \code{info.tail = TRUE}, the function uses at least \code{extend.tail}
#'     as the length of the tail.
#'
#' @param twopt an object of class poly.est.two.pts.pairwise
#'     containing the two-point information
#'
#' @param verbose If \code{TRUE}, current progress is shown; if
#'     \code{FALSE}, no output is produced.
#'
#' @param rf.lim limit of recombination fraction convergence.
#'
#' @param tol the desired accuracy during the sequential phase.
#'
#' @param tol.final the desired accuracy for the final map.
#'
#' @param est.given.0.rf logical. If TRUE returns a map forcing all
#'     recombination fractions equal to 0 (1e-5)
#'
#'@param reestimate.single.ph.configuration logical. If \code{TRUE}
#' returns a map without reestimating the map parameters for cases
#' where there are only one possible linkage phase configuration.
#'
#' @param sub.map.size.diff.limit the maximum accepted length
#'     xdifference between the current submap and the previous one. If the
#'     size exceeds the limit, the marker will not be inserted. If
#'     \code{NULL}, the it will insert all markers.
#'
#' @param phase.number.limit the maximum number of linkage phases of the sub-maps. If the
#'     size exceeds the limit, the marker will not be inserted. If
#'     \code{NULL}, the it will insert all markers.
#'
#' @return An object of class 'mappoly.map'
#'
#' @examples
#'  \dontrun{
#'     data(hexafake)
#'     mrk.subset<-make_seq_mappoly(hexafake, 'seq1')
#'     red.mrk<-elim_redundant(mrk.subset)
#'     unique.mrks<-make_seq_mappoly(red.mrk)
#'     counts.web<-cache_counts_twopt(unique.mrks, get.from.web = TRUE)
#'     subset.pairs<-est_pairwise_rf(input.seq = unique.mrks,
#'                                   count.cache = counts.web,
#'                                   n.clusters = 16,
#'                                   verbose=TRUE)
#'     system.time(
#'     subset.map <- est_rf_hmm_sequential(input.seq = unique.mrks,
#'                                         thres.twopt = 5,
#'                                         thres.hmm = 3,
#'                                         extend.tail = 50,
#'                                         tol = 0.1,
#'                                         tol.final = 10e-3,
#'                                         twopt = subset.pairs,
#'                                         verbose = TRUE))
#'
#'      ## Retrieving simulated linkage phase
#'      sim.ph.P.file <- system.file('doc', 'phase_sim_hexa_P.csv', package = 'mappoly')
#'      ph.P <- read.csv2(sim.ph.P.file)
#'      ph.P <- ph_matrix_to_list(ph.P[,-1])
#'      sim.ph.Q.file <- system.file('doc', 'phase_sim_hexa_Q.csv', package = 'mappoly')
#'      ph.Q <- read.csv2(sim.ph.Q.file)
#'      ph.Q <- ph_matrix_to_list(ph.Q[,-1])
#'
#'      ## Estimated linkage phase
#'      ph.P.est <- subset.map$maps[[1]]$seq.ph$P
#'      ph.Q.est <- subset.map$maps[[1]]$seq.ph$Q
#'
#'      ##Notice that two estimated homologous in parent P are different
#'      ##from the simulated ones
#'      compare_haplotypes(m = 6, h1 = ph.P[names(ph.P.est)], h2 = ph.P.est)
#'      compare_haplotypes(m = 6, h1 = ph.Q[names(ph.Q.est)], h2 = ph.Q.est)
#'
#'      try(dev.off(), silent = TRUE)
#'      plot(subset.map, 1)
#'
#'    }
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2018) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _submited_. \url{https://doi.org/10.1101/415232}
#'
#' @export est_rf_hmm_sequential
est_rf_hmm_sequential <- function(input.seq,
                                  thres.twopt,
                                  thres.hmm,
                                  info.tail = TRUE,
                                  extend.tail = NULL,
                                  twopt,
                                  verbose = FALSE,
                                  rf.lim = 0.5,
                                  tol = 0.1,
                                  tol.final = 10e-3,
                                  est.given.0.rf=FALSE,
                                  reestimate.single.ph.configuration = TRUE,
                                  sub.map.size.diff.limit = Inf,
                                  phase.number.limit = Inf) {
  ## checking for correct object
  if (!any(class(input.seq) == "mappoly.sequence"))
    stop(deparse(substitute(input.seq)),
         " is not an object of class 'mappoly.sequence'")
  if (!any(class(twopt) == "poly.est.two.pts.pairwise"))
    stop(deparse(substitute(twopt)),
         " is not an object of class 'poly.est.two.pts.pairwise'")
  if(rf.lim > mf_h(sub.map.size.diff.limit))
    rf.lim <- mf_h(sub.map.size.diff.limit)
  ## Information about the map
  if(verbose) {
    cat("Number of markers :", length(input.seq$seq.num), "\n")
    cat("----------------------------------------\n")
  }
  if(is.null(extend.tail))
    extend.tail<-length(input.seq$seq.num)
  flag <- FALSE
  if(length(input.seq$seq.num) == 2)
  {
    if(pos_twopt_est(input.seq))
      stop("Impossible to estimate recombination fraction due to lack of information")
    maps<-vector("list", 1)
    res.temp <- twopt$pairwise[[paste(input.seq$seq.num, collapse = "-")]]
    dp <- get(twopt$data.name, pos =1)$dosage.p[input.seq$seq.num]
    dq <- get(twopt$data.name, pos =1)$dosage.q[input.seq$seq.num]
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
    return(structure(list(info = list(m = input.seq$m, n.mrk = length(input.seq$seq.num),
                                      data.name = input.seq$data.name, ph.thresh = abs(res.temp[2,1])),
                          maps = maps),
                     class = "mappoly.map"))
  }
  z1<-pos_twopt_est(input.seq)
  if(min(which(!z1 == TRUE)) > 1)
  {
    s<-min(which(!z1 == TRUE)):length(input.seq$seq.num)
    warning("It was not possible to use two-point info for the first ", min(which(!z1 == TRUE))-1, " marker(s)")
    input.seq<-make_seq_mappoly(input.obj = get(input.seq$data.name, pos = 1), input.seq$seq.num[s], data.name = input.seq$data.name)
  }
  max.mrkname.size<-max(nchar(input.seq$seq.mrk.names))
  max.id.size<-max(nchar(input.seq$seq.num))
  ## Submap for the first two markers
  j <- 2
  while(j <= length(input.seq$seq.num)) {
    if(z1[j]){
      j<-j+1
      next()
    }
    seq.num <- input.seq$seq.num[c(j-1, j)]
    all.ph <- input.ph <- NULL
    ## Matrix for submap distances
    sub.maps <- matrix(NA, length(input.seq$seq.num), length(input.seq$seq.num))
    dimnames(sub.maps) <- list(input.seq$seq.mrk.names, input.seq$seq.mrk.names)
    sub.maps[1] <- 0
    ## Make sequence and constructing map with two markers
    seq.cur <- make_seq_mappoly(input.obj = get(input.seq$data.name, pos= 1),
                                arg = seq.num,
                                data.name = input.seq$data.name)
    if (verbose){
      cat("(1):\t", format(as.character(input.seq$seq.num[j-1]), width = max.id.size),
          " - ", input.seq$seq.mrk.names[j-1], sep="")
      cat("\n(2):\t", format(as.character(input.seq$seq.num[j]), width = max.id.size),
          " - ", format(input.seq$seq.mrk.names[j], width = max.mrkname.size), " --> ", sep="")
    }
    if(z1[j-1])
    {
      if(verbose)
      {
        cat("\n-------\n")
        cat("marker (", input.seq$seq.num[j-1], ")",
            input.seq$seq.mrk.names[j-1], " was not included.")
        cat("\n-------\n")
      }
      j <- j + 1
      flag<-TRUE
      next()
    }
    input.ph <- ls_linkage_phases(input.seq = seq.cur,
                                  thres = thres.twopt,
                                  twopt = twopt,
                                  prev.info = input.ph)
    cur.map <- est_rf_hmm(input.seq = seq.cur,
                          input.ph = input.ph,
                          twopt = twopt,
                          tol = tol,
                          est.given.0.rf = est.given.0.rf,
                          rf.lim = rf.lim,
                          verbose = FALSE,
                          reestimate.single.ph.configuration = reestimate.single.ph.configuration)
    ## Checking for quality (map size and LOD Score)
    size <- sapply(cur.map$maps, function(x) round(sum(imf_h(x$seq.rf)),1))
    sub.maps[2,seq.cur$seq.mrk.names]  <- c(0, imf_h(cur.map$maps[[1]]$seq.rf))
    LOD <- get_LOD(cur.map, sorted = FALSE)
    if(verbose){
      #cat("\n-------\n")
      inf<-rbind(size, round(LOD,2))
      dimnames(inf)<-list(c("length", "LOD"), paste0("conf.", 1:ncol(inf)))
      cat(inf[1,])
      #print(inf)
    }
    w <- size < sub.map.size.diff.limit & LOD < thres.hmm
    if(any(w))
    {
      ## selecting configuration that yield higher likelihood
      ## given the map length limit
      conf<-min(which(w))
      break()
    } else {
      if(verbose)
      {
        cat("\n-------\n")
        cat("marker (", input.seq$seq.num[j-1], ")",
            input.seq$seq.mrk.names[j-1], " was not included.")
        cat("\n-------\n")
      }
      j <- j + 1
      flag<-TRUE
    }
  }
  if(j == length(input.seq$seq.num) && flag)
    stop("It was not possible to build a map using the informed threshold")
  ## cat(sum(sub.maps[j,], na.rm = TRUE), "\n")
  cur.map$maps <- cur.map$maps[conf]
  #print(cur.map)
  all.ph$config.to.test<-list(cur.map$maps[[1]]$seq.ph)
  if(length(input.seq$seq.num) == j)
  {
    final.seq.num<-as.numeric(names(all.ph$config.to.test[[1]]$P))
    final.ph <- input.ph
    final.ph$config.to.test <- all.ph$config.to.test
    final.ph$rec.frac<-rep(0.01, length(final.ph$config.to.test[[1]]$P)-1)
    final.ph$seq.num<-final.seq.num
    final.seq <- make_seq_mappoly(input.obj = get(input.seq$data.name, pos=1),
                                  arg = final.seq.num,
                                  data.name = input.seq$data.name)
    res <- est_rf_hmm(input.seq = final.seq,
                      input.ph = final.ph,
                      twopt = twopt,
                      rf.lim = rf.lim,
                      tol = tol.final,
                      est.given.0.rf = est.given.0.rf,
                      verbose = verbose)
    return(cur.map)
  }
  all.ph.temp<-all.ph
  i1 <- j
  ## Submap for the remaining markers
  for(i in (j+1):length(input.seq$seq.num))
  {
    #print(i)
    if (verbose)
      #cat("\n")
      #cat("----------------------------\n")
      #cat("----------------------------\n")
      cat("\n(", i, "):\t",
          format(as.character(input.seq$seq.num[i]), width = max.id.size),
          " - ",
          format(input.seq$seq.mrk.names[i], width = max.mrkname.size),
          " --> ",
          sep="")
    seq.num <- NULL
    if(info.tail)
      seq.num <- get_full_info_tail(cur.map)$maps[[1]]$seq.num
    if(length(seq.num) < extend.tail)
      seq.num <- tail(cur.map$maps[[1]]$seq.num, extend.tail)
    seq.cur <- make_seq_mappoly(input.obj = get(input.seq$data.name),
                                arg = seq.num,
                                data.name = input.seq$data.name)
    ## eliminating phase cofigurations based on two point
    all.ph.temp$config.to.test[[1]]$P<-all.ph$config.to.test[[1]]$P[as.character(seq.num)]
    all.ph.temp$config.to.test[[1]]$Q<-all.ph$config.to.test[[1]]$Q[as.character(seq.num)]
    input.ph <- ls_linkage_phases(input.seq = seq.cur,
                                  thres = thres.twopt,
                                  twopt = twopt,
                                  mrk.to.add = input.seq$seq.num[i],
                                  prev.info = all.ph.temp)

    if(length(input.ph$config.to.test) > phase.number.limit) {
      if(verbose)
        cat("not included (too many linkage phases)...")
      next()
    }
    seq.num <- c(seq.num, input.seq$seq.num[i])
    seq.cur <- make_seq_mappoly(input.obj = get(input.seq$data.name, pos=1),
                                arg = seq.num,
                                data.name = input.seq$data.name)
    if(verbose) cat(length(input.ph$config.to.test), "ph(s) / ")
    cur.map.temp <- est_rf_hmm(input.seq = seq.cur,
                               input.ph = input.ph,
                               twopt = twopt,
                               rf.lim = rf.lim,
                               tol = tol,
                               est.given.0.rf = est.given.0.rf,
                               verbose = FALSE,
                               reestimate.single.ph.configuration = reestimate.single.ph.configuration)
    ## Checking for quality (map size and LOD Score)
    size <- sapply(cur.map.temp$maps, function(x) round(sum(imf_h(x$seq.rf)),1))
    last.rf <- imf_h(rev(cur.map.temp$maps[[1]]$seq.rf)[1])
    sub.maps[i,seq.cur$seq.mrk.names]  <- c(0, imf_h(cur.map.temp$maps[[1]]$seq.rf))
    size.dif<-sum(sub.maps[i,], na.rm = TRUE)-sum(sub.maps[i1,], na.rm = TRUE)
    LOD <- get_LOD(cur.map.temp, sorted = FALSE)
    if(verbose){
      #cat("\n----------------------------\n")
      inf<-rbind(size, round(LOD,2))
      dimnames(inf)<-list(c("length", "LOD"), paste0("conf.", 1:ncol(inf)))
      #print(inf)
      cat(formatC(round(size.dif, 1), 2, format = "f"),"(", round(last.rf,1),") tail: ", length(seq.num), sep = "")
    }
    #w <- size < sub.map.size.diff.limit & LOD < thres.hmm
    ###############################################
    ###############################################
    ## TESTING ANOTHER APPROACH
    ## Instead of using the length of the submap, use the
    ## difference between the current and the previous map
    w <- size.dif < sub.map.size.diff.limit & LOD < thres.hmm
    if(any(w))
    {
      ## selecting configuration that yield higher likelihood
      ## given the map length limit
      conf<-min(which(w))
      i1 <- i
    } else {
      if(verbose)
        cat(" --> not included...")
      next()
    }
    cur.map<-cur.map.temp
    cur.map$maps <- cur.map$maps[conf]
    #print(cur.map)
    all.ph$config.to.test[[1]]$P<-c(all.ph$config.to.test[[1]]$P,
                                    rev(cur.map$maps[[1]]$seq.ph$P)[1])
    all.ph$config.to.test[[1]]$Q<-c(all.ph$config.to.test[[1]]$Q,
                                    rev(cur.map$maps[[1]]$seq.ph$Q)[1])

  }
  final.seq.num<-as.numeric(names(all.ph$config.to.test[[1]]$P))
  final.ph <- input.ph
  final.ph$config.to.test <- all.ph$config.to.test
  final.ph$rec.frac<-rep(0.01, length(final.ph$config.to.test[[1]]$P)-1)
  final.ph$seq.num<-final.seq.num
  final.seq <- make_seq_mappoly(input.obj = get(input.seq$data.name, pos=1),
                                arg = final.seq.num,
                                data.name = input.seq$data.name)
  res <- est_rf_hmm(input.seq = final.seq,
                    input.ph = final.ph,
                    twopt = twopt,
                    rf.lim = rf.lim,
                    tol = tol.final,
                    est.given.0.rf = est.given.0.rf,
                    verbose = verbose)
  res$sub.maps<-sub.maps
  return(res)
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
      format(apply(M, 1, function(y) cat(c("\t", y[1], "\t", y[2:(x$info$m + 1)], rep(" ", 4), y[(x$info$m + 2):(x$info$m * 2 + 3)], "\n"), collapse = "")))
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
#' @export
plot.mappoly.map <- function(x, col.cte = 2, config = "best", ...) {
  grid::grid.newpage()
  m <- x$info$m
  vp1 <- viewport(x = 0, y = 0.85, width = 1, height = 0.15, just = c("left", "bottom"))
  vp2 <- viewport(x = 0, y = 0.6, width = 1, height = 0.15, just = c("left", "bottom"))
  vp3 <- viewport(x = 0.04, y = 0.4, width = 0.92, height = 0.2, just = c("left", "bottom"))
  pushViewport(vp1)
  draw_homologous(m, y.pos = 0.1)
  if (config == "best")
    config <- which.min(abs(get_LOD(x, sorted=FALSE)))
  if (!is.numeric(config))
    stop("config must be numeic")
  if (length(x$maps) < config) {
    warning("config should be equal or less than", length(x$maps), "\nplotting the best configuration")
    config <- which.min(abs(get_LOD(x, sorted = FALSE)))
  }
  P <- x$maps[[config]]$seq.ph$P
  if(length(col.cte)==1)
    col.cte<-rep(col.cte, length(P))
  y <- seq(from = 0.05, to = 0.95, length.out = length(P))
  for (i in 1:length(P))
    draw_alleles(m, y[i], P[[i]], col.cte = col.cte[i], y.pos = 0.1)
  upViewport()
  pushViewport(vp2)
  draw_homologous(m)
  Q <- x$maps[[config]]$seq.ph$Q
  y <- seq(from = 0.05, to = 0.95, length.out = length(Q))
  for (i in 1:length(Q))
    draw_alleles(m, y[i], Q[[i]], col.cte = col.cte[i])
  upViewport()
  pushViewport(vp3)
  d <- cumsum(c(0, imf_h(x$maps[[config]]$seq.rf)))
  x1 <- d/max(d)
  grid.roundrect(x = -0.01, y = 0.5, width = 1.02, height = 0.07, r = unit(2, "mm"), gp = gpar(fill = "white", col = "black", lwd = 0.7), just = c("left",
                                                                                                                                                   "center"))
  x2 <- seq(from = 0.01125, to = 0.99, length.out = length(x1))
  for (i in 1:length(x1)) {
    grid.lines(x = x1[i], y = c(0.45, 0.55), gp = gpar(lwd = 1))
    grid.lines(x = c(x1[i], x2[i]), y = c(0.55, 1.19), gp = gpar(lwd = 0.5, lty = 1, col = "gray"))
    grid.text(label = x$maps[[config]]$seq.num[i], x = x2[i], y = 1.24,
              rot = 90, just = "right", gp = gpar(cex = 0.5))
  }
  dmax <- max(d)  #-max(d)%%10#+10
  d1 <- seq(0, dmax, length.out = 20)
  x3 <- seq(from = 0, to = dmax/max(d), length.out = length(d1))
  grid.lines(x = c(-0.01, 1.02), y = c(0, 0), gp = gpar(lwd = 0.5))
  for (i in 1:length(x3)) {
    grid.lines(x = x3[i], y = c(0, -0.03), gp = gpar(lwd = 2))
    grid.text(label = round(d1[i], 1), x = x3[i], y = -0.045, rot = 90, just = "right", gp = gpar(cex = 0.8))
  }
  grid.text(label = "centimorgans (cM)", x = 0.5, y = -0.5)
  upViewport()
}

#' draw homologous
#' @param void interfunction to be documented
#' @keywords internal
draw_homologous <- function(m, y.pos = 0.5) {
  st <- y.pos - m/40
  s <- seq(from = st, by = 0.1, length.out = m)
  for (i in s) grid.roundrect(x = 0.04, y = i, width = 0.92, height = 0.09, r = unit(1, "mm"), gp = gpar(fill = "gray", col = "darkgray", lwd = 0.1),
                              just = c("left", "center"))
}

#' draw alleles
#' @param void interfunction to be documented
#' @keywords internal
draw_alleles <- function(m, x, v, col.cte = 2, y.pos = 0.5) {
  st <- y.pos - m/40
  s <- seq(from = st, by = 0.1, length.out = m)
  u <- is.na(match(1:m, v)) + 1
  gp <- list(gpar(fill = col.cte, col = col.cte, lineend = "square", linejoin = 	"bevel"),
             gpar(fill = "white", col = "white", lineend = "square", linejoin = 	"bevel"),
             just = c("center", "center"))
  for (i in 1:length(s)) grid.rect(x = x, y = s[i], width = 0.006, height = 0.075, gp = gp[[u[i]]])
}

#' Get the tail of a marker sequence up to the point where the markers
#' provide no additional infomation.
#'
#' @param void interfunction to be documented
#' @keywords internal
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
