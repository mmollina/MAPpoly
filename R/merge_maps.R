#' Merge two maps
#' 
#' Estimates the linkage phase and recombination fraction between two pre-built maps 
#' and creates a new map by merging them.
#' 
#' \code{merge_maps} uses two-point information, under a given LOD threshold, to reduce the 
#' linkage phase search space. The remaining linkage phases are tested using the genotype 
#' probabilities
#' 
#' @param map.list a list of objects of class \code{mappoly.map} to be merged.
#' 
#' @param twopt an object of class \code{poly.est.two.pts.pairwise}
#'     containing the two-point information
#'     
#' @param  thres.twopt the threshold used to determine if the linkage
#'     phases compared via two-point analysis should be considered 
#'     for the search space reduction (default = 3)
#'      
#' @param thres.hmm the threshold used to determine which linkage 
#'     phase configurations should be retuned when merging two maps.
#'     If "best" (default), returns only the best linkage phase 
#'     configuration. NOTE: if merging multiple maps, it always uses 
#'     the "best" linkage phase configuration at each block insertion.
#'     
#' @param genoprob.list a list of objects of class \code{mappoly.genoprob} 
#'     containing the genotype probabilities for the maps to be merged. 
#'     If \code{NULL} (default), the probabilities are computed.
#'                     
#' @param tol the desired accuracy (default = 10e-04)
#'                     
#' @return An object of class 'mappoly.map' with the following structure:
#' \item{m}{the ploidy level}
#' \item{mrk.names}{the names of markers present in the sequence}
#' \item{data.name}{name of the dataset of class \code{mappoly.data}}
#' \item{ph.thres}{the LOD threshold used to define the linkage phase configurations to test}
#' \item{maps}{a list containing the sequence of markers, their recombination fractions,
#' the linkage phase configuration for all markers in both parents P and Q and the 
#' map's joint likelihood}
#' 
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' 
#' @examples
#' \dontrun{
#' #### Tetraploid example #####
#' map1<-get_submap(solcap.dose.map[[1]], 1:5)
#' map2<-get_submap(solcap.dose.map[[1]], 6:15)
#' map3<-get_submap(solcap.dose.map[[1]], 16:30)
#' full.map<-get_submap(solcap.dose.map[[1]], 1:30)
#' s<-make_seq_mappoly(tetra.solcap, full.map$maps[[1]]$seq.num)
#' counts<-cache_counts_twopt(input.seq = s, cached = TRUE)
#' twopt <- est_pairwise_rf(input.seq = s, count.cache = counts)
#' merged.maps<-merge_maps(map.list = list(map1, map2, map3), 
#'                         twopt = twopt,
#'                         thres.twopt = 3)
#' plot(merged.maps$map, mrk.names = TRUE)                       
#' plot(full.map, mrk.names = TRUE)                       
#' best.phase <- merged.maps$map$maps[[1]]$seq.ph
#' names.id<-names(best.phase$P)
#' compare_haplotypes(m = 4, best.phase$P[names.id], 
#'                    full.map$maps[[1]]$seq.ph$P[names.id]) 
#' compare_haplotypes(m = 4, best.phase$Q[names.id], 
#'                    full.map$maps[[1]]$seq.ph$Q[names.id])
#'                     
#' #### Hexaploid example #####
#' map1<-get_submap(maps.hexafake[[1]], 1:5)
#' map2<-get_submap(maps.hexafake[[1]], 6:15)
#' map3<-get_submap(maps.hexafake[[1]], 16:30)
#' full.map<-get_submap(maps.hexafake[[1]], 1:30)
#' s<-make_seq_mappoly(hexafake, full.map$maps[[1]]$seq.num)
#' counts<-cache_counts_twopt(input.seq = s, cached = TRUE)
#' twopt <- est_pairwise_rf(input.seq = s, count.cache = counts)
#' merged.maps<-merge_maps(map.list = list(map1, map2, map3), 
#'                         twopt = twopt,
#'                         thres.twopt = 3)
#' plot(merged.maps$map, mrk.names = TRUE)                       
#' plot(full.map, mrk.names = TRUE)                       
#' best.phase <- merged.maps$map$maps[[1]]$seq.ph
#' names.id<-names(best.phase$P)
#' compare_haplotypes(m = 6, best.phase$P[names.id], 
#'                    full.map$maps[[1]]$seq.ph$P[names.id]) 
#' compare_haplotypes(m = 6, best.phase$Q[names.id], 
#'                    full.map$maps[[1]]$seq.ph$Q[names.id])
#'}
#' 
#' @export
merge_maps<-function(map.list, 
                     twopt,
                     thres.twopt = 10,
                     genoprob.list = NULL,
                     thres.hmm = "best",
                     tol = 10e-5){
  ## Checking class of arguments
  if(any(sapply(map.list, class) != "mappoly.map"))
    stop("at least one object in ", deparse(substitute(map.list)), " is not an object of class 'mappoly.map'")
  if (length(unique(sapply(map.list, function(x) x$info$data.name))) != 1)
    stop("MAPpoly won't merge maps from different datasets.")
  if (!any(class(twopt) == "poly.est.two.pts.pairwise")){
    stop(deparse(substitute(twopt)), " is not an object of class 'poly.est.two.pts.pairwise'")    
  }
  ## Check twopt consistency
  s.temp<-make_seq_mappoly(get(map.list[[1]]$info$data.name), 
                           unlist(sapply(map.list, function(x) x$info$mrk.names)), 
                           map.list[[1]]$info$data.name)
  check <- check_pairwise(s.temp, twopt)
  if (any(check != 0)){
    cat("There is no information for pairs: \n")
    print(check)
    stop()
  }
  ## For merging multiple maps 
  if(length(map.list) > 2)
    thres.hmm = "best"
  ## choosing best linkage phase configuration
  i.lpc <- sapply(map.list, function(x) which.min(get_LOD(x, sorted = FALSE)))
  ## Calculating genoprob if needed
  if (is.null(genoprob.list) | 
      length(genoprob.list) != length(map.list) |
      !identical(sapply(genoprob.list, function(x) names(x$map)), sapply(map.list, function(x) x$info$mrk.names)))
  {
    ##message("Calculating genoprob.list")
    genoprob.list <- vector("list", length(map.list))
    for(i in 1:length(genoprob.list)){
      r <- map.list[[i]]$maps[[1]]$seq.rf
      r[r<1e-5] <- 1e-5
      map.list[[i]]$maps[[1]]$seq.rf <- r
      suppressMessages(genoprob.list[[i]] <- calc_genoprob(map.list[[i]], phase.config = i.lpc[[i]], verbose = FALSE))
    }
  }
  ## checking ploidy level consistency
  m <- unique(sapply(map.list, function(x) x$info$m))
  if(length(m) != 1)
    stop("MAPpoly won't merge maps with different ploidy levels.")
  ## number of genotipic states
  ngam <- choose(m, m/2)
  ## Number of genotypes in the offspring
  ngen <- ngam * ngam
  ## number of markers
  nmrk <- sapply(map.list, function(x) length(x$info$mrk.names))
  ## number of individuals
  nind <- dim(genoprob.list[[1]]$probs)[3]
  ## the threshold for visiting states: 1/ngen
  thresh.cut.path <- 1/ngen
  ## Hash table: homolog combination --> states to visit in both parents
  A<-as.matrix(expand.grid(0:(ngam-1), 
                           0:(ngam-1))[,2:1])
  rownames(A) <- dimnames(genoprob.list[[1]]$probs)[[1]]
  if(length(map.list) == 2){
    ## h: states to visit in both parents
    ## e: probability distribution 
    e.first <- h.first <- vector("list", nind)
    for(i in 1:nind){
      a <- genoprob.list[[1]]$probs[,dim(genoprob.list[[1]]$probs)[2],i]  
      e.first[[i]] <- a[a > thresh.cut.path]
      h.first[[i]] <- A[names(e.first[[i]]), , drop = FALSE]
    }
    ## Second map
    tpt.temp <- make_pairs_mappoly(input.twopt = twopt, s.temp)
    rf.matrix <- rf_list_to_matrix(input.twopt = tpt.temp,
                                   thresh.LOD.ph =  thres.twopt, 
                                   thresh.LOD.rf = thres.twopt, 
                                   shared.alleles = TRUE, 
                                   verbose = FALSE)
    w<-generate_all_link_phases_elim_equivalent_haplo(block1 = map.list[[1]]$maps[[i.lpc[1]]], 
                                                      block2 = map.list[[2]]$maps[[i.lpc[2]]],
                                                      rf.matrix = rf.matrix,
                                                      m = m, 
                                                      max.inc = 0)
    ## get test.maps and conditional probabilities
    test.maps<-p<-vector("list", length(w))
    rem <- logical(length(w))
    for(i in 1:length(w))
    {
      test.maps[[i]] <- map.list[[2]]
      test.maps[[i]]$maps[[1]]$seq.ph <- w[[i]]
      r <- test.maps[[i]]$maps[[1]]$seq.rf
      r[r<1e-5] <- 1e-5
      test.maps[[i]]$maps[[1]]$seq.rf <- r
      suppressMessages(p[[i]] <- calc_genoprob(test.maps[[i]], verbose = FALSE))
    }
    ## h: states to visit in both parents
    ## e: probability distribution 
    h.second<-e.second<-vector("list", length(w))
    for(j in 1:length(w)){
      etemp<-htemp<-vector("list", nind)
      for(i in 1:nind){
        a <- p[[j]]$probs[,1,i]  
        etemp[[i]] <- a[a > thresh.cut.path]
        htemp[[i]] <- A[names(etemp[[i]]), , drop = FALSE]
      }
      h.second[[j]] <- htemp
      e.second[[j]] <- etemp
    }
    configs<-vector("list", length(test.maps))
    names(configs) <- paste0("Phase_config.", 1:length(test.maps))
    res<-matrix(NA, nrow = length(w), ncol = 2, dimnames = list(names(configs), c("log_like", "rf")))
    ## HMM
    for(i in 1:length(w)){
      #cat("testing", i, "of", length(h), "\n")
      #cat(".")
      h.test<-c(list(h.first), h.second[i])
      e.test<-c(list(e.first), e.second[i])
      restemp<-est_haplo_hmm(m = m, 
                             n.mrk = length(h.test), 
                             n.ind = nind, 
                             haplo = h.test, 
                             emit = e.test, 
                             rf_vec = rep(0.01, length(h.test)-1), 
                             verbose = FALSE, 
                             use_H0 = FALSE, 
                             tol = tol) 
      res[i,]<-unlist(restemp)
    }
    #cat("\n")
    for(i in 1:length(test.maps)){
      P<-c(map.list[[1]]$maps[[i.lpc[1]]]$seq.ph$P, test.maps[[i]]$maps[[1]]$seq.ph$P)
      Q<-c(map.list[[1]]$maps[[i.lpc[1]]]$seq.ph$Q, test.maps[[i]]$maps[[1]]$seq.ph$Q)
      names(P)<-names(Q)<-c(map.list[[1]]$maps[[i.lpc[1]]]$seq.num, test.maps[[i]]$maps[[1]]$seq.num)
      configs[[i]]<-list(P = P, Q = Q)
    }
    res<-res[order(res[,"log_like"], decreasing = TRUE),,drop = FALSE]
    ## Updating map
    output.map<-map.list[[1]]
    seq.num<-as.numeric(names(configs[[1]]$P))
    output.map$info$mrk.names <- c(map.list[[1]]$info$mrk.names, map.list[[2]]$info$mrk.names) 
    output.map$info$n.mrk <- length(output.map$info$mrk.names)
    for(i in 1:nrow(res))
    {
      seq.rf <- c(map.list[[1]]$maps[[i.lpc[1]]]$seq.rf, res[i, "rf"], map.list[[2]]$maps[[i.lpc[2]]]$seq.rf)
      output.map$maps[[i]] <- list(seq.num = seq.num, 
                                   seq.rf = seq.rf, 
                                   seq.ph = configs[[rownames(res)[i]]],
                                   loglike = res[i, "log_like"])
    }
    if(thres.hmm == "best")
      thres.hmm <- 10e-10
    output.map <- filter_map_at_hmm_thres(output.map, thres.hmm)
    return(output.map)
  } else { 
    ##sequential phasing
    ## FIXME: instead of computing genoprob for each round of block insetion$
    ## Inherit from previous rounds
    out.map <- map.list[[1]]
    for(i in 2:length(map.list)){
      out.map<-merge_maps(map.list = list(out.map, map.list[[i]]), 
                          twopt = twopt, 
                          tol = tol, 
                          thres.twopt = thres.twopt)  
    }
    out.map$info$seq.num <- unlist(sapply(map.list, function(x) x$info$seq.num))
    out.map$info$seq.dose.p <- unlist(sapply(map.list, function(x) x$info$seq.dose.p))
    out.map$info$seq.dose.q <- unlist(sapply(map.list, function(x) x$info$seq.dose.q))
    out.map$info$sequence <- unlist(sapply(map.list, function(x) x$info$sequence))
    out.map$info$sequence.pos <- unlist(sapply(map.list, function(x) x$info$sequence.pos))
    out.map$info$chisq.pval <- unlist(sapply(map.list, function(x) x$info$chisq.pval))
    ##spliting to reestimate
    map.list2 <- vector("list", length(map.list))
    for(i in 1:length(map.list)){
      suppressMessages(map.list2[[i]] <- get_submap(out.map, mrk.pos = match(map.list[[i]]$info$mrk.names, out.map$info$mrk.names), reestimate.rf = FALSE))
    }
    genoprob.list <- vector("list", length(map.list2))
    for(i in 1:length(genoprob.list)){
      r <- map.list2[[i]]$maps[[1]]$seq.rf
      r[r<1e-5] <- 1e-5
      map.list2[[i]]$maps[[1]]$seq.rf <- r
      suppressMessages(genoprob.list[[i]] <- calc_genoprob(map.list2[[i]], verbose = FALSE))
    }
    temp.map<-est_map_haplo_given_genoprob(map.list2, genoprob.list, tol = tol)
    out.map$maps[[1]]$loglike <- temp.map$map[[1]] 
    w <- NULL
    for(i in 1:(length(map.list)-1)){
      w <- c(w, map.list[[i]]$maps[[1]]$seq.rf, temp.map$map[[2]][i])
    }
    out.map$maps[[1]]$seq.rf <- c(w, map.list[[length(map.list)]]$maps[[1]]$seq.rf)
    temp.map$genoprob$map <- cumsum(c(0, imf_h(temp.map$genoprob$map)))
    dimnames(temp.map$genoprob$probs)[[2]] <- names(temp.map$genoprob$map) <- paste0("M_", 1:length(temp.map$genoprob$map))
    return(out.map)
  }
}
