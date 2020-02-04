#' Merge two maps
#' 
#' Estimates the linkage phase and recombination fraction between two pre-built maps 
#' and creates a new map by merging them.
#' 
#' \code{merge_maps} uses two-point information, under a given LOD threshold, to reduce the 
#' linkage phase search space. The remaining linkage phases are tested using the genotype 
#' probabilities
#' 
#' @param input.map1 an object of class \code{mappoly.map}, which comprehends the first map to be merged
#' 
#' @param input.map2 an object of class \code{mappoly.map}, which comprehends the second map to be merged 
#' 
#' @param twopt an object of class \code{poly.est.two.pts.pairwise}
#'     containing the two-point information
#'     
#' @param  thres.twopt the threshold used to determine if the linkage
#'     phases compared via two-point analysis should be considered 
#'     for the search space reduction (default = 3)
#'      
#' @param thres.hmm the threshold used to determine which linkage 
#'     phase configurations should be retuned (default = 10)
#'     
#' @param genoprob.map1 an object of class \code{mappoly.genoprob} containing the 
#' genotype probabilities for the first map
#' 
#' @param phase.config.map1 which phase configuration should be used for \code{input.map1}. "best" (default) 
#'                     will choose the configuration associated with the maximum likelihood
#'                     
#' @param phase.config.map2 which phase configuration should be used for \code{input.map2}. "best" (default) 
#'                     will choose the configuration associated with the maximum likelihood
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
#' full.map<-get_submap(solcap.dose.map[[1]], 1:15)
#' s<-make_seq_mappoly(tetra.solcap, full.map$maps[[1]]$seq.num)
#' counts<-cache_counts_twopt(input.seq = s, get.from.web = TRUE)
#' twopt <- est_pairwise_rf(input.seq = s, count.cache = counts)
#' genoprob.map1 <- calc_genoprob(map1)
#' merged.maps<-merge_maps(input.map1 = map1, 
#'                         input.map2 = map2, 
#'                         twopt = twopt,
#'                         thres.twopt = 3,
#'                         genoprob.map1 = genoprob.map1)
#' 
#' plot(merged.maps, mrk.names = TRUE)                       
#' plot(full.map, mrk.names = TRUE)                       
#' best.phase <- merged.maps$maps[[1]]$seq.ph
#' names.id<-names(best.phase$P)
#' plot_compare_haplotypes(m = 6,
#'                         hom.allele.p1 = best.phase$P[names.id],
#'                         hom.allele.q1 = best.phase$Q[names.id],
#'                         hom.allele.p2 = full.map$maps[[1]]$seq.ph$P[names.id],
#'                         hom.allele.q2 = full.map$maps[[1]]$seq.ph$Q[names.id])  
#'                         
#' #### Hexaploid example #####
#' map1<-get_submap(maps.hexafake[[1]], 1:5)
#' map2<-get_submap(maps.hexafake[[1]], 6:15)
#' full.map<-get_submap(maps.hexafake[[1]], 1:15)
#' s<-make_seq_mappoly(hexafake, full.map$maps[[1]]$seq.num)
#' counts<-cache_counts_twopt(input.seq = s, get.from.web = TRUE)
#' twopt <- est_pairwise_rf(input.seq = s, count.cache = counts)
#' genoprob.map1 <- calc_genoprob(map1)
#' merged.maps<-merge_maps(input.map1 = map1, 
#'                         input.map2 = map2, 
#'                         twopt = twopt,
#'                         thres.twopt = 3,
#'                         genoprob.map1 = genoprob.map1)
#' 
#' plot(merged.maps, mrk.names = TRUE)                       
#' plot(full.map, mrk.names = TRUE)                       
#' best.phase <- merged.maps$maps[[1]]$seq.ph
#' names.id<-names(best.phase$P)
#' plot_compare_haplotypes(m = 6,
#'                         hom.allele.p1 = best.phase$P[names.id],
#'                         hom.allele.q1 = best.phase$Q[names.id],
#'                         hom.allele.p2 = full.map$maps[[1]]$seq.ph$P[names.id],
#'                         hom.allele.q2 = full.map$maps[[1]]$seq.ph$Q[names.id])  
#'}
#' 
#' @export
merge_maps<-function(input.map1, 
                     input.map2, 
                     twopt,
                     thres.twopt = 3,
                     thres.hmm = 10,
                     genoprob.map1 = NULL,
                     phase.config.map1 = "best",
                     phase.config.map2 = "best",
                     tol = 10e-4){
  ## Checking class of arguments
  if(!inherits(input.map1, "mappoly.map")) {
    stop(deparse(substitute(input.map1)), " is not an object of class 'mappoly.map'")
  }
  if(!inherits(input.map2, "mappoly.map")) {
    stop(deparse(substitute(input.map2)), " is not an object of class 'mappoly.map'")
  }
  if(input.map1$info$data.name != input.map2$info$data.name)
    stop("MAPpoly won't merge maps from different datasets.")
  if (!any(class(twopt) == "poly.est.two.pts.pairwise")){
    stop(deparse(substitute(twopt)), " is not an object of class 'poly.est.two.pts.pairwise'")    
  }
  ## Check twopt consistency
  s.temp<-make_seq_mappoly(get(input.map1$info$data.name), 
                           c(input.map1$info$mrk.names, input.map2$info$mrk.names), 
                           data.name = input.map1$info$data.name)
  check <- check_pairwise(s.temp, twopt)
  if (any(check != 0)){
    cat("There is no information for pairs: \n")
    print(check)
    stop()
  }
  
  ## choosing the linkage phase configuration
  ## map 1
  LOD.conf1 <- get_LOD(input.map1, sorted = FALSE)
  if(phase.config.map1 == "best") {
    i.lpc1 <- which.min(LOD.conf1)
  }  else if (phase.config.map1 > length(LOD.conf1)) {
    stop("invalid linkage phase configuration")
  } else i.lpc1 <- phase.config.map1
  ## map 2
  LOD.conf2 <- get_LOD(input.map2, sorted = FALSE)
  if(phase.config.map2 == "best") {
    i.lpc2 <- which.min(LOD.conf2)
  }  else if (phase.config.map2 > length(LOD.conf2)) {
    stop("invalid linkage phase configuration")
  } else i.lpc2 <- phase.config.map2
  
  ## FIXME: include a check for the appropriate phase of the provided genoprob.map's
  if(is.null(genoprob.map1)){
    message("Calculating genoprob.map1")
    genoprob.map1 <- calc_genoprob(input.map1, phase.config = i.lpc1)
  }
  if(!inherits(genoprob.map1 , "mappoly.genoprob")) {
    stop("'", deparse(substitute(genoprob.map1)), "' is not an object of class 'mappoly.genoprob'")
  }
  if(!identical(names(genoprob.map1$map), input.map1$info$mrk.names)){
    warning("'", deparse(substitute(genoprob.map1)), "' is inconsistent with 'input.map1'.\n  Recalculating genoprob.map1.")
    genoprob.map1 <- calc_genoprob(input.map1)
  }
  ## ploidy
  m1 <- input.map1$info$m
  m2 <- input.map2$info$m
  if(m1!=m2)
    stop("MAPpoly won't merge maps with different ploidy levels.")
  m <- m1
  ## number of genotipic states
  ngam <- choose(m, m/2)
  ## Number of genotypes in the offspring
  ngen <- ngam^2
  ## number of markers
  nmrk1 <- input.map1$info$n.mrk
  nmrk2 <- input.map2$info$n.mrk
  ## number of individuals
  nind <- dim(genoprob.map1$probs)[3]
  ## the threshold for visiting states: 1/ngen
  thresh.cut.path <- 1/ngen
  
  ## Hash table: homolog combination --> states to visit in both parents
  A<-as.matrix(expand.grid(0:(ngam-1), 
                           0:(ngam-1))[,2:1])
  rownames(A) <- dimnames(genoprob.map1$probs)[[1]]
  ## First map
  ## h: states to visit in both parents
  ## e: probability distribution 
  e.first <- h.first <- vector("list", nind)
  for(i in 1:nind){
    a <- genoprob.map1$probs[,dim(genoprob.map1$probs)[2],i]  
    e.first[[i]] <- a[a > thresh.cut.path]
    h.first[[i]] <- A[names(e.first[[i]]), , drop = FALSE]
  }
  ## Second map
  tpt.temp <- make_pairs_mappoly(input.twopt = twopt, s.temp)
  rf.matrix <- rf_list_to_matrix(input.twopt = tpt.temp,
                                 thresh.LOD.ph =  thres.twopt, 
                                 thresh.LOD.rf = thres.twopt, 
                                 shared.alleles = TRUE)
  w<-generate_all_link_phases_elim_equivalent_haplo(block1 = input.map1$maps[[i.lpc1]], 
                                                    block2 = input.map2$maps[[i.lpc2]],
                                                    rf.matrix = rf.matrix,
                                                    m = m, 
                                                    max.inc = 0)
  ## get test.maps and conditional probabilities
  test.maps<-p<-vector("list", length(w))
  for(i in 1:length(w))
  {
    input.map2$maps[[1]]$seq.ph <- w[[i]]
    test.maps[[i]] <- input.map2
    p[[i]] <- calc_genoprob(test.maps[[i]], verbose = FALSE)
  }
  ## h: states to visit in both parents
  ## e: probability distrobution 
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
    cat(".")
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
  cat("\n")
  for(i in 1:length(test.maps)){
    P<-c(input.map1$maps[[i.lpc1]]$seq.ph$P, test.maps[[i]]$maps[[1]]$seq.ph$P)
    Q<-c(input.map1$maps[[i.lpc1]]$seq.ph$Q, test.maps[[i]]$maps[[1]]$seq.ph$Q)
    names(P)<-names(Q)<-c(input.map1$maps[[i.lpc1]]$seq.num, test.maps[[i]]$maps[[1]]$seq.num)
    configs[[i]]<-list(P = P, Q = Q)
  }
  res<-res[order(res[,"log_like"], decreasing = TRUE),,drop = FALSE]
  ## Updating map
  output.map<-input.map1
  seq.num<-as.numeric(names(configs[[1]]$P))
  output.map$info$mrk.names <- c(input.map1$info$mrk.names, input.map2$info$mrk.names) 
  output.map$info$n.mrk <- length(output.map$info$mrk.names)
  for(i in 1:nrow(res))
  {
    seq.rf <- c(input.map1$maps[[i.lpc1]]$seq.rf, res[i, "rf"], input.map2$maps[[i.lpc2]]$seq.rf)
    output.map$maps[[i]] <- list(seq.num = seq.num, 
                                 seq.rf = seq.rf, 
                                 seq.ph = configs[[rownames(res)[i]]],
                                 loglike = res[i, "log_like"])
  }
  output.map <- filter_map_at_hmm_thres(output.map, thres.hmm)
  return(output.map)
}
