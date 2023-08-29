#' Remove markers from a map
#' 
#' This function creates a new map by removing markers from an existing one.
#'
#' @param input.map an object of class \code{mappoly.map}
#' @param mrk a vector containing markers to be removed from the input map, 
#'            identified by their names or positions
#' @return an object of class \code{mappoly.map}
#'
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' 
#' @examples
#'  sub.map <- get_submap(maps.hexafake[[1]], 1:50, reestimate.rf = FALSE)
#'  plot(sub.map, mrk.names = TRUE)
#'  mrk.to.remove <- c("M_1", "M_23", "M_34")
#'  red.map <- drop_marker(sub.map, mrk.to.remove)
#'  plot(red.map, mrk.names = TRUE)
#' 
#' @export
drop_marker <- function(input.map, mrk, verbose = TRUE)
{
  ## Checking class of arguments
  if (!inherits(input.map, "mappoly.map")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  }
  ## Checking markers to be removed
  if(is.numeric(mrk) & all(mrk >= 0)){
    if(any(mrk > input.map$info$n.mrk))
      stop("At least one marker position is not contained in the map.")
  } else if(is.character(mrk)){
    if(any(!mrk%in%input.map$info$mrk.names)){
      stop("At least one marker position is not contained in the map.")
    } else {
      mrk <- which(input.map$info$mrk.names%in%mrk)
    }
  }
  ## Getting new map
  suppressMessages(output.map <- get_submap(input.map = input.map,
                                            mrk.pos = c(1:input.map$info$n.mrk)[-mrk], 
                                            phase.config = 1, 
                                            reestimate.rf = FALSE))
  if(length(input.map$maps) > 1){
    for(i in 2:length(input.map$maps)){
      suppressMessages(temp.map <- get_submap(input.map = input.map,
                                              mrk.pos = c(1:input.map$info$n.mrk)[-mrk], 
                                              phase.config = 1, 
                                              reestimate.rf = FALSE))
      output.map$maps[[i]] <- temp.map$maps[[1]]
    }
  }
  if (verbose) message("
    INFO:
    -----------------------------------------
    The recombination fractions provided were
    obtained using the marker positions in the 
    input map; For accurate values, plese 
    reestimate the map using functions 'reest_rf', 
    'est_full_hmm_with_global_error' or 
    'est_full_hmm_with_prior_prob'")
  return(output.map)
}

#' Add a single marker to a map
#' 
#' Creates a new map by adding a marker in a given position in a pre-built map. 
#'
#' \code{add_marker} splits the input map into two sub-maps to the left and the 
#' right of the given position. Using the genotype probabilities, it computes 
#' the log-likelihood of all possible linkage phases under a two-point threshold
#' inherited from function \code{\link[mappoly]{rf_list_to_matrix}}. 
#'
#' @param input.map an object of class \code{mappoly.map}
#'  
#' @param mrk the name of the marker to be inserted
#' 
#' @param pos the name of the marker after which the new marker should be added.
#'            One also can inform the numeric position (between markers) were the 
#'            new marker should be added. To insert a marker at the beginning of a 
#'            map, use \code{pos = 0}
#'            
#' @param rf.matrix an object of class \code{mappoly.rf.matrix} containing the recombination
#'                  fractions and the number of homologues sharing alleles between pairwise 
#'                  markers on \code{input.map}. It is important that \code{shared.alleles = TRUE}
#'                  in function \code{\link[mappoly]{rf_list_to_matrix}} when computing \code{rf.matrix}.
#' 
#' @param genoprob an object of class \code{mappoly.genoprob} containing the genotype probabilities
#' for all marker positions on \code{input.map}
#' 
#' @param phase.config which phase configuration should be used. "best" (default) 
#'                     will choose the maximum likelihood configuration
#'                     
#' @param tol the desired accuracy (default = 10e-04)
#' 
#' @param extend.tail the length of the chain's tail that should
#'     be used to calculate the likelihood of the map. If \code{NULL} (default), 
#'     the function uses all markers positioned. 
#' 
#' @param r.test for internal use only
#'
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
#' 
#' @return A list of class \code{mappoly.map} with two elements: 
#' 
#' i) info:  a list containing information about the map, regardless of the linkage phase configuration:
#' \item{ploidy}{the ploidy level}
#' \item{n.mrk}{number of markers}
#' \item{seq.num}{a vector containing the (ordered) indices of markers in the map, 
#'                according to the input file}
#' \item{mrk.names}{the names of markers in the map}
#' \item{seq.dose.p1}{a vector containing the dosage in parent 1 for all markers in the map}
#' \item{seq.dose.p2}{a vector containing the dosage in parent 2 for all markers in the map}
#' \item{chrom}{a vector indicating the sequence (usually chromosome) each marker belongs 
#'                 as informed in the input file. If not available, 
#'                 \code{chrom = NULL}}
#' \item{genome.pos}{physical position (usually in megabase) of the markers into the sequence}
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
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' 
#' @examples
#' \donttest{
#' sub.map <- get_submap(maps.hexafake[[1]], 1:20, reestimate.rf = FALSE)
#' plot(sub.map, mrk.names = TRUE)
#' s <- make_seq_mappoly(hexafake, sub.map$info$mrk.names)
#' tpt <- est_pairwise_rf(s)
#' rf.matrix <- rf_list_to_matrix(input.twopt = tpt,
#'                                thresh.LOD.ph = 3, 
#'                                thresh.LOD.rf = 3,
#'                                shared.alleles = TRUE)
#' ###### Removing marker "M_1" (first) #######
#' mrk.to.remove <- "M_1"
#' input.map <- drop_marker(sub.map, mrk.to.remove)
#' plot(input.map, mrk.names = TRUE)
#' ## Computing conditional probabilities using the resulting map
#' genoprob <- calc_genoprob(input.map)
#' res.add.M_1 <- add_marker(input.map = input.map,
#'                         mrk = "M_1",
#'                         pos = 0,
#'                         rf.matrix = rf.matrix,
#'                         genoprob = genoprob,
#'                         tol = 10e-4)  
#'  plot(res.add.M_1, mrk.names = TRUE)                       
#'  best.phase <- res.add.M_1$maps[[1]]$seq.ph
#'  names.id <- names(best.phase$P)
#'  plot_compare_haplotypes(ploidy = 6,
#'                          hom.allele.p1 = best.phase$P[names.id],
#'                          hom.allele.q1 = best.phase$Q[names.id],
#'                          hom.allele.p2 = sub.map$maps[[1]]$seq.ph$P[names.id],
#'                          hom.allele.q2 = sub.map$maps[[1]]$seq.ph$Q[names.id])
#'                          
#' ###### Removing marker "M_10" (middle or last) #######
#' mrk.to.remove <- "M_10"
#' input.map <- drop_marker(sub.map, mrk.to.remove)
#' plot(input.map, mrk.names = TRUE)
#' # Computing conditional probabilities using the resulting map
#' genoprob <- calc_genoprob(input.map)
#' res.add.M_10 <- add_marker(input.map = input.map,
#'                         mrk = "M_10",
#'                         pos = "M_9",
#'                         rf.matrix = rf.matrix,
#'                         genoprob = genoprob,
#'                         tol = 10e-4)  
#'  plot(res.add.M_10, mrk.names = TRUE)                       
#'  best.phase <- res.add.M_10$maps[[1]]$seq.ph
#'  names.id <- names(best.phase$P)
#'  plot_compare_haplotypes(ploidy = 6,
#'                          hom.allele.p1 = best.phase$P[names.id],
#'                          hom.allele.q1 = best.phase$Q[names.id],
#'                          hom.allele.p2 = sub.map$maps[[1]]$seq.ph$P[names.id],
#'                          hom.allele.q2 = sub.map$maps[[1]]$seq.ph$Q[names.id]) 
#' }
#' @export
add_marker <- function(input.map,  
                       mrk, 
                       pos, 
                       rf.matrix, 
                       genoprob = NULL, 
                       phase.config = "best", 
                       tol = 10e-4,
                       extend.tail = NULL,
                       r.test = NULL, 
                       verbose = TRUE){
  ## Checking class of arguments
  if(!inherits(input.map, "mappoly.map")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  }
  if(!inherits(mrk, "character")) {
    stop("Please, provide the marker name")
  }
  if(!inherits(rf.matrix, "mappoly.rf.matrix")) {
    stop(deparse(substitute(rf.matrix)), " is not an object of class 'mappoly.rf.matrix'")
  }
  if(is.null(rf.matrix$ShP)){
    stop(deparse(substitute(rf.matrix)), " should contain the number of homologues sharing alleles.
    Use 'shared.alleles = TRUE' in function 'rf_list_to_matrix'")
  }
  if(!all(c(input.map$info$mrk.names, mrk)%in%colnames(rf.matrix$rec.mat))){
    stop(deparse(substitute(rf.matrix)), " does not contain all necessary information about 'input.map' and 'mrk'.")
  }
  
  ## choosing the linkage phase configuration
  LOD.conf <- get_LOD(input.map, sorted = FALSE)
  if(phase.config  ==  "best") {
    i.lpc <- which.min(LOD.conf)
  }  
  else if (phase.config > length(LOD.conf)) {
    stop("invalid linkage phase configuration")
  } 
  else i.lpc <- phase.config
  
  ## Checking genoprob
  if(is.null(genoprob)){
    if (verbose) message("Calculating genoprob.")
    genoprob <- calc_genoprob(input.map, phase.config = i.lpc, verbose = FALSE)
  }
  if(!inherits(genoprob, "mappoly.genoprob")) {
    stop("'", deparse(substitute(genoprob)), "' is not an object of class 'mappoly.genoprob'")
  }
  if(!identical(names(genoprob$map), input.map$info$mrk.names)){
    warning("'", deparse(substitute(genoprob)), "' is inconsistent with 'input.map'.\n  Recalculating genoprob.")
    genoprob <- calc_genoprob(input.map, phase.config = i.lpc, verbose = FALSE)
  }
  
  ## ploidy
  ploidy <- input.map$info$ploidy
  ## Number of genotypes in the offspring
  n.gen <- dim(genoprob$probs)[1]
  ## number of markers
  n.mrk <- dim(genoprob$probs)[2]
  ## number of individuals
  n.ind <- dim(genoprob$probs)[3]
  ## number of genotypic states
  ngam <- choose(ploidy, ploidy/2)
  ## the threshold for visiting states: 1/n.gen
  thresh.cut.path <- 1/n.gen
  ## Indexing position were the marker should be inserted
  if(is.numeric(pos)){
    if(!(pos  >= 0 & pos <= (n.mrk+1))){
      stop(deparse(substitute(pos)), " out of bounds.")
    } 
  } 
  else if(is.character(pos)){
    if(!pos%in%input.map$info$mrk.names){
      stop("The reference marker is not contained in the map.")
    } else {
      pos <- which(input.map$info$mrk.names%in%pos)
    }
  }
  
  ## Hash table: homolog combination --> states to visit in both parents
  A <- as.matrix(expand.grid(c(1:ngam) - 1, c(1:ngam) - 1)[,2:1])
  rownames(A) <- dimnames(genoprob$probs)[[1]]
  
  ## Indexing marker to be inserted
  if(!mrk%in%colnames(rf.matrix$rec.mat)){
    stop(deparse(substitute(mrk)), " is not in 'rf.matrix'")
  }
  dat <- get(input.map$info$data.name, pos = 1)
  #mrk.id <- match(mrk, dat$mrk.names)
  mrk.id <- mrk
  ## Adding marker: beginning of sequence
  if(pos  ==  0){
    ## h: states to visit in both parents
    ## e: probability distribution (ignored in this version) 
    e.right <- h.right <- vector("list", n.ind)
    for(i in 1:n.ind){
      a.right <-  genoprob$probs[,pos+1,i]
      e.right[[i]] <- a.right[a.right > thresh.cut.path]
      h.right[[i]] <- A[names(e.right[[i]]), , drop = FALSE]
    }
    if(is.null(r.test)){
      cur.map <- rev_map(input.map)
      cur.map <- get_full_info_tail(cur.map, extend.tail)
      r.test <- generate_all_link_phases_elim_equivalent_haplo(block1 = c(cur.map$maps[[i.lpc]], 
                                                                          mrk.names = list(cur.map$info$mrk.names)), 
                                                               block2 = mrk.id, 
                                                               rf.matrix = rf.matrix, 
                                                               ploidy = ploidy, max.inc = 0)
    }
  } else if(pos > 0 & pos < n.mrk){   ## Adding marker: middle positions
    ## h: states to visit in both parents
    ## e: probability distribution (ignored in this version) 
    e.left <- h.left <- e.right <- h.right <- vector("list", n.ind)
    for(i in 1:n.ind){
      a.left <- genoprob$probs[,pos,i]  
      a.right <-  genoprob$probs[,pos+1,i]
      e.left[[i]] <- a.left[a.left > thresh.cut.path]
      h.left[[i]] <- A[names(e.left[[i]]), , drop = FALSE]
      e.right[[i]] <- a.right[a.right > thresh.cut.path]
      h.right[[i]] <- A[names(e.right[[i]]), , drop = FALSE]
    }
    ## Info to the left
    map.left <- get_submap(input.map = input.map,
                           mrk.pos = 1:pos,
                           phase.config = i.lpc, 
                           reestimate.rf = FALSE, 
                           verbose = FALSE)
    cur.map <- get_full_info_tail(map.left, extend.tail)
    r.left <- generate_all_link_phases_elim_equivalent_haplo(block1 = c(cur.map$maps[[i.lpc]], 
                                                                        mrk.names = list(cur.map$info$mrk.names)),
                                                             block2 = mrk.id, 
                                                             rf.matrix = rf.matrix, 
                                                             ploidy = ploidy, max.inc = 0)
    ## Info to the right
    map.right <- get_submap(input.map = input.map,
                            mrk.pos = (pos+1):input.map$info$n.mrk,
                            phase.config = i.lpc, 
                            reestimate.rf = FALSE, verbose = FALSE)
    cur.map <- rev_map(map.right)
    cur.map <- get_full_info_tail(input.obj = cur.map, extend = extend.tail)
    r.right <- generate_all_link_phases_elim_equivalent_haplo(block1 = c(cur.map$maps[[i.lpc]], 
                                                                         mrk.names = list(cur.map$info$mrk.names)), 
                                                              block2 = mrk.id, 
                                                              rf.matrix = rf.matrix, 
                                                              ploidy = ploidy, max.inc = 0)
    ## Using left and right maps 
    r.test <- c(r.left, r.right)
    r.test <- unique(r.test)
    
  } else if(pos  ==  n.mrk){
    ## h: states to visit in both parents
    ## e: probability distribution (ignored in this version) 
    e.left <- h.left <- vector("list", n.ind)
    for(i in 1:n.ind){
      a.left <- genoprob$probs[,pos,i]  
      e.left[[i]] <- a.left[a.left > thresh.cut.path]
      h.left[[i]] <- A[names(e.left[[i]]), , drop = FALSE]
    }
    cur.map <- get_full_info_tail(input.map, extend.tail)
    r.test <- generate_all_link_phases_elim_equivalent_haplo(block1 = c(cur.map$maps[[i.lpc]], 
                                                                        mrk.names = list(cur.map$info$mrk.names)), 
                                                             block2 = mrk.id, 
                                                             rf.matrix = rf.matrix, 
                                                             ploidy = ploidy, max.inc = 0)
  } 
  else stop("should not get here!")
  ## gathering maps to test and conditional probabilities
  test.maps <- mrk.genoprobs <- vector("list", length(r.test))
  for(i in 1:length(test.maps))
  {
    ## This sub-map is just to create a framework
    hap.temp <- get_submap(input.map, 
                           c(1,1), 
                           reestimate.rf = FALSE, verbose = FALSE)
    hap.temp <- filter_map_at_hmm_thres(hap.temp, thres.hmm = 10e-10)
    hap.temp$maps[[1]]$seq.num <- rep(match(mrk.id, dat$mrk.names), 2)
    hap.temp$maps[[1]]$seq.ph <- list(P = c(r.test[[i]]$P, r.test[[i]]$P),
                                      Q = c(r.test[[i]]$Q, r.test[[i]]$Q))
    hap.temp$maps[[1]]$seq.rf <- 10e-6
    hap.temp$info$mrk.names <- rep(mrk,2)
    test.maps[[i]] <- hap.temp
    ## States to visit for inserted biallelic SNP 
    mrk.genoprobs[[i]] <- calc_genoprob(test.maps[[i]], verbose = FALSE)
  }
  ## Inserted marker  
  ## h: states to visit in both parents
  ## e: probability distribution 
  h <- e <- vector("list", length(r.test))
  for(j in 1:length(r.test)){
    etemp <- htemp <- vector("list", n.ind)
    for(i in 1:n.ind){
      a <- mrk.genoprobs[[j]]$probs[,1,i]  
      etemp[[i]] <- a[a > thresh.cut.path]
      htemp[[i]] <- A[names(etemp[[i]]), , drop = FALSE]
    }
    h[[j]] <- htemp
    e[[j]] <- etemp
  }
  configs <- vector("list", length(test.maps))
  names(configs) <- paste0("Phase_config.", 1:length(test.maps))
  res <- matrix(NA, nrow = length(h), ncol = 3, dimnames = list(names(configs), c("log_like", "rf1", "rf2")))
  ## Computing three-point multiallelic map using HMM
  for(i in 1:length(h)){
    if (verbose) {
      if(i%%50 == 0) cat("\n")
      cat(".")
    }
    if(pos  ==  0){
      h.test <- c(h[i], list(h.right))
      names(h.test) <- c("hap1", "hap2")
      e.test <- c(e[i], list(e.right))
      restemp <- est_haplo_hmm(ploidy = ploidy, 
                               n.mrk = length(h.test), 
                               n.ind = n.ind, 
                               haplo = h.test, 
                               #emit = e.test, 
                               rf_vec = rep(0.01, length(h.test)-1), 
                               verbose = FALSE, 
                               use_H0 = FALSE, 
                               tol = tol) 
      temp <- unlist(restemp)
      res[i,1:length(temp)] <- temp
      P <- c(test.maps[[i]]$maps[[1]]$seq.ph$P[1],
             input.map$maps[[1]]$seq.ph$P)
      Q <- c(test.maps[[i]]$maps[[1]]$seq.ph$Q[1], 
             input.map$maps[[1]]$seq.ph$Q)
      names(P) <- names(Q) <- c(test.maps[[i]]$maps[[1]]$seq.num[1], 
                                input.map$maps[[1]]$seq.num)
      configs[[i]] <- list(P = P, Q = Q)
    } 
    else if(pos > 0 & pos < n.mrk){
      h.test <- c(list(h.left), h[i], list(h.right))
      names(h.test) <- c("hap1", "hap2", "hap3")
      e.test <- c(list(e.left), e[i], list(e.right))
      restemp <-est_haplo_hmm(ploidy = ploidy, 
                              n.mrk = length(h.test), 
                              n.ind = n.ind, 
                              haplo = h.test, 
                              #emit = e.test, 
                              rf_vec = rep(0.01, length(h.test)-1), 
                              verbose = FALSE, 
                              use_H0 = FALSE, 
                              tol = tol) 
      temp <- unlist(restemp)
      res[i,1:length(temp)] <- temp
      P <- c(map.left$maps[[1]]$seq.ph$P, 
             test.maps[[i]]$maps[[1]]$seq.ph$P[1],
             map.right$maps[[1]]$seq.ph$P)
      Q <- c(map.left$maps[[1]]$seq.ph$Q, 
             test.maps[[i]]$maps[[1]]$seq.ph$Q[1], 
             map.right$maps[[1]]$seq.ph$Q)
      names(P) <- names(Q) <- c(map.left$maps[[1]]$seq.num, 
                                test.maps[[i]]$maps[[1]]$seq.num[1], 
                                map.right$maps[[1]]$seq.num)
      configs[[i]] <- list(P = P, Q = Q)
    } 
    else if(pos  ==  n.mrk){
      h.test <- c(list(h.left), h[i])
      names(h.test) <- c("hap1", "hap2")
      e.test <- c(list(e.left), e[i])
      restemp <- est_haplo_hmm(ploidy = ploidy, 
                               n.mrk = length(h.test), 
                               n.ind = n.ind, 
                               haplo = h.test, 
                               #emit = e.test, 
                               rf_vec = rep(0.01, length(h.test)-1), 
                               verbose = FALSE, 
                               use_H0 = FALSE, 
                               tol = tol) 
      temp <- unlist(restemp)
      res[i,1:length(temp)] <- temp
      P <- c(input.map$maps[[1]]$seq.ph$P, 
             test.maps[[i]]$maps[[1]]$seq.ph$P[1])
      Q <- c(input.map$maps[[1]]$seq.ph$Q, 
             test.maps[[i]]$maps[[1]]$seq.ph$Q[1])
      names(P) <- names(Q) <- c(input.map$maps[[1]]$seq.num, 
                                test.maps[[i]]$maps[[1]]$seq.num[1])
      configs[[i]] <- list(P = P, Q = Q)
    }
  }
  if (verbose) cat("\n")
  res <- res[order(res[,"log_like"], decreasing = TRUE),,drop = FALSE]
  
  
  ## Updating map
  output.map <- input.map
  seq.num <- as.numeric(names(configs[[1]]$P))
  output.map$info$seq.num <- seq.num
  output.map$info$mrk.names <- dat$mrk.names[seq.num]
  output.map$info$n.mrk <- length(output.map$info$mrk.names)
  output.map$maps <- vector("list", nrow(res))
  for(i in 1:nrow(res))
  {
    ## Updating recombination fractions (approximated)
    if(pos  ==  0){
      seq.rf <- as.numeric(c(res[i, "rf1"], input.map$maps[[i.lpc]]$seq.rf))
    } else if(pos > 0 & pos < n.mrk){
      seq.rf <- as.numeric(c(head(input.map$maps[[i.lpc]]$seq.rf, n = pos - 1),
                             res[i, c("rf1", "rf2")], 
                             tail(input.map$maps[[i.lpc]]$seq.rf, n = input.map$info$n.mrk - pos - 1)))
    } else if(pos  ==  n.mrk){
      seq.rf <- as.numeric(c(input.map$maps[[i.lpc]]$seq.rf, res[i, "rf1"]))
    }
    output.map$maps[[i]] <- list(seq.num = seq.num, 
                                 seq.rf = seq.rf, 
                                 seq.ph = configs[[rownames(res)[i]]],
                                 loglike = res[i, "log_like"])
  }
  output.map$info$seq.dose.p1 <- dat$dosage.p1[output.map$info$mrk.names]
  output.map$info$seq.dose.p2 <- dat$dosage.p2[output.map$info$mrk.names]
  output.map$info$chrom <- dat$chrom[output.map$info$mrk.names]
  output.map$info$genome.pos <- dat$genome.pos[output.map$info$mrk.names]
  output.map$info$seq.ref <-  dat$seq.ref[output.map$info$mrk.names]
  output.map$info$seq.alt <-  dat$seq.alt[output.map$info$mrk.names]
  output.map$info$chisq.pval <- dat$chisq.pval[output.map$info$mrk.names]
  return(output.map)
}
