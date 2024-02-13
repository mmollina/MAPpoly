#' Prior probability for genotyping error
#'
#' If \code{restricted = TRUE}, it restricts the prior to the 
#' possible classes under Mendelian, non double-reduced segregation 
#' given dosage of the parents
#' @keywords internal
genotyping_global_error <- function(x, ploidy, restricted = TRUE,  error = 0.01, th.prob = 0.95)
{
  if(restricted){
    x1 <- x[1:(ploidy+1)]
    if(sum(x1 > th.prob)  ==  1){
      x2 <- x[ploidy + 2:3]
      id <- segreg_poly(ploidy, dP = x2[1], dQ = x2[2]) > 0
      x3 <- x1[id]
      o <- which.max(x3)
      x3[o] <- 1-error
      x3[-o] <- error/(sum(id)-1)
      x1[match(names(x3), names(x1))] <- x3
      return(x1)
    }
    return(x1)
  } else {
    x1 <- x[1:(ploidy+1)]
    if(sum(x1 > th.prob) == 1){
      o <- which.max(x1)
      x1[o] <- 1-error
      x1[-o] <- error/(length(x1)-1)
      return(x1)
    }
    return(x1)
  }
}

#' Re-estimate genetic map given a global genotyping error
#'
#' This function considers a global error when re-estimating
#' a genetic map using Hidden Markov models. Since this function 
#' uses the whole transition space in the HMM, its computation 
#' can take a while, especially for hexaploid maps. 
#'
#' @param input.map an object of class \code{mappoly.map}
#' 
#' @param error the assumed global error rate (default = NULL)
#' 
#' @param tol the desired accuracy (default = 10e-04)
#' 
#' @param restricted if \code{TRUE} (default), restricts the prior to the 
#'                   possible classes under Mendelian, non double-reduced segregation 
#'                   given dosage of the parents
#'                    
#' @param th.prob the threshold for using global error or genotype 
#'     probability distribution if present in the dataset (default = 0.95)
#'      
#' @param verbose if \code{TRUE}, current progress is shown; if
#'     \code{FALSE} (default), no output is produced
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
#' @examples
#'     submap <- get_submap(solcap.dose.map[[1]], mrk.pos = 1:20, verbose = FALSE)
#'     err.submap <- est_full_hmm_with_global_error(submap, 
#'                                                  error = 0.01, 
#'                                                  tol = 10e-4, 
#'                                                  verbose = TRUE)
#'     err.submap
#'     plot_map_list(list(dose = submap, err = err.submap), 
#'                   title = "estimation procedure")
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
#' @export est_full_hmm_with_global_error
est_full_hmm_with_global_error <- function(input.map, error = NULL, tol = 10e-4, 
                                           restricted = TRUE, 
                                           th.prob = 0.95, 
                                           verbose = FALSE)
{
  if (!inherits(input.map, "mappoly.map")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  }
  output.seq <- input.map
  mrk.names <- get(input.map$info$data.name, pos = 1)$mrk.names[input.map$maps[[1]]$seq.num]
  ## 
  if(!mappoly::is.prob.data(get(input.map$info$data.name, pos = 1))){
    geno.temp <- get(input.map$info$data.name, pos = 1)$geno.dose[mrk.names,]
    ind.names <- get(input.map$info$data.name, pos = 1)$ind.names
    gen <- vector("list", length(ind.names))
    names(gen) <- ind.names
    mrk <- ind <- NULL
    dp <- get(input.map$info$data.name, pos = 1)$dosage.p1[input.map$maps[[1]]$seq.num]
    dq <- get(input.map$info$data.name, pos = 1)$dosage.p2[input.map$maps[[1]]$seq.num]
    names(dp) <- names(dq) <- mrk.names
    d.pq <- data.frame(dp = dp, 
                     dq = dq)
    d.pq$mrk <- mrk.names
    for(i in names(gen))
    {
      a <- matrix(0, nrow(geno.temp), input.map$info$ploidy+1, dimnames = list(mrk.names, 0:input.map$info$ploidy))
      for(j in rownames(a)){
        if(geno.temp[j,i]  ==  input.map$info$ploidy+1){
          a[j,] <- segreg_poly(ploidy = input.map$info$ploidy, dP = dp[j], dQ = dq[j])
        } else {
          a[j,geno.temp[j,i]+1] <- 1          
        }
      }
      a <- as.data.frame(a)
      a$mrk <- rownames(a)
      a.temp <- t(merge(a, d.pq, sort = FALSE)[,-c(1)])
      if(!is.null(error))
        a.temp <- apply(a.temp, 2, genotyping_global_error, ploidy = input.map$info$ploidy, 
                      restricted = restricted, error = error, th.prob = th.prob)
      else
        a.temp <- a.temp[1:(input.map$info$ploidy+1), ]
      colnames(a.temp) <- a[,1]
      gen[[i]] <- a.temp
    }
  } 
  else {##
    geno.temp <- subset(get(input.map$info$data.name, pos = 1)$geno, mrk%in%mrk.names)
    ind.names <- get(input.map$info$data.name, pos = 1)$ind.names
    gen <- vector("list", length(ind.names))
    names(gen) <- ind.names
    mrk <- ind <- NULL
    d.pq <- data.frame(dp = get(input.map$info$data.name, pos = 1)$dosage.p1[input.map$maps[[1]]$seq.num], 
                     dq = get(input.map$info$data.name, pos = 1)$dosage.p2[input.map$maps[[1]]$seq.num])
    d.pq$mrk <- mrk.names
    for(i in names(gen))
    {
      a <- subset(geno.temp, ind%in%i)
      a <- a[match(mrk.names, a$mrk),]
      a.temp <- t(merge(a, d.pq, sort = FALSE)[,-c(1:2)])
      if(!is.null(error))
        a.temp <- apply(a.temp, 2, genotyping_global_error, ploidy = input.map$info$ploidy, 
                      restricted = restricted, error = error, th.prob = th.prob)
      else
        a.temp <- a.temp[1:(input.map$info$ploidy+1), ]
      colnames(a.temp) <- a[,1]
      gen[[i]] <- a.temp
    }
  }
  if (verbose) cat("
 ----------------------------------------------
 INFO: running HMM using full transition space:
       this operation may take a while.
-----------------------------------------------\n")
  for(i in 1:length(input.map$maps))
  {
    YP <- input.map$maps[[i]]$seq.ph$P
    YQ <- input.map$maps[[i]]$seq.ph$Q
    map <- poly_hmm_est(ploidy = as.numeric(input.map$info$ploidy),
                      n.mrk = as.numeric(input.map$info$n.mrk),
                      n.ind = as.numeric(length(gen)),
                      p = as.numeric(unlist(YP)),
                      dp = as.numeric(cumsum(c(0, sapply(YP, function(x) sum(length(x)))))),
                      q = as.numeric(unlist(YQ)),
                      dq = as.numeric(cumsum(c(0, sapply(YQ, function(x) sum(length(x)))))),
                      g = as.double(unlist(gen)),
                      rf = as.double(input.map$maps[[i]]$seq.rf),
                      verbose = verbose,
                      tol = tol)
    output.seq$maps[[i]]$seq.rf <- map$rf
    output.seq$maps[[i]]$loglike <- map$loglike
  }
  return(output.seq)
}

#' Re-estimate genetic map using dosage prior probability distribution
#'
#' This function considers dosage prior distribution when re-estimating
#' a genetic map using Hidden Markov models
#'
#' @param input.map an object of class \code{mappoly.map}
#' 
#' @param dat.prob an object of class \code{mappoly.data} containing the
#'                 probability distribution of the genotypes
#'                 
#' @param phase.config which phase configuration should be used. "best" (default) 
#'                     will choose the maximum likelihood configuration
#'                     
#' @param tol the desired accuracy (default = 10e-04)
#' 
#' @param verbose if \code{TRUE}, current progress is shown; if
#'     \code{FALSE} (default), no output is produced
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
#' @examples
#'     submap <- get_submap(solcap.dose.map[[1]], mrk.pos = 1:20, verbose = FALSE)
#'     prob.submap <- est_full_hmm_with_prior_prob(submap,
#'                                                 dat.prob = tetra.solcap.geno.dist,
#'                                                 tol = 10e-4, 
#'                                                 verbose = TRUE)
#'     prob.submap
#'     plot_map_list(list(dose = submap, prob = prob.submap), 
#'                   title = "estimation procedure")
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
#' @export est_full_hmm_with_prior_prob
#'
est_full_hmm_with_prior_prob <- function(input.map, dat.prob = NULL, phase.config = "best", 
                                         tol = 10e-4, verbose = FALSE)
{
  if (!inherits(input.map, "mappoly.map")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  }
  ## choosing the linkage phase configuration
  LOD.conf <- get_LOD(input.map, sorted = FALSE)
  if(phase.config  ==  "best") {
    i.lpc <- which.min(LOD.conf)
  } else if (phase.config > length(LOD.conf)) {
    stop("invalid linkage phase configuration")
  } else i.lpc <- phase.config
  if(is.null(dat.prob)){
    if(nrow(get(input.map$info$data.name, pos = 1)$geno) == get(input.map$info$data.name, pos = 1)$n.mrk){
      stop("
            The dataset associated to 'input.map'
            contains no genotypic probability distribution.
            Please provide provide dataset in argument 
           'dat.prob'.
           ")
    } else {
      dat.prob <- get(input.map$info$data.name, pos = 1)
    }
  }
  mrk <- NULL
  original.map.mrk <- input.map$info$mrk.names
  dat.prob.pos <- match(original.map.mrk, dat.prob$mrk.names)
  which.is.na <- which(is.na(dat.prob.pos))
  if(length(which.is.na) > 0)
    stop("Markers", original.map.mrk[which.is.na], "are not present in the 'dat.prob' object")
  temp.map <- input.map
  temp.map$info$seq.num <- temp.map$maps[[i.lpc]]$seq.num <- dat.prob.pos
  names(temp.map$maps[[i.lpc]]$seq.ph$P) <- names(temp.map$maps[[i.lpc]]$seq.ph$Q) <- dat.prob.pos
  if(!all(sort(get(temp.map$info$data.name, pos = 1)$ind.names) %in% sort(get(input.map$info$data.name, pos = 1)$ind.names)))
    stop("The individuals are different in the new and original datasets")
  geno <- subset(dat.prob$geno, mrk%in%original.map.mrk)
  geno.new <- NULL
  for(i in unique(geno$ind))
    geno.new <- rbind(geno.new, geno[geno[,"ind"]  ==  i, ][match(original.map.mrk, geno[,"mrk"]),])
  g <- as.double(t(geno.new[, -c(1:2)]))
  if (verbose) cat("
 ----------------------------------------------
 INFO: running HMM using full transition space:
       this operation may take a while.
-----------------------------------------------\n")
  map.res <- poly_hmm_est(ploidy = as.numeric(temp.map$info$ploidy),
                          n.mrk = as.numeric(temp.map$info$n.mrk),
                          n.ind = get(input.map$info$data.name, pos = 1)$n.ind,
                          p = as.numeric(unlist(temp.map$maps[[i.lpc]]$seq.ph$P)),
                          dp = as.numeric(cumsum(c(0, sapply(temp.map$maps[[i.lpc]]$seq.ph$P, function(x) sum(length(x)))))),
                          q = as.numeric(unlist(temp.map$maps[[i.lpc]]$seq.ph$Q)),
                          dq = as.numeric(cumsum(c(0, sapply(temp.map$maps[[i.lpc]]$seq.ph$Q, function(x) sum(length(x)))))),
                          g = g,
                          rf = temp.map$maps[[i.lpc]]$seq.rf,
                          verbose = verbose,
                          tol = tol)
  if(verbose)
    cat("\n")
  temp.map$maps[[i.lpc]]$seq.rf <- map.res$rf
  temp.map$maps[[i.lpc]]$loglike <- map.res$loglike
  return(temp.map)
}

