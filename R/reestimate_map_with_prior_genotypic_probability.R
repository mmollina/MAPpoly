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
#'     submap <- get_submap(solcap.dose.map[[1]], mrk.pos = 1:20, verbose = FALSE)
#'     prob.submap <- est_full_hmm_with_prior_prob(submap,
#'                                                 dat.prob = tetra.solcap.geno.dist,
#'                                                 tol=10e-4, 
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
#'     \url{https://doi.org/10.1534/g3.119.400378}
#'
#' @export est_full_hmm_with_prior_prob
#'
est_full_hmm_with_prior_prob<-function(input.map, dat.prob = NULL, phase.config = "best", 
                                       tol = 10e-4, verbose = FALSE)
{
  if (!inherits(input.map, "mappoly.map")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  }
  ## choosing the linkage phase configuration
  LOD.conf <- get_LOD(input.map, sorted = FALSE)
  if(phase.config == "best") {
    i.lpc <- which.min(LOD.conf)
  } else if (phase.config > length(LOD.conf)) {
    stop("invalid linkage phase configuration")
  } else i.lpc <- phase.config
  if(is.null(dat.prob)){
    if(nrow(get(input.map$info$data.name, pos=1)$geno)==get(input.map$info$data.name, pos=1)$n.mrk){
      stop("
            The dataset associated to 'input.map'
            contains no genotypic probability distribution.
            Please provide provide dataset in argument 
           'dat.prob'.
           ")
    } else {
      dat.prob <- get(input.map$info$data.name, pos=1)
      }
  }
  mrk<-NULL
  original.map.mrk <- input.map$info$mrk.names
  dat.prob.pos <- match(original.map.mrk, dat.prob$mrk.names)
  which.is.na<-which(is.na(dat.prob.pos))
  if(length(which.is.na) > 0)
    stop("Markers", original.map.mrk[which.is.na], "are not present in the 'dat.prob' object")
  temp.map<-input.map
  temp.map$info$seq.num <- temp.map$maps[[i.lpc]]$seq.num<-dat.prob.pos
  names(temp.map$maps[[i.lpc]]$seq.ph$P)<-names(temp.map$maps[[i.lpc]]$seq.ph$Q)<-dat.prob.pos
  if(!all(sort(get(temp.map$info$data.name, pos = 1)$ind.names) %in% sort(get(input.map$info$data.name, pos = 1)$ind.names)))
    stop("The individuals are different in the new and original datasets")
  geno<-subset(dat.prob$geno, mrk%in%original.map.mrk)
  geno.new<-NULL
  for(i in unique(geno$ind))
    geno.new<-rbind(geno.new, geno[geno[,"ind"] == i, ][match(original.map.mrk, geno[,"mrk"]),])
  g <- as.double(t(geno.new[, -c(1:2)]))
  if (verbose) cat("
 ----------------------------------------------
 INFO: running HMM using full transition space:
       this operation may take a while.
-----------------------------------------------\n")
  map.res<-poly_hmm_est(m = as.numeric(temp.map$info$m),
                        n.mrk = as.numeric(temp.map$info$n.mrk),
                        n.ind = get(input.map$info$data.name, pos=1)$n.ind,
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
  temp.map$maps[[i.lpc]]$seq.rf<-map.res$rf
  temp.map$maps[[i.lpc]]$loglike<-map.res$loglike
  return(temp.map)
}
