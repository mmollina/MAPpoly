#' Reestimate genetic map using dosage prior distribution
#'
#' This function considers dosage prior distribution when reestimating
#' a genetic map using Hidden Markov models
#'
#' @param input.map an object of class \code{mappoly.map}
#' 
#' @param phase.config which phase configuration should be used. "best" (default) 
#'                     will choose the maximum likelihood configuration
#'                     
#' @param tol the desired accuracy (default = 10e-04)
#' 
#' @param verbose if \code{TRUE}, current progress is shown; if
#'     \code{FALSE} (default), no output is produced
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
#' @examples
#'   \dontrun{
#'   
#'   solcap.p<-vector("list", 12)
#'   names(solcap.p)<-names(solcap.dose.map)
#'   for(i in 1:12)
#'     solcap.p[[i]] <- est_full_hmm_with_prior_dist(solcap.dose.map[[i]], verbose = FALSE)
#'  w<-NULL
#'  for(i in 1:12)
#'    w<-c(w, c(solcap.dose.map[i], 
#'              solcap.p[i]))
#'              
#'  names(w) <- apply(expand.grid(c("dose", "prior"), paste0("LG_", 1:12), 
#'                              stringsAsFactors = FALSE)[,2:1], 1, paste, 
#'                  collapse = "_")
#'                  
#'  op <- par(cex.axis = .7)
#'  plot_map_list(w, horiz = FALSE, col = rep(gg_color_hue(2), 12))
#'  par(op)
#'  legend("bottomright", legend = c("Dosage based", "Prior"), pch=15, col = rep(gg_color_hue(2)))
#'   
#' }
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
#' @export est_full_hmm_with_prior_dist
#'
est_full_hmm_with_prior_dist<-function(input.map, phase.config = "best", 
                                       tol = 10e-4, verbose = TRUE)
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
    if(nrow(get(input.map$info$data.name, pos=1)$geno)==get(input.map$info$data.name, pos=1)$n.mrk) 
      stop("
            The data set associated to 'input.map'
            contains no genotypic probability distribution.
            Please provide provide data set in argument 
           'dat.dist'.
           ")
    dat.dist <- get(input.map$info$data.name, pos=1)
    mrk<-NULL
    original.map.mrk<-get(input.map$info$data.name, pos=1)$mrk.names[input.map$maps[[i.lpc]]$seq.num]
    dat.dist.pos<-match(original.map.mrk, dat.dist$mrk.names)
    which.is.na<-which(is.na(dat.dist.pos))
    if(length(which.is.na) > 0)
      stop("Markers", original.map.mrk[which.is.na], "are not present in the 'dat.dist' object")
    temp.map<-input.map
    temp.map$maps[[i.lpc]]$seq.num<-dat.dist.pos
    names(temp.map$maps[[i.lpc]]$seq.ph$P)<-names(temp.map$maps[[i.lpc]]$seq.ph$Q)<-dat.dist.pos
    if(!all(sort(get(temp.map$info$data.name, pos = 1)$ind.names) %in% sort(get(input.map$info$data.name, pos = 1)$ind.names)))
      stop("The individuals in the new data set are not contained in the original data set")
    geno<-subset(get(temp.map$info$data.name, pos = 1)$geno, mrk%in%original.map.mrk)
    geno.new<-NULL
    for(i in unique(geno$ind))
      geno.new<-rbind(geno.new, geno[geno[,"ind"] == i, ][match(original.map.mrk, geno[,"mrk"]),])
    g <- as.double(t(geno.new[, -c(1:2)]))
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
    temp.map$maps[[i.lpc]]$seq.rf<-map.res$rf
    temp.map$maps[[i.lpc]]$loglike<-map.res$loglike
    return(temp.map)
}
