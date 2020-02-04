#' Reestimate the recombination fractions in a genetic map
#' 
#' This function reestimates the recombination fractions between all markers in a given map.
#'
#' @param input.map An object of class \code{mappoly.map}
#' 
#' @param input.mat An object of class \code{mappoly.rf.matrix}
#' 
#' @param tol tolerance for determining convergence (default = 10e-03)
#' 
#' @param phase.config which phase configuration should be used. "best" (default) 
#'                     will choose the maximum likelihood configuration
#'                     
#' @param method indicates whether to use \code{'hmm'} (Hidden Markov Models) 
#' or \code{'ols'} (Ordinary Least Squares) to reestimate the recombination fractions
#'     
#' @param weight if \code{TRUE} (default), it uses the LOD scores to perform a weighted
#'    regression when the Ordinary Least Saquares is chosen
#'    
#' @param verbose if \code{TRUE} (default), current progress is shown; if
#'     \code{FALSE}, no output is produced
#'     
#' @return a numeric vector of size \code{m} indicating which
#'     homologous in h2 represents the homologous in h1. If there is
#'     no correspondence, i.e. different homologous, it returns NA for
#'     that homologous
#'     
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \url{https://doi.org/10.1534/g3.119.400378}     
#'     
#' @export reest_rf
#' 
reest_rf<-function(input.map, input.mat = NULL, tol = 10e-3,  phase.config = "all",
                    method = c("hmm", "ols"), weight = TRUE, verbose = TRUE)
{
  if (!inherits(input.map, "mappoly.map")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  }
  method<-match.arg(method)
  ## choosing the linkage phase configuration
  LOD.conf <- get_LOD(input.map, sorted = FALSE)
  if(phase.config == "best") {
    i.lpc <- which.min(LOD.conf)
  } else if(phase.config == "all"){
    i.lpc <- seq_along(LOD.conf) } else if (phase.config > length(LOD.conf)) {
    stop("invalid linkage phase configuration")
  } else i.lpc <- phase.config
  if(is.null(input.mat) & method == "ols")
    stop("'input.mat' is expected when 'method = ols'")
  output.map<-input.map
  if(method == "ols")
  {
    for(j in i.lpc){
      id<-get(input.map$info$data.name, pos =1)$mrk.names[input.map$maps[[i.lpc]]$seq.num]
      y<-as.numeric((imf_h(as.dist(input.mat$rec.mat[id, id]))))
      w<-as.numeric((imf_h(as.dist(input.mat$lod.mat[id, id]))))
      v<-t(combn(id,2))
      z<-cumsum(imf_h(c(0,input.map$maps[[i.lpc]]$seq.rf)))
      names(z)<-id
      x<-numeric(nrow(v))
      names(x)<-names(y)<-apply(v, 1, paste0, collapse="-")
      for(i in 1:nrow(v))
        x[i]<-z[v[i,2]]-z[v[i,1]]
      if(weight)
        model <- lm(y ~ x-1, weights=w)
      else
        model <- lm(y ~ x-1)
      new <- data.frame(x = z)
      u<-predict(model, new)
      output.map$maps[[j]]$seq.rf<-mf_h(diff(u))
      output.map$maps[[j]]$loglike<-0
    }
  }
  else if(method == "hmm")
  {
    for(j in i.lpc){
    s<-make_seq_mappoly(input.obj = get(input.map$info$data.name, pos =1),
                        arg = input.map$maps[[j]]$seq.num,
                        data.name = input.map$info$data.name)
    mtemp<-est_rf_hmm_single(input.seq = s,
                             input.ph.single = input.map$maps[[j]]$seq.ph,
                             rf.temp = input.map$maps[[j]]$seq.rf,
                             tol = tol, verbose = verbose)
    output.map$maps[[j]]<-mtemp
    }
  }
  return(output.map)
}
