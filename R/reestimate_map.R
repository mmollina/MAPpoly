#' Re-estimate the recombination fractions in a genetic map
#' 
#' This function re-estimates the recombination fractions between all markers in a given map.
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
#' @param method indicates whether to use \code{'hmm'} (Hidden Markov Models), 
#'  \code{'ols'} (Ordinary Least Squares) or \code{'wMDS_to_1D_pc'} (weighted MDS 
#'  followed by fitting a one dimensional principal curve) to re-estimate the 
#'  recombination fractions.
#'     
#' @param weight if \code{TRUE} (default), it uses the LOD scores to perform a weighted
#'    regression when the Ordinary Least Squares is chosen
#'    
#' @param verbose if \code{TRUE} (default), current progress is shown; if
#'     \code{FALSE}, no output is produced
#'     
#' @param high.prec logical. If \code{TRUE} uses high precision 
#' (long double) numbers in the HMM procedure implemented in C++,
#' which can take a long time to perform (default = FALSE)  
#'    
#' @param max.rf.to.break.EM for internal use only.   
#'     
#' @param input.mds An object of class \code{mappoly.map}
#'      
#' @return An updated object of class \code{mappoly.pcmap} whose 
#'         order was used in the \code{input.map}
#'     
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \url{https://doi.org/10.1534/g3.119.400378}     
#'     
#'     Stam P (1993) Construction of integrated genetic-linkage maps 
#'     by means of a new computer package: Joinmap. _Plant J_ 3:739–744
#'     \url{https://doi.org/10.1111/j.1365-313X.1993.00739.x}
#'     
#' @export reest_rf
#' 
reest_rf<-function(input.map, input.mat = NULL, tol = 10e-3,  phase.config = "all",
                   method = c("hmm", "ols", "wMDS_to_1D_pc"), weight = TRUE, verbose = TRUE, 
                   high.prec = FALSE, max.rf.to.break.EM = 0.5, input.mds = NULL)
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
    }
    output.map <- loglike_hmm(output.map)
  }
  else if(method == "hmm")
  {
    s<-make_seq_mappoly(input.obj = get(input.map$info$data.name, pos =1),
                        arg = input.map$info$mrk.names,
                        data.name = input.map$info$data.name)
    for(j in i.lpc){
      mtemp<-est_rf_hmm_single(input.seq = s,
                               input.ph.single = input.map$maps[[j]]$seq.ph,
                               rf.temp = input.map$maps[[j]]$seq.rf,
                               tol = tol, verbose = verbose, 
                               high.prec = high.prec,
                               max.rf.to.break.EM = max.rf.to.break.EM)
      output.map$maps[[j]]<-mtemp
    }
  } else if(method == "wMDS_to_1D_pc")
  {
    if(is.null(input.mds))
      stop("Provide 'mds.mappoly' object.")
    if (!inherits(input.mds, "mappoly.pcmap")) {
      stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.pcmap'")
    }
    pos <- input.mds$locimap[match(input.map$info$mrk.names, input.mds$locimap$locus), "position"]
    for(j in i.lpc){
      output.map$maps[[j]]$seq.rf <- mf_h(diff(pos))
    }
    output.map <- loglike_hmm(output.map)
  }
  return(output.map)
}
