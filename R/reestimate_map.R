#' Reestimate the recombination fractions in a genetic map
#'
#' @param input.map An object of class \code{mappoly.map}
#' @param input.mat An object of class \code{mappoly.rf.matrix}
#' @param tol tolerance for determining convergence.
#' @param phase.config should be a string \code{'best'} or the position of the
#'     configuration to be plotted. If \code{'best'}, plot the configuration
#'     with the highest likelihood.
#' @param method indicates whether to use Hidden Markov Models or Ordinary
#'     Least Squares to reestimate the recombination fraction
#' @param weight if \code{TRUE}, it uses the LOD scores to perform a weighted
#'    regression when the Ordinary Least Saquares is chosen
#' @param verbose if \code{TRUE}, current progress is shown; if
#'     \code{FALSE}, no output is produced.
#' @return a numeric vector of size \code{m} indicating which
#'     homologous in h2 represents the homologous in h1. If there is
#'     no correspondence, i.e. different homologous, it returns NA for
#'     that homologous.
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2018) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _submited_. \url{https://doi.org/10.1101/415232}     
#'     
#' @export reest_map
#' 
reest_map<-function(input.map, input.mat = NULL, tol = 10e-3,  phase.config = "best",
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
  } else if (phase.config > length(LOD.conf)) {
    stop("invalid linkage phase configuration")
  } else i.lpc <- phase.config
  if(is.null(input.mat) & method == "ols")
    stop("'input.mat' is expected when 'method = ols'")
  if(method == "ols")
  {
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
    output.map<-input.map
    output.map$maps<-output.map$maps[i.lpc]
    output.map$maps[[1]]$seq.rf<-mf_h(diff(u))
    output.map$maps[[1]]$loglike<-0
    return(output.map)
  }
  else if(method == "hmm")
  {
    s<-make_seq_mappoly(input.obj = get(input.map$info$data.name, pos =1),
                        arg = input.map$maps[[i.lpc]]$seq.num,
                        data.name = input.map$info$data.name)
    mtemp<-est_rf_hmm_single(input.seq = s,
                             input.ph.single = input.map$maps[[i.lpc]]$seq.ph,
                             rf.temp = input.map$maps[[i.lpc]]$seq.rf,
                             tol = tol, verbose = verbose)
    output.map<-input.map
    output.map$maps<-output.map$maps[i.lpc]
    output.map$maps[[1]]<-mtemp
    return(output.map)
  }
}
