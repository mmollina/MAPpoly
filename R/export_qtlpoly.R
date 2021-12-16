#' Export to QTLpoly 
#' 
#' Compute homolog probabilities for all individuals in the full-sib
#' population given a map and conditional genotype probabilities, and exports
#' the results to be used for QTL mapping in the QTLpoly package.
#'
#' @param input.genoprobs an object of class \code{mappoly.genoprob}
#'
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
#' 
#'@examples
#'    \donttest{
#'      ## tetraploid example
#'      w1 <- calc_genoprob(solcap.dose.map[[1]])
#'      h.prob <- export_qtlpoly(w1)
#'   }
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari M., Olukolu B. A.,  Pereira G. da S., 
#'     Khan A., Gemenet D., Yencho G. C., Zeng Z-B. (2020), 
#'     Unraveling the Hexaploid Sweetpotato Inheritance 
#'     Using Ultra-Dense Multilocus Mapping, 
#'     _G3: Genes, Genomes, Genetics_. 
#'     \doi{10.1534/g3.119.400620} 
#'     
#' @export export_qtlpoly
#' 
export_qtlpoly <- function(input.genoprobs, verbose = TRUE){
  if(inherits(input.genoprobs, "mappoly.genoprob")) 
    input.genoprobs <- list(input.genoprobs)
  if(!inherits(input.genoprobs, "list"))
    stop(deparse(substitute(input.genoprobs)), 
         " is not an object of class 'mappoly.genoprob' neither a list containing 'mappoly.genoprob' objects.")
  if (any(!sapply(input.genoprobs, inherits, "mappoly.genoprob"))) 
    stop(deparse(substitute(input.genoprobs)), 
         " is not an object of class 'mappoly.genoprob' neither a list containing 'mappoly.genoprob' objects.")
  res <- NULL
  for(j in 1:length(input.genoprobs)){
    if (verbose) cat("\nLinkage group ", j, "...")
    stt.names <- dimnames(input.genoprobs[[j]]$probs)[[1]] ## state names
    mrk.names <- dimnames(input.genoprobs[[j]]$probs)[[2]] ## mrk names
    ind.names <- dimnames(input.genoprobs[[j]]$probs)[[3]] ## individual names
    v <- c(2,4,6,8,10,12)
    names(v) <- choose(v,v/2)^2
    ploidy <- v[as.character(length(stt.names))]
    hom.prob <- array(NA, dim = c(ploidy*2, length(mrk.names), length(ind.names)))
    dimnames(hom.prob) <- list(letters[1:(2*ploidy)], mrk.names, ind.names)
    for(i in letters[1:(2*ploidy)])
      hom.prob[i,,] <- apply(input.genoprobs[[j]]$probs[grep(stt.names, pattern = i),,], c(2,3), function(x) round(sum(x, na.rm = TRUE),4))
    map <- data.frame(map.position = input.genoprobs[[j]]$map, marker = names(input.genoprobs[[j]]$map))
    out <- list(LG = j, map = map, probs = hom.prob)
    res <- c(res, out)
  }
  if (verbose) cat("\n")
  structure(list(info = list(ploidy = ploidy, n.ind = length(ind.names)) , homoprob = res), class = "mappoly.qtlpoly")
}
