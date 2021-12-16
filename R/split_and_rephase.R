#' Divides map in sub-maps and re-phase them
#'
#' The function splits the input map in sub-maps 
#' given a distance threshold of neighboring markers 
#' and evaluates alternative phases between the sub-maps.
#'
#' @param input.map an object of class \code{mappoly.map}
#' 
#' @param twopt an object of class \code{mappoly.twopt}
#'     containing the two-point information for the markers contained 
#'     in \code{input.map}
#'     
#' @param gap.threshold distance threshold of neighboring markers 
#'                      where the map should be spitted. The default 
#'                      value is 5 cM
#' 
#' @param size.rem.cluster the size of the marker cluster (in number of markers) 
#'                         from which the cluster should be removed. The default 
#'                         value is 1
#' 
#' @param phase.config which phase configuration should be used. "best" (default) 
#'                     will choose the maximum likelihood phase configuration 
#'                     
#' @param tol.final the desired accuracy for the final map (default = 10e-04)   
#'
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
#' 
#' @return An object of class \code{mappoly.map}
#' 
#' @examples
#'  map <- get_submap(solcap.dose.map[[1]], 1:20, verbose = FALSE)
#'  tpt <- est_pairwise_rf(make_seq_mappoly(map))
#'  new.map <- split_and_rephase(map, tpt, 1, 1)
#'  map
#'  new.map
#'  plot_map_list(list(old.map = map, new.map = new.map), col = "ggstyle")
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
#' @export split_and_rephase
split_and_rephase <- function(input.map,
                              twopt,
                              gap.threshold = 5, 
                              size.rem.cluster = 1,
                              phase.config = "best",
                              tol.final = 10e-4,
                              verbose = TRUE){
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
  id <- which(imf_h(input.map$maps[[i.lpc]]$seq.rf) > gap.threshold)
  if(length(id) == 0){
    if(verbose) cat("no submaps found\n")
    return(input.map)
  } 
  id <- cbind(c(1, id+1), c(id, input.map$info$n.mrk))
  temp.map <- input.map
  ## Selecting map segments larger then the specified threshold
  segments <- id[apply(id, 1, diff) > size.rem.cluster - 1, , drop = FALSE]
  if(length(segments) == 0) stop("all markers were eliminated\n")
  ## Dividing map in sub-maps
  temp.maps <- vector("list", nrow(segments))
  if (verbose) {
    ns <- nrow(segments)
    if(ns == 1){
      cat("one submap found ...\n")
      map <- get_submap(input.map, c(segments[1, 1]:segments[1, 2]), tol.final = tol.final, verbose = FALSE)
      return(filter_map_at_hmm_thres(map, 10e-4))
    } 
    else cat(ns, "submaps found ...\n")
  }
  for(i in 1:length(temp.maps)){
    temp.id <- c(segments[i, 1]:segments[i, 2])
    if(length(temp.id) > 1)
      temp.maps[[i]] <- get_submap(input.map, temp.id, reestimate.rf = FALSE, verbose = FALSE)
    else
      temp.maps[[i]] <- input.map$info$mrk.names[temp.id]    
  }
  newmap <- temp.maps[[1]]
  for(i in 2:length(temp.maps)){
    if (verbose) cat("adding block", i, "of", length(temp.maps), "\n")
    if(is.character(newmap) & is.character(temp.maps[[i]])){
      s <- make_seq_mappoly(get(input.map$info$data.name, pos = 1), 
                            c(newmap, temp.maps[[i]]), 
                            input.map$info$data.name)
      tpt <- est_pairwise_rf(s, verbose = FALSE)
      newmap <- est_rf_hmm_sequential(s,tpt, verbose = FALSE)
    } 
    else if(is.character(newmap) & !is.character(temp.maps[[i]])){
      newmap <- add_marker(input.map = temp.maps[[i]], mrk = newmap, pos = 0,
                           rf.matrix = rf_list_to_matrix(twopt,
                                                         thresh.LOD.ph = 5, 
                                                         thresh.LOD.rf = 5,
                                                         shared.alleles = TRUE),
                           verbose = FALSE)
    } 
    else if(!is.character(newmap) & is.character(temp.maps[[i]])){
      newmap <- add_marker(newmap, temp.maps[[i]], 
                           newmap$info$n.mrk, 
                           rf.matrix = rf_list_to_matrix(twopt,
                                                         thresh.LOD.ph = 5, 
                                                         thresh.LOD.rf = 5,
                                                         shared.alleles = TRUE),
                           verbose = FALSE)
    } 
    else {
      newmap <- merge_maps(list(newmap, temp.maps[[i]]), 
                           twopt = twopt, 
                           thres.twopt = 10, 
                           thres.hmm = 50)
    }
    newmap <- filter_map_at_hmm_thres(newmap, thres.hmm = 0.01)
  }
  newmap$info$seq.num <- newmap$maps[[1]]$seq.num
  newmap$info$seq.ref <- temp.map$info$seq.ref[newmap$info$mrk.names]
  newmap$info$seq.alt <- temp.map$info$seq.alt[newmap$info$mrk.names]
  new.map <- reest_rf(input.map = newmap, tol = tol.final, verbose = FALSE)
  return(new.map)
}
