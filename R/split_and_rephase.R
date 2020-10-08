#' Divides map in sub-maps and re-phase them
#'
#' The function splits the input map in sub-maps 
#' given a distance threshold of neighboring markers 
#' and evaluates alternative phases between the sub-maps.
#'
#' @param input.map an object of class \code{mappoly.map}
#' 
#' @param twopt an object of class \code{poly.est.two.pts.pairwise}
#'     containing the two-point information for the markers contained 
#'     in \code{input.map}
#'     
#' @param gap.threshold distance threshold of neighboring markers 
#'                      where the map should be spitted. The default 
#'                      value is 5 cM
#' 
#' @param remove.single Should isolated markers be removed?
#' 
#' @param phase.config which phase configuration should be used. "best" (default) 
#'                     will choose the phase configuration associated with the
#'                     maximum likelihood
#' @param tol.final the desired accuracy for the final map (default = 10e-04)   
#'
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
#' 
#' @return An object of class \code{mappoly.map}
#' 
#' @examples
#'  map <- get_submap(maps.hexafake[[1]], 1:20, reestimate.rf = FALSE, reestimate.phase = FALSE)
#'  tpt <- est_pairwise_rf(make_seq_mappoly(map))
#'  new.map <- split_and_rephase(map, tpt, 5)
#'  map
#'  new.map
#'  plot_map_list(list(old.map = map, new.map = new.map))
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
#' @export split_and_rephase
#'
#' @importFrom utils capture.output
#' 
split_and_rephase<-function(input.map,
                            twopt,
                            gap.threshold = 5, 
                            remove.single = TRUE,
                            phase.config = "best",
                            tol.final = 10e-4,
                            verbose = TRUE){
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
  id <- which(imf_h(input.map$maps[[i.lpc]]$seq.rf) > gap.threshold)
  if(length(id)==0) return(input.map)
  id <- cbind(c(1, id+1), c(id, input.map$info$n.mrk))
  ## Removing single markers at the beginning of the group
  i <- 1
  while(diff(id[i, ]) == 0){
    i <- i + 1
  }
  invisible(capture.output(suppressMessages(
    input.map <- get_submap(input.map, i:input.map$info$n.mrk, reestimate.rf = FALSE))))
  ## Dividing map in submaps
  id <- which(imf_h(input.map$maps[[1]]$seq.rf) > gap.threshold)
  id <- cbind(c(1, id+1), c(id, input.map$info$n.mrk))
  temp.maps <- vector("list", nrow(id))
  if (verbose) cat(nrow(id), "submaps found ...\n")
  for(i in 1:length(temp.maps)){
    temp.id <- c(id[i, 1]:id[i, 2])
    if(length(temp.id) > 1)
      invisible(capture.output(suppressMessages(
        temp.maps[[i]] <- get_submap(input.map, temp.id, reestimate.rf = FALSE))))
    else
      temp.maps[[i]] <- input.map$info$mrk.names[temp.id]    
  }
  newmap<-temp.maps[[1]]
  for(i in 2:length(temp.maps)){
    if (verbose) cat("Adding block", i, "of", length(temp.maps), "\n")
    if(!is.character(temp.maps[[i]])){
      invisible(capture.output(suppressMessages(
        newmap<-merge_maps(list(newmap, temp.maps[[i]]), 
                           twopt = twopt, 
                           thres.twopt = 10, 
                           thres.hmm = 50))))
      newmap <- filter_map_at_hmm_thres(newmap, thres.hmm = 0.01)
    } else {
      if(remove.single) next()
      invisible(capture.output(suppressMessages(
        maptemp <- add_marker(newmap, temp.maps[[i]], 
                              newmap$info$n.mrk, 
                              rf.matrix = rf_list_to_matrix(twopt,
                                                            thresh.LOD.ph = 5, 
                                                            thresh.LOD.rf = 5,
                                                            shared.alleles = TRUE),
                              verbose = verbose))))
      newmap <- filter_map_at_hmm_thres(maptemp, thres.hmm = 0.01)
    }
  }
  new.map <- reest_rf(input.map = newmap, tol = tol.final, verbose = FALSE)
  return(new.map)
}
