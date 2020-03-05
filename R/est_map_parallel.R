#' Function to build a map in parallel with submaps
#'
#' This function estimates the map for a sequence using batches (submaps)
#' and parallel computing.
#' 
#' @param void
#'
#' @keywords internal
#'
#' @author Gabriel Gesteira, \email{gabrielgesteira@usp.br}
#' @export
est_map_parallel = function(data, markers, partial_tpt, n.batches = 4, n.cores = 4, platform = 'auto', start.set = 5, thres.twopt = 10, thres.hmm = 10, submap.size.diff = 3, thres.twopt2 = 3, thres.hmm2 = 10, phase.n.lim = 20, tol = 1e-3) {
    ## Getting and checking platform
    if (platform == 'auto'){
        platform = .Platform$OS.type
    }
    if (!any(platform %in% c('unix','windows'))){
        stop("This function is not supported in your platform. Current supported platforms are Unix and Windows.")
    }
    ## Checking n.batches
    if (n.batches ==1) stop("For only one batch, use the function 'est_rf_hmm_sequential' instead.")
    ## Subsetting markers into n.batches
    size_batches = vector(mode = "integer",length = n.batches)
    n = length(markers)%/%n.batches
    rem = length(markers)%%n.batches
    for (i in 1:length(size_batches)){
        n1 = n
        if (rem > 0) {
            n1 = n + 1
            rem = rem - 1
        }
        size_batches[i] = n1
    }
    ## Creating list to handle submaps
    submaps = list()
    for (i in 1:length(size_batches)){
        submaps[[i]] = list()
        submaps[[i]][[1]] = make_seq_mappoly(data, arg = markers[(sum(cumsum(size_batches)[i-1])+1):cumsum(size_batches)[i]], data.name = partial_tpt$data.name)
        submaps[[i]][[2]] = make_pairs_mappoly(partial_tpt, input.seq = submaps[[i]][[1]])
    }
    ## Running hmm_sequential for each submap
    if (platform == 'windows') cl = parallel::makeCluster(n.cores)
    else cl = parallel::makeCluster(n.cores, type = 'FORK')
    #parallel::clusterEvalQ(cl, require(mappoly))
    parallel::clusterExport(cl, varlist = c('est_rf_hmm_sequential'))
    #start.set2,thres.twopt2,thres.hmm2,submap.size.diff2,phase.n.lim2,tol2 / 'start.set','thres.twopt','thres.hmm','submap.size.diff','phase.n.lim','tol'
    submaps2 = parallel::parLapply(cl,submaps,function(x){
    y = est_rf_hmm_sequential(input.seq = x[[1]],
                              start.set = start.set,
                              thres.twopt = thres.twopt,
                              thres.hmm = thres.hmm,
                              info.tail = TRUE,
                              twopt = x[[2]],
                              sub.map.size.diff.limit = submap.size.diff,
                              phase.number.limit = phase.n.lim,
                              reestimate.single.ph.configuration = TRUE,
                              tol = tol,
                              tol.final = tol/10)
    return(y)
})
    parallel::stopCluster(cl)
    ## Running iteratively until one submap is left
    submaps_merged = submaps2
    while((length(submaps_merged) > 1) && (class(submaps_merged[[length(submaps_merged)]]) == "mappoly.map")){
        submaps_merged2 = submap_merge(submaps_merged, 
                                       tpt = partial_tpt, 
                                       n.cores = n.cores,
                                       platform = platform,
                                       thres.twopt2, 
                                       thres.hmm2,
                                       tol = tol)
        submaps_merged = submaps_merged2
    } 
    ## Reestimating rfs
    submaps_merged_reest = reest_rf(input.map = submaps_merged[[1]],
                                    input.mat = partial_tpt,
                                    phase.config = 'best',
                                    method = 'hmm', 
                                    verbose = F)
    return(submaps_merged)
}

#' Function to merge submaps two-by-two
#'
#' @param submaps a list of \code{mappoly.map} objects
#' @param tpt a two points object
#' @param n.cores number of cores
#' @param platform automatically detected
#' @param thres.twopt2 two-point LOD threshold to eliminate phases
#' @param thres.hmm2 LOD threshold to eliminate map phases
#' @param tol precision for hmm convergence
#'
#' @keywords internal
#'
#' @author Gabriel Gesteira, \email{gabrielgesteira@usp.br}
#'
submap_merge = function(submaps, tpt, n.cores, platform, thres.twopt2, thres.hmm2, tol)
{
    ## Getting submaps two-by-two
    combs = length(submaps)%/%2
    remaining = length(submaps)%%2 # It is supposed to be always 1 or 0  
    combinations = list()
    if (combs > 1) {
        j=1
        for (i in 1:combs){
            combinations[[i]] = list()
            combinations[[i]][[1]] = submaps[[j]]
            combinations[[i]][[2]] = submaps[[(j+1)]]
            j = j+2
        }
    } else {
        combinations[[1]] = list()
        combinations[[1]][[1]] = submaps[[1]]
        combinations[[1]][[2]] = submaps[[2]]
    }  
    ## Merging each submap (two by two)
    if (length(combinations) > 1){
        if (platform == 'windows') cl = parallel::makeCluster(n.cores)
        else cl = parallel::makeCluster(n.cores, type = 'FORK')
          parallel::clusterEvalQ(cl, require(mappoly))
        submaps3 = parallel::parLapply(cl,combinations,function(x,tpt2=tpt,thres.twopt22=thres.twopt2,thres.hmm22=thres.hmm2,tol2=tol){
    y = merge_maps(x[[1]],
                   x[[2]],
                   tpt2,
                   thres.twopt = thres.twopt22,
                   thres.hmm = thres.hmm22,
                   genoprob.map1 = NULL,
                   phase.config.map1 = "best",
                   phase.config.map2 = "best",
                   tol = tol2)
    return(y)
})
          parallel::stopCluster(cl)
        }
    else {
        submaps3 = list()
        submaps3[[1]] = merge_maps(combinations[[1]][[1]], 
                                   combinations[[1]][[2]],
                                   tpt,
                                   thres.twopt = thres.twopt2,
                                   thres.hmm = thres.hmm2, 
                                   genoprob.map1 = NULL, 
                                   phase.config.map1 = "best",
                                   phase.config.map2 = "best", 
                                   tol = 0.001)
    }    
    if (remaining > 0){
        submaps3[[(length(submaps3)+1)]] = submaps[[length(submaps)]]
    }
    return(submaps3)
}

##' Wrapper for function est_rf_hmm_sequential
##' 
##' @param void
##' 
##' @keywords internal
##' 
func.sequential = function(x,start.set,thres.twopt,thres.hmm,submap.size.diff,phase.n.lim,tol){
    y = est_rf_hmm_sequential(input.seq = x$seq,
                              start.set = start.set,
                              thres.twopt = thres.twopt,
                              thres.hmm = thres.hmm,
                              info.tail = TRUE,
                              twopt = x$tpt,
                              sub.map.size.diff.limit = submap.size.diff,
                              phase.number.limit = phase.n.lim,
                              reestimate.single.ph.configuration = TRUE,
                              tol = tol,
                              tol.final = tol/10)
    return(y)
}

##' Wrapper for function merge_maps
##' 
##' @param void
##' 
##' @keywords internal
##' 
func.merge = function(x,tpt,thres.twopt2,thres.hmm2,tol){
    y = merge_maps(x[[1]],
                   x[[2]],
                   tpt,
                   thres.twopt = thres.twopt2,
                   thres.hmm = thres.hmm2,
                   genoprob.map1 = NULL,
                   phase.config.map1 = "best",
                   phase.config.map2 = "best",
                   tol = tol)
    return(y)
}

##' Function to check gaps
##'
##' This function was created to facilitate gap checking and weird marker removal.
##' It is intended to be run in interactive mode.
##'
##' @keywords interval
##' 
##' @author Gabriel Gesteira, \email{gabrielgesteira@usp.br}
##'
## check_gaps = function(map, tpt, max.gap = 0.1, window = 10){
##     gaps = which(map$maps[[1]]$seq.rf >= max.gap)
##     if (length(gaps) >= 1){
##         mat.rf = rf_list_to_matrix(tpt)
##         mrk.to.rem = vector()
##         for (i in 1:length(gaps)){
##             keep = 1
##             mrk.up.pos = gaps[i] - round(window/2,0)
##             mrk.low.pos = gaps[i] + round(window/2,0)
##             range = map$info$mrk.names[mrk.up.pos:mrk.low.pos]
##             cat('\nA gap of ',map$maps[[1]]$seq.rf[gaps[i]],' was found between markers ', map$info$mrk.names[gaps[i]], ' and ', map$info$mrk.names[gaps[i]+1], '\nPlease look at the plot and choose a marker to remove.\n')
##             plot(mat.rf, ord = range, index = TRUE)
##             while (keep == 1){
##                 rm.mrk = readline('Enter the name of the marker to be removed (press ENTER to skip): ')
##                 if (rm.mrk == '') keep = 0
##                 else if (!rm.mrk %in% map$info$mrk.names){
##                     keep = 1
##                     cat('Please inform a valid marker name, or press ENTER to skip.\n')
##                 } else {
##                     keep = 0
##                     mrk.to.rem = c(mrk.to.rem, rm.mrk)
##                     rm.mrk = ''
##                 }
##             }
##         }
##         if (length(mrk.to.rem) >= 1){
##             cat('Removing selected markers: ', mrk.to.rem, '\n')
##             posits = which(!map$info$mrk.names %in% mrk.to.rem)
##             updated_map = get_submap(input.map = map,
##                                      mrk.pos = posits)
##             return(updated_map)
##         } else cat('\nNo markers were removed from the map.\n')
##     } else cat('\nNo gaps were found considering the informed value.\n')
## }


##' Check phase conformity
##'
##' This function is intended to check if all homologues in two given maps
##' are the same for both parents. The function reduces the sequence
##' to the common markers only before making the comparison.
##' 
##' @keywords internal
##'
##' @author Gabriel Gesteira
##' 
## check_phases = function(map1, map2, verbose = TRUE){
##     ## Defining arguments
##     ploidy = map1$info$m
##     mrks.1 = map1$info$mrk.names
##     mrks.2 = map2$info$mrk.names
##     common.mrk = intersect(mrks.1,mrks.2)
##     pos.mrk.1 = which(mrks.1 %in% common.mrk)
##     pos.mrk.2 = which(mrks.2 %in% common.mrk)
##         if (verbose) {
##         cat('\nMap 1 has ', length(mrks.1), ' markers.\n')
##         cat('Map 2 has ', length(mrks.2), ' markers.\n')
##         cat('There are ', length(common.mrk), 'common markers.\n')
##         cat('The following markers from Map 1 were discarded: ', mrks.1[which(!mrks.1 %in% common.mrk)], '\n')
##         cat('The following markers from Map 2 were discarded: ', mrks.2[which(!mrks.2 %in% common.mrk)], '\n')
##     }
##     ## Subsetting phases
##     phases.1.P = map1$maps[[1]]$seq.ph$P[pos.mrk.1]
##     phases.2.P = map2$maps[[1]]$seq.ph$P[pos.mrk.2]
##     phases.1.Q = map1$maps[[1]]$seq.ph$Q[pos.mrk.1]
##     phases.2.Q = map2$maps[[1]]$seq.ph$Q[pos.mrk.2]

##     ## Creating matrices to handle homologues
##     phases.1.P.hom = phases.2.P.hom = phases.1.Q.hom = phases.2.Q.hom = matrix(NA, length(common.mrk), ploidy)

##     ## Filling matrices
##     for (i in 1:length(common.mrk)){
##         phases.1.P.hom[i,phases.1.P[[i]]] = 1
##         phases.1.Q.hom[i,phases.1.Q[[i]]] = 1
##         phases.2.P.hom[i,phases.2.P[[i]]] = 1
##         phases.2.Q.hom[i,phases.2.Q[[i]]] = 1
##     }
##     phases.1.P.hom[which(is.na(phases.1.P.hom))] = 0
##     phases.1.Q.hom[which(is.na(phases.1.Q.hom))] = 0
##     phases.2.P.hom[which(is.na(phases.2.P.hom))] = 0
##     phases.2.Q.hom[which(is.na(phases.2.Q.hom))] = 0

##     ## Checking if homologues match
##     P = apply(phases.1.P.hom == phases.2.P.hom, 2, any)
##     Q = apply(phases.1.Q.hom == phases.2.Q.hom, 2, any)

##     return(rbind(P,Q))
## }
