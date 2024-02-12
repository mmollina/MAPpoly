#' Design linkage map framework in two steps: i) estimating the recombination fraction with 
#' HMM approach for each parent separately using only markers segregating individually 
#' (e.g. map 1 - P1:3 x P2:0, P1: 2x4; map 2 - P1:0 x P2:3, P1:4 x P2:2); ii) merging both 
#' maps and re-estimate recombination fractions.
#' 
#' @param input.seq object of class \code{mappoly.sequence}
#' @param twopt object of class \code{mappoly.twopt}
#' @param start.set number of markers to start the phasing procedure (default = 4)
#' @param thres.twopt the LOD threshold used to determine if the linkage phases compared via two-point 
#' analysis should be considered for the search space reduction (default = 5)
#' @param thres.hmm the LOD threshold used to determine if the linkage phases compared via hmm analysis 
#' should be evaluated in the next round of marker inclusion (default = 50)
#' @param extend.tail the length of the chain's tail that should be used to calculate the likelihood of 
#' the map. If NULL (default), the function uses all markers positioned. Even if info.tail = TRUE, 
#' it uses at least extend.tail as the tail length
#' @param inflation.lim.p1 the maximum accepted length difference between the current and the previous 
#' parent 1 sub-map defined by arguments info.tail and extend.tail. If the size exceeds this limit, the marker will 
#' not be inserted. If NULL(default), then it will insert all markers.
#' @param inflation.lim.p2 same as `inflation.lim.p1` but for parent 2 sub-map.
#' @param phase.number.limit the maximum number of linkage phases of the sub-maps defined by arguments info.tail 
#' and extend.tail. Default is 20. If the size exceeds this limit, the marker will not be inserted. If Inf, 
#' then it will insert all markers.
#' @param tol the desired accuracy during the sequential phase of each parental map (default = 10e-02)
#' @param tol.final the desired accuracy for the final parental map (default = 10e-04)
#' @param verbose If TRUE (default), current progress is shown; if FALSE, no output is produced
#' @param method indicates whether to use 'hmm' (Hidden Markov Models), 'ols' (Ordinary Least Squares) 
#' to re-estimate the recombination fractions while merging the parental maps (default:hmm)
#' 
#' @return list containing three \code{mappoly.map} objects:1) map built with markers with segregation information from parent 1; 
#' 2) map built with markers with segregation information from parent 2; 3) maps in 1 and 2 merged
#' 
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu} with documentation and minor modifications by Cristiane Taniguti \email{chtaniguti@tamu.edu}
#' 
#' 
#' @export
framework_map <- function(input.seq, 
                          twopt, 
                          start.set = 10, 
                          thres.twopt = 10, 
                          thres.hmm = 30, 
                          extend.tail = 30,
                          inflation.lim.p1 = 5,
                          inflation.lim.p2 = 5,
                          phase.number.limit = 10,
                          tol = 10e-3,
                          tol.final = 10e-4,
                          verbose = TRUE,
                          method = "hmm"){
  
  if (!inherits(input.seq, "mappoly.sequence")) {
    stop(deparse(substitute(input.seq)), " is not an object of class 'mappoly.sequence'")
  }
  
  if (!inherits(twopt, "mappoly.twopt")) {
    stop(deparse(substitute(twopt)), " is not an object of class 'mappoly.twopt'")
  }
  ##### Map for P1 ####
  s.p1 <- make_seq_mappoly(input.seq, info.parent = "p1")
  if(length(s.p1$seq.num) > 0){
    tpt.p1 <- make_pairs_mappoly(twopt, s.p1)
    map.p1 <- est_rf_hmm_sequential(input.seq = s.p1,
                                    twopt = tpt.p1,
                                    start.set = start.set,
                                    thres.twopt = thres.twopt, 
                                    thres.hmm = thres.hmm,
                                    extend.tail = extend.tail,
                                    sub.map.size.diff.limit = inflation.lim.p1, 
                                    phase.number.limit = phase.number.limit,
                                    reestimate.single.ph.configuration = TRUE,
                                    tol = tol,
                                    tol.final = tol.final, 
                                    verbose = verbose)
  } else {
    warning("No linkage information for parent 1")
    map.p1 <- NULL
    map.p1.p2 <- NULL
  }
  ##### Map for P2 ####
  s.p2 <- make_seq_mappoly(input.seq, info.parent = "p2")
  if(length(s.p2$seq.num) > 0){
    tpt.p2 <- make_pairs_mappoly(twopt, s.p2)
    map.p2 <- est_rf_hmm_sequential(input.seq = s.p2,
                                    twopt = tpt.p2,
                                    start.set = start.set,
                                    thres.twopt = thres.twopt, 
                                    thres.hmm = thres.hmm,
                                    extend.tail = extend.tail,
                                    sub.map.size.diff.limit = inflation.lim.p2, 
                                    phase.number.limit = phase.number.limit,
                                    reestimate.single.ph.configuration = TRUE,
                                    tol = tol,
                                    tol.final = tol.final, 
                                    verbose = verbose)
  } else {
    warning("No linkage information for parent 2")
    map.p2 <- NULL
    map.p1.p2 <- NULL
  }
  #### Merging P1 and P2 ####
  
  if(length(s.p1$seq.num) > 0 & length(s.p2$seq.num) > 0){
    rf.mat <- rf_list_to_matrix(twopt, 
                                thresh.LOD.ph = thres.twopt,
                                shared.alleles = TRUE, 
                                verbose = verbose)
    
    map.p1.p2 <- merge_parental_maps(map.p1 = map.p1, 
                                     map.p2 = map.p2, 
                                     full.seq = input.seq, 
                                     full.mat = rf.mat, 
                                     method = method,
                                     verbose = verbose)
  } 
  
  init.map.list <- list(map.p1 = map.p1, 
                        map.p2 = map.p2,
                        map.p1.p2 = map.p1.p2)
  
  return(init.map.list)
}

#' Add markers that are informative in both parents using HMM approach and evaluating difference 
#' in LOD and gap size
#' 
#' @param input.map.list list containing three \code{mappoly.map} objects:1) map built with markers with segregation information from parent 1; 
#' 2) map built with markers with segregation information from parent 2; 3) maps in 1 and 2 merged
#' @param input.seq object of class \code{mappoly.sequence} containing all markers for specific group
#' @param twopt object of class \code{mappoly.twopt}
#' @param thres.twopt the LOD threshold used to determine if the linkage phases compared via two-point 
#' analysis should be considered for the search space reduction (default = 5)
#' @param init.LOD the LOD threshold used to determine if the marker will be included or not after hmm analysis  (default = 30)
#' @param verbose If TRUE (default), current progress is shown; if FALSE, no output is produced
#' @param method indicates whether to use 'hmm' (Hidden Markov Models), 'ols' (Ordinary Least Squares) or 'wMDS_to_1D_pc' 
#' (weighted MDS followed by fitting a one dimensional principal curve) to re-estimate the recombination fractions after adding markers
#' @param input.mds An object of class \code{mappoly.map}
#' @param max.rounds integer defining number of times to try to fit the remaining markers in the sequence 
#' @param gap.threshold threshold for gap size
#' @param size.rem.cluster threshold for number of markers that must contain in a segment after a gap is removed to keep this segment in the sequence 
#'  
#' @return object of class \code{mappoly.map2}
#'   
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu} with documentation and minor modifications by Cristiane Taniguti \email{chtaniguti@tamu.edu}
#' 
#' @export
update_framework_map <- function(input.map.list,
                                 input.seq, 
                                 twopt,
                                 thres.twopt = 10,
                                 init.LOD = 30,
                                 verbose = TRUE,
                                 method = "hmm",
                                 input.mds = NULL,
                                 max.rounds = 50,
                                 size.rem.cluster = 2,
                                 gap.threshold = 4)
{
  
  if (!all(sapply(input.map.list, function(x) inherits(x, "mappoly.map")))) {
    stop(deparse(substitute(twopt)), " is not an object of class 'mappoly.map'")
  }
  if (!inherits(input.seq, "mappoly.sequence")) {
    stop(deparse(substitute(input.seq)), " is not an object of class 'mappoly.sequence'")
  }
  
  #### Inserting remaining markers ####
  if(is.null(input.map.list$map.p1.p2)) {
    warning("Map only have linkage information for one of the parents.")
    cur.map <- input.map.list[[which(sapply(input.map.list, function(x) !is.null(x)))]]
  } else {
    cur.map <- input.map.list$map.p1.p2
  }
  cur.seq <- input.seq
  LOD <- init.LOD
  ll <- map.list <- vector("list", max.rounds)
  la <- numeric(max.rounds)
  i <- 1
  dat <- get(cur.map$info$data.name, pos = 1)
  mrk.to.include <- setdiff(cur.seq$seq.mrk.names, cur.map$info$mrk.names)
  
  rf.mat <- rf_list_to_matrix(twopt, 
                              thresh.LOD.ph = thres.twopt,
                              shared.alleles = TRUE, verbose = verbose)
  
  rem.mrk <- NULL
  while(length(mrk.to.include) > 0 & i <= length(la)){
    cur.gen <- calc_genoprob(cur.map, verbose = verbose)
    cur.res <- add_md_markers(input.map = cur.map,
                              mrk.to.include = mrk.to.include,
                              input.seq = cur.seq,
                              input.matrix = rf.mat, 
                              input.genoprob = cur.gen, 
                              input.data = dat, 
                              input.mds = input.mds,
                              thresh = LOD, 
                              method = "hmm", 
                              verbose = verbose)
    
    a <- rem_mrk_clusters(input.map = cur.res$map, 
                          size.rem.cluster = size.rem.cluster, 
                          gap.threshold = gap.threshold)
    
    rem.mrk <- c(rem.mrk, setdiff(cur.res$map$info$mrk.names, a$info$mrk.names))
    if(a$info$n.mrk < cur.res$map$info$n.mrk){
      if(method == "hmm"){
        cur.res$map <- reest_rf(a, method = "hmm", verbose = verbose)
      } else if(method == "wMDS_to_1D_pc"){
        cur.res$map <- reest_rf(a, method = "wMDS_to_1D_pc", input.mds = input.mds, verbose=verbose)
      } 
    }
    #### resulting map
    cur.map <- map.list[[i]] <- cur.res$map
    ### resulting log-likelihood vector
    ll[[i]] <- cur.res$ll[cur.res$map$info$mrk.names]
    if(i != 1){
      if(map.list[[i-1]]$info$n.mrk >= map.list[[i]]$info$n.mrk)
      {
        LOD <- max(cur.res$ll, na.rm = TRUE) - 10      
        if(LOD <= 0) LOD <- 10e-5
      }
    }
    la[i] <- LOD
    mrk.to.include <- setdiff(cur.seq$seq.mrk.names, cur.map$info$mrk.names)
    mrk.to.include <- setdiff(mrk.to.include, rem.mrk)
    plot_map_list(map.list[1:i])
    text(rep(0, i), 1:i+.5, labels = unlist(sapply(map.list, function(x) x$info$n.mrk)), adj=-1)
    i <- i + 1
    if(verbose){
      cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
      print(summary_maps(map.list[!sapply(map.list, is.null)], verbose = verbose))
      cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")    
      if(length(mrk.to.include) == 0) cat(paste("No more markers to add. Stopping in iteration:", i))
    }
  }
  #### Post-mapping ####
  id <- which(la!=0)
  structure(list(single =   list(map.p1 = input.map.list$map.p1, 
                                 map.p2 = input.map.list$map.p2,
                                 map.p1.p2 = input.map.list$map.p1.p2),
                 both = list(map.list = map.list[id],
                             lod.thresh = la[id],
                             calc.lod = ll[id])), class = "mappoly.map2")
}


#' Build merged parental maps
#' 
#' @param map.p1 object of class \code{mappoly.map} with parent 1 phased 
#' @param map.p2 object of class \code{mappoly.map} with parent 2 phased
#' @param full.seq object of class \code{mappoly.sequence} containing parent 1 and parent 2 markers
#' @param full.mat object of class \code{mappoly.rf.matrix} containing two-points recombination 
#' fraction estimations for parent 1 and parent 2 markers
#' 
#' @param method indicates whether to use 'hmm' (Hidden Markov Models), 'ols' (Ordinary Least Squares) 
#' to re-estimate the recombination fractions
#' 
#' @importFrom dplyr left_join `%>%`
#' 
#' @return object of class \code{mappoly.map} with both parents information
#' 
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu} with documentation and minor modifications by Cristiane Taniguti \email{chtaniguti@tamu.edu}
#' 
#' @keywords internal
#' 
merge_parental_maps <- function(map.p1, 
                                map.p2, 
                                full.seq, 
                                full.mat, 
                                method = c("ols", "hmm"),
                                verbose = TRUE){
  
  method <- match.arg(method)
  df.map <- data.frame(mrk = full.seq$seq.mrk.names,
                       pos = seq(0, 100, length.out = length(full.seq$seq.mrk.names)), 
                       row.names = full.seq$seq.mrk.names)
  mrk.p1 <- data.frame(mrk = map.p1$info$mrk.names, p1 = TRUE)
  mrk.p2 <- data.frame(mrk = map.p2$info$mrk.names, p2 = TRUE)
  mrk.md <- data.frame(mrk = setdiff(full.seq$seq.mrk.names, 
                                     c(mrk.p1$mrk, mrk.p2$mrk)), 
                       p1p2 = TRUE)
  df.map <- df.map %>% left_join(mrk.p1, by = "mrk") %>% 
    left_join(mrk.md, by = "mrk") %>% left_join(mrk.p2, by = "mrk")
  mrk.p1.p2 <- df.map$mrk[which(df.map$p1 | df.map$p2)]
  seq.num <- full.seq$seq.num[match(mrk.p1.p2, full.seq$seq.mrk.names)]
  seq.ph = list(P = c(map.p1$maps[[1]]$seq.ph$P, 
                      map.p2$maps[[1]]$seq.ph$P)[as.character(seq.num)],
                Q = c(map.p1$maps[[1]]$seq.ph$Q, 
                      map.p2$maps[[1]]$seq.ph$Q)[as.character(seq.num)])
  map.p1.p2 <- .mappoly_map_skeleton(ploidy = full.seq$ploidy, 
                                     n.mrk = length(mrk.p1.p2), 
                                     seq.num = seq.num, 
                                     mrk.names = mrk.p1.p2, 
                                     seq.dose.p1 = sapply(seq.ph$P, function(x) sum(as.logical(x))), 
                                     seq.dose.p2 = sapply(seq.ph$Q, function(x) sum(as.logical(x))),
                                     data.name = full.seq$data.name, 
                                     seq.rf = mf_h(diff(df.map$pos[match(mrk.p1.p2, df.map$mrk)])),
                                     seq.ph = seq.ph)
  if(method == "ols"){
    map.p1.p2 <- reest_rf(map.p1.p2, method = "ols", input.mat = full.mat, verbose = verbose)
  } else if(method == "hmm"){
    map.p1.p2 <- reest_rf(map.p1.p2, method = "hmm", verbose = verbose)
  } else {
    stop("Invalid method.", call. = FALSE)
  }
  map.p1.p2
}

#' Add markers to a pre-existing sequence using HMM analysis and evaluating difference in LOD
#' 
#' @param input.map An object of class \code{mappoly.map}
#' @param mrk.to.include vector for marker names to be included
#' @param input.seq an object of class \code{mappoly.sequence} containing all markers (the ones in the mappoly.map and also the ones to be included)
#' @param input.matrix object of class \code{mappoly.rf.matrix}
#' @param input.genoprob an object of class \code{mappoly.genoprob} obtained with calc_genoprob of the input.map object 
#' @param input.data an object of class \code{mappoly.data}
#' @param input.mds An object of class \code{mappoly.map}
#' @param thresh the LOD threshold used to determine if the marker will be included or not after hmm analysis  (default = 30)
#' @param extend.tail the length of the chain's tail that should be used to calculate the likelihood of 
#' the map. If NULL (default), the function uses all markers positioned. Even if info.tail = TRUE, 
#' it uses at least extend.tail as the tail length
#' @param method indicates whether to use 'hmm' (Hidden Markov Models), 'ols' (Ordinary Least Squares) or 'wMDS_to_1D_pc' 
#' (weighted MDS followed by fitting a one dimensional principal curve) to re-estimate the recombination fractions after adding markers
#' @param verbose If TRUE (default), current progress is shown; if FALSE, no output is produced
#' 
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu} with documentation and minor modifications by Cristiane Taniguti \email{chtaniguti@tamu.edu}
#' 
#' @keywords internal
#' 
add_md_markers <- function(input.map, 
                           mrk.to.include,
                           input.seq,
                           input.matrix, 
                           input.genoprob, 
                           input.data,
                           input.mds = NULL,
                           thresh = 500, 
                           extend.tail = 50,
                           method = c("hmm", "wMDS_to_1D_pc"), 
                           verbose = TRUE){
  method <- match.arg(method)
  if(method == "wMDS_to_1D_pc" & is.null(input.mds))
    stop("you must provide 'input.mds' when selecting method = 'wMDS_to_1D_pc'", call. = FALSE)
  p <- match(input.map$info$seq.num, input.seq$seq.num)
  if(any(diff(p) <= 0))
    stop("map and sequence have different orders", call. = FALSE)
  
  # Positioned markers  
  id.names <- intersect(input.map$info$mrk.names, input.seq$seq.mrk.names)
  # Markers to be positioned
  id2.names <-  setdiff(input.seq$seq.mrk.names, input.map$info$mrk.names)
  id2.names <-  intersect(id2.names, mrk.to.include)
  
  id <- match(id.names, input.seq$seq.mrk.names)  
  id2 <- match(id2.names, input.seq$seq.mrk.names)  
  
  Pl <- Ql <- vector("list", length(id2))
  nm <- character(length(id2))
  ll <- numeric(length(id2))
  names(ll) <- id2.names
  for(i in 1:length(id2)){
    u <- which(id - id2[i] < 0)
    
    if(verbose) cat(crayon::bgMagenta(stringr::str_pad(paste0(round(100*i/length(id2),1), "%"), width = 6)), "----> ")
    if(length(u) == 0 ){
      pos.test <- 0
      if(verbose) cat(crayon::bgRed(id2[i]), "--", id[pos.test + 1], "...|", sep = "")
    } else if(length(u) == length(id)){
      pos.test <- length(id)
      if(verbose) cat("|...",id[pos.test - 1], "--", crayon::bgRed(id2[i]), sep = "")
    } else{
      pos.test <- max(u)    
      if(verbose) cat("|...", id[pos.test], "--", crayon::bgRed(id2[i]), "--", id[pos.test + 1], "...|", sep = "")
    }
    
    temp.map <- add_marker(input.map = input.map, 
                           mrk = id2.names[i], 
                           pos = pos.test,
                           genoprob = input.genoprob, 
                           rf.matrix = input.matrix,
                           extend.tail = extend.tail,
                           verbose = FALSE)
    l <- get_LOD(temp.map)
    if(length(l) > 1){
      ll[i] <- l[2]
      if(verbose) cat("~~~~~",round(l[2], 4),"~~~~~")
      if(l[2] < thresh){
        if(verbose) cat(" ---> skip it!\n")
        next()
      } 
    } else ll[i] <- NA
    a1 <- temp.map$maps[[1]]$seq.ph$P[pos.test + 1]
    a2 <- temp.map$maps[[1]]$seq.ph$Q[pos.test + 1]
    nm[i] <- input.seq$seq.num[id2][i]
    Pl[[i]] <- a1[[1]]
    Ql[[i]] <- a2[[1]]
    if(verbose) cat("\n")
  }
  names(Pl) <- names(Ql) <- nm
  ix <- which(nm == "")
  if(length(ix) > 0){
    Pl <- Pl[-ix]
    Ql <- Ql[-ix]    
  }
  i1<-input.seq$seq.num
  i2<-c(as.numeric(names(Pl)), input.map$info$seq.num)
  seq.num <- intersect(i1, i2)
  seq.ph = list(P = c(input.map$maps[[1]]$seq.ph$P, Pl)[as.character(seq.num)],
                Q = c(input.map$maps[[1]]$seq.ph$Q, Ql)[as.character(seq.num)])
  id<-match(seq.num, input.seq$seq.num)
  map.p1.p2.final <- .mappoly_map_skeleton(ploidy = input.seq$ploidy, 
                                           n.mrk = length(seq.num), 
                                           seq.num = seq.num, 
                                           mrk.names =   input.seq$seq.mrk.names[id], 
                                           seq.dose.p1 = sapply(seq.ph$P, function(x) sum(as.logical(x))), 
                                           seq.dose.p2 = sapply(seq.ph$Q, function(x) sum(as.logical(x))), 
                                           data.name = input.seq$data.name, 
                                           seq.rf = rep(0.01, length(seq.num) - 1),
                                           seq.ph = seq.ph,
                                           chrom = input.seq$chrom[id], 
                                           genome.pos = input.seq$genome.pos[id])
  if(method == "hmm"){
    map.p1.p2.final <- reest_rf(map.p1.p2.final, method = "hmm", verbose = FALSE)
  } else if(method == "wMDS_to_1D_pc"){
    map.p1.p2.final <- reest_rf(map.p1.p2.final, method = "wMDS_to_1D_pc", input.mds = input.mds, verbose = FALSE)
  } 
  list(map = map.p1.p2.final, ll = ll)
}

# remove markers if causing gap > threshold
rem_mrk_clusters <- function(input.map, 
                             size.rem.cluster = 1, 
                             gap.threshold = 3){
  id <- which(imf_h(input.map$maps[[1]]$seq.rf) > gap.threshold)
  id <- cbind(c(1, id+1), c(id, input.map$info$n.mrk))
  ## Selecting map segments larger then the specified threshold
  segments <- id[apply(id, 1, diff) > size.rem.cluster - 1, , drop = FALSE]
  if(dim(segments)[1] == 0) stop("All markers were discarted using the defined gap.threshold.")
  id<-NULL
  for(i in 1:nrow(segments)){
    id <- c(id, segments[i,1]:segments[i,2])  
  }
  output.map <- get_submap(input.map = input.map, mrk.pos = id, 
                           reestimate.rf = FALSE, verbose = FALSE)
  return(output.map)
}


#' Plot object mappoly.map2
#' 
#' @param x object of class \code{mappoly.map2}
#' 
#' @import ggplot2
#' @importFrom ggpubr font
#' 
#' @export
plot_mappoly.map2 <- function(x){
  { 
    Var1 <- value <- Var2 <- position <- iteration <- hmm <- parent <- LODScore <- NULL
    df1 <- t(sapply(x$both$map.list, function(x) summary(diff(extract_map(x)))))
    nit <- nrow(df1)
    df11 <- reshape2::melt(df1)
    df11$hmm <- "partial"
    if(!is.na(names(x$both)[4])){
      df2 <- t(sapply(x$both$map.list.err, function(x) summary(diff(extract_map(x)))))   
      df12 <- reshape2::melt(df2)
      df12$hmm <- "full"
      df <- rbind(df11, df12)
      p1<-ggplot(df, aes(x=as.factor(Var1), y=value, group=Var2)) +
        geom_line(aes(color=Var2)) +
        geom_point(aes(color=Var2)) + 
        ylab("nearest neighbor marker distance") +
        xlab("Iteration") + 
        #annotate(geom="label", x = 1:nit, y = df1[,6] + .5, 
        #         label=sapply(x$both$map.list, function(x) x$info$n.mrk),
        #         color="red") + 
        theme(legend.position = "none")  + facet_wrap(.~hmm)
    } else{
      df <- df11
      p1<-ggplot(df, aes(x=as.factor(Var1), y=value, group=Var2)) +
        geom_line(aes(color=Var2)) +
        geom_point(aes(color=Var2)) + 
        ylab("nearest neighbor marker distance") +
        xlab("Iteration") + 
        annotate(geom="label", x = 1:nit, y = df1[,6] + .5, 
                 label=sapply(x$both$map.list, function(x) x$info$n.mrk),
                 color="red") + 
        theme(legend.position = "none")
    }
  }  
  {
    df31 <- lapply(x$both$map.list, function(x) extract_map(x))
    df41 <- reshape2::melt(df31)
    df41$hmm <- "partial"
    colnames(df41)[1:2] <- c("position", "iteration")
    len.map1 <- sapply(df31, max)  
    if(!is.na(names(x$both)[4])){
      df32 <- lapply(x$both$map.list.err, function(x) extract_map(x))    
      df42 <- reshape2::melt(df32)
      df42$hmm <- "full"
      colnames(df42)[1:2] <- c("position", "iteration")
      len.map2 <- sapply(df32, max)  
      df4 <- rbind(df41, df42)
      #len.map <- c(len.map1, len.map2)
    } else{
      df4 <- df41
    }
    p2<-ggplot(df4, aes(y=-position, x = as.factor(iteration), color = hmm)) +
      geom_point(shape = 95, size = 5) +
      # annotate(geom="label", y = -len.map + 1, x = 1:length(len.map)  + 0.2, 
      #          label=round(len.map, 1),
      #          color="red", size = 3) + 
      xlab("iteration") + facet_wrap(.~hmm) + theme(legend.position = "none")+
      scale_color_manual(values=c('#E69F00', '#56B4E9'))
  }
  
  single_maps <- x$single[which(!sapply(x$single, is.null))]
  df5 <- lapply(single_maps, function(x) extract_map(x))
  df6 <- reshape2::melt(df5)
  colnames(df6) <- c("position", "parent")
  p3<-ggplot(df6, aes(x=position, y=parent, group = as.factor(parent))) +
    geom_point(ggplot2::aes(color = as.factor(parent)), shape = 108, size = 5, show.legend = FALSE) +
    scale_color_brewer(palette="Dark2")
  
  
  df7 <- reshape2::melt(lapply(x$both$calc.lod, na.omit))
  colnames(df7) <- c("LODScore", "iteration")
  p4<-ggplot(df7, aes(x=as.factor(iteration), y=LODScore, group = as.factor(iteration))) +
    geom_boxplot(fill='#A4A4A4', color="black") + xlab("iteration")
  
  f1 <- ggarrange(p1, p3 + font("x.text", size = 9),
                  ncol = 2, nrow = 1, widths = c(2,1))
  f2 <- ggarrange(p2, p4 + font("x.text", size = 9),
                  ncol = 2, nrow = 1, widths = c(2,1))
  ggarrange(f1, f2, ncol = 1, nrow = 2)
}
