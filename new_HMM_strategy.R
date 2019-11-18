## Hexaploid map example
require(mappoly)
map<-maps.hexafake[[1]]
plot(map)
s <- make_seq_mappoly(hexafake, map$maps[[1]]$seq.num)
counts <- get_cache_two_pts_from_web(m = 6)
tpt<-est_pairwise_rf(input.seq = s,
                     count.cache = counts,
                     n.clusters = 16,
                     verbose=TRUE)
mat<-rf_list_to_matrix(tpt)
plot(mat)

blocks<-find_marker_blocks(input.seq = s,
                           search.type = "orig.ord",
                           rf.limit = 1e-04,
                           seq.limit = 1000,
                           ord.limit = 2,
                           reconstruct = TRUE,
                           extend.tail = 10,
                           n.clusters = 1,
                           ph.thres = 5,
                           rf.mat = mat,
                           tol = 1e-02,
                           tol.final = 1e-02,
                           error = NULL,
                           verbose = TRUE,
                           count.cache = counts,
                           ask = TRUE)

get_codom_level<-function(map){
  A<-cbind(ph_list_to_matrix(L = map$maps[[1]]$seq.ph$P, m = map$info$m),
           ph_list_to_matrix(L = map$maps[[1]]$seq.ph$Q, m = map$info$m))
  codom<-length(unique(apply(A, 2, paste, collapse  = "")))
  return(codom)
}


blocks$bins<-blocks$bins[!sapply(blocks$bins, function(x) all(is.na(x)))]
plot(blocks)
blocks.filt<-filter_marker_blocks(blocks, rf.thres = 0.002, ph.thresh = 3)
plot(blocks.filt)
i1<-as.numeric(sapply(blocks.filt$bins, function(x) x$maps[[1]]$seq.num))
plot(x = 1:length(s$seq.num), y = rep(10, length(s$seq.num)), 
     cex = 1, col = "darkorchid", pch = 1, ylim = c(0,10))
next.id.new<-c(1:length(s$seq.num))
next.id.new<-matrix(next.id.new, nrow = 2)
next.id.col<-is.na(match(s$seq.num, i1))
next.id.col<-matrix(next.id.col, nrow = 2)
for(i in which(!apply(next.id.col, 2, all))){
  lines(x = c(next.id.new[1,i], mean(next.id.new[,i])), y = c(10,9))  
  lines(x = c(next.id.new[2,i], mean(next.id.new[,i])), y = c(10,9))  
}
pos.pairs<-apply(next.id.new, 2, mean)
points(x = pos.pairs, 
       y = rep(9, ncol(next.id.new)), 
       cex = sapply(blocks.filt$bins, get_codom_level)/2, 
       col = c(4,2)[apply(next.id.col, 2, all)+1], pch = 19, 
       ylim = c(0,10))
## Function
my_map_func_test<-function(hap1, hap2, 
                           thresh.cut.path = 10e-3, 
                           thresh.twopt = 5, 
                           tol = 10e-4){
  ## get conditional probabilities for the first block (fixed)
  p.fix <- calc_genoprob(hap1, verbose = FALSE)
  ## Hash table: homolog combination --> states to visit in both parents
  A<-as.matrix(expand.grid(0:(sqrt(dim(p.fix$probs)[1])-1), 
                           0:(sqrt(dim(p.fix$probs)[1])-1))[,2:1])
  rownames(A) <- dimnames(p.fix$probs)[[1]]
  ## h: states to visit in both parents
  ## e: probability distribution 
  e.fix <- h.fix <- vector("list", dim(p.fix$probs)[3])
  for(i in 1:dim(p.fix$probs)[3]){
    a <- p.fix$probs[,dim(p.fix$probs)[2],i]  
    e.fix[[i]] <- a[a > thresh.cut.path]
    h.fix[[i]] <- A[names(e.fix[[i]]), , drop = FALSE]
  }
  ## Second block
  stemp<-make_seq_mappoly(hexafake, c(hap1$maps[[1]]$seq.num, hap2$maps[[1]]$seq.num))
  tpt.temp <- make_pairs_mappoly(input.twopt = tpt, stemp)
  M <- mat_share(x1 = hap1,
                 x2 = hap2,
                 twopt = tpt.temp,
                 count.cache = counts,
                 thres = thresh.twopt)
  w<-generate_all_link_phases_elim_equivalent_haplo(block1 = hap1$maps[[1]], 
                                                    block2 = hap2$maps[[1]],
                                                    M = M, 
                                                    m = hap1$info$m, 
                                                    max.inc = 0)
  ## get testmaps and conditional probabilities
  testmaps<-p<-vector("list", length(w))
  for(i in 1:length(w))
  {
    hap2$maps[[1]]$seq.ph <- w[[i]]
    testmaps[[i]] <- hap2
    p[[i]] <- calc_genoprob(testmaps[[i]], verbose = FALSE)
  }
  ## h: states to visit in both parents
  ## e: probability distrobution 
  h<-e<-vector("list", length(w))
  for(j in 1:length(w)){
    etemp<-htemp<-vector("list", dim(p.fix$probs)[3])
    for(i in 1:dim(p.fix$probs)[3]){
      a <- p[[j]]$probs[,1,i]  
      etemp[[i]] <- a[a > thresh.cut.path]
      htemp[[i]] <- A[names(etemp[[i]]), , drop = FALSE]
    }
    h[[j]] <- htemp
    e[[j]] <- etemp
  }
  configs<-vector("list", length(testmaps))
  names(configs) <- paste0("Phase_config.", 1:length(testmaps))
  res<-matrix(NA, nrow = length(h), ncol = 2, dimnames = list(names(configs), c("log_like", "rf")))
  ## HMM
  for(i in 1:length(h)){
    #cat("testing", i, "of", length(h), "\n")
    cat(".")
    h.test<-c(list(h.fix), h[i])
    e.test<-c(list(e.fix), e[i])
    restemp<-est_haplo_hmm(m = hap1$info$m, 
                           nmar = length(h.test), 
                           nind= dim(p.fix$probs)[3], 
                           haplo = h.test, 
                           emit = e.test, 
                           rf_vec = rep(0.01, length(h.test)-1), 
                           verbose = FALSE, 
                           use_H0 = FALSE, 
                           tol = tol) 
    res[i,]<-unlist(restemp)
  }
  cat("\n")
  for(i in 1:length(testmaps)){
    P<-c(hap1$maps[[1]]$seq.ph$P, testmaps[[i]]$maps[[1]]$seq.ph$P)
    Q<-c(hap1$maps[[1]]$seq.ph$Q, testmaps[[i]]$maps[[1]]$seq.ph$Q)
    names(P)<-names(Q)<-c(hap1$maps[[1]]$seq.num, testmaps[[i]]$maps[[1]]$seq.num)
    configs[[i]]<-list(P = P, Q = Q)
  }
  res<-res[order(res[,"log_like"], decreasing = TRUE),,drop = FALSE]
  return(list(stats = res, phase.config  = configs))  
}

#### across LG ####
z <- blocks.filt$bins
iter<-0
while(length(z) > 1){
  iter <- iter + 1
  block.list <-z
  block.at.end <- NULL
  flag<-FALSE
  if(length(block.list)%%2 == 1){
    block.list <- block.list[-length(block.list)]
    block.at.end <- block.list[[length(block.list)]]
    flag<-TRUE
  }
  is.same <- id <- matrix(1:length(block.list), nrow = 2)
  if(iter == 1){
    pos<-pos.pairs[!apply(next.id.col, 2, all)]
    if(flag){
      last.pos <- pos[length(pos)]
      pos <- pos[-length(pos)]      
    }
    pos<-matrix(pos, nrow = 2)
    new.pos <- numeric((length(block.list)/2))
  } else {
    pos <- new.pos
    if(flag){
      last.pos <- tail(pos, 1)
      pos <- pos[-length(pos)] 
    }
    pos<-matrix(pos, nrow = 2)
    new.pos <- numeric((length(block.list)/2))
  }
  z <- vector("list", (length(block.list)/2))
  proc.time <- numeric((length(block.list)/2))
  keep.hmm.based<-logical(length(block.list)/2)
  cat("Merging", length(block.list), "blocks\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  for(i in 1:(length(block.list)/2)){
    cat("Blocks",  id[1,i], "and", id[2,i], ": ")
    hap1 <- block.list[[id[1,i]]]
    hap2 <- block.list[[id[2,i]]]
    pt<-system.time(z.temp<-my_map_func_test(hap1 = hap1, hap2 = hap2, thresh.twopt = 5))
    id.ans<-match(c(hap1$maps[[1]]$seq.num, hap2$maps[[1]]$seq.num), map$maps[[1]]$seq.num)
    suppressMessages(a1<-get_submap(map, id.ans, reestimate.rf = FALSE))
    a2<-a1
    a2$maps[[1]]$seq.ph<-z.temp$phase.config[[rownames(z.temp$stats)[1]]]
    if(nrow(z.temp$stats) > 1){
      keep.hmm.based[i] <- abs(z.temp$stats[2,1]-z.temp$stats[1,1]) > 3    
    } else keep.hmm.based[i]<-TRUE
    proc.time[i]<-pt[3]/nrow(z.temp$stats)
    z[[i]] <- a2
    P.compare <- compare_haplotypes(m = 6, h1 = a1$maps[[1]]$seq.ph$P, h2 = a2$maps[[1]]$seq.ph$P)
    Q.compare <- compare_haplotypes(m = 6, h1 = a1$maps[[1]]$seq.ph$Q, h2 = a2$maps[[1]]$seq.ph$Q)
    is.same[,i]<-c(P.compare$is.same.haplo, Q.compare$is.same.haplo)
    lines(x = c(pos[1,i], mean(pos[,i])), y = c(10-iter,9-iter))  
    lines(x = c(pos[2,i], mean(pos[,i])), y = c(10-iter,9-iter))  
    points(x = mean(pos[,i]), 
           y = 9 - iter, 
           cex = get_codom_level(a2)/2, 
           col = c(2,4)[keep.hmm.based[i]+1], 
           pch = 19)
    new.pos[i]<-mean(pos[,i])
  }
  if(flag){
    cat("Block",  id[2,i]+1, "at the end of the chain...")
    hap1 <- z[[i]]
    hap2 <- block.at.end
    pt<-system.time(z.temp<-my_map_func_test(hap1 = hap1, hap2 = hap2, thresh.twopt = 5))
    id.ans<-match(c(hap1$maps[[1]]$seq.num, hap2$maps[[1]]$seq.num), map$maps[[1]]$seq.num)
    suppressMessages(a1<-get_submap(map, id.ans, reestimate.rf = FALSE))
    a2<-a1
    a2$maps[[1]]$seq.ph<-z.temp$phase.config[[rownames(z.temp$stats)[1]]]
    if(nrow(z.temp$stats) > 1){
      keep.hmm.based[i] <- abs(z.temp$stats[2,1]-z.temp$stats[1,1]) > 3    
    } else keep.hmm.based[i]<-TRUE
    proc.time[i]<-pt[3]/nrow(z.temp$stats)
    z[[i]] <- a2
    P.compare <- compare_haplotypes(m = 6, h1 = a1$maps[[1]]$seq.ph$P, h2 = a2$maps[[1]]$seq.ph$P)
    Q.compare <- compare_haplotypes(m = 6, h1 = a1$maps[[1]]$seq.ph$Q, h2 = a2$maps[[1]]$seq.ph$Q)
    is.same[,ncol(is.same)]<-c(P.compare$is.same.haplo, Q.compare$is.same.haplo)
    lines(x = c(last.pos, mean(pos[,i])), y = c(10-iter,9-iter))
    points(x = mean(pos[,i]), 
           y = 9 - iter, 
           cex = get_codom_level(a2)/2, 
           col = c(2,4)[keep.hmm.based[i]+1], 
           pch = 19)
    new.pos[i]<-mean(pos[,i])
  }
  if(length(z)!=1){
    (is.same<-is.same[,keep.hmm.based])
    z<-z[keep.hmm.based]
    new.pos<-new.pos[keep.hmm.based]    
  }
}
x11()

plot(z[[1]])
res<-reest_rf(z[[1]])
res
plot(res)

id.ans<-match(z[[1]]$maps[[1]]$seq.num, map$maps[[1]]$seq.num)
suppressMessages(a1<-get_submap(map, id.ans, reestimate.rf = FALSE))
(P.compare <- compare_haplotypes(m = 6, h1 = a1$maps[[1]]$seq.ph$P, h2 = z[[1]]$maps[[1]]$seq.ph$P))
(Q.compare <- compare_haplotypes(m = 6, h1 = a1$maps[[1]]$seq.ph$Q, h2 = z[[1]]$maps[[1]]$seq.ph$Q))
