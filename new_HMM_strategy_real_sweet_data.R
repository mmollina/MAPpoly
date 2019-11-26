## Hexaploid map example
require(mappoly)

#load Ch12 info 
#load("~/repos/MAPpoly/objects_to_test_in_ch12.RData")

map<-maps.final[[12]]
#plot(map)
s <- make_seq_mappoly(BT.trifida.triloba.298, map$maps[[1]]$seq.num[1:100])
counts <- get_cache_two_pts_from_web(m = s$m)
tpt<-est_pairwise_rf(input.seq = s,
                     count.cache = counts,
                     n.clusters = 16,
                     verbose=TRUE)
mat<-rf_list_to_matrix(tpt)
plot(mat)
init.block.size <- 3

system.time(blocks<-find_marker_blocks(input.seq = s,
                           search.type = "orig.ord",
                           rf.limit = 1e-04,
                           seq.limit = 50000,
                           ord.limit = init.block.size,
                           reconstruct = TRUE,
                           extend.tail = 10,
                           n.clusters = 16,
                           ph.thres = 5,
                           rf.mat = mat,
                           tol = 1e-02,
                           tol.final = 1e-02,
                           error = NULL,
                           verbose = TRUE,
                           count.cache = counts,
                           ask = TRUE))


blocks$bins<-blocks$bins[!sapply(blocks$bins, function(x) all(is.na(x)))]
plot(blocks)
blocks.filt<-filter_marker_blocks(blocks, rf.thres = 0.002, ph.thresh = 3)
plot(blocks.filt)
i1<-as.numeric(sapply(blocks.filt$bins, function(x) x$maps[[1]]$seq.num))

plot(x = 1:length(s$seq.num), y = rep(10, length(s$seq.num)), 
     cex = 1, col = "darkorchid", pch = 1, ylim = c(0,10))

length(blocks.filt$bins)

next.id.new<-c(1:length(s$seq.num))



next.id.new<-matrix(next.id.new, nrow = init.block.size)



next.id.new<-next.id.new[,-ncol(next.id.new)]
next.id.col<-is.na(match(s$seq.num, i1))
next.id.col<-matrix(next.id.col, nrow = init.block.size)
next.id.col<-next.id.col[,-ncol(next.id.col)]
dim(next.id.col)
dim(next.id.new)

for(i in which(!apply(next.id.col, 2, all))){
  lines(x = c(next.id.new[1,i], mean(next.id.new[,i])), y = c(10,9))  
  lines(x = c(next.id.new[2,i], mean(next.id.new[,i])), y = c(10,9))  
  lines(x = c(next.id.new[3,i], mean(next.id.new[,i])), y = c(10,9))  
}
pos.pairs<-apply(next.id.new, 2, mean)
points(x = pos.pairs, 
       y = rep(9, ncol(next.id.new)), 
       cex = sapply(blocks.filt$bins, get_codom_level)/2, 
       col = c(4,2)[apply(next.id.col, 2, all)+1], pch = 19, 
       ylim = c(0,10))


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
    ## IMPLEMENT: if both blocks are biallelic, use two-point again....
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


plot(z[[1]])
res<-reest_rf(z[[1]])
res
plot(res)

id.ans<-match(z[[1]]$maps[[1]]$seq.num, map$maps[[1]]$seq.num)
suppressMessages(a1<-get_submap(map, id.ans, reestimate.rf = FALSE))
(P.compare <- compare_haplotypes(m = 6, h1 = a1$maps[[1]]$seq.ph$P, h2 = z[[1]]$maps[[1]]$seq.ph$P))
(Q.compare <- compare_haplotypes(m = 6, h1 = a1$maps[[1]]$seq.ph$Q, h2 = z[[1]]$maps[[1]]$seq.ph$Q))
plot_compare_haplotypes(6, 
                        a1$maps[[1]]$seq.ph$P, 
                        a1$maps[[1]]$seq.ph$Q,
                        z[[1]]$maps[[1]]$seq.ph$P,
                        z[[1]]$maps[[1]]$seq.ph$Q)







