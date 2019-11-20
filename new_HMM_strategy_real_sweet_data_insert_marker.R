## Hexaploid map example
require(mappoly)

#load Ch12 info 
load("~/repos/MAPpoly/objects_to_test_in_ch12.RData")
source("~/repos/MAPpoly/util_func.R")
map<-maps.final[[12]]
#plot(map)
s <- make_seq_mappoly(BT.trifida.triloba.298, map$maps[[1]]$seq.num[1:200])
counts <- get_cache_two_pts_from_web(m = s$m)
tpt<-est_pairwise_rf(input.seq = s,
                     count.cache = counts,
                     n.clusters = 16,
                     verbose=TRUE)

map1<-get_submap(map, mrk.pos = 1:200, reestimate.rf = FALSE)
window <- 20


for(i in 1:)
i<-1

length(map1$maps[[1]]$seq.num)




length(map$maps[[1]]$seq.num)
hap<-get_submap(map, 1:10)

hap2<-get_submap(map, 16:25)
full.hap<-my_map_func_test(hap1, hap2, thresh.cut.path = 0.001, thresh.twopt = 5)
full.hap$stats
mrk <- map$maps[[1]]$seq.num[11:15]
thresh.cut.path = 10e-3 
thresh.twopt = 5 
tol = 10e-4


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
stemp<-make_seq_mappoly(get(hap1$info$data.name), c(hap1$maps[[1]]$seq.num, mrk.test), data.name = hap1$info$data.name)
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


