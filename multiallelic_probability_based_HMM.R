## This script tests the multiallelic probability based HMM
## Loading MAPpoly
require(mappoly)

## Hexaploid map example
map<-maps.hexafake[[1]]
plot(map)

## Two-point distance matrix based on the HMM-based pre-built map
pos<-cumsum(c(0, imf_h(map$maps[[1]]$seq.rf)))
M<-matrix(0, map$info$n.mrk, map$info$n.mrk)
for(i in 1:map$info$n.mrk){
  for(j in 1:map$info$n.mrk){
    M[i,j] <- abs(diff(pos[c(j,i)]))
  }
}
image(M)

## Cluster markers in groups with < 0.1 cM
a <- hclust(as.dist(M), method = "complete")
plot(a)
ct <- cutree(a, h = 0.1)

## get submaps and conditional probabilities
submaps<-p<-NULL
cte<-1
for(i in unique(ct))
{
  print(which(ct == i))
  if(sum(ct == i) > 2){
    submaps[[cte]] <- get_submap(map, which(ct == i), reestimate.rf = FALSE)
    p[[cte]] <- calc_genoprob(submaps[[cte]], verbose = FALSE)    
    cte <- cte + 1
  }
}

id.mark<-match(unlist(sapply(submaps, function(x) x$maps[[1]]$seq.num)), map$maps[[1]]$seq.num)
## Map of the selected markers
map.subset<-get_submap(map, id.mark, reestimate.rf = FALSE)
plot(map.subset)

## Hash table: homolog combination --> states to visit in both parents
A<-as.matrix(expand.grid(0:19, 0:19)[,2:1])
rownames(A) <- dimnames(p[[1]]$probs)[[1]]
thresh <- 0.001
n.mrk<-length(p)

## h: states to visit in both parents
## e: probability distrobution 
h<-e<-vector("list", n.mrk)
for(j in 1:n.mrk){
  etemp<-htemp<-vector("list", 300)
  for(i in 1:300){
    a<-p[[j]]$probs[,ceiling(length(p[[j]]$map)/2),i]  
    etemp[[i]]<-a[a>thresh]
    htemp[[i]]<-A[names(etemp[[i]]), , drop = FALSE]
  }
  h[[j]] <- htemp
  e[[j]] <- etemp
}

## HMM
res<-est_haplo_hmm(m = 6, 
                   nmar = length(h), 
                   nind= 300, 
                   haplo = h, 
                   emit = e, 
                   rf_vec = rep(0.01, length(h)-1), 
                   verbose = TRUE, 
                   use_H0 = FALSE, 
                   tol = 10e-4) 

plot(cumsum(imf_h(res[[2]])))

rest<-which(is.na(match(map$maps[[1]]$seq.num, id.mark)))

submap.test<-get_submap(map, 34:35)
print(submap.test, detailed = TRUE)
p.test<-calc_genoprob(submap.test)
e.test<-h.test<-vector("list", 300)
for(i in 1:300){
    a<-p.test$probs[,ceiling(length(p.test$map)/2),i]  
    e.test[[i]]<-a[a>thresh]
    h.test[[i]]<-A[names(e.test[[i]]), , drop = FALSE]
}



e.in<-h.in<-NULL
##
h.in[[1]]<-c(list(h.test), h[1:74])
e.in[[1]]<-c(list(e.test), e[1:74])
##
h.in[[2]]<-c(h[1], list(h.test), h[2:74])
e.in[[2]]<-c(e[1], list(e.test), e[2:74])
##
h.in[[3]]<-c(h[1:2], list(h.test), h[3:74])
e.in[[3]]<-c(e[1:2], list(e.test), e[3:74])
##
h.in[[4]]<-c(h[1:3], list(h.test), h[4:74])
e.in[[4]]<-c(e[1:3], list(e.test), e[4:74])
##
h.in[[5]]<-c(h[1:4], list(h.test), h[5:74])
e.in[[5]]<-c(e[1:4], list(e.test), e[5:74])
##
h.in[[6]]<-c(h[1:5], list(h.test), h[6:74])
e.in[[6]]<-c(e[1:5], list(e.test), e[6:74])
##
h.in[[7]]<-c(h[1:6], list(h.test), h[7:74])
e.in[[7]]<-c(e[1:6], list(e.test), e[7:74])
##
h.in[[8]]<-c(h[1:7], list(h.test), h[8:74])
e.in[[8]]<-c(e[1:7], list(e.test), e[8:74])
##
h.in[[9]]<-c(h[1:8], list(h.test), h[9:74])
e.in[[9]]<-c(e[1:8], list(e.test), e[9:74])
##
h.in[[10]]<-c(h[1:9], list(h.test), h[10:74])
e.in[[10]]<-c(e[1:9], list(e.test), e[10:74])

save.image(file = "test.RData")

res<-NULL
for(i in 1:length(h.in)){
  res[[i]]<-est_haplo_hmm(m = 6, 
                      nmar = length(h.in[[i]]), 
                      nind= 300, 
                      haplo = h.in[[i]], 
                      emit = e.in[[i]], 
                      rf_vec = rep(0.01, length(h.in[[i]])-1), 
                      verbose = TRUE, 
                      use_H0 = FALSE, 
                      tol = 10e-3) 
}

plot(cumsum(c(0, imf_h(res[[1]][[2]]))))
text(x = 60, y = 0, labels = round(res[[1]][[1]]))
for(i in 2:10)
{
  points(cumsum(c(0, imf_h(res[[i]][[2]]))), col = i)
  text(x = 60, y = seq(10, 80, length.out = 9)[i], labels = round(res[[i]][[1]]), col = i)
}

which.max(sapply(res, function(x) x[[1]]))
sapply(submaps, function(x) x$maps[[1]]$seq.num)[1:8]
submap.test$maps[[1]]$seq.num


