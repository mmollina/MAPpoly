# Install and restart current MAPpoly branch
rm(list = ls())
require(mappoly)

# So far, this dataset is private, however it will be released 
# as soon as https://doi.org/10.1101/689638 is accepted

# Loading ch12., I. trifida reference
dat.ch12 <- read_geno_dist(file.in = "~/repos/BT_map/src/genotype_calling/trifida/BT_trifida_12", prob.thres = 0.8)
match(dat.ch12$ind.names, unique(dat.ch12$geno$ind))
print(dat.ch12, detailed = TRUE)
plot(dat.ch12)
dat.ch12.1 <- filter_missing(input.data = dat.ch12, type = "marker", filter.thres = .4, inter = FALSE)
match(dat.ch12.1$ind.names, unique(dat.ch12.1$geno$ind))
dat.ch12.2 <- filter_missing(input.data = dat.ch12.1, type = "individual", filter.thres = .33, inter = FALSE)
match(dat.ch12.2$ind.names, unique(dat.ch12.2$geno$ind))

print(dat.ch12.2, detailed = TRUE)
plot(dat.ch12.2)
seq.chi.ch12.2 <- filter_segregation(input.data = dat.ch12.2, chisq.pval.thres = 5e-6, inter = FALSE)
seq.chi.ch12.2 <- make_seq_mappoly(dat.ch12.2, seq.chi.ch12.2$keep)
print(seq.chi.ch12.2, detailed = TRUE)
plot(seq.chi.ch12.2)
red.mrk <- elim_redundant(seq.chi.ch12.2)
plot(red.mrk)
seq.ch12.2 <- make_seq_mappoly(red.mrk)
counts <- cache_counts_twopt(seq.ch12.2, get.from.web = TRUE)
twopt.pairs <- est_pairwise_rf(input.seq = seq.ch12.2, 
                               n.clusters = 14,
                               count.cache = counts)
mat<-rf_list_to_matrix(twopt.pairs, n.clusters = 8)
id<-rownames(get_genomic_order(seq.ch12.2))
fields::image.plot(mat$rec.mat[id, id], col = rev(viridis::inferno(20)))
w<-rf_snp_filter(input.twopt = twopt.pairs, 
                 thresh.LOD.ph = 5, 
                 thresh.LOD.rf = 5, 
                 thresh.perc = 0.1)
plot(w)
mat1<-make_mat_mappoly(input.mat = mat, input.seq = w)
id<-rownames(get_genomic_order(w))
fields::image.plot(mat1$rec.mat[id, id], col = rev(viridis::inferno(20)))
x1 <- make_seq_mappoly(input.obj = dat.ch12.2, w$seq.mrk.names[1:400])
plot(x1$sequence.pos/1e6)
tpt1 <- make_pairs_mappoly(input.twopt = twopt.pairs, input.seq = x1)
subset.map <- est_rf_hmm_sequential(input.seq = x1,
                                    start.set = 5,
                                    thres.twopt = 10,
                                    thres.hmm = 10,
                                    extend.tail = 10,
                                    tol = 0.1,
                                    tol.final = 10e-3,
                                    twopt = tpt1,
                                    phase.number.limit = 20, 
                                    sub.map.size.diff.limit = 2,
                                    verbose = TRUE,
                                    high.prec = FALSE)
plot(subset.map)
M1<-est_full_hmm_with_global_error(input.map = subset.map, error = 0.05, tol = 10e-4, th.prob = 0.95, verbose = TRUE)  
#M1<-est_full_hmm_with_prior_dist(input.map = subset.map, tol = 10e-3, verbose = TRUE)  
plot(M1)
save.image("test.RData")

vec <- round(imf_h(M1$maps[[1]]$seq.rf),2)
res <- rep(1, length(vec) + 1)
ct <- 1
threshold <- 0 
for(i in 1:length(vec)){
  if(vec[i] > threshold){
    ct <- ct+1
    res[i+1] <- ct
  } else{
    res[i+1] <- ct
  }
}



id <- which(table(res) > 4)
mat.rec <- mat.lod <- matrix(NA, length(id), length(id))
for(i in 1:(length(id)-1)){
  for(j in (i+1):length(id)){
    cat("\nComputing rf between marker blocks", i , "and", j, "...")
    suppressMessages(map1 <- get_submap(M1, which(res==id[i]), reestimate.rf = FALSE))
    suppressMessages(map2 <- get_submap(M1, which(res==id[j]), reestimate.rf = FALSE))
    twopt.sub <- make_pairs_mappoly(tpt1, 
                                    make_seq_mappoly(dat.ch12.2, 
                                                     c(map1$maps[[1]]$seq.num, 
                                                       map2$maps[[1]]$seq.num)))
    M<-mat_share(map1,
                 map2,
                 twopt.sub,
                 count.cache = counts.web,
                 thres = 3)
    rf<-est_rf_marker_blocks(block1 = map1,
                             block2 = map2,
                             ph1 = "best",
                             ph2 = "best",
                             M = M,
                             max.inc = 0,
                             block1.tail = NULL,
                             tol = 0.01)
    mat.rec[i,j]<-rf$rf.stats[1,"rf"]
    mat.lod[i,j]<-rf$rf.stats[1,"rf_LOD"]
  }
}
mat.rec[mat.rec > .1]<-NA
fields::image.plot(imf_h(mat.rec), col = rev(viridis::inferno(20)))

#### Find Blocks
##
x2<-make_seq_mappoly(dat.ch12.2, id)
blocks_seq<-find_marker_blocks(input.seq = x2,
                               search.type = 'seq',
                               ask = TRUE,
                               seq.limit = 5000,
                               reconstruct = TRUE,
                               ph.thres = 5,
                               n.clusters = 12, 
                               count.cache = counts,
                               tol = 10e-2,
                               tol.final = 10e-3,
                               error = 0.2)
blocks_seq$bins<-blocks_seq$bins[!is.na(blocks_seq$bins)]
plot(blocks_seq)
blo<-filter_marker_blocks(blocks_seq, mf_h(.1))
plot(blo)
x<-blo
m <- x$bins[[1]]$info$m
n.bins <- length(x$bins)
n.mrk <- sapply(x$bins, function(x) x$info$n.mrk)
tot.n.mrk <- sum(n.mrk)
hap.len <- table(n.mrk)
single.bins <- lapply(x$bins, function(x) x$maps[[which.max(sapply(x$maps, function(x) x$loglike))]])
P.hap <- lapply(single.bins, function(x) ph_list_to_matrix(x$seq.ph$P, m))
Q.hap <- lapply(single.bins, function(x) ph_list_to_matrix(x$seq.ph$Q, m))
cod <- numeric(n.bins)
for (i in 1:n.bins) cod[i] <- length(unique(c(apply(P.hap[[i]], 2, paste, collapse = ""), apply(Q.hap[[i]], 2, paste, collapse = ""))))
bins.rf <- sapply(single.bins, function(x) x$seq.rf)
blo$bins<-blo$bins[cod > 2]
length(blo$bins)
plot(blo)
blo

plot(sapply(blo$bins, function(x) mean(dat.ch12.2$sequence.pos[x$maps[[1]]$seq.num])))

(id <- 1:length(blo$bins))
mat.rec <- mat.lod <- matrix(NA, length(id), length(id))
rf.list<-vector("list", choose(length(id), 2))
ct<-1
for(i in 1:(length(id)-1)){
  for(j in (i+1):length(id)){
    cat("\nComputing rf between marker blocks", i , "and", j, "...")
    map1<-blo$bins[[i]]
    #plot(map1)
    map2<-blo$bins[[j]]
    #plot(map2)
    twopt.sub <- make_pairs_mappoly(twopt.pairs, 
                                    make_seq_mappoly(dat.ch12.2, 
                                                     c(map1$maps[[1]]$seq.num, 
                                                       map2$maps[[1]]$seq.num)))
    M<-mat_share(map1,
                 map2,
                 twopt.sub,
                 count.cache = counts.web,
                 thres = 5)
    rf<-est_rf_marker_blocks(block1 = map1,
                             block2 = map2,
                             ph1 = "best",
                             ph2 = "best",
                             M = M,
                             max.inc = 0,
                             block1.tail = NULL,
                             tol = 0.01)
    rf.list[[ct]] <- rf
    ct<-ct+1
    if(all(is.na(rf$rf.stats[1]))) next()
    mat.rec[j,i]<-mat.rec[i,j]<-rf$rf.stats[1,"rf"]
    mat.lod[j,i]<-mat.lod[i,j]<-rf$rf.stats[1,"rf_LOD"]
    fields::image.plot(mat.rec, col = rev(viridis::inferno(20)))
  }
}
### 
load(file = "~/repos/MAPpoly/test.RData")
mat.rec <- mat.lod <- matrix(NA, length(id), length(id))
ct<-0
for(i in 1:(length(id)-1)){
  for(j in (i+1):length(id)){
    ct<-ct+1
    rf<-rf.list[[ct]]
    if(all(is.na(rf$rf.stats[1]))) next()
    mat.rec[j,i]<-mat.rec[i,j]<-rf$rf.stats[1,"rf"]
    mat.lod[j,i]<-mat.lod[i,j]<-abs(rf$rf.stats[2,"ph_LOD"])
  }
}
load(file = "~/repos/MAPpoly/test.RData")
mat.rec[mat.lod < 3]<-NA
mat.rec[mat.rec > .25]<-NA
fields::image.plot(imf_h(mat.rec), col = rev(viridis::inferno(20)))
#save.image(file = "~/repos/MAPpoly/test.RData")






bla<-generate_all_link_phases_elim_equivalent_haplo(map1$maps[[1]], map2$maps[[1]], M, m, max.inc)


new.map0<-get_submap(subset.map, match(c(map1$maps[[1]]$seq.num, map2$maps[[1]]$seq.num),  subset.map$maps[[1]]$seq.num))
new.map<-new.map0                    
new.map$maps[[1]]$seq.ph <- list(P = c(rf_map1_map2$phM1[[1]]$config.to.test[[1]]$P,
                                       rf_map1_map2$phM2[[1]]$config.to.test[[1]]$P),
                                 Q = c(rf_map1_map2$phM1[[1]]$config.to.test[[1]]$Q,
                                       rf_map1_map2$phM2[[1]]$config.to.test[[1]]$Q))
new.map$maps[[1]]$seq.rf[map1$info$n.mrk]<-rf_map1_map2$rf.stats[1, "rf"]
plot_compare_haplotypes(m = 6, 
                        hom.allele.p1 = new.map0$maps[[1]]$seq.ph$P, 
                        hom.allele.p2 = new.map$maps[[1]]$seq.ph$P, 
                        hom.allele.q1 = new.map0$maps[[1]]$seq.ph$Q,
                        hom.allele.q2 = new.map$maps[[1]]$seq.ph$Q)
