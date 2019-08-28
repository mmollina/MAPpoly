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
plot(mat, ord = rownames(get_genomic_order(seq.ch12.2)))
w<-rf_snp_filter(input.twopt = twopt.pairs, 
                 thresh.LOD.ph = 5, 
                 thresh.LOD.rf = 5, 
                 thresh.perc = 0.1)
plot(w)
mat1<-make_mat_mappoly(input.mat = mat, input.seq = w)
plot(mat1, ord = rownames(get_genomic_order(w)))

x1 <- make_seq_mappoly(input.obj = dat.ch12., w$seq.mrk.names[1:10])
tpt1 <- make_pairs_mappoly(input.twopt = twopt.pairs, input.seq = x1)


subset.map <- est_rf_hmm_sequential(input.seq = x1,
                                    start.set = 3,
                                    thres.twopt = 5,
                                    thres.hmm = 10,
                                    extend.tail = 10,
                                    tol = 0.1,
                                    tol.final = 10e-3,
                                    twopt = tpt1,
                                    verbose = TRUE,
                                    high.prec = FALSE)
bla<-est_full_hmm_with_global_error(input.map = subset.map, error = 0.1, tol = 10e-3, th.prob = 0.95)  

s1 <- make_seq_mappoly(hexafake, subset.map$maps[[1]]$seq.num[1:50])                                     
map1 <- get_submap(subset.map, 1:50)

s2 <- make_seq_mappoly(hexafake, subset.map$maps[[1]]$seq.num[51:100])                                     
map2 <- get_submap(subset.map, 51:100)

twopt.sub <- make_pairs_mappoly(subset.pairs, 
                                make_seq_mappoly(hexafake, 
                                                 c(map1$maps[[1]]$seq.num, 
                                                   map2$maps[[1]]$seq.num)))           
M<-mat_share(map1,
             map2,
             twopt.sub,
             count.cache = counts.web,
             thres = 3)

system.time(
bla1<-est_rf_marker_blocks(block1 = map1,
                           block2 = map2,
                           ph1 = "best",
                           ph2 = "best",
                           M = M,
                           max.inc = 0,
                           block1.tail = NULL,
                           tol = 0.01))
system.time(
bla2<-est_rf_marker_blocks(block1 = map1,
                           block2 = map2,
                           ph1 = "best",
                           ph2 = "best",
                           M = M,
                           max.inc = 0,
                           block1.tail = NULL,
                           tol = 0.01))
identical(bla1, bla2)

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
