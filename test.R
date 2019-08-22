data(hexafake)
mrk.subset<-make_seq_mappoly(hexafake, 1:200)
red.mrk<-elim_redundant(mrk.subset)
unique.mrks<-make_seq_mappoly(red.mrk)
counts.web<-cache_counts_twopt(unique.mrks, get.from.web = TRUE)
subset.pairs<-est_pairwise_rf(input.seq = unique.mrks,
                              count.cache = counts.web)
subset.map <- est_rf_hmm_sequential(input.seq = unique.mrks,
                                    thres.twopt = 5,
                                    thres.hmm = 10,
                                    extend.tail = 10,
                                    tol = 0.1,
                                    tol.final = 10e-3,
                                    twopt = subset.pairs,
                                    verbose = TRUE,
                                    high.prec = FALSE)

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
