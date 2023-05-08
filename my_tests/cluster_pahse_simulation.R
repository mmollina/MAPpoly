#### Functions ####
require(mappoly)
require(tidyverse)
rm(list = ls())
dev.off()
#### Simulation ####
ploidy <- 6
n.mrk <- 500
n.ind <- 200
map.length <- 100
seed <- 5887
h.temp <- sim_homologous(ploidy = ploidy, 
                         n.mrk = n.mrk, 
                         max.d = ploidy, 
                         max.ph = ploidy, 
                         seed = seed)
dat <- poly_cross_simulate(ploidy = ploidy, 
                           rf.vec = mf_h(diff(seq(0, map.length, length.out = n.mrk))), 
                           n.mrk = n.mrk,
                           n.ind = n.ind, 
                           h.temp, 
                           seed = seed)
h1 <- h.temp$hom.allele.p
h2 <- h.temp$hom.allele.q
names(h1) <- names(h2) <- seq_along(dat$mrk.names)

s <- make_seq_mappoly(dat, "all")
twopt <- est_pairwise_rf(s, ncpus = 24)
m <- rf_list_to_matrix(twopt)
#o <- mds_mappoly(m)
#full.seq <- make_seq_mappoly(o)
full.seq <- make_seq_mappoly(get_genomic_order(s))

##### Map for P1 ####
s.p1 <- make_seq_mappoly(s, info.parent = "p1")
full.seq.p1 <- s.p1$seq.mrk.names[s.p1$seq.mrk.names%in%full.seq$seq.mrk.names]
full.seq.p1 <- make_seq_mappoly(dat, full.seq.p1)
plot(m, ord = full.seq.p1, main.text = "Chromofull.seqme 1")
map.p1 <- est_rf_hmm_sequential(input.seq = full.seq.p1,
                                twopt = twopt,
                                start.set = 5,
                                thres.twopt = 10, 
                                thres.hmm = 10,
                                extend.tail = 10,
                                info.tail = TRUE, 
                                sub.map.size.diff.limit = 5, 
                                phase.number.limit = 20,
                                reestimate.single.ph.configuration = TRUE,
                                tol = 10e-2,
                                tol.final = 10e-4)
plot(map.p1)
h.p1 <- map.p1$maps[[1]]$seq.ph$P
compare_haplotypes(ploidy, h1[names(h.p1)], h.p1)
g1 <- calc_genoprob_single_parent(map.p1,
                                  step = 0,
                                  info.parent = 1,
                                  uninfo.parent = 2,
                                  global.err = 0.0,
                                  phase.config = "best", 
                                  verbose = TRUE)
image(t(g1$probs[,,3]))
##### Map for P2 ####
s.p2 <- make_seq_mappoly(s, info.parent = "p2")
full.seq.p2 <- s.p2$seq.mrk.names[s.p2$seq.mrk.names%in%full.seq$seq.mrk.names]
full.seq.p2 <- make_seq_mappoly(dat, full.seq.p2)
plot(m, ord = full.seq.p2)
map.p2 <- est_rf_hmm_sequential(input.seq = full.seq.p2,
                                twopt = twopt,
                                start.set = 5,
                                thres.twopt = 10, 
                                thres.hmm = 10,
                                extend.tail = 10,
                                info.tail = TRUE, 
                                sub.map.size.diff.limit = 5, 
                                phase.number.limit = 20,
                                reestimate.single.ph.configuration = TRUE,
                                tol = 10e-2,
                                tol.final = 10e-4)
plot(map.p2, mrk.names = TRUE)
h.p2 <- map.p2$maps[[1]]$seq.ph$Q
compare_haplotypes(ploidy, h2[names(h.p2)], h.p2)
g2 <- calc_genoprob_single_parent(map.p2,
                                  step = 0,
                                  info.parent = 2,
                                  uninfo.parent = 1,
                                  global.err = 0.0,
                                  phase.config = "best", 
                                  verbose = TRUE)
image(t(g2$probs[,,10]))

##### Merging P1 and P2 ####
map.p1
map.p2
full.seq 
full.mat <- rf_list_to_matrix(twopt, 1, 1)
map.p1.p2 <- merge_parental_maps(map.p1, 
                                 map.p2, 
                                 full.seq, full.mat, 
                                 method = "ols")
plot(map.p1.p2)
#####################round1##############################
map.1 <- map.p1.p2
gen.1 <- calc_genoprob(map.p1.p2)
mat.1 <- rf_list_to_matrix(twopt, 5, 5, shared.alleles = TRUE)
mrk.to.include <- setdiff(full.seq$seq.mrk.names, map.1$info$mrk.names)
input.seq <- full.seq
x.1 <- add_md_markers(input.map = map.1, 
                      input.seq = full.seq,
                      mrk.to.include = mrk.to.include,
                      input.matrix = mat.1, 
                      input.genoprob = gen.1, 
                      input.data = dat, 
                      thresh = 100, method = "ols")

plot(x.1, mrk.names = T,cex = .7)
plot_genome_vs_map(x.1)
h.p1 <- x.1$maps[[1]]$seq.ph$P
compare_haplotypes(ploidy, h1[names(h.p1)], h.p1)
h.p2 <- x.1$maps[[1]]$seq.ph$Q
compare_haplotypes(ploidy, h2[names(h.p2)], h.p2)
plot_compare_haplotypes(ploidy = x.1$info$ploidy, 
                        hom.allele.p1 = h1[names(h.p1)], 
                        hom.allele.q1 = h2[names(h.p2)], 
                        hom.allele.p2 = h.p1, 
                        hom.allele.q2 = h.p2)

#####################round2##############################
map.2 <- x.1
gen.2 <- calc_genoprob(map.2)
mat.2 <- rf_list_to_matrix(twopt, 5, 5, shared.alleles = TRUE)
mrk.to.include <- setdiff(full.seq$seq.mrk.names, map.2$info$mrk.names)
input.seq <- full.seq
x.2 <- add_md_markers(input.map = map.2, 
                      input.seq = full.seq,
                      mrk.to.include = mrk.to.include,
                      input.matrix = mat.2, 
                      input.genoprob = gen.2, 
                      input.data = dat, 
                      thresh = 50, method = "hmm")
x.2
plot(x.2, mrk.names = T,cex = .7)
plot_genome_vs_map(x.2)
h.p1 <- x.2$maps[[1]]$seq.ph$P
compare_haplotypes(ploidy, h1[names(h.p1)], h.p1)
h.p2 <- x.2$maps[[1]]$seq.ph$Q
compare_haplotypes(ploidy, h2[names(h.p2)], h.p2)
plot_compare_haplotypes(ploidy = x.2$info$ploidy, 
                        hom.allele.p1 = h1[names(h.p1)], 
                        hom.allele.q1 = h2[names(h.p2)], 
                        hom.allele.p2 = h.p1, 
                        hom.allele.q2 = h.p2)


#####################round3##############################
map.3 <- x.2
gen.3 <- calc_genoprob(map.3)
mat.3 <- rf_list_to_matrix(twopt, 5, 5, shared.alleles = TRUE)
mrk.to.include <- setdiff(full.seq$seq.mrk.names, map.3$info$mrk.names)
input.seq <- full.seq
x.3 <- add_md_markers(input.map = map.3, 
                      input.seq = full.seq,
                      mrk.to.include = mrk.to.include,
                      input.matrix = mat.3, 
                      input.genoprob = gen.3, 
                      input.data = dat, 
                      thresh = 10, method = "hmm")

plot(x.3, mrk.names = T,cex = .7)
plot_genome_vs_map(x.3)
h.p1 <- x.3$maps[[1]]$seq.ph$P
compare_haplotypes(ploidy, h1[names(h.p1)], h.p1)
h.p2 <- x.3$maps[[1]]$seq.ph$Q
compare_haplotypes(ploidy, h2[names(h.p2)], h.p2)
plot_compare_haplotypes(ploidy = x.3$info$ploidy, 
                        hom.allele.p1 = h1[names(h.p1)], 
                        hom.allele.q1 = h2[names(h.p2)], 
                        hom.allele.p2 = h.p1, 
                        hom.allele.q2 = h.p2)
