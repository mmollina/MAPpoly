require(mappoly)
data("hexafake")
print(hexafake, detailed = TRUE)

##all markers
all.mrk<-make_seq_mappoly(input.obj = sweetpotato, arg = 'all')
print(all.mrk)
filt.mrk<-elim_redundant(all.mrk)
plot(filt.mrk)
print(filt.mrk)
counts<-cache_counts_twopt(input.seq = all.mrk,
                           get.from.web = TRUE)

all.pairs<-est_pairwise_rf(input.seq = all.mrk, count.cache = counts,
                           n.clusters = 16, batch.size = 1000000)

all.mat <- rf_list_to_matrix(all.pairs, n.clusters = 16)
plot(all.mat)
link.lgs<-group_mappoly(input.mat = all.mat,
                        input.seq = filt.lgs.1.2.3.seq,
                        expected.groups = 3,
                        comp.mat = TRUE)
link.lgs
lg1<-make_seq_mappoly(link.lgs, arg = 1)
m.lg1<-make_mat_mappoly(input.mat = all.mat, input.seq = lg1)
lg1.pairs<-make_pairs_mappoly(all.pairs, lg1)
lg1.filt<-rf_snp_filter(lg1.pairs, 7, 7, 0.15, thresh.perc = 0.1)
m.lg1.filt<-make_mat_mappoly(input.mat = all.mat, input.seq = lg1.filt)
plot(m.lg1.filt)
ord.lg1.mds<-mds_mappoly(m.lg1, n = c(10, 601, 600, 496,487))
plot(ord.lg1.mds)
seq.mds.lg1<-make_seq_mappoly(input.obj = ord.lg1.mds)

phases<-ls_linkage_phases(input.seq = seq.mds.lg1, thres = 3, twopt = all.pairs)



map.lg1.hmm<-est_rf_hmm_sequential(input.seq = seq.mds.lg1,
                                   thres.twopt = 10,
                                   thres.hmm = 10,
                                   info.tail = TRUE,
                                   extend.tail = 100,
                                   twopt = all.pairs,
                                   verbose = TRUE,
                                   reestimate.single.ph.configuration = FALSE,
                                   tol = 0.1,
                                   tol.final = 10e-2,
                                   phase.number.limit = 20)
map.lg1.hmm
dev.off()
plot(map.lg1.hmm)

x<-names(map.lg1.hmm$maps[[1]]$seq.ph$P)
phP<-ph_matrix_to_list(read.csv2(file = "inst/doc/phase_sim_hexa_P.csv")[,-1])
phQ<-ph_matrix_to_list(read.csv2(file = "inst/doc/phase_sim_hexa_Q.csv")[,-1])

phP.map<-map.lg1.hmm$maps[[1]]$seq.ph$P
compare_haplotypes(6, h1 = phP.map, h2 = phP[names(phP.map)])
phQ.map<-map.lg1.hmm$maps[[1]]$seq.ph$Q
compare_haplotypes(6, h1 = phQ.map, h2 = phQ[names(phQ.map)])

map.hmm<-reest_map(input.map = map.lg1.hmm, input.mat = all.mat,
                   tol = 10e-3,  phase.config = "best",
                    method = "hmm", weight = TRUE, verbose = TRUE)

map.ols<-reest_map(input.map = map.lg1.hmm, input.mat = all.mat,
                   tol = 10e-3,  phase.config = "best",
                   method = "ols", weight = TRUE, verbose = TRUE)

genoprob<-calc_genoprob(input.map = map.ols,  phase.config = "best", verbose = TRUE)

