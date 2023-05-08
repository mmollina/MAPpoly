#### Functions ####
require(mappoly)
source("my_tests/temp_utils.R")
rm(list = ls())
source("~/repos/official_repos/MAPpoly/my_tests/temp_utils.R")
# require(tidyverse)
# BT_trifida <- readRDS("~/BT_trifida.rds")
# a <- str_split_fixed(BT_trifida$mrk.names, "_|S", 3)
# dat.in <- data.frame(mrk.names = BT_trifida$mrk.names,
#                      dose.p1 = BT_trifida$dosage.p,
#                      dose.p2 = BT_trifida$dosage.q,
#                      chrom = as.numeric(a[,2]),
#                      geno.pos = as.numeric(a[,3]),
#                      BT_trifida$geno.dose)
# BT_trifida <- table_to_mappoly(dat.in, ploidy = 6)
# saveRDS(BT_trifida, "~/repos/BT_trifida_dose.rds")
##### Initial computations ######
BT_trifida_dose <- readRDS("~/repos/BT_trifida_dose.rds")
BT_trifida_dose
BT_trifida_dose <- filter_individuals(BT_trifida_dose)

ploidy = 6
filt.mrk <- filter_segregation(BT_trifida_dose, 
                               chisq.pval.thres = 0.05/28651, 
                               inter = FALSE)
full.seq <- make_seq_mappoly(filt.mrk)
plot(full.seq)
full.seq.ch1 <- make_seq_mappoly(full.seq, names(which(full.seq$chrom == 1)))
plot(full.seq.ch1)
full.tpt.ch1 <- est_pairwise_rf(full.seq.ch1, ncpus = 24)
full.mat.ch1 <- rf_list_to_matrix(full.tpt.ch1, shared.alleles = TRUE)
filt.full.seq.ch1 <- rf_snp_filter(full.tpt.ch1, probs = c(0.025,0.975))
plot(filt.full.seq.ch1)
filt.full.mat.ch1 <- make_mat_mappoly(full.mat.ch1, filt.full.seq.ch1)
mds.ord.ch1 <- mds_mappoly(filt.full.mat.ch1)
full.mds.seq.ch1 <- make_seq_mappoly(mds.ord.ch1)
plot(filt.full.mat.ch1, ord = full.mds.seq.ch1, fact = 10)

##### Map for P1 ####
s.p1 <- make_seq_mappoly(full.mds.seq.ch1, info.parent = "p1")
plot(s.p1)
tpt.p1 <- make_pairs_mappoly(full.tpt.ch1, s.p1)
plot(filt.full.mat.ch1, ord = s.p1)
map.p1 <- est_rf_hmm_sequential(input.seq = s.p1,
                                twopt = tpt.p1,
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
g1 <- calc_genoprob_single_parent(map.p1,
                                  step = 0,
                                  info.parent = 1,
                                  uninfo.parent = 2,
                                  global.err = 0.0,
                                  phase.config = "best", 
                                  verbose = TRUE)
image(t(g1$probs[,,3]))


##### Map for P2 ####
s.p2 <- make_seq_mappoly(full.mds.seq.ch1, info.parent = "p2")
plot(s.p2)
tpt.p2 <- make_pairs_mappoly(full.tpt.ch1, s.p2)
plot(filt.full.mat.ch1, ord = s.p2)
map.p2 <- est_rf_hmm_sequential(input.seq = s.p2,
                                twopt = tpt.p2,
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
plot(map.p2)
g2 <- calc_genoprob_single_parent(map.p2,
                                  step = 0,
                                  info.parent = 2,
                                  uninfo.parent = 1,
                                  global.err = 0.0,
                                  phase.config = "best", 
                                  verbose = TRUE)
image(t(g2$probs[,,3]))


##### Merging P1 and P2 ####
map.p1
map.p2
full.seq <- full.seq.ch1
full.mat <- full.mat.ch1
map.p1.p2 <- merge_parental_maps(map.p1, 
                                 map.p2, 
                                 full.seq, full.mat, 
                                 method = "ols")
plot(map.p1.p2)
#####################round1##############################
map.1 <- map.p1.p2
gen.1 <- calc_genoprob(map.p1.p2)
mat.1 <- rf_list_to_matrix(full.tpt.ch1, 5, 5, shared.alleles = TRUE)
mrk.to.include <- setdiff(filt.full.seq.ch1$seq.mrk.names, map.1$info$mrk.names)
input.seq <- filt.full.seq.ch1
x.1 <- add_md_markers(input.map = map.1, 
                      input.seq = filt.full.seq.ch1,
                      mrk.to.include = mrk.to.include,
                      input.matrix = mat.1, 
                      input.genoprob = gen.1, 
                      input.data = dat, 
                      thresh = 100, 
                      method = "ols")


input.map = map.1
input.seq = full.mds.seq.ch1
mrk.to.include = mrk.to.include
input.matrix = mat.1
input.genoprob = gen.1
input.data = BT_trifida_dose
thresh = 100
method = "ols"
verbose = T


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
