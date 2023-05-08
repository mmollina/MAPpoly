#### Loading packages and functions ####
.libPaths(c( "/home3/mmollin/R/rocker-rstudio/4.1", .libPaths()))
require(mappoly)
require(tidyverse)
require(ggpubr)
rm(list = ls())
source("~/repos/official_repos/MAPpoly/my_tests/temp_utils.R")

setwd("~/repos/current_work/rose")
#### Reading dataset ####
#dat_BExMG <- dat.fs[[1]]
#dat_SWxBE <- dat.fs[[2]]
dat.fs <- readRDS(file="data/multiple_fulsib_rose.rds")
dat <- dat.fs[[1]]
s <- make_seq_mappoly(dat, "seq2")
tpt <- est_pairwise_rf(s, ncpus = 32)
mat <- rf_list_to_matrix(tpt, ncpus = 10)
sf <- rf_snp_filter(tpt, 
                    thresh.LOD.ph = 5,
                    thresh.LOD.rf = 5,
                    probs = c(0.05,0.95), 
                    ncpus = 10)
mat2 <- make_mat_mappoly(mat, sf)
mds.o <- mds_mappoly(mat2)
#seq.init <- make_seq_mappoly(get_genomic_order(sf))
seq.init <- make_seq_mappoly(mds.o)
pdf("matrix.pdf", width = 10, height = 10)
plot(mat, ord = seq.init, fact = 10)
dev.off()
t1<-system.time(init.map.list <- build_map_1(seq.init, 
                                         tpt, 
                                         start.set = 5,
                                         inflation.lim.p1 = 10, 
                                         inflation.lim.p2 = 10, 
                                         thres.twopt = 7,
                                         thres.hmm = 10,
                                         phase.number.limit =10))

print(paste(round(t1[3]/60, 4), "minutes"))

pdf("map_single_parents.pdf", width = 10, height = 10)
plot_map_list(init.map.list)
dev.off()

init.map.list.err <- map_error_model(init.map.list, error = 0.05, tol = 10e-3)
# plot_map_list(init.map.list.err)
# init.map.list.err$map.p1.p2 <- rem_mrk_clusters(init.map.list.err$map.p1.p2, size.rem.cluster = 1, 1)
# init.map.list.err <- map_error_model(init.map.list.err, error = 0.05, tol = 10e-3)
pdf("map_single_parents_err.pdf", width = 10, height = 10)
plot_map_list(init.map.list.err)
dev.off()

t2<-system.time(map.list <- build_map_2(init.map.list.err, 
                                    seq.init, 
                                    tpt, 
                                    thres.twopt = 10, 
                                    init.LOD = 100))

print(paste(round(t2[3]/60, 4), "minutes"))

pdf("map_list.pdf", width = 10, height = 10)
plot(map.list)
dev.off()

i <- which(duplicated(sapply(map.list$both$map.list, function(x) x$info$n.mrk)))
map.list.reduced <- map.list
map.list.reduced$both$map.list <- map.list.reduced$both$map.list[-i]
map.list.reduced$both$lod.thresh <- map.list.reduced$both$lod.thresh[-i]
map.list.reduced$both$calc.lod <- map.list.reduced$both$calc.lod[-i]

pdf("map_list_reduced.pdf", width = 10, height = 10)
plot(map.list.reduced)
dev.off()

final.map <- map_error_model(map.list.reduced, error = 0.15, tol = 10e-4)
pdf("final_map.pdf", width = 20, height = 20)
plot(final.map)
dev.off()
summary_maps(final.map$both$map.list.err)
plot(final.map$both$map.list.err[[1]])
save.image(file = "Ch2_map.rda")

