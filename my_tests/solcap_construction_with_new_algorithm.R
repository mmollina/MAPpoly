require("mappoly")
require("tidyverse")
require("ggpubr")
source("~/repos/official_repos/MAPpoly/my_tests/temp_utils.R")
ch <- 2
input.seq <- make_seq_mappoly(tetra.solcap, paste0("seq",ch))
twopt <- est_pairwise_rf(input.seq)
init.map.list <- build_map_1(input.seq, 
                             twopt, 
                             start.set = 5,
                             inflation.lim.p1 = 15, 
                             inflation.lim.p2 = 10)
plot_map_list(init.map.list)
res <- build_map_2(init.map.list, 
                   input.seq, 
                   twopt, 
                   thres.twopt = 5, 
                   init.LOD = 30)
plot(res)
iter <- 6
map <- est_full_hmm_with_global_error(res$both$map.list[[iter]], error = 0.05, verbose = T)
plot(map, mrk.names = T, cex = .4)

ph.new.p1 <- res$both$map.list[[iter]]$maps[[1]]$seq.ph$P
names(ph.new.p1) <- res$both$map.list[[iter]]$info$mrk.names
ph.old.p1 <- solcap.err.map[[ch]]$maps[[1]]$seq.ph$P
names(ph.old.p1) <- solcap.err.map[[ch]]$info$mrk.names
ph.new.p2 <- res$both$map.list[[iter]]$maps[[1]]$seq.ph$Q
names(ph.new.p2) <- res$both$map.list[[iter]]$info$mrk.names
ph.old.p2 <- solcap.err.map[[ch]]$maps[[1]]$seq.ph$Q
names(ph.old.p2) <- solcap.err.map[[ch]]$info$mrk.names
id <- intersect(names(ph.old.p1), names(ph.new.p1))
compare_haplotypes(4, ph.old.p1[id], ph.new.p1[id])
compare_haplotypes(4, ph.old.p2[id], ph.new.p2[id])
plot_compare_haplotypes(4, ph.old.p1[id], ph.old.p2[id], ph.new.p1[id], ph.new.p2[id])
compare_maps(res$both$map.list[[2]], solcap.err.map[[ch]])

