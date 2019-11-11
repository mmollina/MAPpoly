require(mappoly)
h.temp <- sim_homologous(m=4, n.mrk=100, max.d=2, max.ph=2, seed=123)
rf.vec <- c(rep(c(rep(0,9),0.05),9),rep(0,9))
dat <- poly_cross_simulate(m = 4, rf.vec = rf.vec, n.mrk=100, n.ind=1000, hom.allele = h.temp, seed = 123)
plot(dat)
s <- make_seq_mappoly(dat, "all")
counts <- cache_counts_twopt(s, get.from.web = TRUE)
tpt<-est_pairwise_rf(input.seq = s,
                     count.cache = counts,
                     n.clusters = 1,
                     verbose=TRUE)
plot(rf_list_to_matrix(tpt))
map <- est_rf_hmm_sequential(input.seq = s,
                                    thres.twopt = 5,
                                    thres.hmm = 10,
                                    extend.tail = 10,
                                    tol = 0.1,
                                    tol.final = 10e-3,
                                    twopt = tpt,
                                    verbose = TRUE)
map <- reest_rf(input.map = map, tol = 10e-5, verbose = TRUE)
print(map, detailed = TRUE)
plot(map)

z<-list(get_submap(input.map = map, 1:10, reestimate.rf = FALSE),
        get_submap(input.map = map, 11:20, reestimate.rf = FALSE),
        get_submap(input.map = map, 21:30, reestimate.rf = FALSE),
        get_submap(input.map = map, 31:40, reestimate.rf = FALSE),
        get_submap(input.map = map, 41:50, reestimate.rf = FALSE),
        get_submap(input.map = map, 51:60, reestimate.rf = FALSE),
        get_submap(input.map = map, 61:70, reestimate.rf = FALSE),
        get_submap(input.map = map, 71:80, reestimate.rf = FALSE),
        get_submap(input.map = map, 81:90, reestimate.rf = FALSE),
        get_submap(input.map = map, 91:100, reestimate.rf = FALSE))
print(z[[2]], detailed = TRUE)

mat.rec <- mat.lod <- matrix(NA, length(z), length(z))
for(i in 1:(length(z)-1)){
  for(j in (i+1):length(z)){
    cat("\nComputing rf between marker blocks", i , "and", j, "...")
    map1 <- z[[i]]
    map2 <- z[[j]]
    twopt.sub <- make_pairs_mappoly(tpt, 
                                    make_seq_mappoly(dat, 
                                                     c(map1$maps[[1]]$seq.num, 
                                                       map2$maps[[1]]$seq.num)))
    #plot(rf_list_to_matrix(twopt.sub))
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
    fields::image.plot(imf_h(mat.rec), col = rev(viridis::inferno(20)))
  }
}
mat.rec[mat.rec > .1]<-NA



