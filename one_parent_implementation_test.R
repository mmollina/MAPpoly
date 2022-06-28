require(mappoly)
require(fields)
ploidy <- 4
n.mrk <- 100
map.len <- 100
h.temp <- sim_homologous(ploidy = ploidy, n.mrk = n.mrk, max.d = ploidy/2, max.ph = 0, seed = 123)
dat <- poly_cross_simulate(ploidy = ploidy, rf.vec = mf_h(map.len/n.mrk), n.mrk = n.mrk, n.ind = 200, h.temp, seed = 123)

## Recombination fraction matrix
s <- make_seq_mappoly(dat, "all")
rf <- est_pairwise_rf2(s)
M <- rf$pairwise$rf
Lr <- rf$pairwise$LOD_rf
Lp <- rf$pairwise$LOD_ph
M[Lp < 5 & Lr < 5] <- NA
image.plot(M, col = rev(tim.colors(64)))

## Map reconstruction given parental phase
s.all <- make_seq_mappoly(dat, "all")
s.all
mapP1 <- est_rf_hmm_single_one_parent(input.seq = s.all, 
                                      input.ph.single = list(P1 = h.temp$hom.allele.p,
                                                             P2 = h.temp$hom.allele.q), 
                                      info.parent = 1, 
                                      uninfo.parent = 2,
                                      highprec = TRUE,  
                                      verbose = TRUE)
mapP1
plot(mapP1)
mapP2 <- est_rf_hmm_single_one_parent(input.seq = s.all, 
                                      input.ph.single = list(P1 = h.temp$hom.allele.p,
                                                             P2 = h.temp$hom.allele.q), 
                                      info.parent = 2, 
                                      uninfo.parent = 1,  
                                      verbose = TRUE)
mapP2
plot(mapP2)

map5 <- get_submap(mapP1, mrk.pos = 1:5, reestimate.rf = F)
plot(map5)
input.ph.single <- map5$maps[[1]]$seq.ph
input.seq <- make_seq_mappoly(map5)


