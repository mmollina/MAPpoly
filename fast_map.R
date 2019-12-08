require(mappoly)
load("~/repos/BT_map/src/denovo_mapping/all_LGs/ch12.RData")
source("~/repos/MAPpoly/my_func.R")
Tf.mrk<-grep("Tf", BT.trifida.triloba$mrk.names[map$maps[[1]]$seq.num], value = TRUE)
s12<-make_seq_mappoly(BT.trifida.triloba, Tf.mrk[1:40])
s12
counts<-cache_counts_twopt(s12, get.from.web = TRUE)
tpt12<-est_pairwise_rf(input.seq = s12,
                       count.cache = counts,
                       n.clusters = 1,
                       verbose=TRUE)
m12<-rf_list_to_matrix(tpt12, thresh.LOD.ph = 10, thresh.LOD.rf = 10, shared.alleles = T)
##Phasing and reestimating
bla<-system.time(map.lg12<-est_rf_hmm_sequential(input.seq = s12,
                                            start.set = 4,
                                            thres.twopt = 10,
                                            thres.hmm = 10,
                                            extend.tail = 10,
                                            twopt = tpt12,
                                            verbose = TRUE,
                                            tol = 10e-2,
                                            tol.final = 10e-3,
                                            phase.number.limit = 20,
                                            sub.map.size.diff.limit =  2,
                                            info.tail = TRUE,
                                            reestimate.single.ph.configuration = TRUE))

##Phasing and reestimating
ble<-system.time(
{s12.sub<-make_seq_mappoly(BT.trifida.triloba, s12$seq.num[1:4])
  map.init<-est_rf_hmm_sequential(input.seq = s12.sub,
                                  start.set = 4,
                                  thres.twopt = 10,
                                  thres.hmm = 10,
                                  extend.tail = 10,
                                  twopt = tpt12,
                                  verbose = TRUE,
                                  tol = 10e-2,
                                  tol.final = 10e-3,
                                  phase.number.limit = 20,
                                  sub.map.size.diff.limit =  2,
                                  info.tail = TRUE,
                                  reestimate.single.ph.configuration = TRUE)
  map.cur<-map.init
  probs.cur <- calc_genoprob(map.cur, verbose = FALSE)
  for(i in 5:10){
    cat(i, ": ")
    res<-my_map_func_test(hap1 = map.cur, 
                          hap2 = as.integer(s12$seq.num[i]), 
                          probs.hap1 = probs.hap1,
                          thresh.cut.path = 1/400, 
                          thresh.twopt = 5, 
                          tol = 10e-2,
                          rf.matrix = m12)
    
    ## include selection of the best phases
    ## e.g.
    # Phase_config.27 251.541470 0.009247245
    # Phase_config.12 249.251101 0.009425017
    # Phase_config.23 248.711777 0.009641096
    # Phase_config.26 247.318847 0.009367241
    # Phase_config.16 233.793064 0.010025366
    # Phase_config.15 229.401411 0.010221266
    # Phase_config.67 211.308975 0.013204765
    # res$stats[1,1]-res$stats[,1]

    map.cur$maps[[1]]$seq.num <- c(map.cur$maps[[1]]$seq.num, s12$seq.num[i])
    map.cur$maps[[1]]$seq.ph <- res$phase.config[[rownames(res$stats)[1]]]
    map.cur$maps[[1]]$seq.rf <- c(map.cur$maps[[1]]$seq.rf, res$stats[1,2])
    map.cur$maps[[1]]$loglike <- res$stats[1,1]
    map.cur$info$n.mrk <- length(map.cur$maps[[1]]$seq.num)
    probs.cur <- res$probs[[1]]
  }
  final.map <- reest_rf(input.map = map.cur, tol = 10e-3)
})
id<-names(map.cur$maps[[1]]$seq.ph$P)[1:10]
id
plot_compare_haplotypes(m = 6,
                        hom.allele.p1 = map.cur$maps[[1]]$seq.ph$P[id], hom.allele.q1 = map.cur$maps[[1]]$seq.ph$Q[id], 
                        hom.allele.p2 = map.lg12$maps[[1]]$seq.ph$P[id], hom.allele.q2 = map.lg12$maps[[1]]$seq.ph$Q[id])
bla;ble
