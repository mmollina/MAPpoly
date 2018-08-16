require(mappoly)
myfunc<-function(i)
{
  data(hexafake)
  mrk.subset<-make_seq_mappoly(hexafake, paste0("seq",i))
  red.mrk<-elim_redundant(mrk.subset)
  unique.mrks<-make_seq_mappoly(red.mrk)
  counts.web<-cache_counts_twopt(unique.mrks, get.from.web = TRUE)
  subset.pairs<-est_pairwise_rf(input.seq = unique.mrks,
                                count.cache = counts.web,
                                n.clusters = 1,
                                verbose=TRUE)
  subset.map <- est_rf_hmm_sequential(input.seq = unique.mrks,
                                      thres.twopt = 10,
                                      thres.hmm = 10,
                                      extend.tail = 50,
                                      tol = 0.1,
                                      tol.final = 10e-3,
                                      twopt = subset.pairs,
                                      verbose = TRUE, 
                                      info.tail = TRUE, 
                                      reestimate.single.ph.configuration = FALSE)
  retunr(subset.map)
}
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
clusterEvalQ(cl = cl, require(mappoly))
maps<-foreach(i=1:3) %dopar% myfunc(i)
probs<-foreach(i=1:3) %dopar% calc_genoprob(maps[[i]])
stopCluster(cl)
save(maps, file = "fake_maps.RData")
