load("~/repos/MAPpoly_vignettes/vignette_1/maps.rda")

solcap.file <- system.file('extdata', 'tetra_solcap_geno_dist.bz2', package = 'mappoly')
tetra.solcap.geno.dist <- read_geno_dist(file.in  = solcap.file, prob.thres = 0.95)

maps1<-MAPs
for(i in 1:12){
  maps1[[i]]$info$data.name<-"dat.dist.mpl"
  maps1[[i]]$maps[[1]]$seq.num <- match(dat.dose.filt$mrk.names[maps1[[i]]$maps[[1]]$seq.num], dat.dist.mpl$mrk.names)
  names(maps1[[i]]$maps[[1]]$seq.ph$P) <- names(maps1[[i]]$maps[[1]]$seq.ph$Q) <- maps1[[i]]$maps[[1]]$seq.num
}

maps2<-MAP.err
for(i in 1:12){
  maps2[[i]]$info$data.name<-"dat.dist.mpl"
  maps2[[i]]$maps[[1]]$seq.num <- match(dat.dose.filt$mrk.names[maps2[[i]]$maps[[1]]$seq.num], dat.dist.mpl$mrk.names)
  names(maps2[[i]]$maps[[1]]$seq.ph$P) <- names(maps2[[i]]$maps[[1]]$seq.ph$Q) <- maps2[[i]]$maps[[1]]$seq.num
}

w2<-lapply(maps1, calc_genoprob)
h.prob.solcap<-calc_homoprob(w2)
print(h.prob.solcap)
plot(h.prob.solcap, ind = "ind_10")

w2.err<-lapply(maps2, calc_genoprob_error)
h.prob.solcap.err<-calc_homoprob(w2.err)
print(h.prob.solcap.err)
plot(h.prob.solcap.err, ind = "ind_10")

maps3<-vector("list", 12)
for(i in 1:12)
maps3[[i]]<-est_full_hmm_with_prior_dist(input.map = maps1[[i]])


for(i in 1:12){
  maps1[[i]]$info$data.name<-"tetra.solcap.geno.dist"  
  maps2[[i]]$info$data.name<-"tetra.solcap.geno.dist"
  maps3[[i]]$info$data.name<-"tetra.solcap.geno.dist"
  
}

est_full_hmm_with_global_error(maps1[[1]], verbose = TRUE, tol = .01)
est_full_hmm_with_prior_dist(maps1[[1]], tol = .01)



solcap.dose.map<-maps1
solcap.err.map<-maps2
solcap.prior.map<-maps3  


my.phase.func<-function(X){
  x<-est_rf_hmm_sequential(input.seq = X$lg,
                           start.set = 10,
                           thres.twopt = 10, 
                           thres.hmm = 10,
                           extend.tail = 200,
                           info.tail = TRUE, 
                           twopt = X$tpt,
                           sub.map.size.diff.limit = 8, 
                           phase.number.limit = 10,
                           reestimate.single.ph.configuration = TRUE,
                           tol = 10e-3,
                           tol.final = 10e-4)
  return(x)
}
system.time({
  cl <- parallel::makeCluster(12)
  parallel::clusterEvalQ(cl, require(mappoly))
  parallel::clusterExport(cl, "dat.dose.filt")
  MAPs.mds <- parallel::parLapply(cl,LGS.mds,my.phase.func)
  parallel::stopCluster(cl)
})

maps1<-MAPs
for(i in 1:12){
  maps1[[i]]$info$data.name<-"tetra.solcap.geno.dist"
}


solcap.mds.map<-maps1
plot_map_list(solcap.mds.map)

bla<-calc_genoprob(xdos[[1]], verbose = TRUE)




