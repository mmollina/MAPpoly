#### Functions ####
require(mappoly)
require(tidyverse)
rm(list = ls())
source("my_tests/temp_utils.R")




##### Oredering #####
ch <- 3
dat <- readRDS("~/repos/current_work/rose/data/dat_BExMG.rds")
plot(dat)
lg <- readRDS(file = "~/repos/current_work/rose/data/lg_BExMG.rds")
full.seq.ch3 <- lg[[ch]]
plot(full.seq.ch3)
full.seq.ch3 <- make_seq_mappoly(filter_segregation(full.seq.ch3, chisq.pval.thres = 0.1, inter = FALSE))
full.tpt.ch3 <- est_pairwise_rf(full.seq.ch3, ncpus = 32)
plot(full.seq.ch3)
full.mat.ch3 <- rf_list_to_matrix(full.tpt.ch3, shared.alleles = TRUE, ncpus = 10)
filt.full.seq.ch3 <- rf_snp_filter(full.tpt.ch3, 
                                   thresh.LOD.ph = 5,
                                   thresh.LOD.rf = 5,
                                   probs = c(0.025,0.975), 
                                   ncpus = 10)
filt.full.mat.ch3 <- make_mat_mappoly(full.mat.ch3, filt.full.seq.ch3)
mds.ord.ch3 <- mds_mappoly(filt.full.mat.ch3)
full.mds.seq.ch3 <- make_seq_mappoly(mds.ord.ch3)
plot(filt.full.mat.ch3, ord = full.mds.seq.ch3, fact = 10)
##### Map for P1 ####
s.p1 <- make_seq_mappoly(full.mds.seq.ch3, info.parent = "p1")
plot(s.p1)
tpt.p1 <- make_pairs_mappoly(full.tpt.ch3, s.p1)
plot(filt.full.mat.ch3, ord = s.p1)
map.p1 <- est_rf_hmm_sequential(input.seq = s.p1,
                                twopt = tpt.p1,
                                start.set = 30,
                                thres.twopt = 10, 
                                thres.hmm = 10,
                                extend.tail = 10,
                                info.tail = TRUE, 
                                sub.map.size.diff.limit = 10, 
                                phase.number.limit = 20,
                                reestimate.single.ph.configuration = TRUE,
                                tol = 10e-2,
                                tol.final = 10e-4)
s.map.p1 <- make_seq_mappoly(map.p1)
plot(full.mat.ch3, ord = s.map.p1)
plot(map.p1)
g1 <- calc_genoprob_single_parent(map.p1,
                                  step = 0,
                                  info.parent = 1,
                                  uninfo.parent = 2,
                                  global.err = 0.0,
                                  phase.config = "best", 
                                  verbose = TRUE)
image(t(g1$probs[,,1]))
##### Map for P2 ####
s.p2 <- make_seq_mappoly(full.mds.seq.ch3, info.parent = "p2")
plot(s.p2)
tpt.p2 <- make_pairs_mappoly(full.tpt.ch3, s.p2)
plot(filt.full.mat.ch3, ord = s.p2)
map.p2 <- est_rf_hmm_sequential(input.seq = s.p2,
                                twopt = tpt.p2,
                                start.set = 5,
                                thres.twopt = 10, 
                                thres.hmm = 50,
                                extend.tail = 10,
                                info.tail = TRUE, 
                                sub.map.size.diff.limit = 10, 
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
full.seq <- full.mds.seq.ch3
full.mat <- full.mat.ch3
save.image("~/repos/current_work/rose/new_algorithm.rda")
map.p1.p2 <- merge_parental_maps(map.p1, 
                                 map.p2, 
                                 full.seq, full.mat, 
                                 method = "hmm")
plot(map.p1.p2)
#####################round1##############################

cur.mat <- rf_list_to_matrix(full.tpt.ch3, 10, 10, shared.alleles = TRUE)
cur.map <- map.p1.p2
cur.seq <- full.mds.seq.ch3
LOD <- 100
ll <- map.list <- vector("list", 500)
la <- numeric(500)
i <- 1
mrk.to.include <- setdiff(cur.seq$seq.mrk.names, cur.map$info$mrk.names)
while(length(mrk.to.include) > 0){
  cur.gen <- calc_genoprob(cur.map)
  cur.res <- add_md_markers(input.map = cur.map,
                        mrk.to.include = mrk.to.include,
                        input.seq = cur.seq,
                        input.matrix = cur.mat, 
                        input.genoprob = cur.gen, 
                        input.data = dat, 
                        input.mds = mds.ord.ch3,
                        thresh = LOD, 
                        method = "hmm", 
                        verbose = TRUE)
  cur.map <- map.list[[i]] <- cur.res$map
  ll[[i]] <- cur.res$ll
  plot(cur.res$ll, col = mp_pallet3(20)[i])
  abline(h = LOD, lty = 2)
  if(i != 1){
    if(map.list[[i-1]]$info$n.mrk == map.list[[i]]$info$n.mrk)
    {
      LOD <- max(cur.res$ll, na.rm = TRUE) - 10      
      if(LOD <= 0) LOD <- 10e-5
    }
  }
  la[i] <- LOD
  mrk.to.include <- setdiff(cur.seq$seq.mrk.names, cur.map$info$mrk.names)
  i <- i + 1
  cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  print(summary_maps(map.list[!sapply(map.list, is.null)], verbose = FALSE))
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")    
}


id <- which(la!=0)
la <- la[id]
ll <- ll[id]
map.list <- map.list[id]

i<-1
df<- NULL
for(i in 1:length(id)){
  x <- na.omit(ll[[i]][ll[[i]] > la[i]])
  df <- rbind(df, data.frame(snp = names(x), LOD = x, iter = as.factor(i)))
}
p<-ggplot(df, aes(x=iter, y=LOD, color=iter, group = iter)) +
  geom_boxplot()
plot_map_list(map.list)
plot(sapply(map.list, function(x) x$info$n.mrk), type = "b", col = 2)
abline(v = 8, col = 4)

a <- est_full_hmm_with_global_error(map.list[[8]], error = 0.05, verbose = T)
plot(a)

