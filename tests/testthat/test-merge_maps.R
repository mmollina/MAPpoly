context("Merge maps")
test_that("merging maps correctly", {
  ##### Tetraploid
  map1<-get_submap(solcap.dose.map[[1]], 1:5, reestimate.rf = FALSE, verbose = FALSE)
  map2<-get_submap(solcap.dose.map[[1]], 6:15, reestimate.rf = FALSE, verbose = FALSE)
  map3<-get_submap(solcap.dose.map[[1]], 16:30, reestimate.rf = FALSE, verbose = FALSE)
  full.map<-get_submap(solcap.dose.map[[1]], 1:30, reestimate.rf = FALSE, verbose = FALSE)
  s<-make_seq_mappoly(tetra.solcap, full.map$maps[[1]]$seq.num)
  twopt <- est_pairwise_rf(input.seq = s)
  merged.maps<-merge_maps(map.list = list(map1, map2, map3), 
                          twopt = twopt,
                          thres.twopt = 3)
  expect_equal(round(mean(unlist(merged.maps$maps[[1]]$seq.ph)),6), 2.495868)
  best.phase <- merged.maps$maps[[1]]$seq.ph
  names.id<-names(best.phase$P)
  x1 <- compare_haplotypes(m = 4, best.phase$P[names.id], 
                           full.map$maps[[1]]$seq.ph$P[names.id]) 
  x2 <- compare_haplotypes(m = 4, best.phase$Q[names.id], 
                           full.map$maps[[1]]$seq.ph$Q[names.id])
  expect_true(x1$is.same.haplo)
  expect_true(x2$is.same.haplo)
})
