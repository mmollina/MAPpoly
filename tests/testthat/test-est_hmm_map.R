context("Estimate HMM map")
test_that("map contructed correctly", {
  ##### Tetraploid
  s<-make_seq_mappoly(tetra.solcap, 1:2)
  map <- est_rf_hmm(input.seq = make_seq_mappoly(tetra.solcap, 1:2), 
                    twopt =est_pairwise_rf(s),
                    thres = 2, tol = 10e-4)
  map <- loglike_hmm(map)
  map <- est_full_hmm_with_global_error(map, error = .05)
  map <- est_full_hmm_with_prior_prob(map, dat.prob = tetra.solcap.geno.dist)
  s<-make_seq_mappoly(tetra.solcap, 3:7)
  expect_is(s, "mappoly.sequence")
  expect_output(str(s), "List of 14")
  tpt<-est_pairwise_rf(s)
  expect_is(tpt, "poly.est.two.pts.pairwise")
  expect_output(str(tpt), "List of 7")
  expect_output(str(tpt$pairwise), "List of 10")
  expect_is(map, "mappoly.map")
  map <- est_rf_hmm(input.seq = s, twopt = tpt, thres = 2, tol = 10e-4)
  print(map, detailed = TRUE)
  expect_is(map, "mappoly.map")
  map <- loglike_hmm(map)
  expect_output(str(map), "List of 13")
  expect_equivalent(map$maps[[1]]$seq.rf, 
                    c(6.787204e-03, 1.283112e-03, 1.137237e-03, 3.272807e-05))
  expect_equivalent(map$maps[[1]]$seq.ph, 
                    list(P=list('3' = c(1, 2, 3, 4), '4' = c(1, 2, 3), '5' = c(1, 2, 3, 4), '6' = c(1, 2, 3), '7' = 0), 
                         Q=list('3' = c(1, 2, 3), '4' = c(1, 2), '5' = c(1, 2, 3), '6' = c(1, 2), '7' = 4)))
  expect_output(print(map), "This is an object of class 'mappoly.map'\\n    Ploidy level:\\t 4 \\n    No. individuals:\\t 160 \\n    No. markers:\\t 5 \\n    No. linkage phases:\\t 1 \\n\\n    ---------------------------------------------\\n    Number of linkage phase configurations:  1\\n    ---------------------------------------------\\n    Linkage phase configuration:  1\\n       map length:\\t 0.93\\n       log-likelihood:\\t -123.7\\n       LOD:\\t\\t 0\\n    ~~~~~~~~~~~~~~~~~~") 
  expect_length(plot(map), 2)
  expect_length(plot(map, left.lim = .2, right.lim = .9), 2)
  expect_equal(round(mean(plot_map_list(list(map, map))[,3]),6), 0.669977)
  expect_equal(round(mean(plot_map_list(list(map, map), horiz = F, col = "ggstyle")[,3]),6), 0.669977)
  expect_is(plot_genome_vs_map(list(map, map)), "ggplot")
  expect_is(plot_genome_vs_map(list(map, map), same.ch.lg = TRUE), "ggplot")
  expect_equivalent(summary_maps(list(map, map))[,3], c("0.93", "0.93", "1.86"))
})
test_that("sequential map contructed correctly", {
  s<-make_seq_mappoly(tetra.solcap, 1:2)
  map <- est_rf_hmm_sequential(input.seq = s, 
                               twopt = est_pairwise_rf(s), 
                               verbose = FALSE)
  ##### Tetraploid
  s<-make_seq_mappoly(tetra.solcap, 3:7)
  expect_is(s, "mappoly.sequence")
  expect_output(str(s), "List of 14")
  tpt<-est_pairwise_rf(s)
  expect_is(tpt, "poly.est.two.pts.pairwise")
  expect_output(str(tpt), "List of 7")
  expect_output(str(tpt$pairwise), "List of 10")
  map <- est_rf_hmm_sequential(input.seq = s, twopt = tpt, thres.twopt = 2, tol.final = 10e-4, detailed.verbose = TRUE)
  map <- est_rf_hmm_sequential(input.seq = s, twopt = tpt, thres.twopt = 2, tol.final = 10e-4, high.prec = TRUE)
  expect_is(map, "mappoly.map")
  expect_output(str(map), "List of 13")
  expect_equivalent(map$maps[[1]]$seq.rf, 
                    c(6.787204e-03, 1.283112e-03, 1.137237e-03, 3.272807e-05))
  expect_equivalent(map$maps[[1]]$seq.ph, 
                    list(P=list('3' = c(1, 2, 3, 4), '4' = c(1, 2, 3), '5' = c(1, 2, 3, 4), '6' = c(1, 2, 3), '7' = 0), 
                         Q=list('3' = c(1, 2, 3), '4' = c(1, 2), '5' = c(1, 2, 3), '6' = c(1, 2), '7' = 4)))
  expect_output(print(map), "This is an object of class 'mappoly.map'\\n    Ploidy level:\\t 4 \\n    No. individuals:\\t 160 \\n    No. markers:\\t 5 \\n    No. linkage phases:\\t 1 \\n\\n    ---------------------------------------------\\n    Number of linkage phase configurations:  1\\n    ---------------------------------------------\\n    Linkage phase configuration:  1\\n       map length:\\t 0.93\\n       log-likelihood:\\t -123.7\\n       LOD:\\t\\t 0\\n    ~~~~~~~~~~~~~~~~~~") 
})