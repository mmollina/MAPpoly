context("Estimate HMM map")
test_that("map contructed correctly", {
  ##### Tetraploid
  s<-make_seq_mappoly(tetra.solcap, 3:7)
  expect_is(s, "mappoly.sequence")
  expect_output(str(s), "List of 14")
  tpt<-est_pairwise_rf(s)
  expect_is(tpt, "poly.est.two.pts.pairwise")
  expect_output(str(tpt), "List of 7")
  expect_output(str(tpt$pairwise), "List of 10")
  map <- est_rf_hmm(input.seq = s, twopt = tpt, thres = 2, tol = 10e-4)
  expect_is(map, "mappoly.map")
  expect_output(str(map), "List of 13")
  expect_equivalent(map$maps[[1]]$seq.rf, 
                    c(6.787204e-03, 1.283112e-03, 1.137237e-03, 3.272807e-05))
  expect_equivalent(map$maps[[1]]$seq.ph, 
                    list(P=list('3' = c(1, 2, 3, 4), '4' = c(1, 2, 3), '5' = c(1, 2, 3, 4), '6' = c(1, 2, 3), '7' = 0), 
                         Q=list('3' = c(1, 2, 3), '4' = c(1, 2), '5' = c(1, 2, 3), '6' = c(1, 2), '7' = 4)))
  expect_output(print(map), "This is an object of class 'mappoly.map'\\n    Ploidy level:\\t 4 \\n    No. individuals:\\t 160 \\n    No. markers:\\t 5 \\n    No. linkage phases:\\t 1 \\n\\n    ---------------------------------------------\\n    Number of linkage phase configurations:  1\\n    ---------------------------------------------\\n    Linkage phase configuration:  1\\n       map length:\\t 0.93\\n       log-likelihood:\\t -123.7\\n       LOD:\\t\\t 0\\n    ~~~~~~~~~~~~~~~~~~") 
})
test_that("sequential map contructed correctly", {
  ##### Tetraploid
  s<-make_seq_mappoly(tetra.solcap, 3:7)
  expect_is(s, "mappoly.sequence")
  expect_output(str(s), "List of 14")
  tpt<-est_pairwise_rf(s)
  expect_is(tpt, "poly.est.two.pts.pairwise")
  expect_output(str(tpt), "List of 7")
  expect_output(str(tpt$pairwise), "List of 10")
  map <- est_rf_hmm_sequential(input.seq = s, twopt = tpt, thres.twopt = 2, tol.final = 10e-4)
  expect_is(map, "mappoly.map")
  expect_output(str(map), "List of 13")
  expect_equivalent(map$maps[[1]]$seq.rf, 
                    c(6.787204e-03, 1.283112e-03, 1.137237e-03, 3.272807e-05))
  expect_equivalent(map$maps[[1]]$seq.ph, 
                    list(P=list('3' = c(1, 2, 3, 4), '4' = c(1, 2, 3), '5' = c(1, 2, 3, 4), '6' = c(1, 2, 3), '7' = 0), 
                         Q=list('3' = c(1, 2, 3), '4' = c(1, 2), '5' = c(1, 2, 3), '6' = c(1, 2), '7' = 4)))
  expect_output(print(map), "This is an object of class 'mappoly.map'\\n    Ploidy level:\\t 4 \\n    No. individuals:\\t 160 \\n    No. markers:\\t 5 \\n    No. linkage phases:\\t 1 \\n\\n    ---------------------------------------------\\n    Number of linkage phase configurations:  1\\n    ---------------------------------------------\\n    Linkage phase configuration:  1\\n       map length:\\t 0.93\\n       log-likelihood:\\t -123.7\\n       LOD:\\t\\t 0\\n    ~~~~~~~~~~~~~~~~~~") 
})