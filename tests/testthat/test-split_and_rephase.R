context("Split and rephase")
test_that("split and rephase the map correctly", {
  ##### Tetraploid
  map<-get_submap(solcap.err.map[[1]], 1:20)
  tpt<-est_pairwise_rf(make_seq_mappoly(map))
  map2<-split_and_rephase(map, tpt, gap.threshold = 2)
  expect_is(map2, "mappoly.map")
})
