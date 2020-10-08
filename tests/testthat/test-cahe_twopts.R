context("Compute genotype counts")
test_that("compute genotype counts correctly", {
  ##### Tetraploid
  s<-make_seq_mappoly(tetra.solcap, 3:7)
  o<-cache_counts_twopt(s)
  expect_is(o, "cache.info")
  expect_null(print(o))
  o2<-get_cache_two_pts_from_web(m = 4)
  expect_is(o2, "cache.info")
})
