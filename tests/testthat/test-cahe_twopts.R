context("Compute genotype counts")
test_that("compute genotype counts correctly", {
  ##### Tetraploid
  s<-make_seq_mappoly(tetra.solcap, 3:7)
  o<-cache_counts_twopt(s)
  expect_is(o, "cache.info")
  expect_null(print(o))
})
