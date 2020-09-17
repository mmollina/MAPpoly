context("Get submap")
test_that("sub-map extracted correctly", {
  ##### Hexaploid
  sub.map1<-get_submap(solcap.dose.map[[1]], 1:5)
  expect_is(sub.map1, "mappoly.map")
  expect_output(str(sub.map1$info), "List of 11")
  expect_equivalent(sub.map1$maps[[1]]$seq.rf, 
                    c(3.113371e-02, 1.000000e-05, 1.838496e-05, 8.084546e-03))
})