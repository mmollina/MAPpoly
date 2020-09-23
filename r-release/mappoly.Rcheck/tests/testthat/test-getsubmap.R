context("Get submap")
test_that("sub-map extracted correctly", {
  ##### Hexaploid
  sub.map1<-get_submap(maps.hexafake[[1]], 1:10)
  expect_is(sub.map1, "mappoly.map")
  expect_output(str(sub.map1$info), "List of 13")
  expect_equivalent(sub.map1$maps[[1]]$seq.rf, 
                    c(4.296132e-04, 1.995705e-04, 5.958305e-05, 6.214492e-04,
                      1.948445e-04, 6.957085e-05, 4.989336e-05, 5.372199e-05,
                      1.096427e-03))
})