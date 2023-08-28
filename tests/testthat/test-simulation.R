context("Simulate datasets")
test_that("simulate datasets correctly", {
  h.temp<-sim_homologous(ploidy=6, n.mrk=20)
  dat<-cross_simulate(ploidy=6, rf.vec=.05, n.mrk=20,
                                     n.ind=20, h.temp, seed=123)
  expect_is(dat, "mappoly.data")
})
