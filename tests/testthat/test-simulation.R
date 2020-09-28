context("Simulate datasets")
test_that("simulate datasets correctly", {
  h.temp<-sim_homologous(m=6, n.mrk=20, max.d=3, max.ph=3, seed=123)
  dat<-poly_cross_simulate(m=6, rf.vec=.05, n.mrk=20,
                                     n.ind=20, h.temp, seed=123)
  expect_is(dat, "mappoly.data")
})
