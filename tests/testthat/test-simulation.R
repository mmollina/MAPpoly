context("Simulate datasets")
test_that("simulate datasets correctly", {
  h.temp<-sim_homologous(ploidy=6, n.mrk=20)
  dat<-cross_simulate(h.temp, 100, n.ind = 100)
  expect_is(dat, "mappoly.data")
})
