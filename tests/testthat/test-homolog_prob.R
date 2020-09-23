context("Homolog probability")
test_that("computes homolog probabilities correctly", {
  x1 <- get_submap(solcap.dose.map[[1]], 1:20, reestimate.rf = F)
  probs.t1<-calc_genoprob(input.map = x1,
                         verbose = TRUE)
  expect_equivalent(round(var(probs.t1$probs), 6), 0.015261)
  hom.t1 <- calc_homoprob(input.genoprobs = probs.t1)
  expect_equal(round(var(hom.t1$homoprob$probability),6),0.187373)
  expect_is(plot(hom.t1), "plotly")
  expect_is(plot(hom.t1, stack = T), "plotly")
  expect_is(plot(hom.t1, use.plotly = FALSE), "ggplot")
})