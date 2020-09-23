context("Preferential Pairing")
test_that("computes pairing probabilities correctly", {
  x1 <- get_submap(solcap.dose.map[[1]], 1:20, reestimate.rf = F)
  probs.t1<-calc_genoprob(input.map = x1,
                         verbose = TRUE)
  expect_equivalent(round(var(probs.t1$probs), 6), 0.015261)
  pref.t1 <- calc_prefpair_profiles(input.genoprobs = probs.t1)
  expect_equivalent(round(apply(pref.t1$prefpair.psi.pval, 2, mean),6), 
                    c(0.822695, 0.876216, 1.000000, 7.917049))
  expect_is(plot(pref.t1), "ggplot")
})