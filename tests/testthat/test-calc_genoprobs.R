context("Conditional probabilities")
test_that("computes genotype probabilities correctly", {
  x1 <- get_submap(solcap.dose.map[[1]], 1:20, reestimate.rf = F)
  probs.t1<-calc_genoprob(input.map = x1,
                         verbose = TRUE)
  expect_equivalent(round(var(probs.t1$probs), 6), 0.015261)
  probs.t2<-calc_genoprob_error(input.map = x1,
                               verbose = TRUE, 
                               error = 0.2)
  expect_equivalent(round(var(probs.t2$probs), 6), 0.013076)
  expect_error(calc_genoprob_dist(input.map = x1,
                                  verbose = TRUE))
  x3 <- get_submap(solcap.prior.map[[1]], 1:20, reestimate.rf = F)
  probs.t3<-calc_genoprob_dist(input.map = x3,
                               verbose = TRUE)
  expect_equivalent(round(var(probs.t3$probs), 6), 0.015315)
})