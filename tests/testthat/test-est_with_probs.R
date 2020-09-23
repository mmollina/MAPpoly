context("Estimate map with probabilities")
test_that("map estimated correctly", {
  map <- get_submap(solcap.prior.map[[1]], 1:10, reestimate.rf = FALSE)
  mape <- est_full_hmm_with_global_error(map)
  mapp <- est_full_hmm_with_prior_prob(map)
  expect_equivalent(map$maps[[1]]$seq.rf, 
                    c(0.030928345, 0.000100000, 0.000100000, 0.006987686, 0.006978596,
                      0.009314603, 0.006089064, 0.014675723, 0.010374833))
  expect_equivalent(mape$maps[[1]]$seq.rf, 
                    c(0.022411061, 0.000010000, 0.000010000, 0.007022810, 0.007013837,
                      0.008713443, 0.005086579, 0.015316109, 0.003529484))
  expect_equivalent(mapp$maps[[1]]$seq.rf, 
                    c(0.022411061, 0.000010000, 0.000010000, 0.007022810, 0.007013837,
                      0.008713443, 0.005086579, 0.015316109, 0.003529484))
  map <- get_submap(solcap.dose.map[[1]], 1:10, reestimate.rf = FALSE)
  mape <- est_full_hmm_with_global_error(map)
  expect_equivalent(map$maps[[1]]$seq.rf, 
                    c(0.030964644, 0.000100000, 0.000100000, 0.006987554, 0.006978459,
                      0.009315223, 0.006088936, 0.014674755, 0.010377167))
  expect_error(est_full_hmm_with_prior_prob(map))
})