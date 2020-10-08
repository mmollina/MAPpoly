context("Estimate map with probabilities")
test_that("map estimated correctly", {
  map <- get_submap(solcap.prior.map[[1]], c(1,3,5,7,9))
  mape <- est_full_hmm_with_global_error(map, to = 10e-4)
  mapp <- est_full_hmm_with_prior_prob(map, to = 10e-4)
  expect_equivalent(map$maps[[1]]$seq.rf, 
                    c(0.037783070, 0.001634628, 0.015620733, 0.016022936))
  expect_equivalent(mape$maps[[1]]$seq.rf, 
                    c(0.02923187, 0.00162774, 0.01649081, 0.01506356))
  expect_equivalent(mapp$maps[[1]]$seq.rf, 
                    c(0.02923187, 0.00162774, 0.01649081, 0.01506356))
  map <- get_submap(solcap.dose.map[[1]], 1:5)
  expect_error(est_full_hmm_with_prior_prob(map))
})