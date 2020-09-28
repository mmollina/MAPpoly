context("Filter functions")
test_that("test filter functions", {
  a1<-sample_data(tetra.solcap, n = 20, type = "marker")
  expect_equal(check_data_sanity(a1), 0)
  a2<-sample_data(tetra.solcap.geno.dist, n = 20, type = "marker")
  expect_equal(check_data_sanity(a2), 0)
  a3<-filter_missing(input.data = a1, type = "marker", 
                    filter.thres = 0.1, inter = FALSE)
  expect_equal(check_data_sanity(a3), 0)
  a4<-filter_missing(input.data = a2, type = "marker", 
                     filter.thres = 0.1, inter = FALSE)
  expect_equal(check_data_sanity(a4), 0)
  a5<-filter_missing(input.data = a1, type = "individual", 
                     filter.thres = 0.1, inter = FALSE)
  expect_equal(check_data_sanity(a5), 0)
  a6<-filter_missing(input.data = a2, type = "marker", 
                     filter.thres = 0.1, inter = FALSE)
  expect_equal(check_data_sanity(a5), 0)
  a1$geno.dose[10, 10] <- 8
  a7<-filter_non_conforming_classes(input.data = a1)
  expect_equal(a7$geno.dose[10, 10],5)
  w1<-filter_segregation(tetra.solcap, chisq.pval.thres = 0.001, inter = FALSE)
  w2<-make_seq_mappoly(w1)
  expect_equal(length(w2$seq.num), 3659)
  w3<-filter_segregation(tetra.solcap, chisq.pval.thres = 1e-10, inter = FALSE)
  w4<-make_seq_mappoly(w3)
  expect_equal(length(w4$seq.num), length(tetra.solcap$mrk.names))
})