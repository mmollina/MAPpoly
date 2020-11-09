context("Read data")
test_that("read data from VCF file correctly", {
  fl = "https://github.com/mmollina/MAPpoly_vignettes/raw/master/data/sweet_sample_ch3.vcf.gz"
  tempfl <- tempfile(pattern = 'chr1_', fileext = '.vcf.gz')
  download.file(fl, destfile = tempfl)
  dat.dose.vcf = read_vcf(file = tempfl, parent.1 = "PARENT1", parent.2 = "PARENT2")
  expect_equal(check_data_sanity(dat.dose.vcf), 0)
  expect_null(print(dat.dose.vcf, detailed = TRUE))
})
test_that("read data from dosage file correctly", {
  fl1 = "https://raw.githubusercontent.com/mmollina/MAPpoly_vignettes/master/data/SolCAP_dosage"
  tempfl <- tempfile()
  download.file(fl1, destfile = tempfl)
  SolCAP.dose <- read_geno(file.in  = tempfl)
  expect_equal(check_data_sanity(SolCAP.dose), 0)
})
test_that("read data from probability file correctly", {
  ft="https://raw.githubusercontent.com/mmollina/MAPpoly_vignettes/master/data/hexa_sample"
  tempfl <- tempfile()
  download.file(ft, destfile = tempfl)
  SolCAP.dose.prob <- read_geno_prob(file.in  = tempfl)
  expect_equal(check_data_sanity(SolCAP.dose.prob), 0)
})
test_that("read data from CSV file correctly", {
  ft="https://raw.githubusercontent.com/mmollina/MAPpoly_vignettes/master/data/tetra_solcap.csv"
  tempfl <- tempfile()
  download.file(ft, destfile = tempfl)
  SolCAP.dose <- read_geno_csv(file.in  = tempfl, ploidy = 4)
  expect_equal(check_data_sanity(SolCAP.dose), 0)
})
test_that("read data from fitpoly file correctly", {
  fl <- "https://raw.githubusercontent.com/mmollina/MAPpoly_vignettes/master/data/fitpoly.dat"
  tempfl <- tempfile()
  download.file(fl, destfile = tempfl)
  fitpoly.dat <- read_fitpoly(file.in = tempfl, ploidy = 4, 
                               parent1 = "P1", parent2 = "P2", 
                               verbose = TRUE)
  expect_equal(check_data_sanity(fitpoly.dat), 0)
  expect_null(print(fitpoly.dat, detailed = TRUE))
})
