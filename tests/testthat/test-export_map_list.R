context("Export map list")
test_that("export map list correctly", {
  ##### Tetraploid
  expect_null(export_map_list(solcap.dose.map[1:2], ""))
})
