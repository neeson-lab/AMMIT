test_that("build from csv works", {
  HALO_filepath <- testthat::test_path("testdata", "halo_output.csv")
  markers <- paste0("M", 1:7)
  locations <-  c("Nucleus", "Cytoplasm", "Nucleus", "Cytoplasm", "Cytoplasm", "Cytoplasm", "Cytoplasm")
  spe <- AMMIT::build_ammit_spe_HALO(HALO_filepath = HALO_filepath, markers = markers, reference = T, locations = locations, filter_dapi = T)
  data("ammit_spe_dummy")
  expect_equal(spe, ammit_spe_dummy)
})
