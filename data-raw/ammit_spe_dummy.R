## code to prepare `ammit_spe_dummy` dataset goes here

HALO_filepath <- system.file("extdata", "halo_output.csv", package="AMMIT")
markers <- paste0("M", 1:7)
locations <-  c("Nucleus", "Cytoplasm", "Nucleus", "Cytoplasm", "Cytoplasm", "Cytoplasm", "Cytoplasm")
ammit_spe_dummy <- AMMIT::build_ammit_spe_HALO(HALO_filepath = HALO_filepath, markers = markers, reference = T, locations = locations, filter_dapi = T)

usethis::use_data(ammit_spe_dummy, overwrite = TRUE)
