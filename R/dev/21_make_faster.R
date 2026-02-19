devtools::load_all("/home/ftom/01_wd/softDevel/OrthoPaiR")
working_dir <- "/home/ftom/workspace/orthology/rice/output"

hdf5_fn <- list.files(path = working_dir, 
                      pattern = ".h5",
                      full.names = TRUE,
                      recursive = TRUE)
hdf5_fn <- grep("reorg", hdf5_fn, value = TRUE, invert = TRUE)
n_core <- 30
rename <- TRUE
verbose <- TRUE
overwrite <- FALSE

reorgOrthopairs(hdf5_fn = hdf5_fn, 
                rename = rename,
                n_core = n_core,
                overwrite = overwrite, 
                verbose = verbose)
