input_dir <- "/home/ftom/01_wd/orthology/download"
in_files <- list.files(path = input_dir, full.names = TRUE, recursive = TRUE)
in_gff <- grep("\\.gff", in_files, value = TRUE)
in_genome <- grep("\\.genome|\\.fa", in_files, value = TRUE)
in_gff_id <- sub("\\..+", "", basename(in_gff))
in_genome_id <- sub("\\..+", "", basename(in_genome))

################################################################################
# Run OrthoPaiR
devtools::load_all("/home/ftom/01_wd/softDevel/OrthoPaiR")

in_list <- orgInputFiles(name = "NB",
                         genome = "/home/ftom/01_wd/genomeData/rice/cultivar_sativa/nb_combined/version_2023/fasta/nb_genome.fa",
                         gff = "/home/ftom/01_wd/genomeData/rice/cultivar_sativa/nb_combined/version_2023/gff/nb_combined_all.gff",
                         validation = FALSE)

in_list <- orgInputFiles(object = in_list, 
                         name = in_gff_id[1:29],
                         genome = in_genome[1:29],
                         gff = in_gff[1:29],
                         validation = FALSE)

working_dir <- "/home/ftom/workspace/orthology/rice_dev"
miniprot_path <- "/home/ftom/conda/envs/miniprot/bin"
blast_path <- "/home/ftom/07_tools/bin"
n_threads <- 32
overwrite <- TRUE
verbose <- TRUE
use_prot <- FALSE
orthopair <- TRUE
run_miniprot <- FALSE
reorg <- TRUE
makegraph <- TRUE
output_table <- TRUE

source("R/dev/01_functions_init_opr.R")
init_start <- Sys.time()
opr <- init_opr(object = in_list, 
                working_dir = working_dir, 
                overwrite = overwrite,
                n_threads = n_threads)
init_end <- Sys.time()

source("R/dev/02_functions_rbh.R")
Rcpp::sourceCpp("inst/src/rbh.cpp")
rbh_start <- Sys.time()
rbh(object = opr,
    blast_path = blast_path,
    n_threads = n_threads,
    overwrite = overwrite)
rbh_end <- Sys.time()

library(data.table)
source("R/dev/03_functions_orthology_direct.R")
opr_start <- Sys.time()
process_rbh_direct(working_dir = working_dir,
                   n_threads = n_threads,
                   load_all_gff = TRUE,
                   verbose = TRUE)
opr_end <- Sys.time()

################################################################################
input_dir <- "/home/ftom/01_wd/orthology/download"
root_dir <- "/home/ftom/workspace/orthology"
input_dir <- file.path(root_dir, "input/benchmark")
output_dir <- file.path(root_dir, "output/benchmark")
devtools::load_all("/home/ftom/01_wd/softDevel/OrthoPaiR")

in_list <- orgInputFiles(name = "Osat",
                         genome = file.path(input_dir, "nb_genome.fa"),
                         gff = file.path(input_dir, "nb.gff"),
                         validation = FALSE)

# Oryza glabberima
in_list <- orgInputFiles(object = in_list,
                         name = "Ogla",
                         genome = file.path(input_dir, "wk21_genome.fa"),
                         gff = file.path(input_dir, "wk21.gff"),
                         validation = FALSE)

working_dir <- file.path(root_dir, "dev_osat_ogla")
miniprot_path <- "/home/ftom/conda/envs/miniprot/bin"
blast_path <- "/home/ftom/07_tools/bin"
n_threads <- 32
overwrite <- TRUE
verbose <- TRUE
use_prot <- FALSE
orthopair <- TRUE
run_miniprot <- FALSE
reorg <- TRUE
makegraph <- TRUE
output_table <- TRUE

source("R/dev/01_functions_init_opr.R")
init_start <- Sys.time()
opr <- init_opr(object = in_list, 
                working_dir = working_dir, 
                overwrite = overwrite,
                n_threads = n_threads)
init_end <- Sys.time()

source("R/dev/02_functions_rbh.R")
Rcpp::sourceCpp("inst/src/rbh.cpp")
rbh_start <- Sys.time()
rbh(object = opr,
    blast_path = blast_path,
    n_threads = n_threads,
    overwrite = overwrite)
rbh_end <- Sys.time()

library(data.table)
source("R/dev/03_functions_orthology.R")
opr_start <- Sys.time()
out <- process_rbh_direct(working_dir = working_dir,
                   n_threads = n_threads,
                   load_all_gff = TRUE,
                   verbose = TRUE)
opr_end <- Sys.time()
opr_end - opr_start
opr_end - init_start

source("R/dev/04_functions_orthology_dataframe.R")
opr_start <- Sys.time()
out <- process_rbh_dataframe(working_dir = working_dir,
                          n_threads = n_threads,
                          load_all_gff = TRUE,
                          verbose = TRUE)
opr_end <- Sys.time()
opr_end - opr_start


old <- getOrthoPair(object)
old$Osat_Ogla$pair_id <- paste(old$Osat_Ogla$original_query_gene, old$Osat_Ogla$original_subject_gene, sep = "_")
out[[1]]$pair_id <- paste(out[[1]]$original_query_gene, out[[1]]$original_subject_gene, sep = "_")

missing_id <- old$Osat_Ogla$pair_id[!old$Osat_Ogla$pair_id %in% out[[1]]$pair_id]
nb_wk21_busco_list$pair_id <- paste(nb_wk21_busco_list$query, nb_wk21_busco_list$subject, sep = "_")
missing_pair_id <- nb_wk21_busco_list$pair_id[nb_wk21_busco_list$pair_id %in% missing_id]
