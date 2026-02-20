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

library(data.table)
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

source("R/dev/03_functions_orthology.R")
opr_start <- Sys.time()
out <- orthopair(working_dir = working_dir,
                 n_threads = n_threads,
                 load_all_gff = TRUE,
                 verbose = TRUE)
opr_end <- Sys.time()
init_runtime <- init_end - init_start
rbh_runtime <- rbh_end - rbh_start
opr_runtime <- opr_end - opr_start
overall_runtime <- c(init_runtime + rbh_runtime + opr_runtime)
overall_runtime

new <- fread(file = "/home/ftom/workspace/orthology/dev_osat_ogla/orthopair/1001_1002.tsv", sep = "\t")
old <- getOrthoPair(hdf5_fn = "/home/ftom/workspace/orthology/output/benchmark/orthopair/hdf5_out/Osat_Ogla.h5", score = T)
old$Osat_Ogla$pair_id <- paste(old$Osat_Ogla$original_query_gene, old$Osat_Ogla$original_subject_gene, sep = "_")
new$pair_id <- paste(new$original_query_gene, new$original_subject_gene, sep = "_")

missing_id <- old$Osat_Ogla$pair_id[!old$Osat_Ogla$pair_id %in% new$pair_id]
nb_wk21_busco_list$pair_id <- paste(nb_wk21_busco_list$query, nb_wk21_busco_list$subject, sep = "_")
missing_pair_id <- nb_wk21_busco_list$pair_id[nb_wk21_busco_list$pair_id %in% missing_id]
View(old$Osat_Ogla[old$Osat_Ogla$pair_id %in% missing_pair_id, ])

new_rbh <- fread(file = "/home/ftom/workspace/orthology/dev_osat_ogla/rbh/1001_1002.rbh", sep = "\t")
os <- readRDS("/home/ftom/workspace/orthology/dev_osat_ogla/input/1_Osat/gff_df.rds")
og <- readRDS("/home/ftom/workspace/orthology/dev_osat_ogla/input/2_Ogla/gff_df.rds")
hit <- match(new_rbh$V1, os[[1]]$tx_index)
new_rbh$qgene <- os[[1]]$gene_id[hit]
hit <- match(new_rbh$V2, og[[1]]$tx_index)
new_rbh$sgene <- og[[1]]$gene_id[hit]
new_rbh$pair_id <- paste(new_rbh$qgene, new_rbh$sgene, sep = "_")
View(new_rbh[new_rbh$pair_id %in% missing_pair_id, ])

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

library(data.table)
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

source("R/dev/03_functions_orthology.R")
opr_start <- Sys.time()
out <- orthopair(working_dir = working_dir,
                 n_threads = n_threads,
                 load_all_gff = TRUE,
                 verbose = TRUE)
opr_end <- Sys.time()
init_runtime <- init_end - init_start
rbh_runtime <- rbh_end - rbh_start
opr_runtime <- opr_end - opr_start
overall_runtime <- c(init_runtime + rbh_runtime + opr_runtime)
overall_runtime
