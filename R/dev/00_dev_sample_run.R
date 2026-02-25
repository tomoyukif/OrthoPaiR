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

og <- subset(out[[1]], select = c("original_query_gene", "original_subject_gene"))
names(og)[1:2] <- c("query", "subject")
eval_busco_list(og = og,
                busco_list = nb_wk21_busco_list,
                del = del,
                query = "query",
                subject = "subject")
# ortholog_pairs    recall precision specificity  accuracy  f1_score   tp
# 1          46782 0.9913876 0.9980732   0.9985866 0.9955285 0.9947192 3108
# fn   tn fp
# 1 27 4239  6

nb_only_gene_id <- c("Os01g0258600",
                     "Os01g0700900",
                     "Os01g0701400",
                     "Os06g0165600",
                     "LOC_Os03g44710",
                     "Os07g0153600",
                     "Os11g0225300",
                     "Os11g0694600",
                     "Os11g0569733"
)
og[og$query %in% nb_only_gene_id, ]
wk_only_gene_id <- c("Ogla-WK21_12G0555300",
                     "Ogla-WK21_12G0555900")
og[og$subject %in% wk_only_gene_id, ]

og$pair_id <- paste(og$query, og$subject, sep = "_")
nb_wk21_busco_list$pair_id <- paste(nb_wk21_busco_list$query, nb_wk21_busco_list$subject, sep = "_")
missing_pair_id <- nb_wk21_busco_list$pair_id[!nb_wk21_busco_list$pair_id %in% og$pair_id]

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
root_dir <- "/home/ftom/workspace/orthology"
input_dir <- file.path(root_dir, "input/benchmark_plant")
devtools::load_all("/home/ftom/01_wd/softDevel/OrthoPaiR")

in_list <- orgInputFiles(name = "Osat",
                         genome = file.path(input_dir, "osat_genome.fa"),
                         gff = file.path(input_dir, "osat.gff"),
                         validation = FALSE)

in_list <- orgInputFiles(object = in_list,
                         name = "Hvul",
                         genome = file.path(input_dir, "hvul_genome.fa"),
                         gff = file.path(input_dir, "hvul.gff"),
                         validation = FALSE)

in_list <- orgInputFiles(object = in_list,
                         name = "Atha",
                         genome = file.path(input_dir, "atha_genome.fa"),
                         gff = file.path(input_dir, "atha.gff"),
                         validation = FALSE)

in_list <- orgInputFiles(object = in_list,
                         name = "Mpol",
                         genome = file.path(input_dir, "mpol_genome.fa"),
                         gff = file.path(input_dir, "mpol.gff"),
                         validation = FALSE)

working_dir <- file.path(root_dir, "dev_plant")
miniprot_path <- "/home/ftom/conda/envs/miniprot/bin"
blast_path <- "/home/ftom/07_tools/bin"
n_threads <- 32
overwrite <- TRUE
verbose <- TRUE

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

source("/home/ftom/workspace/orthology/script/00_functions/02_eval_pairing.R")
osat_hvul_busco_list <- read.csv(file.path(root_dir, "input/benchmark_plant/benchmark_list/osat_hvul_buscolist.csv"))
osat_atha_busco_list <- read.csv(file.path(root_dir, "input/benchmark_plant/benchmark_list/osat_atha_buscolist.csv"))
osat_mpol_busco_list <- read.csv(file.path(root_dir, "input/benchmark_plant/benchmark_list/osat_mpol_buscolist.csv"))
osat_hvul_del_list <- read.csv(file.path(root_dir, "input/benchmark_plant/benchmark_list/osat_hvul_falsepositive.csv"))
osat_atha_del_list <- read.csv(file.path(root_dir, "input/benchmark_plant/benchmark_list/osat_atha_falsepositive.csv"))
osat_mpol_del_list <- read.csv(file.path(root_dir, "input/benchmark_plant/benchmark_list/osat_mpol_falsepositive.csv"))
og <- subset(out[[1]], select = c("original_query_gene", "original_subject_gene"))
names(og)[1:2] <- c("query", "subject")
eval_busco_list(og = og,
                busco_list = osat_hvul_busco_list,
                del = osat_hvul_del_list,
                query = "query",
                subject = "subject")
og <- subset(out[[2]], select = c("original_query_gene", "original_subject_gene"))
names(og)[1:2] <- c("query", "subject")
eval_busco_list(og = og,
                busco_list = osat_atha_busco_list,
                del = osat_atha_del_list,
                query = "query",
                subject = "subject")
og <- subset(out[[3]], select = c("original_query_gene", "original_subject_gene"))
names(og)[1:2] <- c("query", "subject")
eval_busco_list(og = og,
                busco_list = osat_mpol_busco_list,
                del = osat_mpol_del_list,
                query = "query",
                subject = "subject")

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
