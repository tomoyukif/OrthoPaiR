root_dir <- "/home/ftom/workspace/orthology"
working_dir <- file.path(root_dir, "dev_plant")
input_dir <- file.path(root_dir, "input/benchmark_plant")
devtools::load_all("/home/ftom/01_wd/softDevel/OrthoPaiR")

source("/home/ftom/workspace/orthology/script/00_functions/02_eval_pairing.R")
osat_hvul_busco_list <- read.csv(file.path(root_dir, "input/benchmark_plant/benchmark_list/osat_hvul_buscolist.csv"))
osat_hvul_del_list <- read.csv(file.path(root_dir, "input/benchmark_plant/benchmark_list/osat_hvul_falsepositive.csv"))

in_list <- orgInputFiles(name = "Osat",
                         genome = file.path(input_dir, "osat_genome.fa"),
                         gff = file.path(input_dir, "osat.gff"),
                         validation = FALSE)
in_list <- orgInputFiles(object = in_list,
                         name = "Hvul",
                         genome = file.path(input_dir, "hvul_genome.fa"),
                         gff = file.path(input_dir, "hvul.gff"),
                         validation = FALSE)
library(data.table)
source("R/dev/01_functions_init_opr.R")
opr <- init_opr(object = in_list, 
                working_dir = working_dir, 
                overwrite = T,
                n_threads = 30)

blast_path <- "/home/ftom/07_tools/bin"
source("R/dev/02_functions_rbh.R")
Rcpp::sourceCpp("inst/src/rbh.cpp")
rbh(object = opr,
    blast_path = blast_path,
    n_threads = 30,
    overwrite = T)

source("R/dev/03_functions_orthology.R")
out <- orthopair(working_dir = working_dir,
                 n_threads = 30,
                 load_all_gff = TRUE,
                 verbose = TRUE)
# ortholog_pairs    recall precision specificity  accuracy  f1_score  tp  fn    tn   fp
# 1          22118 0.8827508 0.2663265   0.9288917 0.9275808 0.4091978 783 104 28177 2157

og <- subset(out[[1]], select = c("original_genome1_gene", "original_genome2_gene"))
names(og)[1:2] <- c("query", "subject")
eval_busco_list(og = og,
                busco_list = osat_hvul_busco_list,
                del = osat_hvul_del_list,
                query = "query",
                subject = "subject")
# ortholog_pairs    recall precision specificity  accuracy  f1_score  tp fn    tn   fp
# 1          22509 0.9560316 0.2956764   0.9334081 0.9340508 0.4516644 848 39 28314 2020

blast_out <- fread(file = "~/workspace/orthology/dev_plant/blast/1_blast.out", sep = "\t")
df1 <- blast_out[1002000000L > blast_out$V1 & 1002000000L < blast_out$V2, c(1:3, 9, 4:7)]
df2 <- blast_out[1002000000L < blast_out$V1 & 1002000000L > blast_out$V2, c(1:3, 9, 4:7)]
names(df1) <- c("qseqid", "sseqid", "pident", "qcovs_q2s", 
                "qstart", "qend", "sstart", "send")
names(df2) <- c("sseqid", "qseqid", "pident", "qcovs_s2q", 
                "sstart", "send", "qstart", "qend")
rbh <- inner_join(df1, df2,
                  by = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "pident"),
                  relationship = "many-to-many")
rbh <- .orgBLASTout(rbh = rbh)
rbh$ci_q2s <- rbh$pident * rbh$qcovs_q2s * 1e-4
rbh$ci_s2q <- rbh$pident * rbh$qcovs_s2q * 1e-4
rbh$mutual_ci <- rbh$ci_q2s * rbh$ci_s2q
rbh <- rbh[order(-rbh$mutual_ci), 1:8]
# gff1 <- readRDS("/home/ftom/workspace/orthology/dev_plant/input/1_Osat/gff_df.rds")
# gff2 <- readRDS("/home/ftom/workspace/orthology/dev_plant/input/2_Hvul/gff_df.rds")
# hit <- match(rbh$qseqid, gff1[[1]]$tx_index)
# rbh$qgeneid <- gff1[[1]]$gene_index[hit]
# hit <- match(rbh$sseqid, gff2[[1]]$tx_index)
# rbh$sgeneid <- gff2[[1]]$gene_index[hit]
# rbh$pair_id <- paste(rbh$qgeneid, rbh$sgeneid, sep = "0")
# rbh <- subset(rbh, subset = !duplicated(pair_id))
# rbh <- lapply(rbh, as.numeric)
# rbh <- as.data.frame(rbh[, 1:8])
# rbh <- fread(file = "~/workspace/orthology/dev_plant/rbh/1001_1002.rbh",
#              header = FALSE,
#              sep = "\t",
#              select = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L),
#              col.names = c("query_tx", "subject_tx", "pident", "q2s_qcovs", 
#                            "s2q_qcovs", "q2s_ci", "s2q_ci", "mutual_ci", "len"),
#              na.strings = "",
#              colClasses = c("integer", "integer", "numeric", "integer", 
#                             "integer", "numeric", "numeric", "numeric", "integer"),
#              verbose = FALSE)
setDT(rbh)
names(rbh) <- c("query_tx", "subject_tx", "pident", "q2s_qcovs", 
                "s2q_qcovs", "q2s_ci", "s2q_ci", "mutual_ci")
rbh <- rbh[!is.na(query_tx) & !is.na(subject_tx) & query_tx != subject_tx]
query_genome_pre <- as.integer(substr(as.character(rbh$query_tx), 1L, genome_width))
subject_genome_pre <- as.integer(substr(as.character(rbh$subject_tx), 1L, genome_width))
if (!load_all_gff || is.null(gff_data_all)) {
    pair_genomes <- unique(c(query_genome_pre, subject_genome_pre))
    gff_data_pair <- lapply(as.character(pair_genomes), function(gid) {
        path <- gff_lookup[[gid]]
        if (is.null(path) || !file.exists(path)) return(NULL)
        return(readRDS(path))
    })
    gff_data_pair <- gff_data_pair[!sapply(gff_data_pair, is.null)]
} else {
    pair_genomes <- unique(c(query_genome_pre, subject_genome_pre))
    gff_data_pair <- lapply(as.character(pair_genomes), function(gid) {
        path <- gff_lookup[[gid]]
        if (is.null(path)) return(NULL)
        idx <- which(sapply(gff_data_all, function(g) {
            if (is.null(g)) return(FALSE)
            el <- if (is.list(g) && length(g) >= 1L) g[[1L]] else g
            if (!is.data.frame(el) || nrow(el) == 0L) return(FALSE)
            tx_idx <- as.character(el$tx_index[1L])
            if (nchar(tx_idx) >= genome_width) {
                gid_check <- as.integer(substr(tx_idx, 1L, genome_width))
                return(gid_check == as.integer(gid))
            }
            return(FALSE)
        }))
        if (length(idx) > 0L) return(gff_data_all[[idx[1L]]])
        return(NULL)
    })
    gff_data_pair <- gff_data_pair[!sapply(gff_data_pair, is.null)]
}
query_genome_id <- as.character(min(pair_genomes))
subject_genome_id <- as.character(max(pair_genomes))
query_gff_full <- NULL
subject_gff_full <- NULL
query_cds_full <- NULL
subject_cds_full <- NULL
for (el in gff_data_pair) {
    gff_df <- if (is.list(el) && length(el) >= 1L) el[[1L]] else el
    cds_df_el <- if (is.list(el) && length(el) >= 2L) el[[2L]] else NULL
    if (is.null(gff_df) || !is.data.frame(gff_df) || nrow(gff_df) == 0L) next
    tx_idx <- as.character(gff_df$tx_index[1L])
    if (nchar(tx_idx) >= genome_width) {
        gid <- as.character(as.integer(substr(tx_idx, 1L, genome_width)))
        if (gid == query_genome_id) {
            query_gff_full <- gff_df
            query_cds_full <- cds_df_el
        } else if (gid == subject_genome_id) {
            subject_gff_full <- gff_df
            subject_cds_full <- cds_df_el
        }
    }
}
query_gff_full <- query_gff_full[order(query_gff_full$seqnames, query_gff_full$start), ]
subject_gff_full <- subject_gff_full[order(subject_gff_full$seqnames, subject_gff_full$start), ]
query_gff_full$gene_index <- as.integer(query_gff_full$gene_index)
query_gff_full$tx_index <- as.integer(query_gff_full$tx_index)
subject_gff_full$gene_index <- as.integer(subject_gff_full$gene_index)
subject_gff_full$tx_index <- as.integer(subject_gff_full$tx_index)
g2g_graph <- .linkGene2Genome_dataframe(genome1_gff_df = query_gff_full,
                                        genome2_gff_df = subject_gff_full,
                                        genome1_cds_df = query_cds_full,
                                        genome2_cds_df = subject_cds_full)
query_tx_match <- match(rbh$query_tx, query_gff_full$tx_index)
subject_tx_match <- match(rbh$subject_tx, subject_gff_full$tx_index)
rbh_df <- data.frame(
    genome1_tx = as.integer(rbh$query_tx),
    genome2_tx = as.integer(rbh$subject_tx),
    genome1_gene = as.integer(query_gff_full$gene_index[query_tx_match]),
    genome2_gene = as.integer(subject_gff_full$gene_index[subject_tx_match]),
    ci_q2s = rbh$q2s_ci,
    ci_s2q = rbh$s2q_ci,
    pident = rbh$pident,
    mutual_ci = rbh$mutual_ci,
    stringsAsFactors = FALSE
)
rbh_df <- rbh_df[!is.na(rbh_df$genome1_gene) & !is.na(rbh_df$genome2_gene), ]
rbh_df$pair_id <- paste(rbh_df$genome1_gene, rbh_df$genome2_gene, sep = "_")
rbh_df <- subset(rbh_df, subset = !duplicated(rbh_df$pair_id))
anchor <- .findAnchors(rbh = rbh_df, g2g_graph = g2g_graph)
t2a_graph <- .link2Anchor(g2g_graph = g2g_graph, anchor = anchor)
orthopair <- .findSyntenicOrtho(rbh = rbh_df,
                                anchor = anchor,
                                g2g_graph = g2g_graph,
                                t2a_graph = t2a_graph,
                                rbh_threshold = NULL)
orthopair <- .findSyntenyBlocks(orthopair = orthopair)
orthopair <- .pickBestPair(orthopair = orthopair)
orthopair <- .classifyOrthoPair(orthopair = orthopair)
orthopair <- .filterOrthopair(orthopair = orthopair, g2g_graph = g2g_graph)
orthopair <- .classifyOrthoPair(orthopair = orthopair)
orthopair <- .reformatOrthoPair(orthopair = orthopair, g2g_graph = g2g_graph)

og <- subset(orthopair, select = c("original_genome1_gene", "original_genome2_gene"))
names(og)[1:2] <- c("query", "subject")
eval_busco_list(og = og,
                busco_list = osat_hvul_busco_list,
                del = osat_hvul_del_list,
                query = "query",
                subject = "subject")
# ortholog_pairs    recall precision specificity  accuracy  f1_score  tp fn    tn   fp
# 1          22362 0.9503946 0.2964135   0.9340344 0.9344992 0.4518896 843 44 28333 2001


source("R/dev/02_functions_rbh.R")
Rcpp::sourceCpp("inst/src/rbh.cpp")
dir.create(rbh_dir, showWarnings = FALSE, recursive = TRUE)
tmp_dir <- "~/workspace/orthology/dev_plant/rbh/tmp"
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
sorttmp_dir <- "~/workspace/orthology/dev_plant/rbh/sort_tmp"
dir.create(sorttmp_dir, showWarnings = FALSE, recursive = TRUE)
rbh_extract(blast_out_list = "/home/ftom/workspace/orthology/dev_plant/blast/1_blast.out",
            out_rbh_fn = "/home/ftom/workspace/orthology/dev_plant/rbh/rbh.debug.tsv",
            strategy = "auto",
            n_threads = 30,
            tmpdir = tmp_dir,
            sort_tmp = sorttmp_dir,
            final_sort_by_mutual_ci = TRUE,
            keep_intermediate = TRUE)

.splitRBHbyGenomePair(rbh_fn = "/home/ftom/workspace/orthology/dev_plant/rbh/rbh.debug.tsv",
                      rbh_dir = "/home/ftom/workspace/orthology/dev_plant/rbh_debug",
                      genome_width = 4L)
