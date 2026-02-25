os <- readDNAStringSet("/home/ftom/workspace/orthology/output/benchmark/orthopair/input/Osat/cds.fa")
names(os) <- 1001000000 + as.numeric(names(os))
og <- readDNAStringSet("/home/ftom/workspace/orthology/output/benchmark/orthopair/input/Ogla/cds.fa")
names(og) <- 1002000000 + as.numeric(names(og))
writeXStringSet(c(os, og), "/home/ftom/workspace/orthology/dev_osat_ogla/test/all_cds_renamed.fa")

blast_args <- paste("-in", "/home/ftom/workspace/orthology/dev_osat_ogla/test/all_cds_renamed.fa", 
                    "-dbtype nucl", 
                    "-out", "/home/ftom/workspace/orthology/dev_osat_ogla/test/all_cds_renamed.blastdb")
system2(command = file.path(blast_path, "makeblastdb"),
        args = blast_args, 
        stdout = FALSE)

blast_args <- paste(paste("-query", "/home/ftom/workspace/orthology/dev_osat_ogla/test/all_cds_renamed.fa"),
                    paste("-db", "/home/ftom/workspace/orthology/dev_osat_ogla/test/all_cds_renamed.blastdb"),
                    "-task blastn -max_target_seqs 10000",
                    "-evalue 1e-4 -strand plus",
                    "-outfmt '6 qseqid sseqid pident qcovs qstart qend sstart send qlen slen'",
                    paste("-num_threads", 30),
                    "-out", "/home/ftom/workspace/orthology/dev_osat_ogla/test/all_cds_renamed.blast.out")
system2(command = file.path(blast_path, "blastn"),
        args = blast_args)

blast_out <- fread(file = "~/workspace/orthology/dev_osat_ogla/blast/1_blast.out", sep = "\t")
df1 <- blast_out[1002000000L > blast_out$V1 & 1002000000L < blast_out$V2, 1:8]
# df1$V1 <- df1$V1 - 1001000000L
# df1$V2 <- df1$V2 - 1002000000L
df2 <- blast_out[1002000000L < blast_out$V1 & 1002000000L > blast_out$V2, 1:8]
# df2$V1 <- df2$V1 - 1002000000L
# df2$V2 <- df2$V2 - 1001000000L
names(df1) <- c("qseqid", "sseqid", "pident", "qcovs_q2s", 
                "qstart", "qend", "sstart", "send")
names(df2) <- c("sseqid", "qseqid", "pident", "qcovs_s2q", 
                "sstart", "send", "qstart", "qend")
rbh <- inner_join(df1, df2, 
                  by = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "pident"),
                  relationship = "many-to-many")

rbh <- fread("~/workspace/orthology/dev_osat_ogla/rbh/1001_1002.rbh", sep = "\t")
# rbh$V1 <- rbh$V1 - 1001000000L
# rbh$V2 <- rbh$V2 - 1002000000L
# names(rbh) <- c("qseqid", "sseqid", "pident", "qcovs_q2s", "qcovs_s2q", 
#                 "ci_q2s", "ci_s2q", "mutual_ci")
# rbh$pair_id <- paste(rbh$qseqid, rbh$sseqid, sep = "_")
rbh <- .orgBLASTout(rbh = rbh)
rbh$ci_q2s <- rbh$pident * rbh$qcovs_q2s * 1e-4
rbh$ci_s2q <- rbh$pident * rbh$qcovs_s2q * 1e-4
rbh$mutual_ci <- rbh$ci_q2s * rbh$ci_s2q
rbh <- rbh[order(-rbh$mutual_ci), ]
fwrite(rbh, "~/workspace/orthology/dev_osat_ogla/rbh/1001_1002.rbh", sep = "\t", row.names = F, col.names = F)

old_df1 <- fread("/home/ftom/workspace/orthology/dev_osat_ogla/test/old_blast_out1.txt", sep = "\t")
old_df2 <- fread("/home/ftom/workspace/orthology/dev_osat_ogla/test/old_blast_out2.txt", sep = "\t")
names(old_df1) <- c("qseqid", "sseqid", "pident", "qcovs_q2s", 
                    "qstart", "qend", "sstart", "send")
names(old_df2) <- c("sseqid", "qseqid", "pident", "qcovs_s2q", 
                    "sstart", "send", "qstart", "qend")
old_rbh <- inner_join(old_df1, old_df2, 
                      by = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "pident"),
                      relationship = "many-to-many")
old_rbh$pair_id <- paste(old_rbh$qseqid, old_rbh$sseqid, sep = "_")
old_rbh <- .orgBLASTout(rbh = old_rbh)
old_rbh$ci_q2s <- old_rbh$pident * old_rbh$qcovs_q2s * 1e-4
old_rbh$ci_s2q <- old_rbh$pident * old_rbh$qcovs_s2q * 1e-4
old_rbh$mutual_ci <- old_rbh$ci_q2s * old_rbh$ci_s2q
old_rbh <- old_rbh[order(-old_rbh$mutual_ci), ]
old_rbh$qseqid <- old_rbh$qseqid + 1001000000
old_rbh$sseqid <- old_rbh$sseqid + 1002000000
fwrite(old_rbh, "~/workspace/orthology/dev_osat_ogla/rbh/1001_1002.rbh", sep = "\t", row.names = F, col.names = F)

gff_list <- list(query_gff = readRDS("~/workspace/orthology/output/benchmark/orthopair/input/Osat/gff_df.rds"),
                 subject_gff = readRDS("~/workspace/orthology/output/benchmark/orthopair/input/Ogla/gff_df.rds"))
q_tx_i <- gff_list$query_gff$type %in% c("transcript", "mRNA")
query_df <- gff_list$query_gff[q_tx_i, ]
q_cds_i <- gff_list$query_gff$type %in% "CDS"
query_gff <- gff_list$query_gff[q_cds_i, ]
hit <- match(unlist((query_gff$Parent)), gff_list$query_gff$ID)
query_gff$tx_index <- gff_list$query_gff$tx_index[hit]
s_tx_i <- gff_list$subject_gff$type %in% c("transcript", "mRNA")
subject_df <- gff_list$subject_gff[s_tx_i, ]
s_cds_i <- gff_list$subject_gff$type %in% "CDS"
subject_gff <- gff_list$subject_gff[s_cds_i, ]
hit <- match(unlist((subject_gff$Parent)), gff_list$subject_gff$ID)
subject_gff$tx_index <- gff_list$subject_gff$tx_index[hit]
g2g_graph <- list(query_df = query_df,
                  subject_df = subject_df,
                  query_gff = query_gff,
                  subject_gff = subject_gff)

hit <- match(rbh$qseqid, gff_list$query_gff$tx_index)
rbh$qgeneid <- gff_list$query_gff$gene_index[hit]
hit <- match(rbh$sseqid, gff_list$subject_gff$tx_index)
rbh$sgeneid <- gff_list$subject_gff$gene_index[hit]
rbh$pair_id <- paste(rbh$qgeneid, rbh$sgeneid, sep = "0")
rbh <- subset(rbh, subset = !duplicated(pair_id))
rbh <- lapply(rbh, as.numeric)
rbh <- as.data.frame(rbh)
anchor <- .findAnchors(rbh = rbh, g2g_graph = g2g_graph)
t2a_graph <- .link2Anchor(anchor = anchor, g2g_graph = g2g_graph)
orthopair <- .findSyntenicOrtho(rbh = rbh,
                                anchor = anchor,
                                g2g_graph = g2g_graph, 
                                t2a_graph = t2a_graph)
orthopair <- .findSyntenyBlocks(orthopair = orthopair)
orthopair <- .pickBestPair(orthopair = orthopair)
orthopair <- .classifyOrthoPair(orthopair = orthopair)
orthopair <- .filterOrthopair(orthopair = orthopair, g2g_graph = g2g_graph)
orthopair <- .classifyOrthoPair(orthopair = orthopair)
orthopair <- .reformatOrthoPair(orthopair = orthopair, 
                                g2g_graph = g2g_graph)
collapsed_id <- .collapseOverlappingGene(orthopair = orthopair, 
                                         g2g_graph = g2g_graph)
orthopair$query_collapse <- collapsed_id$query_collapse
orthopair$subject_collapse <- collapsed_id$subject_collapse


hit <- match(old_rbh$qseqid, gff_list$query_gff$tx_index)
old_rbh$qgeneid <- gff_list$query_gff$gene_index[hit]
hit <- match(old_rbh$sseqid, gff_list$subject_gff$tx_index)
old_rbh$sgeneid <- gff_list$subject_gff$gene_index[hit]
old_rbh$pair_id <- paste(old_rbh$qgeneid, old_rbh$sgeneid, sep = "0")
old_rbh <- subset(old_rbh, subset = !duplicated(pair_id))
old_rbh <- lapply(old_rbh, as.numeric)
old_rbh <- as.data.frame(old_rbh)
old_anchor <- .findAnchors(rbh = old_rbh, g2g_graph = g2g_graph)
old_t2a_graph <- .link2Anchor(anchor = old_anchor, g2g_graph = g2g_graph)
old_orthopair <- .findSyntenicOrtho(rbh = old_rbh,
                                    anchor = old_anchor,
                                    g2g_graph = g2g_graph, 
                                    t2a_graph = old_t2a_graph)
old_orthopair <- .findSyntenyBlocks(orthopair = old_orthopair)
old_orthopair <- .pickBestPair(orthopair = old_orthopair)
old_orthopair <- .classifyOrthoPair(orthopair = old_orthopair)
old_orthopair <- .filterOrthopair(orthopair = old_orthopair, 
                                  g2g_graph = g2g_graph)
old_orthopair <- .classifyOrthoPair(orthopair = old_orthopair)
old_orthopair <- .reformatOrthoPair(orthopair = old_orthopair, 
                                    g2g_graph = g2g_graph)
collapsed_id <- .collapseOverlappingGene(orthopair = old_orthopair, 
                                         g2g_graph = g2g_graph)
old_orthopair$query_collapse <- collapsed_id$query_collapse
old_orthopair$subject_collapse <- collapsed_id$subject_collapse


og <- subset(out[[1]], select = c("original_query_gene", "original_subject_gene"))
names(og)[1:2] <- c("query", "subject")
eval_busco_list(og = og,
                busco_list = nb_wk21_busco_list,
                del = del,
                query = "query",
                subject = "subject")

old_og <- subset(old_orthopair, select = c("original_query_gene", "original_subject_gene"))
names(old_og)[1:2] <- c("query", "subject")
eval_busco_list(og = old_og,
                busco_list = nb_wk21_busco_list,
                del = del,
                query = "query",
                subject = "subject")
nb_wk21_busco_list$pair_id <- paste(nb_wk21_busco_list$query, nb_wk21_busco_list$subject, sep = "_")
old_og$pair_id <- paste(old_og$query, old_og$subject, sep = "_")
og$pair_id <- paste(og$query, og$subject, sep = "_")
old_og_hit <- old_og$pair_id[old_og$pair_id %in% nb_wk21_busco_list$pair_id]
og_hit <- og$pair_id[og$pair_id %in% nb_wk21_busco_list$pair_id]
old_og_hit <- old_og_hit[!old_og_hit %in% og_hit]
[1] "Os01g0231600_Ogla-WK21_01G0154200" "Os07g0168000_Ogla-WK21_07G0093600"
[3] "Os11g0579800_Ogla-WK21_11G0499600" "Os04g0577300_Ogla-WK21_04G0625800"
[5] "Os09g0347900_Ogla-WK21_09G0246200"

og_missing <- nb_wk21_busco_list$pair_id[!nb_wk21_busco_list$pair_id %in% og$pair_id]

out <- out[[1]]
out[out$original_query_gene == "Os01g0231600", ]
out[out$original_query_gene == "Os11g0579800", ]
out[out$original_query_gene == "Os09g0347900", ]
out[out$original_query_gene == "Os07g0168000", ]
out[out$original_query_gene == "Os04g0577300", ]
out[out$original_subject_gene == "Ogla-WK21_01G0154200", ]
out[out$original_subject_gene == "Ogla-WK21_11G0499600", ]
out[out$original_subject_gene == "Ogla-WK21_09G0246200", ]
out[out$original_subject_gene == "Ogla-WK21_07G0093600", ]
out[out$original_subject_gene == "Ogla-WK21_04G0625800", ]
rbh <- fread("/home/ftom/workspace/orthology/dev_osat_ogla/rbh/1001_1002.rbh", sep = "\t")
qg <- readRDS("/home/ftom/workspace/orthology/dev_osat_ogla/input/1_Osat/gff_df.rds")
hit <- match(rbh$V1, qg[[1]]$tx_index)
rbh$qgeneid <- qg[[1]]$gene_id[hit]
sg <- readRDS("/home/ftom/workspace/orthology/dev_osat_ogla/input/2_Ogla/gff_df.rds")
hit <- match(rbh$V2, sg[[1]]$tx_index)
rbh$sgeneid <- sg[[1]]$gene_id[hit]
rbh$pair_id <- paste(rbh$qgeneid, rbh$sgeneid, sep = "_")

rbh[rbh$qgeneid == "Os01g0231600" & rbh$sgeneid == "Ogla-WK21_01G0154200", ]
rbh[rbh$qgeneid == "Os11g0579800" & rbh$sgeneid == "Ogla-WK21_11G0499600", ]
rbh[rbh$qgeneid == "Os09g0347900" & rbh$sgeneid == "Ogla-WK21_09G0246200", ]
rbh[rbh$qgeneid == "Os07g0168000" & rbh$sgeneid == "Ogla-WK21_07G0093600", ]
rbh[rbh$qgeneid == "Os04g0577300" & rbh$sgeneid == "Ogla-WK21_04G0625800", ]
rbh[rbh$pair_id %in% og_missing, ]

hit <- match(old_rbh$qseqid, qg[[1]]$tx_index)
old_rbh$qgeneid <- qg[[1]]$gene_id[hit]
hit <- match(old_rbh$sseqid, sg[[1]]$tx_index)
old_rbh$sgeneid <- sg[[1]]$gene_id[hit]
old_rbh[old_rbh$qgeneid == "Os01g0231600" & old_rbh$sgeneid == "Ogla-WK21_01G0154200", ]
old_rbh[old_rbh$qgeneid == "Os11g0579800" & old_rbh$sgeneid == "Ogla-WK21_11G0499600", ]
old_rbh[old_rbh$qgeneid == "Os09g0347900" & old_rbh$sgeneid == "Ogla-WK21_09G0246200", ]
old_rbh[old_rbh$qgeneid == "Os07g0168000" & old_rbh$sgeneid == "Ogla-WK21_07G0093600", ]
old_rbh[old_rbh$qgeneid == "Os04g0577300" & old_rbh$sgeneid == "Ogla-WK21_04G0625800", ]

g2g_graph <- list(query_df = qg[[1]],
                  subject_df = sg[[1]],
                  query_gff = qg[[1]],
                  subject_gff = sg[[1]])
names(rbh) <- names(old_rbh)
hit <- match(rbh$qseqid, g2g_graph$query_gff$tx_index)
rbh$qgeneid <- g2g_graph$query_gff$gene_index[hit]
hit <- match(rbh$sseqid, g2g_graph$subject_gff$tx_index)
rbh$sgeneid <- g2g_graph$subject_gff$gene_index[hit]
rbh$pair_id <- paste(rbh$qgeneid, rbh$sgeneid, sep = "0")
rbh <- subset(rbh, subset = !duplicated(pair_id))
rbh <- lapply(rbh, as.numeric)
rbh <- as.data.frame(rbh)
anchor <- .findAnchors(rbh = rbh, g2g_graph = g2g_graph)
t2a_graph <- .link2Anchor(anchor = anchor, g2g_graph = g2g_graph)
orthopair <- .findSyntenicOrtho(rbh = rbh,
                                anchor = anchor,
                                g2g_graph = g2g_graph, 
                                t2a_graph = t2a_graph)

hit <- match(old_rbh$qseqid, g2g_graph$query_gff$tx_index)
old_rbh$qgeneid <- g2g_graph$query_gff$gene_index[hit]
hit <- match(old_rbh$sseqid, g2g_graph$subject_gff$tx_index)
old_rbh$sgeneid <- g2g_graph$subject_gff$gene_index[hit]
old_rbh$pair_id <- paste(old_rbh$qgeneid, old_rbh$sgeneid, sep = "0")
old_rbh <- subset(old_rbh, subset = !duplicated(pair_id))
old_rbh <- lapply(old_rbh, as.numeric)
old_rbh <- as.data.frame(old_rbh)
old_anchor <- .findAnchors(rbh = old_rbh, g2g_graph = g2g_graph)
t2a_graph <- .link2Anchor(anchor = old_anchor, g2g_graph = g2g_graph)
old_orthopair <- .findSyntenicOrtho(rbh = old_rbh,
                                    anchor = old_anchor,
                                    g2g_graph = g2g_graph, 
                                    t2a_graph = t2a_graph)
old_orthopair <- .findSyntenyBlocks(orthopair = old_orthopair)
old_orthopair <- .pickBestPair(orthopair = old_orthopair)
old_orthopair <- .classifyOrthoPair(orthopair = old_orthopair)
old_orthopair <- .filterOrthopair(orthopair = old_orthopair, 
                                  g2g_graph = g2g_graph)
old_orthopair <- .classifyOrthoPair(orthopair = old_orthopair)
old_orthopair <- .reformatOrthoPair(orthopair = old_orthopair, 
                                    g2g_graph = g2g_graph)