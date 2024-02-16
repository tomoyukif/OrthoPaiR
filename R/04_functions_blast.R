#' Define a function to execute BLAST search
#' @export
#' @importFrom rBLAST makeblastdb blast
#' @importFrom Biostrings readDNAStringSet
#'
rbh <- function(object,
                db1 = NULL,
                db2 = NULL,
                n_threads = 1,
                n_batch = NULL,
                makedb = TRUE,
                max_target_seqs = 100000,
                pident = 90,
                qcovs = 0,
                evalue = 1e-10
){
    stopifnot(inherits(x = object, "SynogDB"))

    fa1 <- .makeFASTA(cds_fn = object$query_cds)
    fa2 <- .makeFASTA(cds_fn = object$subject_cds)

    if(makedb | {is.null(db1) | is.null(db2)}){
        makeblastdb(file = fa1, dbtype = "nucl")
        makeblastdb(file = fa2, dbtype = "nucl")
        db1 <- fa1
        db2 <- fa2
    }

    fa1 <- readDNAStringSet(fa1)
    fa2 <- readDNAStringSet(fa2)
    db1 <- blast(db = db1, type = "blastn")
    db2 <- blast(db = db2, type = "blastn")

    options(scipen = 10^6)
    blast_out1 <- .blast_search(fa = fa1,
                                db = db2,
                                n_threads = n_threads,
                                n_batch = n_batch,
                                max_target_seqs = max_target_seqs,
                                evalue = 1e-3,
                                task = "-task blastn")
    blast_out2 <- .blast_search(fa = fa2,
                                db = db1,
                                n_threads = n_threads,
                                n_batch = n_batch,
                                max_target_seqs = max_target_seqs,
                                evalue = 1e-3,
                                task = "-task blastn")

    .h5creategroup(object$h5,"blast")
    .h5overwrite(obj = blast_out1, file = object$h5, "blast/blast_q2s")
    .h5overwrite(obj = blast_out2, file = object$h5, "blast/blast_s2q")

    # RBBH
    blast_best1 <- .blast_filter(blast_out = blast_out1,
                                 pident = 0,
                                 qcovs = 0,
                                 evalue = 1e-3,
                                 best = TRUE)

    blast_best2 <- .blast_filter(blast_out = blast_out2,
                                 pident = 0,
                                 qcovs = 0,
                                 evalue = 1e-3,
                                 best = TRUE)

    # RBH
    blast_out1 <- .blast_filter(blast_out = blast_out1,
                                pident = pident,
                                qcovs = qcovs,
                                evalue = evalue,
                                best = FALSE)

    blast_out2 <- .blast_filter(blast_out = blast_out2,
                                pident = pident,
                                qcovs = qcovs,
                                evalue = evalue,
                                best = FALSE)

    rbbh_out <- .find_reciprocal(df1 = blast_best1, df2 = blast_best2)
    rbh_out <- .find_reciprocal(df1 = blast_out1, df2 = blast_out2)
    options(scipen = 0)
    .h5overwrite(obj = rbbh_out, file = object$h5, "blast/rbbh")
    .h5overwrite(obj = rbh_out, file = object$h5, "blast/rbh")
}

#' @importFrom Biostrings readDNAStringSet writeXStringSet
.makeFASTA <- function(cds_fn){
    if(length(cds_fn) != 1){
        for(i in seq_along(cds_fn)){
            if(i == 1){
                cds <- readDNAStringSet(cds_fn[i])
            } else {
                cds <- c(cds, readDNAStringSet(cds_fn[i]))
            }
        }
        cds_fn <- sub("_cds.fa", "_merge_cds.fa", cds_fn[i])
        writeXStringSet(cds, cds_fn)
    }
    return(cds_fn)
}

.blast_search <- function(fa,
                          db,
                          n_threads,
                          n_batch,
                          max_target_seqs,
                          evalue,
                          task){
    blast_args <- paste(task,
                        "-best_hit_overhang 0.1",
                        "-best_hit_score_edge 0.1",
                        paste("-max_target_seqs",
                              max_target_seqs),
                        paste("-evalue", evalue),
                        paste("-num_threads", n_threads))

    if(is.null(n_batch)){
        out <- predict(db, fa,
                       silent = TRUE,
                       BLAST_args = blast_args,
                       custom_format = "qseqid sseqid pident evalue qcovs")

    } else {
        batch <- split(seq_along(fa), cut(seq_along(fa), n_batch))
        out <- NULL
        for(i in batch){
            tmp <- predict(db, fa[i],
                           silent = TRUE,
                           BLAST_args = blast_args,
                           custom_format = "qseqid sseqid pident evalue qcovs")
            out <- rbind(out, tmp)
        }
    }
    return(out)
}

.blast_filter <- function(blast_out,
                          pident = 0,
                          qcovs = 0,
                          evalue = Inf,
                          best){
    if(best){
        n_hit <- unlist(tapply(blast_out$qseqid, blast_out$qseqid, length))
        single_hit <- blast_out$qseqid %in% names(n_hit)[n_hit == 1]
        single_hit <- blast_out[single_hit, ]
        mult_hit <- blast_out$qseqid %in% names(n_hit)[n_hit > 1]
        mult_hit <- blast_out[mult_hit, ]

        filter <- tapply(mult_hit$evalue, mult_hit$qseqid, min)
        hit <- match(mult_hit$qseqid, names(filter))
        filter <- filter[hit]
        mult_hit <- subset(mult_hit, subset = evalue == filter)

        qcov_pident <- mult_hit$pident * mult_hit$qcovs * 1e-2

        filter <- tapply(qcov_pident, mult_hit$qseqid, max)
        hit <- match(mult_hit$qseqid, names(filter))
        filter <- filter[hit]
        mult_hit <- subset(mult_hit, subset = qcov_pident == filter)

        blast_out <- rbind(single_hit, mult_hit)

    } else {
        filter <- blast_out$pident >= pident &
            blast_out$qcovs >= qcovs &
            blast_out$evalue <= evalue
        blast_out <- subset(blast_out, subset = filter)
    }
    return(blast_out)
}

.find_reciprocal <- function(df1, df2){
    id1 <- paste(df1$qseqid, df1$sseqid, sep = "_")
    id2 <- paste(df2$sseqid, df2$qseqid, sep = "_")
    rhit <- id1 %in% id2
    out <- subset(df1, subset = rhit, select = c(qseqid, sseqid))
    out <- unique(out)

    out_id <- paste(out$qseqid, out$sseqid, sep = "_")
    out$q2s_pident <- df1$pident[match(out_id, id1)]
    out$q2s_qcovs <- df1$qcovs[match(out_id, id1)]
    out$q2s_evalue <- df1$evalue[match(out_id, id1)]
    out$s2q_pident <- df2$pident[match(out_id, id2)]
    out$s2q_qcovs <- df2$qcovs[match(out_id, id2)]
    out$s2q_evalue <- df2$evalue[match(out_id, id2)]

    return(out)
}
