#' Execute Reciprocal BLAST Hits (RBH) search
#'
#' This function performs a Reciprocal BLAST Hits (RBH) search between query and subject genomes.
#'
#' @param object A OrthoPairDB object.
#' @param db1 Path to the BLAST database for query genome (default is NULL).
#' @param db2 Path to the BLAST database for subject genome (default is NULL).
#' @param n_threads Number of threads to use for BLAST (default is 1).
#' @param n_batch Number of sequences to process in each batch (default is NULL).
#' @param makedb Logical indicating whether to create BLAST databases (default is TRUE).
#' @param max_target_seqs Maximum number of target sequences (default is 100000).
#' @param pident Minimum percentage identity for BLAST hits (default is 90).
#' @param qcovs Minimum query coverage for BLAST hits (default is 0).
#' @param evalue E-value threshold for BLAST (default is 1e-10).
#' @export
#' @importFrom rBLAST makeblastdb blast
#' @importFrom Biostrings readDNAStringSet
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
    # Check if the input object is of class "OrthoPairDB"
    stopifnot(inherits(x = object, "OrthoPairDB"))

    # Open the HDF5 file
    h5 <- H5Fopen(object$h5)
    # Ensure the HDF5 file is closed when the function exits
    on.exit(H5Fclose(h5))

    # Generate FASTA files from the CDS sequences
    fa1 <- .makeFASTA(cds_fn = as.vector(h5$files$query_cds))
    fa2 <- .makeFASTA(cds_fn = as.vector(h5$files$subject_cds))

    # Create BLAST databases if necessary
    if(makedb | is.null(db1) | is.null(db2)){
        makeblastdb(file = fa1, dbtype = "nucl")
        makeblastdb(file = fa2, dbtype = "nucl")
        db1 <- fa1
        db2 <- fa2
    }

    # Read the FASTA files and create BLAST database connections
    fa1 <- readDNAStringSet(fa1)
    fa2 <- readDNAStringSet(fa2)
    db1 <- blast(db = db1, type = "blastn")
    db2 <- blast(db = db2, type = "blastn")

    # Set scientific notation options
    options(scipen = 10^6)

    # Perform BLAST searches in both directions
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

    # Create HDF5 groups and save BLAST outputs
    .h5creategroup(object$h5,"blast")
    .h5overwrite(obj = blast_out1, file = object$h5, "blast/blast_q2s")
    .h5overwrite(obj = blast_out2, file = object$h5, "blast/blast_s2q")

    # Filter BLAST results for best hits
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

    # Filter BLAST results based on user-defined criteria
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

    # Find Reciprocal Best BLAST Hits (RBBH) and Reciprocal BLAST Hits (RBH)
    rbbh_out <- .find_reciprocal(df1 = blast_best1, df2 = blast_best2)
    rbh_out <- .find_reciprocal(df1 = blast_out1, df2 = blast_out2)
    options(scipen = 0)

    # Save RBBH and RBH results to HDF5 file
    .h5overwrite(obj = rbbh_out, file = object$h5, "blast/rbbh")
    .h5overwrite(obj = rbh_out, file = object$h5, "blast/rbh")
    .h5overwrite(obj = as.character(Sys.time()), file = object$h5, "timestamp/blast")
}

#' Merge CDS files into a single FASTA file
#'
#' This function reads multiple CDS files and merges them into a single FASTA file.
#'
#' @param cds_fn A character vector containing paths to CDS files.
#' @return The path to the merged FASTA file.
#' @importFrom Biostrings readDNAStringSet writeXStringSet
.makeFASTA <- function(cds_fn){
    # Check if there is more than one CDS file
    if(length(cds_fn) != 1){
        # Initialize an empty DNAStringSet object
        for(i in seq_along(cds_fn)){
            if(i == 1){
                # Read the first CDS file
                cds <- readDNAStringSet(cds_fn[i])
            } else {
                # Append the subsequent CDS files
                cds <- c(cds, readDNAStringSet(cds_fn[i]))
            }
        }
        # Create a new filename for the merged CDS file
        cds_fn <- sub("_cds.fa", "_merge_cds.fa", cds_fn[i])
        # Write the merged CDS sequences to a new FASTA file
        writeXStringSet(cds, cds_fn)
    }
    # Return the path to the (merged) CDS file
    return(cds_fn)
}

#' Perform BLAST search
#'
#' This function performs a BLAST search using the specified parameters.
#'
#' @param fa A DNAStringSet object containing the query sequences.
#' @param db A BLAST database connection.
#' @param n_threads Number of threads to use for BLAST (default is 1).
#' @param n_batch Number of sequences to process in each batch (default is NULL).
#' @param max_target_seqs Maximum number of target sequences (default is 100000).
#' @param evalue E-value threshold for BLAST (default is 1e-10).
#' @param task BLAST task to perform (default is "-task blastn").
#'
#' @return A data.frame containing the BLAST search results.
#'
.blast_search <- function(fa,
                          db,
                          n_threads,
                          n_batch,
                          max_target_seqs,
                          evalue,
                          task){
    # Construct the BLAST arguments string
    blast_args <- paste(task,
                        "-best_hit_overhang 0.1",
                        "-best_hit_score_edge 0.1",
                        paste("-max_target_seqs", max_target_seqs),
                        paste("-evalue", evalue),
                        paste("-num_threads", n_threads))

    # Check if batch processing is required
    if(is.null(n_batch)){
        # Perform the BLAST search without batching
        out <- predict(db, fa,
                       silent = TRUE,
                       BLAST_args = blast_args,
                       custom_format = "qseqid sseqid pident evalue qcovs")
    } else {
        # Split the query sequences into batches
        batch <- split(seq_along(fa), cut(seq_along(fa), n_batch))
        out <- NULL
        # Perform the BLAST search for each batch and combine the results
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

#' Filter BLAST search results
#'
#' This function filters BLAST search results based on specified criteria for percentage identity,
#' query coverage, and e-value. It can also identify the best hits for each query sequence.
#'
#' @param blast_out A data.frame containing the BLAST search results.
#' @param pident Minimum percentage identity for BLAST hits (default is 0).
#' @param qcovs Minimum query coverage for BLAST hits (default is 0).
#' @param evalue Maximum e-value for BLAST hits (default is Inf).
#' @param best Logical indicating whether to filter for the best hits (default is FALSE).
#'
#' @return A filtered data.frame of BLAST search results.
#'
.blast_filter <- function(blast_out,
                          pident = 0,
                          qcovs = 0,
                          evalue = Inf,
                          best = FALSE){
    if(best){
        # Identify the best hits for each query sequence
        n_hit <- unlist(tapply(blast_out$qseqid, blast_out$qseqid, length))

        # Separate single-hit and multiple-hit sequences
        single_hit <- blast_out$qseqid %in% names(n_hit)[n_hit == 1]
        single_hit <- blast_out[single_hit, ]
        mult_hit <- blast_out$qseqid %in% names(n_hit)[n_hit > 1]
        mult_hit <- blast_out[mult_hit, ]

        # Filter multiple-hit sequences by the lowest e-value
        filter <- tapply(mult_hit$evalue, mult_hit$qseqid, min)
        hit <- match(mult_hit$qseqid, names(filter))
        filter <- filter[hit]
        mult_hit <- subset(mult_hit, subset = evalue == filter)

        # Further filter by the highest product of percentage identity and query coverage
        qcov_pident <- mult_hit$pident * mult_hit$qcovs * 1e-2
        filter <- tapply(qcov_pident, mult_hit$qseqid, max)
        hit <- match(mult_hit$qseqid, names(filter))
        filter <- filter[hit]
        mult_hit <- subset(mult_hit, subset = qcov_pident == filter)

        # Combine single-hit and multiple-hit sequences
        blast_out <- rbind(single_hit, mult_hit)

    } else {
        # Filter based on percentage identity, query coverage, and e-value
        filter <- blast_out$pident >= pident &
            blast_out$qcovs >= qcovs &
            blast_out$evalue <= evalue
        blast_out <- subset(blast_out, subset = filter)
    }
    return(blast_out)
}

#' Find Reciprocal Best Hits (RBH)
#'
#' This function identifies Reciprocal Best Hits (RBH) between two sets of BLAST search results.
#'
#' @param df1 A data.frame containing BLAST search results from query to subject.
#' @param df2 A data.frame containing BLAST search results from subject to query.
#' @return A data.frame containing the reciprocal best hits with relevant BLAST statistics.
.find_reciprocal <- function(df1, df2){
    # Create unique identifiers for BLAST hits in both directions
    id1 <- paste(df1$qseqid, df1$sseqid, sep = "_")
    id2 <- paste(df2$sseqid, df2$qseqid, sep = "_")

    # Identify reciprocal hits
    rhit <- id1 %in% id2
    out <- subset(df1, subset = rhit, select = c(qseqid, sseqid))
    out <- unique(out)

    # Extract relevant BLAST statistics for reciprocal hits
    out_id <- paste(out$qseqid, out$sseqid, sep = "_")
    out$q2s_pident <- df1$pident[match(out_id, id1)]
    out$q2s_qcovs <- df1$qcovs[match(out_id, id1)]
    out$q2s_evalue <- df1$evalue[match(out_id, id1)]
    out$s2q_pident <- df2$pident[match(out_id, id2)]
    out$s2q_qcovs <- df2$qcovs[match(out_id, id2)]
    out$s2q_evalue <- df2$evalue[match(out_id, id2)]

    return(out)
}
