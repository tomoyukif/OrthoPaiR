#' Execute Reciprocal BLAST Hits (RBH) search
#'
#' This function performs a Reciprocal BLAST Hits (RBH) search between query and subject genomes.
#'
#' @param object A OrthoPairDB object.
#' @param db1 Path to the BLAST database for query genome (default is NULL).
#' @param db2 Path to the BLAST database for subject genome (default is NULL).
#' @param n_threads Number of threads to use for BLAST (default is 1).
#' @export
#' @importFrom Biostrings readDNAStringSet
#' @importFrom rhdf5 h5read
rbh <- function(object,
                n_threads = 1,
                use_prot = FALSE,
                blast_path = "",
                diamond_path = ""
){
    # Check if the input object is of class "OrthoPairDB"
    stopifnot(inherits(x = object, "OrthoPairDB"))
    
    # Open the HDF5 file
    # h5 <- H5Fopen(object$h5)
    # Ensure the HDF5 file is closed when the function exits
    # on.exit(H5Fclose(h5))
    files <- h5read(object$h5, "files")
    
    fa1_fn <- files$query_cds
    fa2_fn <- files$subject_cds
    
    # Set scientific notation options
    options(scipen = 9999)
    
    # Perform BLAST searches in both directions
    message("BLAST: query to subject.")
    blast_out1 <- .blast_search(fa = fa1_fn,
                                db = fa2_fn,
                                blast_path = blast_path,
                                n_threads = n_threads)
    
    message("BLAST: subject to query.")
    blast_out2 <- .blast_search(fa = fa2_fn,
                                db = fa1_fn,
                                blast_path = blast_path,
                                n_threads = n_threads)
    
    # Filter BLAST results for best hits
    message("Generating RBH lists.")
    rbh_out <- .getRBH(df1 = blast_out1, df2 = blast_out2)
    # rbh_out <- .blast_filter(rbh_out = rbh_out)
    
    if(use_prot){
        # Generate FASTA files from the CDS sequences
        check1 <- "no_query_prot" %in% files$query_prot
        check2 <- "no_subject_prot" %in% files$subject_prot
        if(!check1 & !check2){
            fa1_fn <- files$query_prot
            fa2_fn <- files$subject_prot
            
            # Perform BLAST searches in both directions
            message("DIAMOND: query to subject.")
            blast_out1 <- .diamond_search(fa = fa1_fn,
                                          db = fa2_fn,
                                          n_threads = n_threads,
                                          diamond_path = diamond_path)
            
            message("DIAMOND: subject to query.")
            blast_out2 <- .diamond_search(fa = fa2_fn,
                                          db = fa1_fn,
                                          n_threads = n_threads,
                                          diamond_path = diamond_path)
            
            # Filter BLAST results for best hits
            message("Update RBH lists.")
            prot_rbh_out <- .getRBH(df1 = blast_out1, df2 = blast_out2)
            rbh_out <- rbind(prot_rbh_out, rbh_out)
            rbh_out <- rbh_out[order(rbh_out$mutual_ci, decreasing = TRUE), ]
            rbh_out <- rbh_out[!duplicated(subset(rbh_out, select = qseqid:sseqid)), ]
        }
    }
    
    rbh_out <- .reorganizeRBH(object = object, rbh_out = rbh_out)
    
    # Save RBH results to HDF5 file
    .h5creategroup(object$h5,"blast")
    .h5overwrite(obj = rbh_out, file = object$h5, "blast/rbh")
    .h5overwrite(obj = as.character(Sys.time()), file = object$h5, "timestamp/blast")
    options(scipen = 0)
}

.getRepTx <- function(gff_cds){
    gff_cds$Parent <- unlist(gff_cds$Parent)
    exon_edges <- paste(start(gff_cds), end(gff_cds), sep = "-")
    exon_str <- tapply(exon_edges, gff_cds$Parent, paste, collapse = "_")
    exon_chr <- tapply(as.character(seqnames(gff_cds)), gff_cds$Parent, unique)
    exon_str <- paste(exon_chr, exon_str, sep = "_")
    rep_tx <- names(exon_chr)[!duplicated(exon_str)]
    return(rep_tx)
}

.getSingleExonTx <- function(gff_cds){
    gff_cds$Parent <- unlist(gff_cds$Parent)
    parent_table <- table(gff_cds$Parent)
    not_single <- names(parent_table[parent_table > 1])
    not_single_gene <- gff_cds$gene_id[gff_cds$Parent %in% not_single]
    single_gene_tx <- gff_cds$Parent[!gff_cds$gene_id %in% not_single_gene]
    single_gene_tx <- unique(single_gene_tx)
    return(single_gene_tx)
}

.reorganizeRBH <- function(object, rbh_out){
    gff_list <- .getGFFlist(object = object)
    hit <- match(rbh_out$qseqid, gff_list$query_gff$tx_index)
    rbh_out$qgeneid <- gff_list$query_gff$gene_index[hit]
    hit <- match(rbh_out$sseqid, gff_list$subject_gff$tx_index)
    rbh_out$sgeneid <- gff_list$subject_gff$gene_index[hit]
    rbh_out$pair_id <- paste(rbh_out$qgeneid, rbh_out$sgeneid, sep = "0")
    rbh_out <- subset(rbh_out, subset = !duplicated(pair_id))
    rbh_out <- lapply(rbh_out, as.numeric)
    rbh_out <- as.data.frame(rbh_out)
    return(rbh_out)
}

#' Merge CDS files into a single FASTA file
#'
#' This function reads multiple CDS files and merges them into a single FASTA file.
#'
#' @param cds_fn A character vector containing paths to CDS files.
#' @return The path to the merged FASTA file.
#' @importFrom Biostrings readDNAStringSet writeXStringSet
.makeFASTA <- function(fasta_fn, type = "cds"){
    # Check if there is more than one CDS file
    if(type == "cds"){
        if(length(fasta_fn) != 1){
            # Initialize an empty DNAStringSet object
            for(i in seq_along(fasta_fn)){
                if(i == 1){
                    # Read the first CDS file
                    out <- readDNAStringSet(fasta_fn[i])
                } else {
                    # Append the subsequent CDS files
                    out <- c(out, readDNAStringSet(fasta_fn[i]))
                }
            }
            # Create a new filename for the merged CDS file
            fasta_fn <- sub("_cds.fa", "_merge_cds.fa", fasta_fn[i])
            # Write the merged CDS sequences to a new FASTA file
            writeXStringSet(out, fasta_fn)
        }
        
    } else {
        if(length(fasta_fn) != 1){
            # Initialize an empty DNAStringSet object
            for(i in seq_along(fasta_fn)){
                if(i == 1){
                    # Read the first CDS file
                    out <- readAAStringSet(fasta_fn[i])
                } else {
                    # Append the subsequent CDS files
                    out <- c(out, readAAStringSet(fasta_fn[i]))
                }
            }
            # Create a new filename for the merged CDS file
            fasta_fn <- sub("_prot.fa", "_merge_prot.fa", fasta_fn[i])
            # Write the merged CDS sequences to a new FASTA file
            writeXStringSet(out, fasta_fn)
        }
    }
    # Return the path to the (merged) CDS file
    return(fasta_fn)
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
#' @param task BLAST task to perform (default is "-task blastn").
#' 
#' @importFrom GenomicRanges reduce GRanges 
#' @importFrom IRanges IRanges width
#' @importFrom GenomeInfoDb seqnames
#' 
#' @return A data.frame containing the BLAST search results.
#'
.blast_search <- function(fa,
                          db,
                          blast_path,
                          n_threads){
    # Construct the BLAST arguments string
    blast_args <- paste(paste("-query", fa),
                        paste("-db", db),
                        "-task blastn -max_target_seqs 100000",
                        "-evalue 1e-4 -strand plus",
                        "-outfmt '6 qseqid sseqid pident qcovs qlen qstart qend sstart send'",
                        paste("-num_threads", n_threads))
    
    check <- try({
        out <- system2(command = file.path(blast_path, "blastn"),
                       args = blast_args, 
                       stdout = TRUE)
    }, silent = TRUE)
    
    if(inherits(check, "try-error")){
        stop("Error in blastn of BLAST.\n", check)
    }
    
    if(length(out) == 0){
        out <- NA
        
    } else {
        out <- do.call("rbind", strsplit(out, split = "\t", fixed = TRUE))
        out <- as.data.frame(out)
        names(out) <- c("qseqid", "sseqid", "pident", "qcovs",
                        "qlen", "qstart", "qend", "sstart", "send")
        out <- subset(out, 
                      subset = !is.na(qseqid) & !is.na(sseqid))
    }
    return(out)
}

#' @importFrom rdiamond diamond_protein_to_protein
.diamond_search <- function(fa, db, n_threads, diamond_path){
    
    # Construct the DIAMOND arguments string
    diamond_args <- paste("blastp",
                          paste("--db", db),
                          paste("--query", fa),
                          "--max-target-seqs 100000 --evalue 1e-4",
                          "--strand plus",
                          "--outfmt 6 qseqid sseqid pident qcovhsp qlen qstart qend sstart send",
                          "--quiet --ultra-sensitive", 
                          paste("--threads", n_threads))
    
    check <- try({
        out <- system2(command = file.path(diamond_path, "diamond"),
                       args = diamond_args, 
                       stdout = TRUE)
    }, silent = TRUE)
    
    if(inherits(check, "try-error")){
        stop("Error in blast of BLAST.\n", check)
    }
    if(length(out) == 0){
        out <- NA
        
    } else {
        out <- do.call("rbind", strsplit(out, split = "\t", fixed = TRUE))
        out <- as.data.frame(out)
        names(out) <- c("qseqid", "sseqid", "pident", "qcovs",
                        "qlen", "qstart", "qend", "sstart", "send")
        out <- subset(out, 
                      subset = !is.na(qseqid) & !is.na(sseqid))
    }
    return(out)
}

#' Filter BLAST search results
#'
#' This function filters BLAST search results based on specified criteria for percentage identity,
#' query coverage, and e-value. It can also identify the best hits for each query sequence.
#'
#' @param blast_out A data.frame containing the BLAST search results.
#'
#' @return A filtered data.frame of BLAST search results.
#'
.blast_filter <- function(rbh_out){
    if(all(is.na(rbh_out[1]))){
        return(rbh_out)
    }
    # Identify the best hits for each query sequence
    q_n_hit <- table(rbh_out$qseqid)
    s_n_hit <- table(rbh_out$sseqid)
    
    # Separate single-hit and multiple-hit sequences
    q_single_hit <- rbh_out$qseqid %in% names(q_n_hit)[q_n_hit == 1]
    s_single_hit <- rbh_out$sseqid %in% names(s_n_hit)[s_n_hit == 1]
    single_hit <- q_single_hit & s_single_hit
    mult_hit <- rbh_out[!single_hit, ]
    single_hit <- rbh_out[single_hit, ]
    
    # Filter multiple-hit sequences by the lowest e-value
    filter <- tapply(mult_hit$mutual_ci, mult_hit$qseqid, max)
    hit <- match(mult_hit$qseqid, names(filter))
    filter <- filter[hit]
    mult_hit <- subset(mult_hit, subset = mutual_ci == filter)
    
    # Combine single-hit and multiple-hit sequences
    blast_out <- rbind(single_hit, mult_hit)
    return(blast_out)
}

#' Find Reciprocal Best Hits (RBH)
#'
#' This function identifies Reciprocal Best Hits (RBH) between two sets of BLAST search results.
#'
#' @param df1 A data.frame containing BLAST search results from query to subject.
#' @param df2 A data.frame containing BLAST search results from subject to query.
#' @return A data.frame containing the reciprocal best hits with relevant BLAST statistics.
#' 
#' @importFrom dplyr inner_join
#' 
.getRBH <- function(df1, df2){
    if(all(is.na(df1)) | all(is.na(df2))){
        return(NA)
    }
    names(df1) <- c("qseqid", "sseqid", "pident", "qcovs_q2s", "qlen", 
                    "qstart", "qend", "sstart", "send")
    names(df2) <- c("sseqid", "qseqid", "pident", "qcovs_s2q", "slen", 
                    "sstart", "send", "qstart", "qend")
    rbh <- inner_join(df1, df2, 
                      c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "pident"),
                      relationship = "many-to-many")
    rbh <- .orgBLASTout(rbh = rbh)
    rbh$ci_q2s <- rbh$pident * rbh$qcovs_q2s * 1e-4
    rbh$ci_s2q <- rbh$pident * rbh$qcovs_s2q * 1e-4
    rbh$mutual_ci <- rbh$ci_q2s * rbh$ci_s2q
    rbh <- rbh[order(rbh$mutual_ci, decreasing = TRUE), ]
    return(rbh)
}


.orgBLASTout <- function(rbh){
    rbh <- lapply(rbh, as.numeric)
    rbh <- as.data.frame(rbh)
    rbh$id <- paste(rbh$qseqid, rbh$sseqid, sep = "_")
    id_dup <- duplicated(rbh$id)
    dup_id <- unique(rbh$id[id_dup])
    hit <- !rbh$id %in% dup_id
    uniq_rbh <- rbh[hit, ]
    dup_rbh <- rbh[!hit, ]
    q_aln <- GRanges(seqnames = dup_rbh$id,
                     ranges = IRanges(start = dup_rbh$qstart,
                                      end = dup_rbh$qend))
    q_ol <- findOverlaps(q_aln, q_aln, type = "within")
    q_ol <- as.matrix(q_ol)
    q_ol <- subset(q_ol, subset = q_ol[, 1] != q_ol[, 2])
    s_aln <- GRanges(seqnames = dup_rbh$id,
                     ranges = IRanges(start = dup_rbh$sstart,
                                      end = dup_rbh$send))
    s_ol <- findOverlaps(s_aln, s_aln, type = "within")
    s_ol <- as.matrix(s_ol)
    s_ol <- subset(s_ol, subset = s_ol[, 1] != s_ol[, 2])
    ol <- as.vector(rbind(q_ol, s_ol))
    ol <- sort(unique(ol))
    dup_rbh_ol <- dup_rbh[ol, ]
    first_hit <- tapply(seq_along(dup_rbh_ol$id), dup_rbh_ol$id, min)
    dup_rbh_ol_first_hit <- dup_rbh_ol[first_hit, ]
    dup_rbh <- rbind(dup_rbh[-ol, ], dup_rbh_ol_first_hit)
    
    aln <- GRanges(seqnames = dup_rbh$id,
                   ranges = IRanges(start = dup_rbh$qstart,
                                    end = dup_rbh$qend))
    n_aln <- tapply(width(aln), as.character(seqnames(aln)), sum)
    n_ident <- dup_rbh$pident * 1e-2 * width(aln)
    n_ident <- tapply(n_ident, dup_rbh$id, sum)
    dup_rbh <- subset(dup_rbh, subset = !duplicated(id))
    hit <- match(dup_rbh$id, names(n_ident))
    dup_rbh$pident <- n_ident[hit] / n_aln[hit] * 100
    rbh <- rbind(subset(dup_rbh, select = c(qseqid:qcovs_q2s, qcovs_s2q)),
                 subset(uniq_rbh, select = c(qseqid:qcovs_q2s, qcovs_s2q)))
    return(rbh)
}

# .getMCI <- function(rbh){
#     # Calculate coverage identity for query to subject (q2s) and subject to query (s2q)
#     ci_q2s <- rbh$pident_q2s * rbh$qcovs_q2s * 1e-4
#     ci_s2q <- rbh$pident_s2q * rbh$qcovs_s2q * 1e-4
#     
#     # Calculate mutual coverage identity as the Euclidean distance of the individual coverage identities
#     mutual_ci <- sqrt(ci_q2s^2 + ci_s2q^2) / sqrt(2)
#     return(list(mutual_ci = mutual_ci, ci_q2s = ci_q2s, ci_s2q = ci_s2q))
# }