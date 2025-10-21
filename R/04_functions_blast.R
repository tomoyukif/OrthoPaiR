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
    
    fa1_fn <- h5read(files$query_h5, "cds_fn")
    fa2_fn <- h5read(files$subject_h5, "cds_fn")
    
    # Avoid changing global options; keep defaults for performance
    
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
    
    if(use_prot){
        # Generate FASTA files from the CDS sequences
        query_prot <- h5read(files$query_h5, "prot_fn")
        subject_prot <- h5read(files$subject_h5, "prot_fn")
        check1 <- "no_query_prot" %in% query_prot
        check2 <- "no_subject_prot" %in% subject_prot
        if(!check1 & !check2){
            fa1_fn <- query_prot
            fa2_fn <- subject_prot
            
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
    # Done
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

#' @importFrom GenomicRanges reduce GRanges 
#' @importFrom IRanges IRanges width
#' @importFrom GenomeInfoDb seqnames
#' 
.blast_search <- function(fa,
                          db,
                          blast_path,
                          n_threads){
    # Construct the BLAST arguments string
    blast_args <- paste(paste("-query", fa),
                        paste("-db", db),
                        "-task blastn -max_target_seqs 200",
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
                          "--max-target-seqs 200 --evalue 1e-4",
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

# Compute RBH using best-hit prefiltering to reduce join size drastically
.getRBH <- function(df1, df2){
    if(all(is.na(df1)) || all(is.na(df2))){
        return(NA)
    }
    names(df1) <- c("qseqid", "sseqid", "pident", "qcovs_q2s", "qlen", 
                    "qstart", "qend", "sstart", "send")
    names(df2) <- c("sseqid", "qseqid", "pident", "qcovs_s2q", "slen", 
                    "sstart", "send", "qstart", "qend")
    # Coerce numeric columns used in scoring
    suppressWarnings({
        df1$pident <- as.numeric(df1$pident)
        df1$qcovs_q2s <- as.numeric(df1$qcovs_q2s)
        df2$pident <- as.numeric(df2$pident)
        df2$qcovs_s2q <- as.numeric(df2$qcovs_s2q)
    })
    # Score = pident * coverage
    df1$score <- df1$pident * df1$qcovs_q2s
    df2$score <- df2$pident * df2$qcovs_s2q
    # Best subject per query (q->s)
    df1 <- df1[order(-df1$score, df1$qseqid), ]
    df1_best <- df1[!duplicated(df1$qseqid), c("qseqid","sseqid","pident","qcovs_q2s")]
    # Best query per subject (s->q)
    df2 <- df2[order(-df2$score, df2$sseqid), ]
    df2_best <- df2[!duplicated(df2$sseqid), c("qseqid","sseqid","pident","qcovs_s2q")]
    # RBH on ids only
    rbh <- merge(df1_best, df2_best, by = c("qseqid","sseqid"), all = FALSE, suffixes = c("_q2s","_s2q"))
    if(nrow(rbh) == 0){
        return(NA)
    }
    # Confidence indices
    rbh$ci_q2s <- rbh$pident_q2s * rbh$qcovs_q2s * 1e-4
    rbh$ci_s2q <- rbh$pident_s2q * rbh$qcovs_s2q * 1e-4
    rbh$mutual_ci <- rbh$ci_q2s * rbh$ci_s2q
    rbh <- rbh[order(-rbh$mutual_ci), ]
    # Rename for downstream compatibility
    names(rbh)[names(rbh) == "pident_q2s"] <- "pident"
    return(rbh)
}


## .orgBLASTout no longer needed with best-hit prefilter; keep minimal passthrough
.orgBLASTout <- function(rbh){
    return(rbh)
}
