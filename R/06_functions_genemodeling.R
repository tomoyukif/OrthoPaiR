#' Map Proteins Between Genomes
#'
#' This function maps proteins between the query and subject genomes using the specified miniprot binary and number of cores.
#'
#' @param object A OrthoPairDB object.
#' @param out_prefix A character string specifying the prefix for the output files.
#' @param miniprot_bin A character string specifying the path to the miniprot binary.
#' @param n_threads An integer specifying the number of cores to use for the mapping.
#' @param len_diff A numeric value specifying the maximum allowable length difference for proteins. Default is 0.2.
#'
#' @return None. The function performs the mapping and writes the results to the output files specified by `out_prefix`.
#'
#' @export
mapProt <- function(object,
                    out_dir,
                    miniprot_bin = "miniprot",
                    conda = "conda",
                    condaenv = NULL,
                    n_threads,
                    len_diff = 0.2){
    stopifnot(inherits(x = object, "OrthoPairDB"))
    
    # Call the mapping engine function with specified parameters
    .mapEngine(object = object,
               out_dir = out_dir,
               miniprot_bin = miniprot_bin,
               conda = conda,
               condaenv = condaenv,
               n_threads = n_threads)
    
    .createFASTA(object = object, out_dir = out_dir, merge = FALSE)
    
    message("Organize and integrate miniprot predicated gene models...")
    
    .orgMiniprot(object = object, out_dir = out_dir, len_diff = len_diff)
    
    .createFASTA(object = object, out_dir = out_dir, merge = TRUE)
    
    .h5overwrite(obj = as.character(Sys.time()), file = object$h5, "timestamp/miniprot")
}

#' Execute Protein Mapping Between Genomes
#'
#' This function maps proteins between the query and subject genomes using the specified miniprot binary and number of cores. It performs both directions of mapping and organizes the results.
#'
#' @importFrom rtracklayer export.gff3
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#' @importFrom S4Vectors mcols mcols<-
.mapEngine <- function(object, subject_prot, query_prot,
                       out_dir, miniprot_bin, conda,
                       condaenv, n_threads,
                       overlap, len_diff){
    # Open the HDF5 file
    h5 <- H5Fopen(object$h5)
    # Ensure the HDF5 file is closed when the function exits
    on.exit(H5Fclose(h5))
    
    # Check if the input object is of class "OrthoPairDB"
    stopifnot(inherits(x = object, "OrthoPairDB"))
    
    dir.create(path = out_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Define output file paths for both directions of mapping
    s2q_out <- file.path(out_dir, "query_miniprot_out")
    q2s_out <- file.path(out_dir, "subject_miniprot_out")
    
    # Run miniprot mapping from subject to query genome
    .miniprot(query_fn = as.vector(h5$files$subject_prot),
              genome_fn = as.vector(h5$files$query_genome),
              out_prefix = s2q_out, miniprot_bin = miniprot_bin,
              conda = conda, condaenv = condaenv,
              n_threads = n_threads)
    
    # Run miniprot mapping from query to subject genome
    .miniprot(query_fn = as.vector(h5$files$query_prot),
              genome_fn = as.vector(h5$files$subject_genome),
              out_prefix = q2s_out, miniprot_bin = miniprot_bin,
              conda = conda, condaenv = condaenv,
              n_threads = n_threads)
    
    # Create HDF5 group for protein mapping results and save GFF files
    .h5creategroup(object$h5,"miniprot")
    .h5overwrite(obj = paste0(q2s_out, ".gff"),
                 file = object$h5, "miniprot/q2s_gff")
    .h5overwrite(obj = paste0(s2q_out, ".gff"),
                 file = object$h5, "miniprot/s2q_gff")
    #
    #     # Import existing GFF files for query and subject genomes
    #     query_gff <- .importAllGFF(object$query_gff)
    #     query_gff <- c(query_gff, import.gff(paste0(s2q_out, ".gff")))
    #     subject_gff <- .importAllGFF(object$subject_gff)
    #     subject_gff <- c(subject_gff, import.gff(paste0(q2s_out, ".gff")))
    #
    #     m_query_gff <- mcols(query_gff)
    #     hit <- names(m_query_gff) %in% c("source", "type", "score", "phase",
    #                                    "ID", "Name", "gene_id", "Parent", "Target")
    #     mcols(query_gff) <- m_query_gff[, hit]
    #     m_subject_gff <- mcols(subject_gff)
    #     hit <- names(m_subject_gff) %in% c("source", "type", "score", "phase",
    #                                    "ID", "Name", "gene_id", "Parent", "Target")
    #     mcols(subject_gff) <- m_subject_gff[, hit]
    #     query_gff$Name <- query_gff$ID
    #     subject_gff$Name <- subject_gff$ID
    #
    #     # Export the final GFF results to output files
    #     export.gff3(query_gff, paste0(s2q_out, ".gff"))
    #     export.gff3(subject_gff, paste0(q2s_out, ".gff"))
}

# # Organize GFF results from subject-to-query and query-to-subject mappings
# s2q_gff <- .orgGFF(gff1 = import.gff(paste0(s2q_out, ".gff")),
#                    gff2 = query_gff,
#                    gff3 = subject_gff,
#                    prefix = "query_",
#                    overlap = overlap,
#                    len_diff = len_diff)
# q2s_gff <- .orgGFF(gff1 = import.gff(paste0(q2s_out, ".gff")),
#                    gff2 = subject_gff,
#                    gff3 = query_gff,
#                    prefix = "subject_",
#                    overlap = overlap,
#                    len_diff = len_diff)
#
# # Merge organized GFF results with original query and subject GFFs
# s2q_gff <- .mergeGFF(gff1 = s2q_gff, gff2 = query_gff)
# q2s_gff <- .mergeGFF(gff1 = q2s_gff, gff2 = subject_gff)
#
# # Fix and filter GFF results, retaining only necessary columns
# s2q_gff <- .fixGFF(gff = s2q_gff)
# q2s_gff <- .fixGFF(gff = q2s_gff)
# m_s2q_gff <- mcols(s2q_gff)
# hit <- names(m_s2q_gff) %in% c("source", "type", "score", "phase",
#                                "ID", "Name", "gene_id", "Parent", "Target")
# mcols(s2q_gff) <- m_s2q_gff[, hit]
# m_q2s_gff <- mcols(q2s_gff)
# hit <- names(m_q2s_gff) %in% c("source", "type", "score", "phase",
#                                "ID", "Name", "gene_id", "Parent", "Target")
# mcols(q2s_gff) <- m_q2s_gff[, hit]
# s2q_gff$Name <- s2q_gff$ID
# q2s_gff$Name <- q2s_gff$ID

.miniprot <- function(query_fn, genome_fn, out_prefix,
                      miniprot_bin,
                      conda, condaenv,
                      n_threads = 1){
    
    if(!is.null(condaenv)){
        log_fn <- paste0(out_prefix, ".log1")
        .condaExe(conda = conda, env = condaenv, command = miniprot_bin,
                  args = paste("-t", n_threads,
                               "-d", paste0(out_prefix, ".mpi"),
                               genome_fn),
                  log_fn = log_fn)
        
        log_fn <- paste0(out_prefix, ".log2")
        error_log_fn <- paste0(out_prefix, "_error.log2")
        .condaExe(conda = conda, env = condaenv, command = miniprot_bin,
                  args = paste("-t", n_threads, "--gff",
                               paste0(out_prefix, ".mpi"),
                               query_fn, ">", paste0(out_prefix, ".gff")),
                  log_fn = log_fn)
        
    } else {
        log_fn <- paste0(out_prefix, ".log1")
        error_log_fn <- paste0(out_prefix, "_error.log1")
        system2(command = miniprot_bin,
                args = paste("-t", n_threads,
                             "-d", paste0(out_prefix, ".mpi"),
                             genome_fn), 
                stderr = log_fn)
        
        log_fn <- paste0(out_prefix, ".log2")
        error_log_fn <- paste0(out_prefix, "_error.log2")
        system2(command = miniprot_bin,
                args = paste("-t", n_threads,
                             "--gff",
                             paste0(out_prefix, ".mpi"),
                             query_fn, ">", paste0(out_prefix, ".gff")), 
                stderr = log_fn)
    }
}

.orgMiniprot <- function(object, out_dir, len_diff){
    # Open the HDF5 file
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    
    out <- .filterGFF(h5 = h5, len_diff = len_diff)
    
    # Fix and filter GFF results, retaining only necessary columns
    out$query <- .fixGFF(gff = out$query)
    out$subject <- .fixGFF(gff = out$subject)
    
    m_query <- mcols(out$query)
    hit <- names(m_query) %in% c("source", "type", "score", "phase",
                                 "ID", "Name", "gene_id", "Parent", "Target")
    mcols(out$query) <- m_query[, hit]
    m_subject <- mcols( out$subject)
    hit <- names(m_subject) %in% c("source", "type", "score", "phase",
                                   "ID", "Name", "gene_id", "Parent", "Target")
    mcols( out$subject) <- m_subject[, hit]
    out$query$Name <- out$query$ID
    out$subject$Name <-  out$subject$ID
    
    # Export the final GFF results to output files
    query_fn <- file.path(out_dir, "miniprot_merge_query.gff")
    subject_fn <- file.path(out_dir, "miniprot_merge_subject.gff")
    .h5overwrite(obj = query_fn, file = object$h5, name = "files/query_gff")
    .h5overwrite(obj = subject_fn, file = object$h5, name = "files/subject_gff")
    
    export.gff3(out$query, query_fn)
    export.gff3(out$subject, subject_fn)
}

.filterGFF <- function(h5, len_diff){
    # Get the GFF lists
    gff_ls <- .getGFFlist(h5 = h5, mp = TRUE)
    
    # Get the CDS lists
    cds_ls <- .getCDSlist(h5 = h5, mp = TRUE)
    
    gff_ls$mp_query_gff <- .filterIdenticalMiniprotTx(original_gff = gff_ls$query_gff,
                                                      mp_gff = gff_ls$mp_query_gff)
    
    gff_ls$mp_subject_gff <- .filterIdenticalMiniprotTx(original_gff = gff_ls$subject_gff,
                                                        mp_gff = gff_ls$mp_subject_gff)
    
    gff_ls$mp_query_gff <- .filterLongMiniprotTx(mp_gff = gff_ls$mp_query_gff,
                                                 original_gff = gff_ls$subject_gff,
                                                 len_diff = len_diff)
    
    gff_ls$mp_subject_gff <- .filterLongMiniprotTx(mp_gff = gff_ls$mp_subject_gff,
                                                   original_gff = gff_ls$query_gff,
                                                   len_diff = len_diff)
    
    mp_query_novel <- .filterMiniprotTxOnNovelLoci(original_gff = gff_ls$query_gff,
                                                   mp_gff = gff_ls$mp_query_gff,
                                                   mp_cds = cds_ls$mp_query_cds)
    
    mp_subject_novel <- .filterMiniprotTxOnNovelLoci(original_gff = gff_ls$subject_gff,
                                                     mp_gff = gff_ls$mp_subject_gff,
                                                     mp_cds = cds_ls$mp_subject_cds)
    gff_ls$query_gff <- suppressWarnings(c(gff_ls$query_gff, mp_query_novel))
    cds_ls$query_cds <- c(cds_ls$query_cds,
                          cds_ls$mp_query_cds[names(cds_ls$mp_query_cds) %in% mp_query_novel$ID])
    valid_longer_query_mp_gff <- .filterMiniprotTxOnGeneLoci(original_gff = gff_ls$query_gff,
                                                             mp_gff = gff_ls$mp_query_gff,
                                                             original_cds = cds_ls$query_cds,
                                                             mp_cds = cds_ls$mp_query_cds)
    
    gff_ls$subject_gff <- suppressWarnings(c(gff_ls$subject_gff, mp_subject_novel))
    cds_ls$subject_cds <- c(cds_ls$subject_cds,
                            cds_ls$mp_subject_cds[names(cds_ls$mp_subject_cds) %in% mp_subject_novel$ID])
    valid_longer_subject_mp_gff <- .filterMiniprotTxOnGeneLoci(original_gff = gff_ls$subject_gff,
                                                               mp_gff = gff_ls$mp_subject_gff,
                                                               original_cds = cds_ls$subject_cds,
                                                               mp_cds = cds_ls$mp_subject_cds)
    
    out <- list(query = suppressWarnings(c(gff_ls$query_gff, valid_longer_query_mp_gff)),
                subject = suppressWarnings(c(gff_ls$subject_gff, valid_longer_subject_mp_gff)))
    
    return(out)
}

.getCDSlist <- function(h5, mp){
    # Import GFF files for query and subject genomes
    if(mp){
        query_cds <- readDNAStringSet(as.vector(h5$files$query_cds))
        mp_query_cds <- readDNAStringSet(as.vector(h5$miniprot$s2q_cds))
        subject_cds <- readDNAStringSet(as.vector(h5$files$subject_cds))
        mp_subject_cds <- readDNAStringSet(as.vector(h5$miniprot$q2s_cds))
        
        # Return the ordered cds data as a list
        out <- list(query_cds = query_cds, subject_cds = subject_cds,
                    mp_query_cds = mp_query_cds, mp_subject_cds = mp_subject_cds)
        
    } else {
        query_cds <- readDNAStringSet(as.vector(h5$files$query_cds))
        subject_cds <- readDNAStringSet(as.vector(h5$files$subject_cds))
        
        # Return the ordered cds data as a list
        out <- list(query_cds = query_cds, subject_cds = subject_cds)
    }
    
    return(out)
}

.filterIdenticalMiniprotTx <- function(original_gff, mp_gff){
    mp_gff_block <- .getCDSblock(gff = mp_gff)
    mp_gff_block_uniq <- mp_gff_block[!duplicated(mp_gff_block)]
    
    # Get CDS blocks from gff2 and remove overlaps from gff1_block_uniq
    original_gff_block <- .getCDSblock(gff = original_gff)
    mp_gff_block_uniq <- mp_gff_block_uniq[!mp_gff_block_uniq %in% original_gff_block]
    
    # Get unique transcripts from gff1
    out <- .getUniqTx(mp_gff = mp_gff, mp_gff_block_uniq = mp_gff_block_uniq)
    return(out)
}

#' @importFrom BiocGenerics start end
.getCDSblock <- function(gff){
    # Identify CDS and transcript indices
    gff_cds_i <- gff$type == "CDS"
    gff_tx_i <- gff$type %in% c("transcript", "mRNA")
    
    # Extract start and end positions of CDS
    gff_cds_start <- start(gff[gff_cds_i])
    gff_cds_end <- end(gff[gff_cds_i])
    
    # Create exon strings combining start and end positions
    gff_cds_exon <- paste(gff_cds_start, gff_cds_end, sep = "-")
    
    # Map CDS to their parent transcripts
    map_to_tx <- match(unlist(gff$Parent[gff_cds_i]), gff$ID[gff_tx_i])
    
    # Identify the first occurrence of each transcript
    first_i <- !duplicated(map_to_tx)
    
    # Get the chromosome names for CDS
    gff_cds_chr <- as.character(seqnames(gff[gff_cds_i]))
    
    # Add chromosome information to the first occurrence of each transcript's CDS block
    gff_cds_exon[first_i] <- paste(gff_cds_chr[first_i], gff_cds_exon[first_i], sep = ":")
    
    # Concatenate CDS blocks for each transcript
    out <- tapply(gff_cds_exon, map_to_tx, paste, collapse = ",")
    
    return(out)
}


#' @importFrom BiocGenerics unlist
.getUniqTx <- function(mp_gff, mp_gff_block_uniq){
    # Identify transcript indices
    mp_gff_tx_i <- mp_gff$type %in% c("transcript", "mRNA")
    
    # Extract unique transcripts using the unique CDS blocks
    out_tx <- mp_gff[mp_gff_tx_i][as.numeric(names(mp_gff_block_uniq))]
    
    # Identify and extract elements associated with the unique transcripts
    out_element <- unlist(mp_gff$Parent[!mp_gff_tx_i]) %in% out_tx$ID
    out_element <- mp_gff[!mp_gff_tx_i][out_element]
    
    # Combine unique transcripts and their associated elements
    out <- c(out_tx, out_element)
    
    return(out)
}

.filterLongMiniprotTx <- function(mp_gff, original_gff, len_diff){
    mp_gff <- .checkCDSlen(mp_gff = mp_gff,
                           original_gff = original_gff,
                           len_diff = len_diff)
    mp_gff <- .checkCDSstretch(mp_gff = mp_gff,
                               original_gff = original_gff,
                               len_diff = len_diff)
    return(mp_gff)
}

.checkCDSlen <- function(mp_gff, original_gff, len_diff){
    # Calculate CDS lengths for the given GFF data
    mp_cds_len <- .getCDSlen(gff = mp_gff)
    index <- as.numeric(names(mp_cds_len))
    mp_tx_i <- mp_gff$type %in% c("transcript", "mRNA")
    names(mp_cds_len) <- sub("\\s.+", "", mp_gff$Target[mp_tx_i][index])
    
    # Calculate CDS lengths for the reference GFF data
    original_cds_len <- .getCDSlen(gff = original_gff)
    index <- as.numeric(names(original_cds_len))
    original_tx_i <- original_gff$type %in% c("transcript", "mRNA")
    names(original_cds_len) <- original_gff$ID[original_tx_i][index]
    
    # Match CDS lengths by transcript IDs
    id_hit <- match(names(mp_cds_len), names(original_cds_len))
    original_cds_len <- original_cds_len[id_hit]
    
    # Determine which transcripts are longer
    longer <- mp_cds_len
    is_original_longer <- longer < original_cds_len
    longer[is_original_longer] <- original_cds_len[is_original_longer]
    
    # Calculate valid transcripts based on length difference criteria
    valid <- abs(mp_cds_len - original_cds_len) / longer <= len_diff
    valid_tx <- mp_gff[mp_tx_i][as.vector(valid)]
    
    # Extract non-transcript elements associated with valid transcripts
    non_tx_gff <- mp_gff[!mp_tx_i]
    out <- c(valid_tx, non_tx_gff[unlist(non_tx_gff$Parent) %in% valid_tx$ID])
    return(out)
}

.checkCDSstretch <- function(mp_gff, original_gff, len_diff){
    # Calculate Tx lengths for the given GFF data
    mp_tx_i <- mp_gff$type %in% c("transcript", "mRNA")
    mp_tx_len <- width(mp_gff[mp_tx_i])
    names(mp_tx_len) <- sub("\\s.+", "", mp_gff$Target[mp_tx_i])
    
    # Calculate Tx lengths for the reference GFF data
    original_cds_i <- original_gff$type %in% "CDS"
    original_cds_gff <- original_gff[original_cds_i]
    original_cds_parents <- unlist(original_cds_gff$Parent)
    original_cds_gff <- original_cds_gff[original_cds_parents %in% names(mp_tx_len)]
    original_cds_parents <- unlist(original_cds_gff$Parent)
    cds_min <- tapply(start(original_cds_gff), original_cds_parents, min)
    cds_max <- tapply(end(original_cds_gff), original_cds_parents, max)
    original_tx_len <- cds_max - cds_min
    
    # Match CDS lengths by transcript IDs
    id_hit <- match(names(mp_tx_len), names(original_tx_len))
    original_tx_len <- original_tx_len[id_hit]
    
    # Determine which transcripts are longer
    longer <- mp_tx_len
    is_original_longer <- longer < original_tx_len
    longer[is_original_longer] <- original_tx_len[is_original_longer]
    
    # Calculate valid transcripts based on length difference criteria
    valid <- abs(mp_tx_len - original_tx_len) / longer <= len_diff
    valid_tx <- mp_gff[mp_tx_i][as.vector(valid)]
    
    # Extract non-transcript elements associated with valid transcripts
    non_tx_gff <- mp_gff[!mp_tx_i]
    out <- c(valid_tx, non_tx_gff[unlist(non_tx_gff$Parent) %in% valid_tx$ID])
    return(out)
}

.getCDSlen <- function(gff){
    # Identify indices for CDS and transcript elements
    gff_cds_i <- gff$type == "CDS"
    gff_tx_i <-  gff$type %in% c("transcript", "mRNA")
    
    # Extract start and end positions for CDS elements
    gff_cds_start <- start(gff[gff_cds_i])
    gff_cds_end <- end(gff[gff_cds_i])
    
    # Calculate the length of each CDS element
    gff_cds_len <- gff_cds_end - gff_cds_start
    
    # Map CDS elements to their corresponding transcripts
    map_to_tx <- match(unlist(gff$Parent[gff_cds_i]), gff$ID[gff_tx_i])
    
    # Sum the lengths of CDS elements for each transcript
    out <- tapply(gff_cds_len, map_to_tx, sum)
    
    return(out)
}

.filterMiniprotTxOnGeneLoci <- function(original_gff,
                                        mp_gff,
                                        original_cds,
                                        mp_cds){
    
    # Find overlapping Miniprot Tx on original Tx
    original_gff_cds_i <- original_gff$type %in% c("CDS")
    mp_gff_cds_i <- mp_gff$type %in% c("CDS")
    mp_gff_tx_i <- mp_gff$type %in% c("transcript", "mRNA")
    ol <- findOverlaps(mp_gff[mp_gff_cds_i], original_gff[original_gff_cds_i])
    ol <- as.data.frame(ol)
    ol$queryHits <- unlist(mp_gff$Parent[mp_gff_cds_i][ol$queryHits])
    ol$subject_gene <- original_gff$gene_id[original_gff_cds_i][ol$subjectHits]
    ol$subjectHits <- unlist(original_gff$Parent[original_gff_cds_i][ol$subjectHits])
    ol <- unique(ol)
    
    # Filter out chimeric Miniprot Tx overlapping more than one gene loci
    check_chimeric <- tapply(X = ol$subject_gene, INDEX = ol$queryHits, FUN = unique)
    check_chimeric <- sapply(check_chimeric, length)
    chimeric_tx <- names(check_chimeric[check_chimeric > 1])
    ol <- subset(ol, subset = !queryHits %in% chimeric_tx)
    
    # Check initial and terminal codons of Tx
    original_tx_init <- .checkInitCodon(cds = original_cds)
    mp_tx_init <- .checkInitCodon(cds = mp_cds)
    original_tx_term <- .checkTermCodon(cds = original_cds)
    mp_tx_term <- .checkTermCodon(cds = mp_cds)
    oiriginal_tx_hit <- match(ol$subjectHits, names(original_tx_init))
    mp_tx_hit <- match(ol$queryHits, names(mp_tx_init))
    ol$mp_valid <- mp_tx_init[mp_tx_hit] & mp_tx_term[mp_tx_hit]
    ol$original_valid <- original_tx_init[oiriginal_tx_hit] & original_tx_term[oiriginal_tx_hit]
    ol <- subset(ol, subset = mp_valid)
    
    # Pick up if Miniprot Tx is only valid
    original_any_valid <- tapply(X = ol$original_valid,
                                 INDEX = ol$subject_gene,
                                 FUN = any)
    mp_any_valid <- tapply(X = ol$mp_valid,
                           INDEX = ol$subject_gene,
                           FUN = any)
    mp_only_valid <- mp_any_valid & !original_any_valid
    mp_only_valid_hit <- match(ol$subject_gene, names(mp_only_valid))
    ol$mp_only_valid <- mp_only_valid[mp_only_valid_hit]
    
    # Check Tx lengths
    mp_tx_len <- width(mp_cds)
    original_tx_len <- width(original_cds)
    mp_tx_hit <- match(ol$queryHits, names(mp_cds))
    ol$mp_tx_len <- mp_tx_len[mp_tx_hit]
    original_tx_hit <- match(ol$subjectHits, names(original_cds))
    ol$original_tx_len <- original_tx_len[original_tx_hit]
    mp_valid_tx_max_len <- tapply(X = ol$mp_tx_len[ol$mp_valid],
                                  INDEX = ol$subject_gene[ol$mp_valid],
                                  FUN = max)
    oiriginal_tx_hit <- match(ol$subject_gene, names(mp_valid_tx_max_len))
    ol$mp_valid_tx_max_len <- mp_valid_tx_max_len[oiriginal_tx_hit]
    original_valid_tx_max_len <- tapply(X = ol$original_tx_len[ol$original_valid],
                                        INDEX = ol$subject_gene[ol$original_valid],
                                        FUN = max)
    oiriginal_tx_hit <- match(ol$subject_gene, names(original_valid_tx_max_len))
    ol$original_valid_tx_max_len <- original_valid_tx_max_len[oiriginal_tx_hit]
    ol$original_valid_tx_max_len[is.na(ol$original_valid_tx_max_len)] <- 0
    ol <- subset(ol, subset = !duplicated(subset(ol, select = c(queryHits, subject_gene))))
    ol$mp_longer_than_original <- ol$mp_valid_tx_max_len > ol$original_valid_tx_max_len
    
    # Pick up longest Miniprot Tx that valid and longer than original Tx at each locus
    ol$mp_max_len_tx <- ol$mp_tx_len == ol$mp_valid_tx_max_len
    valid_longer_mp_tx <- ol$queryHits[ol$mp_longer_than_original & ol$mp_max_len_tx]
    valid_longer_mp_tx_locus <- ol$subject_gene[ol$mp_longer_than_original & ol$mp_max_len_tx]
    
    valid_longer_mp_tx_gff <- mp_gff[mp_gff$ID %in% valid_longer_mp_tx]
    valid_longer_mp_tx_gff$Parent <- lapply(valid_longer_mp_tx_locus, c)
    mp_gff_element_i <- !mp_gff$type %in% c("gene", "transcript", "mRNA")
    valid_longer_mp_tx_element <- mp_gff[mp_gff_element_i][unlist(mp_gff$Parent[mp_gff_element_i]) %in% valid_longer_mp_tx_gff$ID]
    valid_longer_mp_gff <- suppressWarnings(c(valid_longer_mp_tx_gff, valid_longer_mp_tx_element))
    return(valid_longer_mp_gff)
}

.checkInitCodon <- function(cds){
    return(substr(cds, 1, 3) == "ATG")
}

.checkTermCodon <- function(cds){
    len <- width(cds)
    return(substr(cds, len - 2, len) %in% c("TAA", "TGA", "TAG"))
}

.filterMiniprotTxOnNovelLoci <- function(original_gff,
                                         mp_gff,
                                         mp_cds){
    if(is.null(mp_gff)){
        return(NULL)
    }
    
    # Find overlapping Miniprot Tx on original Tx
    original_gff_cds_i <- original_gff$type %in% c("CDS")
    mp_gff_cds_i <- mp_gff$type %in% c("CDS")
    mp_gff_tx_i <- mp_gff$type %in% c("transcript", "mRNA")
    ol <- findOverlaps(mp_gff[mp_gff_cds_i], original_gff[original_gff_cds_i])
    ol <- as.data.frame(ol)
    ol$queryHits <- unlist(mp_gff$Parent[mp_gff_cds_i][ol$queryHits])
    tx_id <- mp_gff$ID[mp_gff_tx_i]
    non_ol_mp_tx <- tx_id[!tx_id %in% ol$queryHits]
    
    if(length(non_ol_mp_tx) == 0){
        return(NULL)
        
    } else {
        non_ol_mp_tx_gff <- mp_gff[mp_gff$ID %in% non_ol_mp_tx]
        mp_gff_element_i <- !mp_gff$type %in% c("gene", "transcript", "mRNA")
        non_ol_mp_tx_element <- mp_gff[mp_gff_element_i][unlist(mp_gff$Parent[mp_gff_element_i]) %in% non_ol_mp_tx_gff$ID]
        non_ol_mp_gff <- c(non_ol_mp_tx_gff, non_ol_mp_tx_element)
    }
    
    # Check if Tx is fully included in another Tx
    non_ol_mp_gff <- .filterTxWithinTx(gff = non_ol_mp_gff)
    
    # Group Tx based on overlaps
    mp_gff_tx_i <- non_ol_mp_gff$type %in% c("transcript", "mRNA")
    rest_mp_tx_gff <- non_ol_mp_gff[mp_gff_tx_i]
    init <- .checkInitCodon(cds = mp_cds) # Check initial codons of Tx
    term <- .checkTermCodon(cds = mp_cds) # Check terminal codons of Tx
    grp <- .groupOverlaps(gff = rest_mp_tx_gff, init = init, term = term)
    
    # Find non-overlapping Tx
    non_ol_tx_id <- rest_mp_tx_gff$ID[grp$members[grp$n_member == 1]]
    if(length(non_ol_tx_id) != 0){
        non_ol_gff <- .orgMiniprotFilteredGFF(gff = mp_gff,
                                              tx_id = non_ol_tx_id)
        out <- non_ol_gff
    }
    grp <- subset(grp, subset = n_member != 1)
    
    # Pick non valid Tx groups
    n_valid_member <- tapply(grp$cds_valid, grp$rep, sum)
    non_valid_rep <- as.numeric(names(n_valid_member)[n_valid_member == 0])
    non_valid_members <- grp$members[grp$rep %in% non_valid_rep]
    non_valid_tx_rep <- rest_mp_tx_gff$ID[non_valid_rep]
    if(length(non_valid_tx_rep) != 0){
        non_valid_gff <- .orgMiniprotFilteredGFF(gff = mp_gff,
                                                 tx_id = non_valid_tx_rep)
        non_valid_member_gff <- .orgMembersGFF(gff = mp_gff,
                                               grp = grp,
                                               rep = non_valid_rep)
        if(is.null(out)){
            out <- c(non_valid_gff, non_valid_member_gff)
            
        } else {
            out <- c(out, non_valid_gff, non_valid_member_gff)
        }
    }
    
    # Regrouping
    grp <- subset(grp, subset = cds_valid)
    valid_rep <- as.numeric(names(n_valid_member)[n_valid_member != 0])
    grp <- subset(grp, subset = rep %in% valid_rep)
    rest_mp_tx_gff <- rest_mp_tx_gff[rest_mp_tx_gff$ID %in% grp$member_id]
    grp <- .groupOverlaps(gff = rest_mp_tx_gff, init = init, term = term)
    n_valid_member <- tapply(grp$cds_valid, grp$rep, sum)
    
    # Find single valid Tx
    single_valid_rep <- as.numeric(names(n_valid_member)[n_valid_member == 1])
    single_valid_members <- grp$members[grp$rep %in% single_valid_rep & grp$cds_valid]
    single_valid_tx <- rest_mp_tx_gff$ID[single_valid_members]
    if(length(single_valid_tx) != 0){
        single_valid_gff <- .orgMiniprotFilteredGFF(gff = mp_gff, tx_id = single_valid_tx)
        if(is.null(out)){
            out <- single_valid_gff
            
        } else {
            out <- c(out, single_valid_gff)
        }
    }
    
    # Pick multiple valid Tx groups
    multiple_valid_rep <- as.numeric(names(n_valid_member)[n_valid_member > 1])
    multiple_valid_tx_rep <- rest_mp_tx_gff$ID[multiple_valid_rep]
    if(length(multiple_valid_tx_rep) != 0){
        multiple_valid_gff <- .orgMiniprotFilteredGFF(gff = mp_gff,
                                                      tx_id = multiple_valid_tx_rep)
        multiple_valid_member_gff <- .orgMembersGFF(gff = mp_gff,
                                                    grp = grp,
                                                    rep = multiple_valid_rep)
        if(is.null(out)){
            out <- c(multiple_valid_gff, multiple_valid_member_gff)
            
        } else {
            out <- c(out, multiple_valid_gff, multiple_valid_member_gff)
        }
    }
    
    if(!is.null(out)){
        out <- .setGeneID(gff = out)
    }
    
    return(out)
}

.filterTxWithinTx <- function(gff){
    gff_cds_i <- gff$type == "CDS"
    within_ol <- findOverlaps(gff[gff_cds_i], gff[gff_cds_i], type = "within")
    within_ol <- as.data.frame(within_ol)
    within_ol <- subset(within_ol, subset = queryHits != subjectHits)
    within_ol$queryHits <- unlist(gff$Parent[gff_cds_i][within_ol$queryHits])
    within_ol$subjectHits <- unlist(gff$Parent[gff_cds_i][within_ol$subjectHits])
    within_ol$pair <- paste(within_ol$queryHits, within_ol$subjectHits, sep = "_")
    within_ol_num <- table(within_ol$pair)
    within_ol_num_id <- sub("_.+", "", names(within_ol_num))
    num_cds <- tapply(X = gff$ID[gff_cds_i],
                      INDEX = unlist(gff$Parent[gff_cds_i]),
                      FUN = length)
    id_hit <- match(within_ol_num_id, names(num_cds))
    within_ol_id <- within_ol_num_id[within_ol_num == num_cds[id_hit]]
    gff_tx_i <- gff$type %in% c("transcript", "mRNA")
    gff_tx <- gff[gff_tx_i][!gff$ID[gff_tx_i] %in% within_ol_id]
    gff_element <- gff[!gff_tx_i][!unlist(gff$Parent[!gff_tx_i]) %in% within_ol_id]
    rest_gff <- c(gff_tx, gff_element)
    return(rest_gff)
}

#' @importFrom igraph graph_from_data_frame V components
.groupOverlaps <- function(gff, init, term){
    ol <- findOverlaps(gff, gff)
    g <- graph_from_data_frame(d = as.data.frame(ol), directed = FALSE)
    grp <- split(V(g)$name, components(g)$membership)
    grp <- lapply(grp, function(x){
        x <- sort(as.numeric(x))
        return(data.frame(rep = x[1], members = x, n_member = length(x)))
    })
    grp <- do.call("rbind", grp)
    
    # Set Tx IDs
    grp$rep_id <- gff$ID[grp$rep]
    grp$member_id <- gff$ID[grp$members]
    
    hit <- match(grp$member_id, names(init))
    grp$cds_valid <- init[hit] & term[hit]
    return(grp)
}

.orgMiniprotFilteredGFF <- function(gff, tx_id){
    gff_tx_i <- gff$type %in% c("transcript", "mRNA")
    gff_element_i <- !gff$type %in% c("gene", "transcript", "mRNA")
    out_tx <-  gff[gff$ID %in% tx_id]
    out_element <- gff[gff_element_i][unlist(gff$Parent[gff_element_i]) %in% out_tx$ID]
    out_gene <- out_tx
    out_gene$type <- "gene"
    out_gene$ID <- paste(out_gene$ID, "gene", sep = ":")
    out_gene$Parent <- lapply(rep("", length(out_gene)), c)
    out_tx$Parent <- lapply(out_gene$ID, c)
    non_ol_gff <- c(out_gene, out_tx, out_element)
}

.orgMembersGFF <- function(gff, grp, rep){
    member_i <- grp$members[grp$rep %in% rep]
    member_i <- member_i[!member_i %in% rep]
    member_tx <- grp$member_id[grp$members %in% member_i]
    member_gff <- gff[gff$ID %in% member_tx]
    hit <- match(member_gff$ID, grp$member_id)
    member_gff$Parent <- lapply(paste0(grp$rep_id[hit], ":gene"), c)
    element_i <- !gff$type %in% c("gene", "transcript", "mRNA")
    element_gff <- gff[element_i][unlist(gff$Parent[element_i]) %in% member_gff$ID]
    out <- c(member_gff, element_gff)
    return(out)
}
################################################################################

.fixGFF <- function(gff){
    entry_type <-  c("gene", "transcript", "mRNA", "five_prime_UTR", "exon",
                     "CDS", "three_prime_UTR")
    gff <- gff[gff$type %in% entry_type]
    gff <- .setGeneID(gff = gff)
    gff <- .fixGFFrange(gff = gff)
    gff <- .fixGFFexon(gff = gff)
    gff <- .fixGFFphase(gff = gff)
    return(gff)
}

.setGeneID <- function(gff){
    gff$gene_id <- NA
    
    gene_i <- gff$type %in% "gene"
    gff$gene_id[gene_i] <- gff$ID[gene_i]
    
    # For transcripts or mRNA, set gene_id to their ID
    tx_i <- gff$type %in% c("transcript", "mRNA")
    tx_p <- unlist(gff$Parent[tx_i])
    if(length(tx_p) > sum(tx_i)){
        tx_p <- sapply(gff$Parent[tx_i], "[", 1)
        gff$Parent[tx_i] <- lapply(tx_p, c)
    }
    hit <- match(tx_p, gff$ID)
    gff$gene_id[tx_i] <- gff$ID[hit]
    
    # For other elements, set gene_id based on their parent transcript/mRNA
    element_i <- !gff$type %in% c("gene", "transcript", "mRNA")
    element_p <- unlist(gff$Parent[element_i])
    if(length(element_p) > sum(element_i)){
        element_p <- sapply(gff$Parent[element_i], "[", 1)
        gff$Parent[element_i] <- lapply(element_p, c)
    }
    hit <- match(element_p, gff$ID[tx_i])
    gff$gene_id[element_i] <- gff$gene_id[tx_i][hit]
    return(gff)
}

#' Fix GFF Range
#'
#' This function adjusts the start and end positions of transcripts and genes to cover the complete range of their member elements.
#'
#' @importFrom rtracklayer start end start<- end<-
.fixGFFrange <- function(gff){
    # Fix the start and end positions of each transcript to
    # cover the whole range of member elements (CDS, exon, and UTRs)
    tx_i <- which(gff$type %in% c("transcript", "mRNA"))
    element_i <- !gff$type %in% c("gene", "transcript", "mRNA")
    
    # Determine minimum start and maximum end positions for member elements
    min_start <- tapply(start(gff[element_i]), unlist(gff$Parent[element_i]), min)
    max_end <- tapply(end(gff[element_i]), unlist(gff$Parent[element_i]), max)
    
    # Match and update transcript start and end positions
    hit <- match(gff$ID[tx_i], names(min_start))
    tx_start <- min_start[hit]
    tx_end <- max_end[hit]
    not_na_start <- !is.na(tx_start)
    not_na_end <- !is.na(tx_end)
    start(gff[tx_i[not_na_start]]) <- tx_start[not_na_start]
    end(gff[tx_i[not_na_end]]) <- tx_end[not_na_end]
    
    # Fix the start and end positions of each gene to
    # cover the whole range of member transcripts
    gene_i <- which(gff$type == "gene")
    tx_i <- which(gff$type %in% c("transcript", "mRNA"))
    
    # Determine minimum start and maximum end positions for transcripts
    min_start <- tapply(start(gff[tx_i]), gff$gene_id[tx_i], min)
    max_end <- tapply(end(gff[tx_i]), gff$gene_id[tx_i], max)
    
    # Match and update gene start and end positions
    hit <- match(gff$gene_id[gene_i], names(min_start))
    gene_start <- min_start[hit]
    gene_end <- max_end[hit]
    start(gff[gene_i]) <- gene_start
    end(gff[gene_i]) <- gene_end
    
    return(gff)
}

.fixGFFexon <- function(gff){
    gff <- gff[gff$type != "exon"]
    gff_exon <- gff[gff$type %in% c("CDS", "five_prime_UTR", "three_prime_UTR")]
    gff_exon$type <- "exon"
    gff_exon$Name <- gff_exon$ID <- paste0(unlist(gff_exon$Parent), ":exon")
    gff_exon$score <- gff_exon$phase <- NA
    
    tx_i <- gff$type %in% c("transcript", "mRNA")
    non_cds <- gff[tx_i][!gff$ID[tx_i] %in% unlist(gff_exon$Parent)]
    if(length(non_cds) > 0){
        non_cds$type <- "exon"
        non_cds$Parent <- lapply(non_cds$ID, c)
        non_cds$Name <- non_cds$ID <- paste0(non_cds$ID, ":exon")
        non_cds$score <- non_cds$phase <- NA
        gff_exon <- c(gff_exon, non_cds)
    }
    
    out <- c(gff, gff_exon)
    out <- out[order(as.numeric(seqnames(out)), start(out), as.numeric(out$type))]
    return(out)
}


#' Fix GFF Phase
#'
#' This function adjusts the phase of CDS features in GFF annotations.
#'
#' @importFrom BiocGenerics start
.fixGFFphase <- function(gff){
    # Extract CDS features
    gff_cds <- gff[gff$type == "CDS"]
    
    # Adjust phase for plus strand
    gff_cds_plus <- .phasePlus(gff_cds = gff_cds)
    
    # Adjust phase for minus strand
    gff_cds_minus <- .phaseMinus(gff_cds = gff_cds)
    
    # Combine non-CDS features with adjusted CDS features
    out <- c(gff[gff$type != "CDS"], gff_cds_plus, gff_cds_minus)
    
    # Set factor levels for feature types
    out$type <- factor(out$type, levels = c("gene", "transcript", "mRNA",
                                            "five_prime_UTR", "exon",
                                            "CDS", "three_prime_UTR"))
    out$type <- droplevels(out$type)
    
    # Order the features by sequence name, start position, and type
    out <- out[order(as.numeric(seqnames(out)), start(out), as.numeric(out$type))]
    return(out)
}
#' Adjust Phase for CDS Features on Plus Strand
#'
#' This function adjusts the phase for CDS features on the plus strand in GFF annotations.
#'
#' @importFrom BiocGenerics start width
#' @importFrom GenomeInfoDb seqnames
.phasePlus <- function(gff_cds){
    # Filter CDS features on the plus strand
    gff_cds_plus <- gff_cds[as.character(strand(gff_cds)) == "+"]
    gff_cds_plus <- gff_cds_plus[order(as.numeric(seqnames(gff_cds_plus)), start(gff_cds_plus))]
    
    # Extract parent IDs for CDS features
    gff_cds_plus_parent <- unlist(gff_cds_plus$Parent)
    
    # Initialize phase for the first CDS in each transcript
    target_i <- which(!duplicated(gff_cds_plus_parent))
    gff_cds_plus$phase[target_i] <- 0
    
    # Calculate the phase for the next CDS feature
    next_phase <- (3 - (width(gff_cds_plus[target_i]) - gff_cds_plus$phase[target_i]) %% 3) %% 3
    names(next_phase) <- gff_cds_plus_parent[target_i]
    gff_cds_plus_parent[target_i] <- "NA"
    
    # Iterate over the remaining CDS features to set their phases
    while(TRUE){
        target_i <- which(!duplicated(gff_cds_plus_parent))[-1]
        if(length(target_i) == 0){
            break
        }
        next_phase <- next_phase[names(next_phase) %in% gff_cds_plus_parent[target_i]]
        gff_cds_plus$phase[target_i] <- next_phase
        next_phase <- (3 - (width(gff_cds_plus[target_i]) - gff_cds_plus$phase[target_i]) %% 3) %% 3
        names(next_phase) <- gff_cds_plus_parent[target_i]
        gff_cds_plus_parent[target_i] <- "NA"
    }
    
    return(gff_cds_plus)
}

#' Adjust Phase for CDS Features on Minus Strand
#'
#' This function adjusts the phase for CDS features on the minus strand in GFF annotations.
#'
#' @importFrom BiocGenerics end width
#' @importFrom GenomeInfoDb seqnames
.phaseMinus <- function(gff_cds){
    # Filter CDS features on the minus strand
    gff_cds_minus <- gff_cds[as.character(strand(gff_cds)) == "-"]
    gff_cds_minus <- gff_cds_minus[order(as.numeric(seqnames(gff_cds_minus)),
                                         end(gff_cds_minus), decreasing = TRUE)]
    
    # Extract parent IDs for CDS features
    gff_cds_minus_parent <- unlist(gff_cds_minus$Parent)
    
    # Initialize phase for the first CDS in each transcript
    target_i <- which(!duplicated(gff_cds_minus_parent))
    gff_cds_minus$phase[target_i] <- 0
    
    # Calculate the phase for the next CDS feature
    next_phase <- (3 - (width(gff_cds_minus[target_i]) - gff_cds_minus$phase[target_i]) %% 3) %% 3
    names(next_phase) <- gff_cds_minus_parent[target_i]
    gff_cds_minus_parent[target_i] <- "NA"
    
    # Iterate over the remaining CDS features to set their phases
    while(TRUE){
        target_i <- which(!duplicated(gff_cds_minus_parent))[-1]
        if(length(target_i) == 0){
            break
        }
        next_phase <- next_phase[names(next_phase) %in% gff_cds_minus_parent[target_i]]
        gff_cds_minus$phase[target_i] <- next_phase
        next_phase <- (3 - (width(gff_cds_minus[target_i]) - gff_cds_minus$phase[target_i]) %% 3) %% 3
        names(next_phase) <- gff_cds_minus_parent[target_i]
        gff_cds_minus_parent[target_i] <- "NA"
    }
    
    return(gff_cds_minus)
}

#' Create FASTA files
#'
#' This function generates FASTA files for query and subject CDS from the OrthoPairDB object.
#'
#' @param object A OrthoPairDB object.
#' @param out_dir Output directory for FASTA files.
#'
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#' @importFrom Biostrings writeXStringSet
.createFASTA <- function(object, out_dir, merge){
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    
    if(merge){
        # Create CDS sequences for query and subject genomes
        q_cds <- .makeCDS(gff = as.vector(h5$files$query_gff),
                          genome = as.vector(h5$files$query_genome))
        s_cds <- .makeCDS(gff = as.vector(h5$files$subject_gff),
                          genome = as.vector(h5$files$subject_genome))
        
        # Define filenames for the CDS FASTA files
        q_cds_fn <- file.path(out_dir, "miniprot_merge_query.cds")
        s_cds_fn <- file.path(out_dir, "miniprot_merge_subject.cds")
        
        # Write the CDS sequences to FASTA files
        writeXStringSet(q_cds, q_cds_fn)
        writeXStringSet(s_cds, s_cds_fn)
        
        # Overwrite the HDF5 file with the new CDS FASTA filenames
        .h5overwrite(obj = q_cds_fn, file = object$h5, "files/query_cds")
        .h5overwrite(obj = s_cds_fn, file = object$h5, "files/subject_cds")
        
    } else {
        # Create CDS sequences for query and subject genomes
        q_cds <- .makeCDS(gff = as.vector(h5$miniprot$s2q_gff),
                          genome = as.vector(h5$files$query_genome))
        s_cds <- .makeCDS(gff = as.vector(h5$miniprot$q2s_gff),
                          genome = as.vector(h5$files$subject_genome))
        
        # Define filenames for the CDS FASTA files
        q_cds_fn <- file.path(out_dir, "miniprot_query.cds")
        s_cds_fn <- file.path(out_dir, "miniprot_subject.cds")
        
        # Write the CDS sequences to FASTA files
        writeXStringSet(q_cds, q_cds_fn)
        writeXStringSet(s_cds, s_cds_fn)
        
        # Overwrite the HDF5 file with the new CDS FASTA filenames
        .h5overwrite(obj = q_cds_fn, file = object$h5, "miniprot/s2q_cds")
        .h5overwrite(obj = s_cds_fn, file = object$h5, "miniprot/q2s_cds")
    }
}

#' @importFrom Biostrings readDNAStringSet
#' @importFrom GenomicFeatures cdsBy extractTranscriptSeqs
#' @importFrom txdbmaker makeTxDbFromGFF
#' @import BSgenome
#'
.makeCDS <- function(gff, genome){
    # Create a TxDb object from the GFF file
    txdb <- makeTxDbFromGFF(file = gff)
    
    # Read the genome file as a DNAStringSet object
    genome <- readDNAStringSet(filepath = genome)
    
    # Extract CDS sequences from the TxDb object
    cds_db <- cdsBy(x = txdb, by = "tx", use.names = TRUE)
    cds <- extractTranscriptSeqs(x = genome, transcripts = cds_db)
    
    # Order CDS sequences by their names
    cds <- cds[order(names(cds))]
    return(cds)
}

#
# .filterMPgenes <- function(df, pattern = "^query_"){
#     df <- df[order(df$gene), ]
#     mp_tx_index <- grepl(pattern, df$tx)
#     eval <- data.frame(id = unique(df$gene))
#
#     mp_valid_max_len <- tapply(df$len[mp_tx_index][df$valid[mp_tx_index]],
#                                df$gene[mp_tx_index][df$valid[mp_tx_index]],
#                                max)
#     hit <- match(eval$id, names(mp_valid_max_len))
#     eval$mp_valid_max_len <- mp_valid_max_len[hit]
#     eval$mp_valid_max_len[is.na(eval$mp_valid_max_len)] <- 0
#
#     non_mp_valid_max_len <- tapply(df$len[!mp_tx_index][df$valid[!mp_tx_index]],
#                                    df$gene[!mp_tx_index][df$valid[!mp_tx_index]],
#                                    max)
#     hit <- match(eval$id, names(non_mp_valid_max_len))
#     eval$non_mp_valid_max_len <- non_mp_valid_max_len[hit]
#     eval$non_mp_valid_max_len[is.na(eval$non_mp_valid_max_len)] <- 0
#
#     max_tx_index <- tapply(df$len[df$valid],
#                            df$gene[df$valid],
#                            which.max)
#     max_tx_index <- sapply(max_tx_index, function(x){
#         if(length(x) == 0){
#             x <- NA
#         }
#         return(x)
#     })
#     hit <- match(eval$id, names(max_tx_index))
#     eval$max_tx_index <- max_tx_index[hit]
#
#     max_tx_is_na <- is.na(eval$max_tx_index)
#     hit <- match(grep(pattern, eval$id[max_tx_is_na], value = TRUE),
#                  df$gene)
#     max_tx_index <- tapply(df$len[hit],
#                            df$gene[hit],
#                            which.max)
#     hit <- match(eval$id[max_tx_is_na], names(max_tx_index))
#     eval$max_tx_index[max_tx_is_na] <- max_tx_index[hit]
#
#     start_index <- match(df$gene, df$gene)
#     hit <- match(df$gene, eval$id)
#     df$index <- eval$max_tx_index[hit] + start_index - 1
#     out <- grep(pattern, df$tx[df$index], value = TRUE)
#     return(out)
# }
#
# .removeInvalidMP <- function(gff, valid_id, pattern){
#     parent <- sapply(gff$Parent, function(x){
#         if(length(x) == 0){
#             x <- NA
#         }
#         return(x)
#     })
#     valid_tx <- gff$ID %in% valid_id
#     valid_elemetns <- parent %in% valid_id
#     non_mp_entries <- !grepl(pattern, gff$ID)
#     valid_entries <- valid_tx | valid_elemetns | non_mp_entries
#     return(gff[valid_entries])
# }
#

################################################################################
#'
#' #' Organize GFF Data
#' #'
#' #' This function organizes and validates GFF data, filtering transcripts and organizing genes, transcripts, and CDS features.
#' #'
#' #' @importFrom BiocGenerics start
#' .orgGFF <- function(gff1, gff2, gff3, prefix, overlap = FALSE, len_diff){
#'     # Adjust levels for transcript types
#'     lv <- levels(gff1$type)
#'     lv[lv == "mRNA"] <- "transcript"
#'     levels(gff1$type) <- lv
#'
#'     # Filter relevant types and sort GFF data
#'     gff1 <- gff1[gff1$type %in% c("gene", "transcript", "five_prime_UTR", "CDS", "three_prime_UTR")]
#'     gff1 <- gff1[order(as.numeric(seqnames(gff1)), start(gff1))]
#'     gff2 <- gff2[order(as.numeric(seqnames(gff2)), start(gff2))]
#'
#'     # Validate GFF data
#'     gff_valid <- .validTX(gff1 = gff1,
#'                           gff2 = gff2,
#'                           gff3 = gff3,
#'                           overlap = overlap,
#'                           len_diff = len_diff)
#'
#'     # Organize transcripts and CDS features
#'     gff_tx <- .orgTX(gff = gff_valid, prefix = prefix)
#'
#'     .solveChimericTX(gff_valid = gff_valid)
#'     .solveChimericTX <- function(gff_valid){
#'         gff_valid_cds_i <- gff_valid$type %in% "CDS"
#'         ol <- findOverlaps(gff_valid[gff_valid_cds_i],
#'                            gff_valid[gff_valid_cds_i])
#'         ol <- as.data.frame(ol)
#'         ol <- subset(ol, subset = queryHits != subjectHits)
#'         ol_sort <- apply(ol, 1, sort)
#'         ol <- subset(ol, subset = !duplicated(t(ol_sort)))
#'         ol$query_tx <- unlist(gff_valid$Parent[gff_valid_cds_i][ol$queryHits])
#'         ol$subject_tx <- unlist(gff_valid$Parent[gff_valid_cds_i][ol$subjectHits])
#'     }
#'
#'     gff_cds <- .orgCDS(gff1 = gff1, gff_tx = gff_tx)
#'     gff_gene <- .orgGene(gff_tx = gff_tx)
#'
#'     # Remove old IDs from transcripts and genes
#'     gff_tx$old_id <- gff_gene$old_id <- NULL
#'
#'     # Combine genes, transcripts, and CDS features
#'     out <- c(gff_gene, gff_tx, gff_cds)
#'
#'     # Sort the final output
#'     out <- out[order(as.numeric(seqnames(out)), start(out), as.numeric(out$type))]
#'
#'     return(out)
#' }
#' #' Validate Transcripts
#' #'
#' #' This function validates transcripts by checking overlaps and length differences with reference GFF data.
#' #'
#' #' @importFrom S4Vectors queryHits
#' #' @importFrom GenomicRanges findOverlaps
#' .validTX <- function(gff1, gff2, gff3, overlap, len_diff){
#'     if(overlap){
#'         # Get unique CDS blocks from gff1
#'         gff1_block <- .getCDSblock(gff = gff1)
#'         gff1_block_uniq <- gff1_block[!duplicated(gff1_block)]
#'
#'         # Get CDS blocks from gff2 and remove overlaps from gff1_block_uniq
#'         gff2_block <- .getCDSblock(gff = gff2)
#'         gff1_block_uniq <- gff1_block_uniq[!gff1_block_uniq %in% gff2_block]
#'
#'         # Get unique transcripts from gff1
#'         out <- .getUniqTx(gff1 = gff1, gff1_block_uniq = gff1_block_uniq)
#'
#'     } else {
#'         # Find overlaps between gff1 and gff2 transcripts
#'         tx_i1 <- gff1$type %in% c("transcript", "mRNA")
#'         tx_i2 <- gff2$type %in% c("transcript", "mRNA")
#'         ol <- findOverlaps(gff1[tx_i1], gff2[tx_i2])
#'         hit <- unique(queryHits(ol))
#'
#'         # Exclude overlapping transcripts
#'         out <- gff1[tx_i1][!seq_len(sum(tx_i1)) %in% hit]
#'     }
#'
#'     # Check length differences with gff3
#'     out <- .checkLength(gff = out, gff3 = gff3, len_diff = len_diff)
#'
#'
#'     return(out)
#' }
#' #' Get CDS Blocks
#' #'
#' #' This function extracts the concatenated CDS blocks for each transcript in the GFF data.
#' #'
#' #' @importFrom BiocGenerics start end
#' .getCDSblock <- function(gff){
#'     # Identify CDS and transcript indices
#'     gff_cds_i <- gff$type == "CDS"
#'     gff_tx_i <- gff$type %in% c("transcript", "mRNA")
#'
#'     # Extract start and end positions of CDS
#'     gff_cds_start <- start(gff[gff_cds_i])
#'     gff_cds_end <- end(gff[gff_cds_i])
#'
#'     # Create exon strings combining start and end positions
#'     gff_cds_exon <- paste(gff_cds_start, gff_cds_end, sep = "-")
#'
#'     # Map CDS to their parent transcripts
#'     map_to_tx <- match(unlist(gff$Parent[gff_cds_i]), gff$ID[gff_tx_i])
#'
#'     # Identify the first occurrence of each transcript
#'     first_i <- !duplicated(map_to_tx)
#'
#'     # Get the chromosome names for CDS
#'     gff_cds_chr <- as.character(seqnames(gff[gff_cds_i]))
#'
#'     # Add chromosome information to the first occurrence of each transcript's CDS block
#'     gff_cds_exon[first_i] <- paste(gff_cds_chr[first_i], gff_cds_exon[first_i], sep = ":")
#'
#'     # Concatenate CDS blocks for each transcript
#'     out <- tapply(gff_cds_exon, map_to_tx, paste, collapse = ",")
#'
#'     return(out)
#' }
#' #' Get Unique Transcripts
#' #'
#' #' This function extracts unique transcripts and their associated elements (e.g., CDS) from a GFF object.
#' #'
#' #' @importFrom BiocGenerics unlist
#' .getUniqTx <- function(gff1, gff1_block_uniq){
#'     # Identify transcript indices
#'     gff1_tx_i <- gff1$type %in% c("transcript", "mRNA")
#'
#'     # Extract unique transcripts using the unique CDS blocks
#'     out_tx <- gff1[gff1_tx_i][as.numeric(names(gff1_block_uniq))]
#'
#'     # Identify and extract elements associated with the unique transcripts
#'     out_element <- unlist(gff1$Parent[!gff1_tx_i]) %in% out_tx$ID
#'     out_element <- gff1[!gff1_tx_i][out_element]
#'
#'     # Combine unique transcripts and their associated elements
#'     out <- c(out_tx, out_element)
#'
#'     return(out)
#' }
#' #' Check Length of Transcripts
#' #'
#' #' This function checks the length of transcripts and filters those that meet the length difference criteria.
#' #'
#' #' @importFrom BiocGenerics unlist
#' .checkLength <- function(gff, gff3, len_diff){
#'     # Calculate CDS lengths for the given GFF data
#'     gff_cds_len <- .getCDSlen(gff = gff)
#'     index <- as.numeric(names(gff_cds_len))
#'     gff_tx_i <- gff$type %in% c("transcript", "mRNA")
#'     names(gff_cds_len) <- sub("\\s.+", "", gff$Target[gff_tx_i][index])
#'
#'     # Calculate CDS lengths for the reference GFF data
#'     gff3_cds_len <- .getCDSlen(gff = gff3)
#'     index <- as.numeric(names(gff3_cds_len))
#'     gff3_tx_i <- gff3$type %in% c("transcript", "mRNA")
#'     names(gff3_cds_len) <- gff3$ID[gff3_tx_i][index]
#'
#'     # Match CDS lengths by transcript IDs
#'     id_hit <- match(names(gff_cds_len), names(gff3_cds_len))
#'     gff3_cds_len <- gff3_cds_len[id_hit]
#'
#'     # Determine which transcripts are longer
#'     longer <- gff_cds_len
#'     is_gff3_longer <- longer < gff3_cds_len
#'     longer[is_gff3_longer] <- gff3_cds_len[is_gff3_longer]
#'
#'     # Calculate valid transcripts based on length difference criteria
#'     valid <- abs(gff_cds_len - gff3_cds_len) / longer <= len_diff
#'     valid_tx <- gff[gff_tx_i][as.vector(valid)]
#'
#'     # Extract non-transcript elements associated with valid transcripts
#'     non_tx_gff <- gff[!gff_tx_i]
#'     out <- c(valid_tx, non_tx_gff[unlist(non_tx_gff$Parent) %in% valid_tx$ID])
#'
#'     return(out)
#' }
#'
#' #' Calculate CDS Lengths
#' #'
#' #' This function calculates the lengths of coding sequences (CDS) for each transcript in a GFF file.
#' #'
#' #' @importFrom BiocGenerics start end
#' .getCDSlen <- function(gff){
#'     # Identify indices for CDS and transcript elements
#'     gff_cds_i <- gff$type == "CDS"
#'     gff_tx_i <-  gff$type %in% c("transcript", "mRNA")
#'
#'     # Extract start and end positions for CDS elements
#'     gff_cds_start <- start(gff[gff_cds_i])
#'     gff_cds_end <- end(gff[gff_cds_i])
#'
#'     # Calculate the length of each CDS element
#'     gff_cds_len <- gff_cds_end - gff_cds_start
#'
#'     # Map CDS elements to their corresponding transcripts
#'     map_to_tx <- match(unlist(gff$Parent[gff_cds_i]), gff$ID[gff_tx_i])
#'
#'     # Sum the lengths of CDS elements for each transcript
#'     out <- tapply(gff_cds_len, map_to_tx, sum)
#'
#'     return(out)
#' }
#'
#' #' Organize Transcripts
#' #'
#' #' This function organizes transcripts in a GFF file by assigning gene IDs and updating transcript IDs.
#' #'
#' #' @importFrom BiocGenerics start
#' #' @importFrom GenomeInfoDb seqnames
#' #' @importFrom GenomicRanges reduce findOverlaps
#' #' @importFrom S4Vectors subjectHits
#' .orgTX <- function(gff, prefix){
#'     # Identify indices for transcript elements
#'     tx_i <- gff$type %in% c("transcript", "mRNA")
#'
#'     # Extract and order transcript elements by chromosome and start position
#'     gff_tx <- gff[tx_i]
#'     gff_tx <- gff_tx[order(as.numeric(seqnames(gff_tx)), start(gff_tx))]
#'
#'     # Remove duplicate transcripts
#'     gff_tx <- unique(gff_tx)
#'
#'     # Reduce transcripts to unique loci
#'     uni_loci <- reduce(gff_tx)
#'
#'     # Map transcripts to unique loci
#'     map_to_loci <- findOverlaps(gff_tx, uni_loci)
#'
#'     # Assign gene IDs based on unique loci
#'     gff_tx$gene_id <- paste0(prefix, "G", sprintf("%05d", subjectHits(map_to_loci)))
#'
#'     # Store old transcript IDs
#'     gff_tx$old_id <- gff_tx$ID
#'
#'     # Order transcripts by gene ID
#'     gff_tx <- gff_tx[order(gff_tx$gene_id)]
#'
#'     # Update transcript IDs based on gene ID
#'     gff_tx$ID <- unlist(tapply(gff_tx$gene_id, gff_tx$gene_id, function(x){
#'         return(paste0(sub("G", "T", x[1]), ".", sprintf("%02d", seq_len(length(x)))))
#'     }))
#'
#'     # Assign parent gene IDs to transcripts
#'     gff_tx$Parent <- lapply(gff_tx$gene_id, c)
#'
#'     return(gff_tx)
#' }
#' #' Organize CDS Elements
#' #'
#' #' This function organizes CDS elements in a GFF file by assigning gene IDs and updating parent IDs.
#' #'
#' .orgCDS <- function(gff1, gff_tx){
#'     # Identify indices for non-transcript elements
#'     cds_i <- !gff1$type %in% c("transcript", "mRNA")
#'
#'     # Extract CDS elements that have parents in the organized transcript data
#'     gff_cds <- gff1[cds_i][unlist(gff1$Parent[cds_i]) %in% gff_tx$old_id]
#'
#'     # Assign gene IDs based on organized transcript data
#'     gff_cds$gene_id <- gff_tx$gene_id[match(unlist(gff_cds$Parent), gff_tx$old_id)]
#'
#'     # Update parent IDs based on organized transcript data
#'     gff_cds$Parent <- lapply(gff_tx$ID[match(unlist(gff_cds$Parent), gff_tx$old_id)], c)
#'
#'     # Update CDS IDs based on parent IDs
#'     gff_cds$ID <- paste0(unlist(gff_cds$Parent), ":CDS")
#'
#'     return(gff_cds)
#' }
#'
#' #' Organize Gene Elements
#' #'
#' #' This function organizes gene elements in a GFF file by creating gene entries from the organized transcript data.
#' #'
#' .orgGene <- function(gff_tx){
#'     # Extract unique gene entries from the transcript data
#'     gff_gene <- gff_tx[!duplicated(gff_tx$gene_id)]
#'
#'     # Set the type of each entry to "gene"
#'     gff_gene$type <- "gene"
#'
#'     # Assign gene IDs to the ID field
#'     gff_gene$ID <- gff_gene$gene_id
#'
#'     # Remove parent information for gene entries
#'     gff_gene$Parent <- lapply(seq_along(gff_gene), function(i) character())
#'
#'     return(gff_gene)
#' }
#' #' Merge GFF Annotations
#' #'
#' #' This function merges GFF annotations by combining two GFF files, removing duplicate gene and transcript entries, and creating a unified annotation.
#' #'
#' #' @importFrom BiocGenerics start
#' #' @importFrom GenomeInfoDb seqnames
#' #' @importFrom GenomicRanges reduce findOverlaps
#' #' @importFrom S4Vectors subjectHits
#' .mergeGFF <- function(gff1, gff2 = NULL, ann_priority){
#'     # Combine the two GFF objects
#'     gff <- c(gff1, gff2)
#'
#'     # Filter for relevant types
#'     gff <- gff[gff$type %in% c("gene", "transcript", "mRNA",
#'                                "five_prime_UTR", "CDS", "three_prime_UTR")]
#'
#'     # Order and reduce gene entries
#'     gene_i <- gff$type == "gene"
#'     gff_gene <- gff[gene_i]
#'     gff_gene <- gff_gene[order(as.numeric(seqnames(gff_gene)), start(gff_gene))]
#'     uni_loci <- reduce(gff_gene)
#'     map_to_loci <- findOverlaps(gff_gene, uni_loci)
#'
#'     # Identify overlapping entries and create ID map
#'     ol_list <- tapply(queryHits(map_to_loci), subjectHits(map_to_loci), c)
#'     ol_list <- ol_list[sapply(ol_list, length) > 1]
#'     id_map <- .getIDmap(gff_gene = gff_gene, ol_list = ol_list)
#'
#'     # Map IDs in the GFF object based on the ID map
#'     gff <- .mapID(gff = gff, id_map = id_map)
#'
#'     # Order the final merged GFF object
#'     gff <- gff[order(as.numeric(seqnames(gff)), start(gff), as.numeric(gff$type))]
#'
#'     return(gff)
#' }
#' #' Create ID Map for Gene Overlaps
#' #'
#' #' This function creates an ID map for overlapping genes by selecting a representative gene ID for each set of overlapping genes.
#' #'
#' .getIDmap <- function(gff_gene, ol_list){
#'     # Create ID map for each set of overlapping genes
#'     out <- lapply(ol_list, function(i){
#'         gene_id <- gff_gene$ID[i]
#'         # Select the representative gene ID that does not start with 'query_' or 'subject_'
#'         set_gene_id <- gene_id[!grepl("^query_|^subject_", gene_id)][1]
#'         id_map <- cbind(gene_id, set_gene_id)
#'         return(id_map)
#'     })
#'     # Combine all ID maps into a single data frame
#'     out <- do.call("rbind", out)
#'     return(out)
#' }
#'
#'
#' .mapID <- function(gff, id_map){
#'     gff_parent <- sapply(gff$Parent, function(x){
#'         if(length(x) == 0){
#'             return(NA)
#'         } else {
#'             return(x)
#'         }
#'     })
#'     hit <- match(gff_parent, id_map[, 1])
#'     gff$Parent[!is.na(hit)] <- lapply(id_map[na.omit(hit), 2], c)
#'     gff$gene_id[!is.na(hit)] <- id_map[na.omit(hit), 2]
#'
#'     rm_gene <- id_map[id_map[, 1] != id_map[, 2], ]
#'     gff <- gff[!gff$ID %in% rm_gene[, 1]]
#'     return(gff)
#' }
#'
#'
#' #' Fix and Standardize GFF Annotations
#' #'
#' #' This function applies a series of fixes and standardizations to GFF annotations.
#' #'
#' .fixGFF <- function(gff){
#'     gff <- .setGeneID(gff = gff)
#'     gff <- .fixGFFrange(gff = gff)
#'     gff <- .fixGFFexon(gff = gff)
#'     gff <- .fixGFFphase(gff = gff)
#'     return(gff)
#' }
#'
#' #' Set Gene IDs for GFF Annotations
#' #'
#' #' This function sets the gene_id attribute for GFF annotations based on their parent-child relationships.
#' #'
#' .setGeneID <- function(gff){
#'     # For transcripts or mRNA, set gene_id to their ID
#'     tx_i <- gff$type %in% c("transcript", "mRNA")
#'     hit <- match(unlist(gff$Parent[tx_i]), gff$ID)
#'     gff$gene_id[tx_i] <- gff$ID[hit]
#'
#'     # For other elements, set gene_id based on their parent transcript/mRNA
#'     element_i <- !gff$type %in% c("gene", "transcript", "mRNA")
#'     hit <- match(unlist(gff$Parent[element_i]), gff$ID[tx_i])
#'     gff$gene_id[element_i] <- gff$gene_id[tx_i][hit]
#'     return(gff)
#' }
#'
#' #' Fix GFF Range
#' #'
#' #' This function adjusts the start and end positions of transcripts and genes to cover the complete range of their member elements.
#' #'
#' #' @importFrom rtracklayer start end start<- end<-
#' .fixGFFrange <- function(gff){
#'     # Fix the start and end positions of each transcript to
#'     # cover the whole range of member elements (CDS, exon, and UTRs)
#'     tx_i <- which(gff$type %in% c("transcript", "mRNA"))
#'     element_i <- !gff$type %in% c("gene", "transcript", "mRNA")
#'
#'     # Determine minimum start and maximum end positions for member elements
#'     min_start <- tapply(start(gff[element_i]), unlist(gff$Parent[element_i]), min)
#'     max_end <- tapply(end(gff[element_i]), unlist(gff$Parent[element_i]), max)
#'
#'     # Match and update transcript start and end positions
#'     hit <- match(gff$ID[tx_i], names(min_start))
#'     tx_start <- min_start[hit]
#'     tx_end <- max_end[hit]
#'     not_na_start <- !is.na(tx_start)
#'     not_na_end <- !is.na(tx_end)
#'     start(gff[tx_i[not_na_start]]) <- tx_start[not_na_start]
#'     end(gff[tx_i[not_na_end]]) <- tx_end[not_na_end]
#'
#'     # Fix the start and end positions of each gene to
#'     # cover the whole range of member transcripts
#'     gene_i <- which(gff$type == "gene")
#'     tx_i <- which(gff$type %in% c("transcript", "mRNA"))
#'
#'     # Determine minimum start and maximum end positions for transcripts
#'     min_start <- tapply(start(gff[tx_i]), gff$gene_id[tx_i], min)
#'     max_end <- tapply(end(gff[tx_i]), gff$gene_id[tx_i], max)
#'
#'     # Match and update gene start and end positions
#'     hit <- match(gff$gene_id[gene_i], names(min_start))
#'     gene_start <- min_start[hit]
#'     gene_end <- max_end[hit]
#'     start(gff[gene_i]) <- gene_start
#'     end(gff[gene_i]) <- gene_end
#'
#'     return(gff)
#' }
#'
#'
#' .fixGFFexon <- function(gff){
#'     gff <- gff[gff$type != "exon"]
#'     gff_exon <- gff[gff$type %in% c("CDS", "five_prime_UTR", "three_prime_UTR")]
#'     gff_exon$type <- "exon"
#'     gff_exon$Name <- gff_exon$ID <- paste0(unlist(gff_exon$Parent), ":exon")
#'     gff_exon$score <- gff_exon$phase <- NA
#'
#'     tx_i <- gff$type %in% c("transcript", "mRNA")
#'     non_cds <- gff[tx_i][!gff$ID[tx_i] %in% unlist(gff_exon$Parent)]
#'     if(length(non_cds) > 0){
#'         non_cds$type <- "exon"
#'         non_cds$Parent <- lapply(non_cds$ID, c)
#'         non_cds$Name <- non_cds$ID <- paste0(non_cds$ID, ":exon")
#'         non_cds$score <- non_cds$phase <- NA
#'         gff_exon <- c(gff_exon, non_cds)
#'     }
#'
#'     out <- c(gff, gff_exon)
#'     out$type <- factor(out$type, levels = c("gene", "transcript", "mRNA",
#'                                             "five_prime_UTR", "exon",
#'                                             "CDS", "three_prime_UTR"))
#'     out$type <- droplevels(out$type)
#'     out <- out[order(as.numeric(seqnames(out)), start(out), as.numeric(out$type))]
#'     return(out)
#' }
#'
#'
#' #' Fix GFF Phase
#' #'
#' #' This function adjusts the phase of CDS features in GFF annotations.
#' #'
#' #' @importFrom BiocGenerics start
#' .fixGFFphase <- function(gff){
#'     # Extract CDS features
#'     gff_cds <- gff[gff$type == "CDS"]
#'
#'     # Adjust phase for plus strand
#'     gff_cds_plus <- .phasePlus(gff_cds = gff_cds)
#'
#'     # Adjust phase for minus strand
#'     gff_cds_minus <- .phaseMinus(gff_cds = gff_cds)
#'
#'     # Combine non-CDS features with adjusted CDS features
#'     out <- c(gff[gff$type != "CDS"], gff_cds_plus, gff_cds_minus)
#'
#'     # Set factor levels for feature types
#'     out$type <- factor(out$type, levels = c("gene", "transcript", "mRNA",
#'                                             "five_prime_UTR", "exon",
#'                                             "CDS", "three_prime_UTR"))
#'     out$type <- droplevels(out$type)
#'
#'     # Order the features by sequence name, start position, and type
#'     out <- out[order(as.numeric(seqnames(out)), start(out), as.numeric(out$type))]
#'     return(out)
#' }
#' #' Adjust Phase for CDS Features on Plus Strand
#' #'
#' #' This function adjusts the phase for CDS features on the plus strand in GFF annotations.
#' #'
#' #' @importFrom BiocGenerics start width
#' #' @importFrom GenomeInfoDb seqnames
#' .phasePlus <- function(gff_cds){
#'     # Filter CDS features on the plus strand
#'     gff_cds_plus <- gff_cds[as.character(strand(gff_cds)) == "+"]
#'     gff_cds_plus <- gff_cds_plus[order(as.numeric(seqnames(gff_cds_plus)), start(gff_cds_plus))]
#'
#'     # Extract parent IDs for CDS features
#'     gff_cds_plus_parent <- unlist(gff_cds_plus$Parent)
#'
#'     # Initialize phase for the first CDS in each transcript
#'     target_i <- which(!duplicated(gff_cds_plus_parent))
#'     gff_cds_plus$phase[target_i] <- 0
#'
#'     # Calculate the phase for the next CDS feature
#'     next_phase <- (3 - (width(gff_cds_plus[target_i]) - gff_cds_plus$phase[target_i]) %% 3) %% 3
#'     names(next_phase) <- gff_cds_plus_parent[target_i]
#'     gff_cds_plus_parent[target_i] <- "NA"
#'
#'     # Iterate over the remaining CDS features to set their phases
#'     while(TRUE){
#'         target_i <- which(!duplicated(gff_cds_plus_parent))[-1]
#'         if(length(target_i) == 0){
#'             break
#'         }
#'         next_phase <- next_phase[names(next_phase) %in% gff_cds_plus_parent[target_i]]
#'         gff_cds_plus$phase[target_i] <- next_phase
#'         next_phase <- (3 - (width(gff_cds_plus[target_i]) - gff_cds_plus$phase[target_i]) %% 3) %% 3
#'         names(next_phase) <- gff_cds_plus_parent[target_i]
#'         gff_cds_plus_parent[target_i] <- "NA"
#'     }
#'
#'     return(gff_cds_plus)
#' }
#'
#' #' Adjust Phase for CDS Features on Minus Strand
#' #'
#' #' This function adjusts the phase for CDS features on the minus strand in GFF annotations.
#' #'
#' #' @importFrom BiocGenerics end width
#' #' @importFrom GenomeInfoDb seqnames
#' .phaseMinus <- function(gff_cds){
#'     # Filter CDS features on the minus strand
#'     gff_cds_minus <- gff_cds[as.character(strand(gff_cds)) == "-"]
#'     gff_cds_minus <- gff_cds_minus[order(as.numeric(seqnames(gff_cds_minus)),
#'                                          end(gff_cds_minus), decreasing = TRUE)]
#'
#'     # Extract parent IDs for CDS features
#'     gff_cds_minus_parent <- unlist(gff_cds_minus$Parent)
#'
#'     # Initialize phase for the first CDS in each transcript
#'     target_i <- which(!duplicated(gff_cds_minus_parent))
#'     gff_cds_minus$phase[target_i] <- 0
#'
#'     # Calculate the phase for the next CDS feature
#'     next_phase <- (3 - (width(gff_cds_minus[target_i]) - gff_cds_minus$phase[target_i]) %% 3) %% 3
#'     names(next_phase) <- gff_cds_minus_parent[target_i]
#'     gff_cds_minus_parent[target_i] <- "NA"
#'
#'     # Iterate over the remaining CDS features to set their phases
#'     while(TRUE){
#'         target_i <- which(!duplicated(gff_cds_minus_parent))[-1]
#'         if(length(target_i) == 0){
#'             break
#'         }
#'         next_phase <- next_phase[names(next_phase) %in% gff_cds_minus_parent[target_i]]
#'         gff_cds_minus$phase[target_i] <- next_phase
#'         next_phase <- (3 - (width(gff_cds_minus[target_i]) - gff_cds_minus$phase[target_i]) %% 3) %% 3
#'         names(next_phase) <- gff_cds_minus_parent[target_i]
#'         gff_cds_minus_parent[target_i] <- "NA"
#'     }
#'
#'     return(gff_cds_minus)
#' }
#'
#' ################################################################################
#'

#'
#' #' Create CDS sequences from GFF and genome files
#' #'
#' #' This function generates CDS sequences from provided GFF and genome files.
#' #'
#' #' @param gff Path to the GFF file.
#' #' @param genome Path to the genome file.
#' #'
#' #' @return A DNAStringSet object containing CDS sequences.
#' #' @importFrom Biostrings readDNAStringSet
#' #' @importFrom GenomicFeatures cdsBy extractTranscriptSeqs
#' #' @importFrom txdbmaker makeTxDbFromGFF
#' #' @import BSgenome
#' #'
#' .makeCDS <- function(gff, genome){
#'     # Create a TxDb object from the GFF file
#'     txdb <- makeTxDbFromGFF(file = gff)
#'
#'     # Read the genome file as a DNAStringSet object
#'     genome <- readDNAStringSet(filepath = genome)
#'
#'     # Extract CDS sequences from the TxDb object
#'     cds_db <- cdsBy(x = txdb, by = "tx", use.names = TRUE)
#'     cds <- extractTranscriptSeqs(x = genome, transcripts = cds_db)
#'
#'     # Order CDS sequences by their names
#'     cds <- cds[order(names(cds))]
#'     return(cds)
#' }
#'
#' #' Update OrthoPairDB files
#' #'
#' #' This function updates the GFF and CDS files in a OrthoPairDB object from the HDF5 file.
#' #'
#' #' @param object A OrthoPairDB object.
#' #'
#' #' @return The updated OrthoPairDB object.
#' #' @importFrom rhdf5 H5Fopen H5Fclose
#' .updateFiles <- function(object){
#'     h5 <- H5Fopen(object$h5)
#'     on.exit(H5Fclose(h5))
#'
#'     # Update GFF and CDS file paths from the HDF5 file
#'     object$query_gff <- as.vector(h5$miniprot$s2q_gff)
#'     object$subject_gff <- as.vector(h5$miniprot$q2s_gff)
#'     object$query_cds <- as.vector(h5$miniprot$s2q_cds)
#'     object$subject_cds <- as.vector(h5$miniprot$q2s_cds)
#'     return(object)
#' }
#'
#' # .makeLCBgr <- function(object, h5){
#' #     # Order the LCB pairs
#' #     pairs <- list(lcb_1to1 = .order(df = h5$sibeliaz$lcb_pairs$lcb_1to1),
#' #                   lcb_non_1to1 = .order(df = h5$sibeliaz$lcb_pairs$lcb_non_1to1))
#' #
#' #     # Convert data frames to GRanges objects
#' #     gr <- .df2gr(pairs = pairs)
#' #
#' #     # Create GRanges objects for gaps in the query and subject genomes
#' #     gap_gr <- list(query = .gapGR(gr = gr$query_1to1,
#' #                                   chrLen = object$genome$query),
#' #                    subject = .gapGR(gr = gr$subject_1to1,
#' #                                     chrLen = object$genome$subject))
#' #     out <- list(gr = gr,
#' #                 gap_gr = gap_gr)
#' #     return(out)
#' # }
#' #
#' # .remakeAnchor <- function(gff_ls, lcb_gr, h5){
#' #     df1 <- .makeMiniprotPairDF(gff1 = gff_ls$query_gff,
#' #                                gff2 = gff_ls$subject_gff)
#' #     df1 <- .filterMPtx(df = df1, gff_ls = gff_ls, query = TRUE)
#' #
#' #     df2 <- .makeMiniprotPairDF(gff1 = gff_ls$subject_gff,
#' #                                gff2 = gff_ls$query_gff)
#' #     df2 <- .filterMPtx(df = df2, gff_ls = gff_ls, query = FALSE)
#' #
#' #     new_anchor <- .getNewAnchors(df1 = df1,
#' #                                  df2 = df2,
#' #                                  h5 = h5,
#' #                                  lcb_gr = lcb_gr,
#' #                                  gff_ls = gff_ls)
#' #
#' #     out <- rbind(subset(h5$anchor, select = c(qseqid, sseqid)),
#' #                  subset(new_anchor, select = c(qseqid, sseqid)))
#' #
#' #     # Overwrite the "anchor" group in the HDF5 file with the filtered orthologs
#' #     .h5overwrite(obj = out, file = object$h5, "anchor")
#' # }
#' #
#' # .filterMPtx <- function(df, gff_ls, query = TRUE){
#' #     if(query){
#' #         df <- subset(df, subset = grepl("query_", gene))
#' #         df$tx_single <- .isSingleExon(tx = df$tx,
#' #                                       gff = gff_ls$query_gff)
#' #         df$target_tx_single <- .isSingleExon(tx = df$target_tx,
#' #                                              gff = gff_ls$subject_gff)
#' #         df <- subset(df, subset = !tx_single & !target_tx_single)
#' #
#' #     } else {
#' #         df <- subset(df, subset = grepl("subject_", gene))
#' #         df$tx_single <- .isSingleExon(tx = df$tx,
#' #                                       gff = gff_ls$subject_gff)
#' #         df$target_tx_single <- .isSingleExon(tx = df$target_tx,
#' #                                              gff = gff_ls$query_gff)
#' #         df <- subset(df, subset = !tx_single & !target_tx_single)
#' #     }
#' #     return(df)
#' # }
#' #
#' # .isSingleExon <- function(tx, gff){
#' #     cds_i <- gff$type == "CDS"
#' #     query_cds_parents <- unlist(gff$Parent[cds_i])
#' #     n_query_cds_parents <- table(query_cds_parents)
#' #     hit <- match(tx, names(n_query_cds_parents))
#' #     out <- n_query_cds_parents[hit] == 1
#' #     return(out)
#' # }
#' #
#' # .getNewAnchors <- function(df1, df2, h5, lcb_gr, gff_ls){
#' #     df1 <- .omitPairedGenes(df = df1,
#' #                             geneid = h5$orthopair_gene$orthopairs$sgeneid,
#' #                             syntenic = h5$orthopair_gene$orthopairs$syntenic)
#' #     df2 <- .omitPairedGenes(df = df2,
#' #                             geneid = h5$orthopair_gene$orthopairs$qgeneid,
#' #                             syntenic = h5$orthopair_gene$orthopairs$syntenic)
#' #
#' #     rbbh <- data.frame(qseqid = c(df1$tx, df2$target_tx),
#' #                        sseqid = c(df1$target_tx, df2$tx))
#' #     rbbh$index <- seq_len(nrow(rbbh))
#' #
#' #     obj <- list(gr = lcb_gr$gr,
#' #                 gap_gr = lcb_gr$gap_gr,
#' #                 rbbh = rbbh,
#' #                 gff = list(query = gff_ls$query_gff,
#' #                            subject = gff_ls$subject_gff))
#' #
#' #     # Filter orthologs in 1-to-1 LCBs
#' #     obj <- .orthoIn1to1lcb(obj = obj)
#' #
#' #     # Filter orthologs in non-1-to-1 LCBs if specified
#' #     obj <- .orthoInNon1to1lcb(obj = obj)
#' #
#' #     # Filter orthologs in 1-to-1 CBI
#' #     obj <- .orthoIn1to1cbi(obj = obj)
#' #
#' #     return(obj$out)
#' # }
#' #
#' # .omitPairedGenes <- function(df, geneid, syntenic){
#' #     hit <- match(df$target_gene, geneid)
#' #     df$paired[!is.na(hit)] <- TRUE
#' #     df$syntenic <- as.logical(syntenic[hit])
#' #     df$orthopair <- df$paired & df$syntenic
#' #     df$orthopair[is.na(df$orthopair)] <- FALSE
#' #     df <- subset(df, subset = !orthopair)
#' #     return(df)
#' # }
#'
#'
#' .makeQueryDF <- function(query_gff, subject_gff, query_cds){
#'     query_cds_init <- .checkInitCodon(cds = query_cds)
#'     query_cds_term <- .checkTermCodon(cds = query_cds)
#'     query_gff_tx_index <- query_gff$type %in% c("transcript", "mRNA")
#'     query_df <- data.frame(tx = query_gff$ID[query_gff_tx_index],
#'                            gene = query_gff$gene_id[query_gff_tx_index])
#'     hit <- match(query_df$tx, names(query_cds_init))
#'     query_df$init <- query_cds_init[hit]
#'     query_df$term <- query_cds_term[hit]
#'     query_df$valid <- query_df$init & query_df$term
#'     query_df$len <- width(query_cds)[hit]
#'     hit <- match(query_df$tx, query_gff$ID)
#'     query_df$target <- sub("\\s.+", "", query_gff$Target[hit])
#'     hit <- match(query_df$target, subject_gff$ID)
#'     query_df$target <- subject_gff$gene_id[hit]
#'     query_df$pair_id <- paste(query_df$gene, query_df$target, sep = "_")
#'     return(query_df)
#' }
#'
#' .makeSubjectDF <- function(subject_gff, query_gff, subject_cds){
#'     subject_cds_init <- .checkInitCodon(cds = subject_cds)
#'     subject_cds_term <- .checkTermCodon(cds = subject_cds)
#'     subject_gff_tx_index <- subject_gff$type %in% c("transcript", "mRNA")
#'     subject_df <- data.frame(tx = subject_gff$ID[subject_gff_tx_index],
#'                              gene = subject_gff$gene_id[subject_gff_tx_index])
#'     hit <- match(subject_df$tx, names(subject_cds_init))
#'     subject_df$init <- subject_cds_init[hit]
#'     subject_df$term <- subject_cds_term[hit]
#'     subject_df$valid <- subject_df$init & subject_df$term
#'     subject_df$len <- width(subject_cds)[hit]
#'     hit <- match(subject_df$tx, subject_gff$ID)
#'     subject_df$target <- sub("\\s.+", "", subject_gff$Target[hit])
#'     hit <- match(subject_df$target, query_gff$ID)
#'     subject_df$target <- query_gff$gene_id[hit]
#'     subject_df$pair_id <- paste(subject_df$gene, subject_df$target, sep = "_")
#'     return(subject_df)
#' }
#'
#' .checkInitCodon <- function(cds){
#'     return(substr(cds, 1, 3) == "ATG")
#' }
#'
#' .checkTermCodon <- function(cds){
#'     len <- width(cds)
#'     return(substr(cds, len - 2, len) %in% c("TAA", "TGA", "TAG"))
#' }
#'
#' .filterMPgenes <- function(df, pattern = "^query_"){
#'     df <- df[order(df$gene), ]
#'     mp_tx_index <- grepl(pattern, df$tx)
#'     eval <- data.frame(id = unique(df$gene))
#'
#'     mp_valid_max_len <- tapply(df$len[mp_tx_index][df$valid[mp_tx_index]],
#'                                df$gene[mp_tx_index][df$valid[mp_tx_index]],
#'                                max)
#'     hit <- match(eval$id, names(mp_valid_max_len))
#'     eval$mp_valid_max_len <- mp_valid_max_len[hit]
#'     eval$mp_valid_max_len[is.na(eval$mp_valid_max_len)] <- 0
#'
#'     non_mp_valid_max_len <- tapply(df$len[!mp_tx_index][df$valid[!mp_tx_index]],
#'                                    df$gene[!mp_tx_index][df$valid[!mp_tx_index]],
#'                                    max)
#'     hit <- match(eval$id, names(non_mp_valid_max_len))
#'     eval$non_mp_valid_max_len <- non_mp_valid_max_len[hit]
#'     eval$non_mp_valid_max_len[is.na(eval$non_mp_valid_max_len)] <- 0
#'
#'     max_tx_index <- tapply(df$len[df$valid],
#'                            df$gene[df$valid],
#'                            which.max)
#'     max_tx_index <- sapply(max_tx_index, function(x){
#'         if(length(x) == 0){
#'             x <- NA
#'         }
#'         return(x)
#'     })
#'     hit <- match(eval$id, names(max_tx_index))
#'     eval$max_tx_index <- max_tx_index[hit]
#'
#'     max_tx_is_na <- is.na(eval$max_tx_index)
#'     hit <- match(grep(pattern, eval$id[max_tx_is_na], value = TRUE),
#'                  df$gene)
#'     max_tx_index <- tapply(df$len[hit],
#'                            df$gene[hit],
#'                            which.max)
#'     hit <- match(eval$id[max_tx_is_na], names(max_tx_index))
#'     eval$max_tx_index[max_tx_is_na] <- max_tx_index[hit]
#'
#'     start_index <- match(df$gene, df$gene)
#'     hit <- match(df$gene, eval$id)
#'     df$index <- eval$max_tx_index[hit] + start_index - 1
#'     out <- grep(pattern, df$tx[df$index], value = TRUE)
#'     return(out)
#' }
#'
#' .removeInvalidMP <- function(gff, valid_id, pattern){
#'     parent <- sapply(gff$Parent, function(x){
#'         if(length(x) == 0){
#'             x <- NA
#'         }
#'         return(x)
#'     })
#'     valid_tx <- gff$ID %in% valid_id
#'     valid_elemetns <- parent %in% valid_id
#'     non_mp_entries <- !grepl(pattern, gff$ID)
#'     valid_entries <- valid_tx | valid_elemetns | non_mp_entries
#'     return(gff[valid_entries])
#' }
#'
