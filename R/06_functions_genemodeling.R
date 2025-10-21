
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
mapProt <- function(in_list,
                    out_dir = "", 
                    miniprot_path = "",
                    overwrite = FALSE,
                    conda_env = NULL,
                    n_threads = 1,
                    len_diff = 0.2){
    if(!inherits(x = in_list, what = "OrthoPairInput")){
        stop("The input object must be a OrthoPairInput class object.",
             call. = FALSE)
    }
    
    miniprot_out_dir <- file.path(out_dir, "miniprot_out")
    dir.create(path = miniprot_out_dir, showWarnings = FALSE, recursive = TRUE)
    
    message("Start running Miniprot...")
    out_files <- .mapEngine(in_list = in_list,
                            miniprot_out_dir = miniprot_out_dir,
                            overwrite = overwrite,
                            miniprot_path = miniprot_path,
                            n_threads = n_threads)
    
    message("Organize and integrate miniprot predicated gene models...")
    
    in_list <- .orgMiniprot(in_list = in_list,
                            out_files = out_files,
                            out_dir = out_dir, 
                            len_diff = len_diff)
    
    in_list <- .createFASTA(in_list = in_list, out_dir = out_dir)
    return(in_list)
}

#' @importFrom rtracklayer export.gff3
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#' @importFrom S4Vectors mcols mcols<-
.mapEngine <- function(in_list,
                       miniprot_out_dir,
                       overwrite,
                       miniprot_path,
                       n_threads){
    if(all(is.na(in_list$genome))){
        stop("No genome information is available.")
    } else {
        omit_genome <- in_list$name[is.na(in_list$genome)]
        if(length(omit_genome) > 0){
            message("Genome sequences are not available for:\n",
                            paste(omit_genome, collapse = "\n"))
        }
    }
    
    if(all(is.na(in_list$prot))){
        stop("No protein information is available.")
        
    } else {
        omit_prot <- in_list$name[is.na(in_list$prot)]
        if(length(omit_prot) > 0){
            message("Protein sequences are not available for:\n",
                            paste(omit_prot, collapse = "\n"))
        }
    }
    
    out_files <- NULL
    for(i in seq_along(in_list$genome)){
        out_files_i <- NULL
        for(j in seq_along(in_list$genome)){
            if(i == j){
                next
            }
            if(is.na(in_list$genome[i]) | is.na(in_list$prot[j])){
                next
            }
            out_fn <- file.path(miniprot_out_dir, 
                                paste0(in_list$name[j], "_mapto_", in_list$name[i]))
            if(file.exists(out_fn) | overwrite){
                # Run miniprot mapping from subject to query genome
                .miniprot(query_fn = in_list$prot[j],
                          genome_fn =  in_list$genome[i],
                          out_prefix = out_fn,
                          miniprot_path = miniprot_path,
                          n_threads = n_threads)
                
                new_cds <- .makeCDS(gff = paste0(out_fn, ".gff"),
                                    genome = in_list$genome[i])
                writeXStringSet(new_cds, paste0(out_fn, ".cds"))
            }
            out_files_i <- c(out_files_i, out_fn)
            names(out_files_i)[length(out_files_i)] <- in_list$name[j]
        }
        out_files <- c(out_files, list(out_files_i))
        names(out_files)[length(out_files)] <- in_list$name[i]
    }
    return(out_files)
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
                      miniprot_path,
                      n_threads = 1){
    
    log_fn <- paste0(out_prefix, ".log1")
    error_log_fn <- paste0(out_prefix, "_error.log1")
    system2(command = file.path(miniprot_path, "miniprot"),
            args = paste("-t", n_threads,
                         "-d", paste0(out_prefix, ".mpi"),
                         genome_fn), 
            stderr = log_fn)
    
    log_fn <- paste0(out_prefix, ".log2")
    error_log_fn <- paste0(out_prefix, "_error.log2")
    system2(command = file.path(miniprot_path, "miniprot"),
            args = paste("-t", n_threads,
                         "--gff",
                         paste0(out_prefix, ".mpi"),
                         query_fn, ">", paste0(out_prefix, ".gff")), 
            stderr = log_fn)
}

.orgMiniprot <- function(in_list, out_files, out_dir, len_diff){
    in_list$update <- rep(FALSE, length(in_list$name))
    for(i in seq_along(out_files)){
        subject_genome <- in_list$name %in% names(out_files)[i]
        subject_gff_fn <- in_list$gff[subject_genome]
        cds_fn <- in_list$cds[subject_genome]
        out_gff <- .filterGFF(in_list = in_list,
                              subject_gff_fn = subject_gff_fn,
                              cds_fn = cds_fn,
                              mp_files = out_files[[i]], 
                              len_diff = len_diff)
        added_tx_id <- out_gff$added_tx_id
        # Fix and filter GFF results, retaining only necessary columns
        out_gff <- .fixGFF(gff = out_gff$out_gff)
        
        m_gff <- mcols(out_gff)
        hit <- names(m_gff) %in% c("source", "type", "score", "phase",
                                   "ID", "Name", "gene_id", "Parent", "Target")
        mcols(out_gff) <- m_gff[, hit]
        out_gff$Name <- out_gff$ID
        
        # Export the final GFF results to output files
        out_fn <- file.path(out_dir, basename(subject_gff_fn))
        export.gff3(out_gff, out_fn)
        in_list$gff[subject_genome] <- out_fn
        in_list$update[subject_genome] <- TRUE
    }
    return(in_list)
}

.filterGFF <- function(in_list, subject_gff_fn, cds_fn, mp_files, len_diff){
    subject_gff <- .getGFFlist(gff_fn = subject_gff_fn)
    subject_cds <- .getCDSlist(cds_fn = cds_fn)
    original_tx_id <- names(subject_cds)
    for(i in seq_along(mp_files)){
        query_genome <- in_list$name %in% names(mp_files)[i]
        query_gff <- .getGFFlist(gff_fn = in_list$gff[query_genome])
        mp_gff <- .getGFFlist(gff_fn = paste0(mp_files[i], ".gff"))
        mp_gff <- .setIDforElements(gff = mp_gff)
        mp_cds <- .getCDSlist(cds_fn = paste0(mp_files[i], ".cds"))
        
        mp_gff <- .filterIdenticalMiniprotTx(original_gff = subject_gff,
                                             mp_gff = mp_gff)
        
        mp_gff <- .filterLongMiniprotTx(mp_gff = mp_gff,
                                        original_gff = query_gff,
                                        len_diff = len_diff)
        
        subject_novel <- .filterMiniprotTxOnNovelLoci(original_gff = subject_gff,
                                                      mp_gff = mp_gff,
                                                      mp_cds = mp_cds)
        
        subject_gff <- suppressWarnings(c(subject_gff, subject_novel))
        subject_cds <- c(subject_cds,
                         mp_cds[names(mp_cds) %in% subject_novel$ID])
        valid_longer_subject_mp_gff <- .filterMiniprotTxOnGeneLoci(original_gff = subject_gff,
                                                                   mp_gff = mp_gff,
                                                                   original_cds = subject_cds,
                                                                   mp_cds = mp_cds)
        
        subject_gff <- suppressWarnings(c(subject_gff, valid_longer_subject_mp_gff))
        subject_cds <- c(subject_cds,
                         mp_cds[names(mp_cds) %in% valid_longer_subject_mp_gff$ID])
    }
    added_tx_id <- names(subject_cds)[names(subject_cds) %in% original_tx_id]
    return(list(out_gff = subject_gff, added_tx_id = added_tx_id))
}

.setIDforElements <- function(gff){
    element_i <- !gff$type %in% c("gene", "transcript", "mRNA")
    gff$ID[element_i] <- paste(unlist(gff$Parent[element_i]),
                               as.character(gff$type[element_i]),
                               sep = ":")
    gff$Name <- gff$ID
    return(gff)
}

.getCDSlist <- function(h5 = NULL, cds_fn = NULL){
    # Import GFF files for query and subject genomes
    if(is.null(cds_fn)){
        query_cds <- readDNAStringSet(as.vector(h5$files$query_cds))
        subject_cds <- readDNAStringSet(as.vector(h5$files$subject_cds))
        
        # Return the ordered cds data as a list
        out <- list(query_cds = query_cds, subject_cds = subject_cds)
        
    } else {
        out <- readDNAStringSet(cds_fn)
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

#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#' @importFrom Biostrings writeXStringSet translate
.createFASTA <- function(in_list, out_dir){
    for(i in seq_along(in_list$name)){
        if(!in_list$update[i]){
            next
        }
        new_cds <- .makeCDS(gff = in_list$gff[i], genome = in_list$genome[i])
        new_cds_fn <- file.path(out_dir, basename(in_list$cds[i]))
        writeXStringSet(new_cds, new_cds_fn)
        in_list$cds[i] <- new_cds_fn
        
        new_prot <- translate(new_cds, if.fuzzy.codon = "solve")
        new_prot <- AAStringSet(sub("\\*.+", "", new_prot))
        new_prot_fn <- file.path(out_dir, basename(in_list$prot[i]))
        writeXStringSet(new_prot, new_prot_fn)
        in_list$prot[i] <- new_prot_fn
    }
    in_list$update <- NULL
    return(in_list)
}

#' @importFrom Biostrings readDNAStringSet
#' @importFrom GenomicFeatures cdsBy extractTranscriptSeqs
#' @importFrom txdbmaker makeTxDbFromGFF
#' @import BSgenome
#'
.makeCDS <- function(gff, genome){
    # Create a TxDb object from the GFF file
    txdb <- suppressMessages({makeTxDbFromGFF(file = gff)})
    
    # Read the genome file as a DNAStringSet object
    genome <- readDNAStringSet(filepath = genome)
    
    # Extract CDS sequences from the TxDb object
    cds_db <- suppressMessages({cdsBy(x = txdb, by = "tx", use.names = TRUE)})
    
    cds <- suppressMessages({extractTranscriptSeqs(x = genome, transcripts = cds_db)})
    
    # Order CDS sequences by their names
    cds <- cds[order(names(cds))]
    return(cds)
}
