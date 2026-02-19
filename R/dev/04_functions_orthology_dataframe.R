#' Process RBH files using R/05_functions_orthology.R algorithm
#'
#' This function uses the same input as process_rbh_direct() but applies
#' the algorithm from R/05_functions_orthology.R (linkGene2Genome() onwards)
#' for speed comparison.
#'
#' @param working_dir Working directory. GFF files are searched in `working_dir/input/<index>_<genome_name>/gff_df.rds`.
#'   Output TSV files are written to `working_dir/orthopair/<genome_pair_id>.tsv`.
#' @param rbh_dir Directory containing split RBH files (<genome_pair_id>.rbh).
#'   If NULL, defaults to `working_dir/rbh`.
#' @param gff_df_paths Character vector of paths to gff_df.rds files (one per genome).
#'   If NULL, automatically searches in `working_dir/input/` folders.
#' @param genome_width Number of leading characters encoding genome ID (default 4)
#' @param pident_threshold_quantile Quantile threshold for filtering anchors (default 0.25)
#' @param mutual_ci_min Minimum mutual_ci for RBH edges (default NULL)
#' @param n_threads Number of cores for parallel processing (default 1L)
#' @param load_all_gff If TRUE, load all GFF data upfront (faster but uses more memory).
#'   If FALSE, load GFF data per pair (slower but memory-efficient for 300+ genomes).
#' @param verbose Print progress messages
#' @return List of data.frames, one per genome pair, with ortholog pairs.
#'   Also writes TSV files to `working_dir/orthopair/`.
#'
#' @importFrom data.table fread fwrite setDT setkey setorder setnames
#' @importFrom parallel mclapply
#' @export
process_rbh_dataframe <- function(working_dir,
                                  rbh_dir = NULL,
                                  gff_df_paths = NULL,
                                  genome_width = 4L,
                                  pident_threshold_quantile = 0.25,
                                  mutual_ci_min = NULL,
                                  n_threads = 1L,
                                  load_all_gff = FALSE,
                                  verbose = TRUE) {
    
    # Set default rbh_dir if not provided
    if (is.null(rbh_dir)) {
        rbh_dir <- file.path(working_dir, "rbh")
    }
    
    # Find all RBH files (exclude same-genome comparisons, e.g. 0001_0001.rbh)
    rbh_files <- list.files(rbh_dir, pattern = "\\.rbh$", full.names = TRUE)
    if (length(rbh_files) > 0L) {
        pair_ids <- sub("\\.rbh$", "", basename(rbh_files))
        parts <- strsplit(pair_ids, "_", fixed = TRUE)
        same_genome <- lengths(parts) == 2L & vapply(parts, function(p) p[1L] == p[2L], logical(1L))
        rbh_files <- rbh_files[!same_genome]
    }
    if (length(rbh_files) == 0L) {
        warning("No RBH files found in: ", rbh_dir, " (or only same-genome comparisons)")
        return(list())
    }
    
    if (verbose) message("[ortho] Found ", length(rbh_files), " RBH pair files")
    
    # Get GFF data paths - search in working_dir/input/<index>_<genome_name>/ folders
    if (is.null(gff_df_paths)) {
        input_dir <- file.path(working_dir, "input")
        if (!dir.exists(input_dir)) {
            stop("Input directory not found: ", input_dir)
        }
        
        # Find all folders matching pattern <index>_<genome_name>
        input_folders <- list.dirs(input_dir, full.names = TRUE, recursive = FALSE)
        # Filter folders that match pattern (have underscore and contain gff_df.rds)
        gff_df_paths <- character(0)
        for (folder in input_folders) {
            gff_fn <- file.path(folder, "gff_df.rds")
            if (file.exists(gff_fn)) {
                gff_df_paths <- c(gff_df_paths, gff_fn)
            }
        }
        
        if (length(gff_df_paths) == 0L) {
            stop("No gff_df.rds files found in ", input_dir, 
                 " subdirectories (expected pattern: <index>_<genome_name>/gff_df.rds)")
        }
        
        if (verbose) message("[ortho] Found ", length(gff_df_paths), " GFF files")
    }
    
    # Create output directory for ortholog pairs
    orthopair_dir <- file.path(working_dir, "orthopair")
    dir.create(orthopair_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Create genome ID to GFF path mapping for efficient lookup
    gff_lookup <- .create_gff_lookup_dataframe(gff_df_paths, genome_width, verbose)
    if (length(gff_lookup) == 0L) {
        stop("Failed to create GFF lookup. Check that GFF files contain valid tx_index values.")
    }
    
    # Load GFF data once if requested (faster but uses more memory)
    gff_data_all <- NULL
    if (load_all_gff) {
        if (verbose) message("[ortho] Loading all GFF data (memory-intensive mode)")
        gff_data_all <- .load_gff_data_dataframe(gff_df_paths, genome_width)
        if (is.null(gff_data_all) || length(gff_data_all) == 0L) {
            stop("Failed to load GFF data. Check GFF file format.")
        }
        if (verbose) message("[ortho] Loaded GFF data from ", length(gff_df_paths), " genomes")
    }
    
    # Process each RBH file in parallel
    if (verbose) message("[ortho] Processing ", length(rbh_files), " genome pairs")
    
    # Worker function
    process_one_rbh_file <- function(rbh_fn) {
        pair_id <- sub("\\.rbh$", "", basename(rbh_fn))
        
        tryCatch({
            result <- .process_one_pair_dataframe(
                rbh_fn = rbh_fn,
                pair_id = pair_id,
                gff_data_all = gff_data_all,
                gff_lookup = gff_lookup,
                genome_width = genome_width,
                pident_threshold_quantile = pident_threshold_quantile,
                mutual_ci_min = mutual_ci_min,
                load_all_gff = load_all_gff,
                verbose = FALSE
            )
            
            # Write result to TSV file
            out_tsv <- file.path(orthopair_dir, paste0(pair_id, ".tsv"))
            if (!is.null(result) && nrow(result) > 0L) {
                fwrite(result, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
            } else {
                # Write empty file with header
                empty_df <- if (!is.null(result)) result[0L, , drop = FALSE] else 
                    data.frame(query_gene = character(0), subject_gene = character(0))
                write.table(empty_df, file = out_tsv, sep = "\t", quote = FALSE, 
                            row.names = FALSE, col.names = TRUE)
            }
            
            return(list(pair_id = pair_id, result = result, error = NULL))
        }, error = function(e) {
            out_tsv <- file.path(orthopair_dir, paste0(pair_id, ".tsv"))
            writeLines("query_gene\tsubject_gene", con = out_tsv)
            return(list(pair_id = pair_id, result = NULL, error = e$message))
        })
    }
    
    # Process in parallel or sequentially
    if (n_threads > 1L && .Platform$OS.type != "windows") {
        if (verbose) message("[ortho] Using parallel processing with ", n_threads, " cores")
        results_list <- parallel::mclapply(rbh_files, process_one_rbh_file, mc.cores = n_threads)
    } else {
        if (verbose) message("[ortho] Using sequential processing")
        results_list <- lapply(rbh_files, process_one_rbh_file)
    }
    
    # Extract results and errors
    results <- lapply(results_list, function(x) x$result)
    errors <- lapply(results_list, function(x) x$error)
    pair_ids <- sapply(results_list, function(x) x$pair_id)
    
    # Report errors
    error_idx <- !sapply(errors, is.null)
    if (any(error_idx)) {
        error_pairs <- pair_ids[error_idx]
        error_msgs <- errors[error_idx]
        warning("[ortho] Errors occurred in ", sum(error_idx), " pairs:\n",
                paste(sprintf("  %s: %s", error_pairs, error_msgs), collapse = "\n"))
    }
    
    # Report NULL results
    null_idx <- sapply(results, is.null)
    if (any(null_idx) && verbose) {
        null_pairs <- pair_ids[null_idx]
        message("[ortho] ", sum(null_idx), " pairs returned NULL results: ", 
                paste(head(null_pairs, 10), collapse = ", "), 
                if (length(null_pairs) > 10) " ..." else "")
    }
    
    # Remove NULL results
    results <- results[!sapply(results, is.null)]
    
    # Summary statistics
    n_total <- length(rbh_files)
    n_with_results <- length(results)
    n_null <- sum(null_idx)
    n_error <- sum(error_idx)
    
    if (verbose) {
        message("[ortho] Processing summary:")
        message("[ortho]   Total pairs: ", n_total)
        message("[ortho]   Pairs with results: ", n_with_results)
        message("[ortho]   Pairs with NULL results: ", n_null)
        message("[ortho]   Pairs with errors: ", n_error)
        message("[ortho] Output TSV files written to: ", orthopair_dir)
    }
    
    return(results)
}

#' Create GFF lookup: genome ID -> GFF file path
.create_gff_lookup_dataframe <- function(gff_df_paths, genome_width, verbose) {
    lookup <- list()
    
    for (path in gff_df_paths) {
        tryCatch({
            gff <- readRDS(path)
            if (is.data.frame(gff) && nrow(gff) > 0L) {
                tx <- gff[gff$type %in% c("mRNA", "transcript"), ]
                if (nrow(tx) > 0L && "tx_index" %in% names(tx)) {
                    tx_idx <- as.character(tx$tx_index[1L])
                    if (nchar(tx_idx) >= genome_width) {
                        genome_id <- as.integer(substr(tx_idx, 1L, genome_width))
                        lookup[[as.character(genome_id)]] <- path
                    }
                }
            }
        }, error = function(e) NULL)
    }
    
    if (verbose && length(lookup) > 0L) {
        message("[ortho] Created GFF lookup for ", length(lookup), " genomes")
    }
    
    return(lookup)
}

#' Load GFF data for all genomes (for linkGene2Genome compatibility)
.load_gff_data_dataframe <- function(gff_df_paths, genome_width) {
    gff_list <- lapply(gff_df_paths, function(path) {
        gff <- readRDS(path)
        if (!is.data.frame(gff)) return(NULL)
        return(gff)
    })
    
    # Filter NULLs
    gff_list <- gff_list[!sapply(gff_list, is.null)]
    if (length(gff_list) == 0L) {
        stop("No valid GFF data loaded")
    }
    
    return(gff_list)
}

#' Process one genome pair RBH file using R/05_functions_orthology.R algorithm
.process_one_pair_dataframe <- function(rbh_fn,
                                        pair_id,
                                        gff_data_all,
                                        gff_lookup,
                                        genome_width,
                                        pident_threshold_quantile,
                                        mutual_ci_min,
                                        load_all_gff,
                                        verbose) {
    
    # Read RBH file
    if (verbose) message("[ortho] Reading RBH file: ", basename(rbh_fn))
    rbh <- fread(
        file = rbh_fn,
        header = FALSE,
        sep = "\t",
        select = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L),
        col.names = c("query_tx", "subject_tx", "pident", "q2s_qcovs", 
                      "s2q_qcovs", "q2s_ci", "s2q_ci", "mutual_ci"),
        na.strings = "",
        colClasses = c("integer", "integer", "numeric", "integer", 
                       "integer", "numeric", "numeric", "numeric"),
        verbose = FALSE
    )
    setDT(rbh)
    
    if (verbose) message("[ortho] Pair ", pair_id, ": Read ", nrow(rbh), " RBH records")
    
    # Filter valid pairs
    rbh <- rbh[!is.na(query_tx) & !is.na(subject_tx) & query_tx != subject_tx]
    if (nrow(rbh) == 0L) {
        if (verbose) message("[ortho] Pair ", pair_id, ": No valid RBH pairs")
        return(NULL)
    }
    
    # Extract genome IDs
    query_genome_pre <- as.integer(substr(as.character(rbh$query_tx), 1L, genome_width))
    subject_genome_pre <- as.integer(substr(as.character(rbh$subject_tx), 1L, genome_width))
    
    # Filter by mutual_ci if specified
    if (!is.null(mutual_ci_min)) {
        rbh <- rbh[mutual_ci >= mutual_ci_min]
        if (nrow(rbh) == 0L) return(NULL)
    }
    
    # Load GFF data for this pair
    if (!load_all_gff || is.null(gff_data_all)) {
        pair_genomes <- unique(c(query_genome_pre, subject_genome_pre))
        gff_data_pair <- lapply(as.character(pair_genomes), function(gid) {
            path <- gff_lookup[[gid]]
            if (is.null(path) || !file.exists(path)) return(NULL)
            return(readRDS(path))
        })
        gff_data_pair <- gff_data_pair[!sapply(gff_data_pair, is.null)]
    } else {
        # Extract GFF data for this pair from pre-loaded data
        pair_genomes <- unique(c(query_genome_pre, subject_genome_pre))
        gff_data_pair <- lapply(as.character(pair_genomes), function(gid) {
            path <- gff_lookup[[gid]]
            if (is.null(path)) return(NULL)
            # Find matching GFF in pre-loaded data
            idx <- which(sapply(gff_data_all, function(g) {
                if (is.null(g) || nrow(g) == 0L) return(FALSE)
                tx <- g[g$type %in% c("mRNA", "transcript"), ]
                if (nrow(tx) == 0L) return(FALSE)
                tx_idx <- as.character(tx$tx_index[1L])
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
    
    if (length(gff_data_pair) == 0L) {
        warning("Could not load GFF data for pair: ", pair_id)
        return(NULL)
    }
    
    # Split GFF data into query and subject based on genome IDs
    query_genome_id <- as.character(min(pair_genomes))
    subject_genome_id <- as.character(max(pair_genomes))
    
    query_gff_full <- NULL
    subject_gff_full <- NULL
    
    for (gff_df in gff_data_pair) {
        if (is.null(gff_df) || nrow(gff_df) == 0L) next
        tx <- gff_df[gff_df$type %in% c("mRNA", "transcript"), ]
        if (nrow(tx) == 0L) next
        tx_idx <- as.character(tx$tx_index[1L])
        if (nchar(tx_idx) >= genome_width) {
            gid <- as.character(as.integer(substr(tx_idx, 1L, genome_width)))
            if (gid == query_genome_id) {
                query_gff_full <- gff_df
            } else if (gid == subject_genome_id) {
                subject_gff_full <- gff_df
            }
        }
    }
    
    if (is.null(query_gff_full) || is.null(subject_gff_full)) {
        warning("Could not identify query/subject GFF for pair: ", pair_id)
        return(NULL)
    }
    
    # Order GFF data (as expected by linkGene2Genome)
    query_gff_full <- query_gff_full[order(query_gff_full$seqnames, query_gff_full$start), ]
    subject_gff_full <- subject_gff_full[order(subject_gff_full$seqnames, subject_gff_full$start), ]
    query_gff_full$gene_index <- as.integer(query_gff_full$gene_index)
    query_gff_full$tx_index <- as.integer(query_gff_full$tx_index)
    subject_gff_full$gene_index <- as.integer(subject_gff_full$gene_index)
    subject_gff_full$tx_index <- as.integer(subject_gff_full$tx_index)
    
    # Create g2g_graph structure (equivalent to linkGene2Genome output)
    g2g_graph <- .linkGene2Genome_dataframe(query_gff = query_gff_full,
                                            subject_gff = subject_gff_full)
    
    # Convert RBH to format expected by R/05_functions_orthology.R
    # Merge with GFF to get gene_index
    query_tx_match <- match(rbh$query_tx, query_gff_full$tx_index)
    subject_tx_match <- match(rbh$subject_tx, subject_gff_full$tx_index)
    
    rbh_df <- data.frame(
        qseqid = as.integer(rbh$query_tx),
        sseqid = as.integer(rbh$subject_tx),
        qgeneid = as.integer(query_gff_full$gene_index[query_tx_match]),
        sgeneid = as.integer(subject_gff_full$gene_index[subject_tx_match]),
        pid = paste(query_gff_full$gene_id[query_tx_match], subject_gff_full$gene_id[subject_tx_match], sep = "_"),
        ci_q2s = rbh$q2s_ci,
        ci_s2q = rbh$s2q_ci,
        pident = rbh$pident,
        mutual_ci = rbh$mutual_ci,
        stringsAsFactors = FALSE
    )
    
    rbh_df <- rbh_df[!is.na(rbh_df$qgeneid) & !is.na(rbh_df$sgeneid), ]
    rbh_df$pair_id <- paste(rbh_df$qgeneid, rbh_df$sgeneid, sep = "_")
    
    if (nrow(rbh_df) == 0L) {
        if (verbose) message("[ortho] Pair ", pair_id, ": No RBH after GFF merge")
        return(NULL)
    }
    
    # Apply R/05_functions_orthology.R algorithm from linkGene2Genome onwards
    # Note: .findAnchors checks length(rbh) == 1, so ensure rbh_df is a data.frame
    anchor <- .findAnchors(rbh = as.data.frame(rbh_df), g2g_graph = g2g_graph)
    
    if (all(is.na(anchor)) || (is.data.frame(anchor) && nrow(anchor) == 0L)) {
        if (verbose) message("[ortho] Pair ", pair_id, ": No anchors found")
        return(NULL)
    }
    
    t2a_graph <- .link2Anchor(g2g_graph = g2g_graph, anchor = anchor)
    # .findSyntenicOrtho expects rbh as data.frame and will rename first 2 columns
    orthopair <- .findSyntenicOrtho(rbh = as.data.frame(rbh_df),
                                    anchor = anchor,
                                    g2g_graph = g2g_graph,
                                    t2a_graph = t2a_graph,
                                    rbh_threshold = NULL)
    
    if (is.null(orthopair) || nrow(orthopair) == 0L) {
        if (verbose) message("[ortho] Pair ", pair_id, ": No syntenic orthologs")
        return(NULL)
    }
    
    orthopair <- .findSyntenyBlocks(orthopair = orthopair)
    orthopair <- .pickBestPair(orthopair = orthopair)
    orthopair <- .classifyOrthoPair(orthopair = orthopair)
    orthopair <- .filterOrthopair(orthopair = orthopair, g2g_graph = g2g_graph)
    orthopair <- .classifyOrthoPair(orthopair = orthopair)
    orthopair <- .reformatOrthoPair(orthopair = orthopair, g2g_graph = g2g_graph)
    
    # Convert to data.table for consistency
    setDT(orthopair)
    
    if (verbose) {
        message("[ortho] Pair ", pair_id, ": ", nrow(orthopair), " ortholog pairs")
    }
    
    return(orthopair)
}

#' Create g2g_graph structure from GFF data (equivalent to linkGene2Genome)
.linkGene2Genome_dataframe <- function(query_gff, subject_gff) {
    # Extract transcript-level data
    q_tx_i <- query_gff$type %in% c("transcript", "mRNA")
    query_df <- query_gff[q_tx_i, ]
    
    # Extract CDS-level data for query
    q_cds_i <- query_gff$type %in% "CDS"
    query_gff_cds <- query_gff[q_cds_i, ]
    if (nrow(query_gff_cds) > 0L && "Parent" %in% names(query_gff_cds)) {
        hit <- match(unlist(query_gff_cds$Parent), query_gff$ID)
        query_gff_cds$tx_index <- query_gff$tx_index[hit]
    }
    
    # Extract transcript-level data for subject
    s_tx_i <- subject_gff$type %in% c("transcript", "mRNA")
    subject_df <- subject_gff[s_tx_i, ]
    
    # Extract CDS-level data for subject
    s_cds_i <- subject_gff$type %in% "CDS"
    subject_gff_cds <- subject_gff[s_cds_i, ]
    if (nrow(subject_gff_cds) > 0L && "Parent" %in% names(subject_gff_cds)) {
        hit <- match(unlist(subject_gff_cds$Parent), subject_gff$ID)
        subject_gff_cds$tx_index <- subject_gff$tx_index[hit]
    }
    
    out <- list(query_df = query_df,
                subject_df = subject_df,
                query_gff = query_gff_cds,
                subject_gff = subject_gff_cds)
    return(out)
}
