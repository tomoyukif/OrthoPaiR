#' Process RBH files directly without graph construction
#'
#' Optimized version that processes RBH files split by genome pairs directly,
#' avoiding graph construction overhead. Uses data.table for speed.
#' Memory-efficient: loads GFF data per genome pair to handle 300+ genomes on 64GB RAM.
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
#'
#' @examples
#' # Test a single pair for debugging:
#' # rbh_files <- list.files(file.path(working_dir, "rbh"), pattern = "\\.rbh$", full.names = TRUE)
#' # test_result <- .process_one_pair_fast(rbh_files[1], ...)
process_rbh_direct <- function(working_dir,
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
  gff_lookup <- .create_gff_lookup(gff_df_paths, genome_width, verbose)
  if (length(gff_lookup) == 0L) {
    stop("Failed to create GFF lookup. Check that GFF files contain valid tx_index values.")
  }
  
  # Load GFF data once if requested (faster but uses more memory)
  gff_data <- NULL
  if (load_all_gff) {
    if (verbose) message("[ortho] Loading all GFF data (memory-intensive mode)")
    gff_data <- .load_gff_data_fast(gff_df_paths, genome_width)
    if (is.null(gff_data) || nrow(gff_data) == 0L) {
      stop("Failed to load GFF data. Check GFF file format.")
    }
    if (verbose) message("[ortho] Loaded ", nrow(gff_data), " transcript records from ", length(gff_df_paths), " genomes")
  }
  
  # Process each RBH file in parallel
  if (verbose) message("[ortho] Processing ", length(rbh_files), " genome pairs")
  
  # Worker function
  process_one_rbh_file <- function(rbh_fn) {
    pair_id <- sub("\\.rbh$", "", basename(rbh_fn))
    
    tryCatch({
      # Process with verbose disabled in parallel workers to avoid message clutter
      # Set verbose=TRUE in the main function call to see detailed progress
      result <- .process_one_pair_fast(
        rbh_fn = rbh_fn,
        pair_id = pair_id,
        gff_data = gff_data,
        gff_lookup = gff_lookup,
        genome_width = genome_width,
        pident_threshold_quantile = pident_threshold_quantile,
        mutual_ci_min = mutual_ci_min,
        load_all_gff = load_all_gff,
        verbose = FALSE  # Suppress verbose in parallel workers (messages don't show in mclapply)
      )
      
      # Write result to TSV file (even if empty, for debugging)
      out_tsv <- file.path(orthopair_dir, paste0(pair_id, ".tsv"))
      if (!is.null(result)) {
        if (nrow(result) > 0L) {
          fwrite(result, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
        } else {
          # Write empty file with header for debugging
          empty_df <- result[0L, , drop = FALSE]
          write.table(empty_df, file = out_tsv, sep = "\t", quote = FALSE, 
                     row.names = FALSE, col.names = TRUE)
        }
      } else {
        # Write empty file even for NULL results
        writeLines("query_gene\tsubject_gene", con = out_tsv)
      }
      
      return(list(pair_id = pair_id, result = result, error = NULL))
    }, error = function(e) {
      # Write empty file on error
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
  n_empty <- n_total - n_with_results - n_null - n_error
  
  if (verbose) {
    message("[ortho] Processing summary:")
    message("[ortho]   Total pairs: ", n_total)
    message("[ortho]   Pairs with results: ", n_with_results)
    message("[ortho]   Pairs with NULL results: ", n_null)
    message("[ortho]   Pairs with errors: ", n_error)
    message("[ortho] Output TSV files written to: ", orthopair_dir)
    
    # If all results are NULL, suggest debugging
    if (n_with_results == 0L && n_total > 0L) {
      message("[ortho] WARNING: All pairs returned NULL results.")
      message("[ortho]   This may indicate:")
      message("[ortho]   1. No anchors found (check pident_threshold_quantile)")
      message("[ortho]   2. GFF data mismatch (check tx_index format)")
      message("[ortho]   3. RBH data format issues")
      message("[ortho]   Try processing a single pair with verbose=TRUE for details")
    }
  }
  
  return(results)
}

#' Create GFF lookup: genome ID -> GFF file path
.create_gff_lookup <- function(gff_df_paths, genome_width, verbose) {
  # Extract genome IDs from file paths or file contents
  lookup <- list()
  
  for (path in gff_df_paths) {
    # Try to extract genome ID from path first (faster)
    path_parts <- strsplit(basename(dirname(path)), "/")[[1]]
    genome_id <- NULL
    
    # Try reading first few lines to get genome ID from tx_index
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

#' Load GFF data efficiently (all genomes)
.load_gff_data_fast <- function(gff_df_paths, genome_width) {
  gff_list <- lapply(gff_df_paths, function(path) {
    gff <- readRDS(path)
    if (!is.data.frame(gff)) return(NULL)
    
    # Extract transcript-level data
    tx <- gff[gff$type %in% c("mRNA", "transcript"), ]
    if (nrow(tx) == 0L) return(NULL)
    
    # Ensure tx_index is integer
    if (is.character(tx$tx_index) || is.factor(tx$tx_index)) {
      tx$tx_index <- suppressWarnings(as.integer(as.character(tx$tx_index)))
    }
    if (is.character(tx$gene_index) || is.factor(tx$gene_index)) {
      tx$gene_index <- suppressWarnings(as.integer(as.character(tx$gene_index)))
    }
    
    # Extract genome ID from tx_index
    tx_char <- as.character(tx$tx_index)
    tx$genome <- as.integer(substr(tx_char, 1L, genome_width))
    
    # Keep only needed columns (include gene_id and transcript_id for orthopair output)
    keep_cols <- c("tx_index", "gene_index", "seqnames", "start", "end", "genome")
    if ("gene_id" %in% names(tx)) keep_cols <- c(keep_cols, "gene_id")
    if ("transcript_id" %in% names(tx)) {
      keep_cols <- c(keep_cols, "transcript_id")
    } else if ("ID" %in% names(tx)) {
      tx$transcript_id <- tx$ID
      keep_cols <- c(keep_cols, "transcript_id")
    }
    tx <- tx[, keep_cols[keep_cols %in% names(tx)], drop = FALSE]
    tx <- tx[!is.na(tx$tx_index) & !is.na(tx$gene_index), ]
    
    setDT(tx)
    setkey(tx, tx_index)
    return(tx)
  })
  
  # Combine all genomes
  gff_combined <- do.call(rbind, gff_list)
  if (is.null(gff_combined) || nrow(gff_combined) == 0L) {
    stop("No valid GFF data loaded")
  }
  setDT(gff_combined)
  setkey(gff_combined, tx_index)
  return(gff_combined)
}

#' Load GFF data for specific genomes (memory-efficient)
.load_gff_for_pair <- function(genome_ids, gff_lookup, genome_width) {
  gff_list <- lapply(as.character(genome_ids), function(gid) {
    path <- gff_lookup[[gid]]
    if (is.null(path) || !file.exists(path)) return(NULL)
    
    gff <- readRDS(path)
    if (!is.data.frame(gff)) return(NULL)
    
    # Extract transcript-level data
    tx <- gff[gff$type %in% c("mRNA", "transcript"), ]
    if (nrow(tx) == 0L) return(NULL)
    
    # Ensure tx_index is integer
    if (is.character(tx$tx_index) || is.factor(tx$tx_index)) {
      tx$tx_index <- suppressWarnings(as.integer(as.character(tx$tx_index)))
    }
    if (is.character(tx$gene_index) || is.factor(tx$gene_index)) {
      tx$gene_index <- suppressWarnings(as.integer(as.character(tx$gene_index)))
    }
    
    # Extract genome ID from tx_index
    tx_char <- as.character(tx$tx_index)
    tx$genome <- as.integer(substr(tx_char, 1L, genome_width))
    
    # Keep only needed columns (include gene_id and transcript_id for orthopair output)
    keep_cols <- c("tx_index", "gene_index", "seqnames", "start", "end", "genome")
    if ("gene_id" %in% names(tx)) keep_cols <- c(keep_cols, "gene_id")
    if ("transcript_id" %in% names(tx)) {
      keep_cols <- c(keep_cols, "transcript_id")
    } else if ("ID" %in% names(tx)) {
      tx$transcript_id <- tx$ID
      keep_cols <- c(keep_cols, "transcript_id")
    }
    tx <- tx[, keep_cols[keep_cols %in% names(tx)], drop = FALSE]
    tx <- tx[!is.na(tx$tx_index) & !is.na(tx$gene_index), ]
    
    setDT(tx)
    setkey(tx, tx_index)
    return(tx)
  })
  
  # Combine genomes for this pair
  gff_combined <- do.call(rbind, gff_list)
  if (is.null(gff_combined) || nrow(gff_combined) == 0L) {
    return(NULL)
  }
  setDT(gff_combined)
  setkey(gff_combined, tx_index)
  return(gff_combined)
}

#' Process one genome pair RBH file
.process_one_pair_fast <- function(rbh_fn,
                                   pair_id,
                                   gff_data,
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
  
  if (verbose) message("[ortho] Pair ", pair_id, ": Read ", nrow(rbh), " RBH records from file")
  
  # Filter valid pairs
  rbh <- rbh[!is.na(query_tx) & !is.na(subject_tx) & query_tx != subject_tx]
  if (nrow(rbh) == 0L) {
    if (verbose) message("[ortho] Pair ", pair_id, ": No valid RBH pairs after filtering self-hits")
    return(NULL)
  }
  
  if (verbose) message("[ortho] Pair ", pair_id, ": ", nrow(rbh), " RBH records after filtering self-hits")
  
  # Extract genome IDs (before filtering, for GFF lookup)
  query_genome_pre <- as.integer(substr(as.character(rbh$query_tx), 1L, genome_width))
  subject_genome_pre <- as.integer(substr(as.character(rbh$subject_tx), 1L, genome_width))
  
  # Filter by mutual_ci if specified
  if (!is.null(mutual_ci_min)) {
    rbh <- rbh[mutual_ci >= mutual_ci_min]
    if (nrow(rbh) == 0L) return(NULL)
  }
  
  # Load GFF data for this pair if not loaded globally
  if (!load_all_gff || is.null(gff_data)) {
    pair_genomes <- unique(c(query_genome_pre, subject_genome_pre))

    if (verbose) message("[ortho] Loading GFF for genomes: ", paste(pair_genomes, collapse = ", "))
    pair_gff_data <- .load_gff_for_pair(pair_genomes, gff_lookup, genome_width)

    if (is.null(pair_gff_data) || nrow(pair_gff_data) == 0L) {
      warning("Could not load GFF data for pair: ", pair_id, 
              " (genomes: ", paste(pair_genomes, collapse = ", "), ")")
      return(NULL)
    }

    if (verbose) message("[ortho] Loaded ", nrow(pair_gff_data), " transcript records for pair: ", pair_id)
    gff_data <- pair_gff_data
  }
  
  # Merge with GFF data to get gene information
  query_gff <- gff_data[gff_data$tx_index %in% rbh$query_tx]
  query_gff <- query_gff[, c("tx_index", "gene_index", "seqnames", "start", "end"), with = FALSE]
  setnames(query_gff, c("gene_index", "seqnames", "start", "end"), 
          c("query_gene", "query_chr", "query_start", "query_end"))
  rbh <- merge(rbh, query_gff, by.x = "query_tx", by.y = "tx_index", all.x = FALSE)
  
  subject_gff <- gff_data[gff_data$tx_index %in% rbh$subject_tx]
  subject_gff <- subject_gff[, c("tx_index", "gene_index", "seqnames", "start", "end"), with = FALSE]
  setnames(subject_gff, c("gene_index", "seqnames", "start", "end"),
          c("subject_gene", "subject_chr", "subject_start", "subject_end"))
  rbh <- merge(rbh, subject_gff, by.x = "subject_tx", by.y = "tx_index", all.x = FALSE)
  
  if (nrow(rbh) == 0L) {
    if (verbose) message("[ortho] No RBH records remaining after GFF merge for pair: ", pair_id)
    return(NULL)
  }
  
  if (verbose) message("[ortho] Pair ", pair_id, ": ", nrow(rbh), " RBH records after GFF merge")
  
  # Add genome IDs (re-extract after filtering)
  rbh$query_genome <- as.integer(substr(as.character(rbh$query_tx), 1L, genome_width))
  rbh$subject_genome <- as.integer(substr(as.character(rbh$subject_tx), 1L, genome_width))
  
  # Find anchors (RBBH) - optimized version
  anchors <- .find_anchors_fast(rbh, pident_threshold_quantile, verbose)
  
  if (is.null(anchors) || nrow(anchors) == 0L) {
    if (verbose) message("[ortho] No anchors found for pair: ", pair_id, " (check pident_threshold_quantile)")
    return(NULL)
  }
  
  if (verbose) message("[ortho] Pair ", pair_id, ": Found ", nrow(anchors), " anchors")
  
  # Link genes to anchors - optimized
  gene_to_anchor <- .link_to_anchors_fast(rbh, anchors, verbose)
  
  # Find syntenic orthologs - optimized
  orthopair <- .find_syntenic_ortho_fast(rbh, anchors, gene_to_anchor, verbose)
  
  if (is.null(orthopair) || nrow(orthopair) == 0L) {
    if (verbose) message("[ortho] No syntenic orthologs for pair: ", pair_id)
    return(NULL)
  }
  
  # Find synteny blocks
  orthopair <- .find_synteny_blocks_fast(orthopair, verbose)
  setorder(orthopair, -mutual_ci)
  orthopair <- orthopair[!duplicated(paste(query_gene, subject_gene, sep = "_"))]
  
  # Classify ortholog pairs
  orthopair <- .classify_orthopair_fast(orthopair, verbose)
  
  # Filter ortholog pairs
  orthopair <- .filter_orthopair_fast(orthopair, verbose)
  
  # Re-classify after filtering
  orthopair <- .classify_orthopair_fast(orthopair, verbose = FALSE)
  
  # Replace query_gene, query_tx, subject_gene, subject_tx with gene_id and transcript_id from GFF
  if ("gene_id" %in% names(gff_data) && "transcript_id" %in% names(gff_data)) {
    gff_id <- gff_data[, .(tx_index, gene_index, gene_id, transcript_id)]
    # Query side: tx_index -> gene_id, transcript_id
    q_id <- gff_id[, .(query_tx = tx_index, query_gene_id = gene_id, query_transcript_id = transcript_id)]
    orthopair <- merge(orthopair, unique(q_id), by = "query_tx", all.x = TRUE)
    # Subject side
    s_id <- gff_id[, .(subject_tx = tx_index, subject_gene_id = gene_id, subject_transcript_id = transcript_id)]
    orthopair <- merge(orthopair, unique(s_id), by = "subject_tx", all.x = TRUE)
    # Replace the four columns
    orthopair[, `:=`(query_gene = query_gene_id, query_tx = query_transcript_id,
                     subject_gene = subject_gene_id, subject_tx = subject_transcript_id)]
    orthopair[, c("query_gene_id", "query_transcript_id", "subject_gene_id", "subject_transcript_id") := NULL]
  }
  
  if (verbose) {
    message("[ortho] Pair ", pair_id, ": ", nrow(orthopair), " ortholog pairs")
  }
  
  return(orthopair)
}

#' Find anchors (RBBH) - optimized version
.find_anchors_fast <- function(rbh, pident_threshold_quantile, verbose) {
  # Filter: only different genomes
  rbh <- rbh[query_genome != subject_genome]
  if (nrow(rbh) == 0L) return(NULL)
  
  # Find RBBH: best hit in each direction per gene
  setorder(rbh, -q2s_ci)
  best_from <- rbh[!duplicated(rbh$query_gene), ]
  
  setorder(rbh, -s2q_ci)
  best_to <- rbh[!duplicated(rbh$subject_gene), ]
  
  # RBBH = intersection
  rbbh <- merge(best_from[, c("query_gene", "subject_gene", "query_genome", "subject_genome"), with = FALSE],
               best_to[, c("query_gene", "subject_gene", "query_genome", "subject_genome"), with = FALSE],
               by = c("query_gene", "subject_gene", "query_genome", "subject_genome"))
  
  # Merge back attributes
  rbbh <- merge(rbbh, rbh, by = c("query_gene", "subject_gene", "query_genome", "subject_genome"))
  
  # If multiple transcript pairs per gene pair, take best by mutual_ci
  setorder(rbbh, -mutual_ci)
  rbbh <- rbbh[!duplicated(rbbh[, c("query_gene", "subject_gene"), with = FALSE]), ]
  anchors <- rbbh
  
  if (nrow(anchors) == 0L) return(NULL)
  
  # Filter by pident threshold
  if (!is.null(pident_threshold_quantile) && nrow(anchors) > 0L) {
    q <- quantile(anchors$pident, c(pident_threshold_quantile, 1 - pident_threshold_quantile), na.rm = TRUE)
    threshold <- q[1] - 1.5 * diff(q)
    anchors <- anchors[pident >= threshold]
  }
  
  # Check for high-copy genes
  anchors <- .check_high_copy_genes_fast(anchors, rbh)
  
  if (nrow(anchors) == 0L) return(NULL)
  
  anchors$anchor_id <- seq_len(nrow(anchors))
  
  return(anchors)
}

#' Check for high-copy genes - optimized
.check_high_copy_genes_fast <- function(anchors, rbh) {
  if (nrow(anchors) == 0L) return(anchors)
  
  # Count RBH per gene
  q_count_dt <- rbh[, .N, by = query_gene]
  q_rbh_count <- q_count_dt$N
  names(q_rbh_count) <- as.character(q_count_dt$query_gene)
  s_count_dt <- rbh[, .N, by = subject_gene]
  s_rbh_count <- s_count_dt$N
  names(s_rbh_count) <- as.character(s_count_dt$subject_gene)
  
  # Find outliers
  q_quant <- quantile(q_rbh_count, c(0.25, 0.75), na.rm = TRUE)
  q_threshold <- q_quant[2] + 1.5 * diff(q_quant)
  s_quant <- quantile(s_rbh_count, c(0.25, 0.75), na.rm = TRUE)
  s_threshold <- s_quant[2] + 1.5 * diff(s_quant)
  
  q_outliers <- names(q_rbh_count)[q_rbh_count > q_threshold]
  s_outliers <- names(s_rbh_count)[s_rbh_count > s_threshold]
  
  # Remove anchors involving outlier genes (unless both are outliers)
  keep <- !(anchors$query_gene %in% q_outliers & !anchors$subject_gene %in% s_outliers) &
          !(anchors$subject_gene %in% s_outliers & !anchors$query_gene %in% q_outliers)
  return(anchors[keep])
}

#' Link genes to anchors - optimized
.link_to_anchors_fast <- function(rbh, anchors, verbose) {
  # Extract unique genes with positions
  query_genes <- unique(rbh[, .(query_gene, query_chr, query_start, query_end)])
  subject_genes <- unique(rbh[, .(subject_gene, subject_chr, subject_start, subject_end)])
  query_anchors <- anchors[, .(query_gene, query_chr)]
  subject_anchors <- anchors[, .(subject_gene, subject_chr)]
  
  # Link query genes to anchors
  query_links <- .link_genes_to_anchors_core(query_genes, query_anchors, "query")
  subject_links <- .link_genes_to_anchors_core(subject_genes, subject_anchors, "subject")
  
  return(list(query_link2anchor = query_links, subject_link2anchor = subject_links))
}

#' Core function to link genes to anchors
.link_genes_to_anchors_core <- function(genes, anchors, type) {
  if (nrow(anchors) == 0L) return(NULL)
  
  # Merge anchor positions
  anchor_pos <- merge(anchors, genes, 
                     by.x = if(type == "query") "query_gene" else "subject_gene",
                     by.y = if(type == "query") "query_gene" else "subject_gene",
                     all.x = TRUE)
  
  chr_col <- if(type == "query") "query_chr.x" else "subject_chr.x"
  anchor_pos <- anchor_pos[!is.na(anchor_pos[[paste0(if(type == "query") "query" else "subject", "_start")]]) &
                          !is.na(anchor_pos[[chr_col]]) &
                          anchor_pos[[chr_col]] == anchor_pos[[paste0(if(type == "query") "query" else "subject", "_chr.y")]], ]
  
  if (nrow(anchor_pos) == 0L) return(NULL)
  
  genes_chr_col <- if(type == "query") "query_chr" else "subject_chr"
  # For each chromosome, find nearest anchor for each gene
  links_list <- lapply(split(genes, genes[[genes_chr_col]]), function(chr_genes) {
    chr <- chr_genes[[genes_chr_col]][1L]
    chr_anchors <- anchor_pos[anchor_pos[[chr_col]] == chr, ]
    if (nrow(chr_anchors) == 0L) return(NULL)
    
    gene_anchor <- data.frame(
      gene = chr_genes[[if(type == "query") "query_gene" else "subject_gene"]],
      anchor = NA_integer_,
      stringsAsFactors = FALSE
    )
    
    start_col <- paste0(if(type == "query") "query" else "subject", "_start")
    anchor_start_col <- paste0(if(type == "query") "query" else "subject", "_start")
    anchor_gene_col <- if(type == "query") "query_gene" else "subject_gene"
    
    for (i in seq_len(nrow(chr_genes))) {
      gene_start <- chr_genes[[start_col]][i]
      dists <- abs(chr_anchors[[anchor_start_col]] - gene_start)
      nearest_idx <- which.min(dists)
      gene_anchor$anchor[i] <- chr_anchors[[anchor_gene_col]][nearest_idx]
    }
    
    return(gene_anchor)
  })
  
  links <- do.call(rbind, links_list)
  return(links)
}

#' Find syntenic orthologs - optimized
.find_syntenic_ortho_fast <- function(rbh, anchors, gene_to_anchor, verbose) {
  # Filter: only keep pairs where both genes link to anchors
  if (!is.null(gene_to_anchor$query_link2anchor)) {
    rbh <- merge(rbh, gene_to_anchor$query_link2anchor,
                by.x = "query_gene", by.y = "gene", all.x = FALSE,
                allow.cartesian = TRUE)
    colnames(rbh)[ncol(rbh)] <- "query_anchor"
  } else {
    rbh$query_anchor <- NA_integer_
  }
  
  if (!is.null(gene_to_anchor$subject_link2anchor)) {
    rbh <- merge(rbh, gene_to_anchor$subject_link2anchor,
                by.x = "subject_gene", by.y = "gene", all.x = FALSE,
                allow.cartesian = TRUE)
    colnames(rbh)[ncol(rbh)] <- "subject_anchor"
  } else {
    rbh$subject_anchor <- NA_integer_
  }
  
  rbh <- rbh[!is.na(query_anchor) & !is.na(subject_anchor)]
  if (nrow(rbh) == 0L) return(NULL)
  
  # Check if anchor pair exists
  anchor_pairs <- paste(anchors$query_gene, anchors$subject_gene, sep = "_")
  rbh$anchor_pair <- paste(rbh$query_anchor, rbh$subject_anchor, sep = "_")
  rbh$is_anchor_pair <- rbh$anchor_pair %in% anchor_pairs
  
  # Filter: only keep syntenic pairs
  rbh <- rbh[is_anchor_pair == TRUE]
  if (nrow(rbh) == 0L) return(NULL)
  
  rbh$anchor_pair <- NULL
  
  return(rbh)
}

#' Find synteny blocks - optimized
.find_synteny_blocks_fast <- function(orthopair, verbose) {
  if (nrow(orthopair) == 0L) return(orthopair)
  
  # Order by query genome, chromosome, and position
  orthopair <- orthopair[order(orthopair$query_genome, orthopair$query_chr, orthopair$query_gene), ]
  
  # Find synteny blocks for query genome
  query_blocks <- .find_synteny_blocks_core_fast(
    orthopair,
    genome_col = "query_genome",
    chr_col = "query_chr",
    gene_col = "query_gene",
    anchor_col = "query_anchor"
  )
  
  # Order by subject genome, chromosome, and position
  orthopair <- orthopair[order(orthopair$subject_genome, orthopair$subject_chr, orthopair$subject_gene), ]
  
  # Find synteny blocks for subject genome
  subject_blocks <- .find_synteny_blocks_core_fast(
    orthopair,
    genome_col = "subject_genome",
    chr_col = "subject_chr",
    gene_col = "subject_gene",
    anchor_col = "subject_anchor"
  )
  
  # Merge back
  orthopair <- orthopair[order(orthopair$query_genome, orthopair$query_chr, orthopair$query_gene), ]
  orthopair$query_synteny_block <- query_blocks[match(
    paste(orthopair$query_genome, orthopair$query_chr, orthopair$query_gene, sep = "_"),
    names(query_blocks)
  )]
  
  orthopair$subject_synteny_block <- subject_blocks[match(
    paste(orthopair$subject_genome, orthopair$subject_chr, orthopair$subject_gene, sep = "_"),
    names(subject_blocks)
  )]
  
  return(orthopair)
}

#' Core function for finding synteny blocks - optimized
.find_synteny_blocks_core_fast <- function(orthopair, genome_col, chr_col, gene_col, anchor_col) {
  block_id <- rep(NA_integer_, nrow(orthopair))
  
  # Process each genome-chromosome combination
  genome_chr <- paste(orthopair[[genome_col]], orthopair[[chr_col]], sep = "_")
  for (gc in unique(genome_chr)) {
    idx <- genome_chr == gc
    if (sum(idx) == 0L) next
    
    anchors <- orthopair[[anchor_col]][idx]
    genes <- orthopair[[gene_col]][idx]
    
    # Find gaps between anchors
    anchor_pos <- which(!is.na(anchors))
    if (length(anchor_pos) == 0L) next
    
    # Create blocks
    current_block <- 1L
    block_id[idx][anchor_pos[1L]] <- current_block
    
    if(length(anchor_pos) > 1){
        for (i in seq(2L, length(anchor_pos))) {
            gap_size <- anchor_pos[i] - anchor_pos[i-1L]
            if (gap_size > 1L) {
                block_id[idx][(anchor_pos[i-1L]+1L):(anchor_pos[i]-1L)] <- current_block
                current_block <- current_block + 1L
            }
            block_id[idx][anchor_pos[i]] <- current_block
        }
    }
    
    # Assign remaining genes to last block
    if (length(anchor_pos) > 0L) {
      last_anchor <- anchor_pos[length(anchor_pos)]
      if (last_anchor < length(block_id[idx])) {
        block_id[idx][(last_anchor+1L):length(block_id[idx])] <- current_block
      }
    }
  }
  
  # Create unique block IDs
  names(block_id) <- paste(orthopair[[genome_col]], orthopair[[chr_col]], 
                          orthopair[[gene_col]], sep = "_")
  
  return(block_id)
}

#' Classify ortholog pairs - optimized (uses igraph only for classification)
.classify_orthopair_fast <- function(orthopair, verbose) {
  if (nrow(orthopair) == 0L) {
    orthopair$class <- character(0L)
    orthopair$SOG <- integer(0L)
    return(orthopair)
  }
  
  # Create graph from ortholog pairs
  edges <- data.frame(
    from = paste0("q_", orthopair$query_gene),
    to = paste0("s_", orthopair$subject_gene),
    stringsAsFactors = FALSE
  )
  
  g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  comp <- igraph::components(g)
  
  # Count genes per component
  comp_df <- data.frame(
    gene = names(comp$membership),
    comp_id = comp$membership,
    stringsAsFactors = FALSE
  )
  comp_df$is_query <- startsWith(comp_df$gene, "q_")
  comp_df$is_subject <- startsWith(comp_df$gene, "s_")
  
  n_query <- tapply(comp_df$is_query, comp_df$comp_id, sum)
  n_subject <- tapply(comp_df$is_subject, comp_df$comp_id, sum)
  
  comp_df$n_query <- n_query[comp_df$comp_id]
  comp_df$n_subject <- n_subject[comp_df$comp_id]
  
  # Classify
  comp_df$class <- NA_character_
  comp_df$class[comp_df$n_query == 1L & comp_df$n_subject == 1L] <- "1to1"
  comp_df$class[comp_df$n_query == 1L & comp_df$n_subject != 1L] <- "1toM"
  comp_df$class[comp_df$n_query != 1L & comp_df$n_subject == 1L] <- "Mto1"
  comp_df$class[comp_df$n_query != 1L & comp_df$n_subject != 1L] <- "MtoM"
  
  # Map back to orthopair
  orthopair$query_gene_label <- paste0("q_", orthopair$query_gene)
  orthopair$subject_gene_label <- paste0("s_", orthopair$subject_gene)
  
  query_match <- match(orthopair$query_gene_label, comp_df$gene)
  orthopair$SOG <- comp_df$comp_id[query_match]
  orthopair$class <- comp_df$class[query_match]
  
  orthopair$query_gene_label <- NULL
  orthopair$subject_gene_label <- NULL
  
  return(orthopair)
}

#' Filter ortholog pairs - optimized
.filter_orthopair_fast <- function(orthopair, verbose) {
  if (nrow(orthopair) == 0L) return(orthopair)
  
  # Filter by CI thresholds (lower quartile - 1.5*IQR)
  sog_best_q2s_dt <- orthopair[, .(max_ci = max(q2s_ci, na.rm = TRUE)), by = SOG]
  sog_best_q2s <- sog_best_q2s_dt$max_ci
  names(sog_best_q2s) <- as.character(sog_best_q2s_dt$SOG)
  
  q <- quantile(sog_best_q2s, c(0.25, 0.75), na.rm = TRUE)
  q2s_threshold <- q[1] - 1.5 * diff(q)
  q2s_filter <- orthopair$q2s_ci >= q2s_threshold
  
  sog_best_s2q_dt <- orthopair[, .(max_ci = max(s2q_ci, na.rm = TRUE)), by = SOG]
  sog_best_s2q <- sog_best_s2q_dt$max_ci
  names(sog_best_s2q) <- as.character(sog_best_s2q_dt$SOG)
  
  q <- quantile(sog_best_s2q, c(0.25, 0.75), na.rm = TRUE)
  s2q_threshold <- q[1] - 1.5 * diff(q)
  s2q_filter <- orthopair$s2q_ci >= s2q_threshold
  
  sog_best_pident_dt <- orthopair[, .(max_pident = max(pident, na.rm = TRUE)), by = SOG]
  sog_best_pident <- sog_best_pident_dt$max_pident
  names(sog_best_pident) <- as.character(sog_best_pident_dt$SOG)
  
  q <- quantile(sog_best_pident, c(0.25, 0.75), na.rm = TRUE)
  pident_threshold <- q[1] - 1.5 * diff(q)
  pident_filter <- orthopair$pident >= pident_threshold
  
  # Keep if any threshold is met
  keep <- q2s_filter | s2q_filter | pident_filter
  orthopair <- orthopair[keep]
  
  # Additional filter for non-anchor pairs
  non_anchor <- !orthopair$is_anchor_pair
  if (any(non_anchor)) {
    sog_best_pident_nonanchor_dt <- orthopair[non_anchor, .(max_pident = max(pident, na.rm = TRUE)), by = SOG]
    sog_best_pident_nonanchor <- sog_best_pident_nonanchor_dt$max_pident
    names(sog_best_pident_nonanchor) <- as.character(sog_best_pident_nonanchor_dt$SOG)
    
    q <- quantile(sog_best_pident_nonanchor, c(0.25, 0.75), na.rm = TRUE)
    random_threshold <- q[1] - 1.5 * diff(q)
    keep_nonanchor <- !non_anchor | orthopair$pident >= random_threshold
    orthopair <- orthopair[keep_nonanchor]
  }
  
  return(orthopair)
}
