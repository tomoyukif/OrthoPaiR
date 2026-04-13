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
#' @return Invisibly `TRUE`. Writes per-pair ortholog TSVs under `working_dir/orthopair/`,
#'   and after all pairs finish, writes `orthopair_pairwise_mutual_ci_stats.tsv`
#'   (mean, sd, median, quantiles of `mutual_ci` per genome pair, from direct ortholog rows only)
#'   and `orthopair_genome_mean_mutual_ci_matrix.tsv` (symmetric matrix of mean `mutual_ci`).
#'
#' @importFrom data.table fread fwrite setDT setkey setorder setnames
#' @importFrom parallel mclapply
#' @importFrom stats median quantile sd
#' @export
orthopair <- function(working_dir,
                      genome_width = 4L,
                      pident_threshold_quantile = 0.25,
                      mutual_ci_min = NULL,
                      n_threads = 1L,
                      load_all_gff = TRUE,
                      verbose = TRUE) {
    rbh_dir <- file.path(working_dir, "rbh")
    
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
        genomes <- unlist(strsplit(pair_id, split = "_"))
        if(genomes[1] == genomes[2]){
            return(NULL)
        }
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
            # Column names are already genome1_* / genome2_* (see .process_one_pair_dataframe)
            
            # Write result to TSV file
            out_tsv <- file.path(orthopair_dir, paste0(pair_id, ".tsv"))
            if (!is.null(result) && nrow(result) > 0L) {
                data.table::fwrite(result, file = out_tsv, sep = "\t",
                                   quote = FALSE, row.names = FALSE)
            } else {
                # Write empty file with header
                if (!is.null(result)) {
                    empty_df <- result[0L, , drop = FALSE]
                } else {
                    empty_df <- data.frame(
                        genome1_gene = character(0),
                        genome2_gene = character(0)
                    )
                }
                write.table(empty_df, file = out_tsv, sep = "\t", quote = FALSE,
                            row.names = FALSE, col.names = TRUE)
            }
            
            return(list(pair_id = pair_id, error = NULL))
        }, error = function(e) {
            out_tsv <- file.path(orthopair_dir, paste0(pair_id, ".tsv"))
            writeLines("genome1_gene\tgenome2_gene", con = out_tsv)
            return(list(pair_id = pair_id, error = e$message))
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
    
    # Summary statistics
    n_total <- length(rbh_files)
    n_success <- sum(!error_idx)
    n_error <- sum(error_idx)
    
    if (verbose) {
        message("[ortho] Processing summary:")
        message("[ortho]   Total pairs: ", n_total)
        message("[ortho]   Pairs with results: ", n_success)
        message("[ortho]   Pairs with errors: ", n_error)
        message("[ortho] Output TSV files written to: ", orthopair_dir)
    }
    
    .summarizeOrthopairMutualCi(orthopair_dir = orthopair_dir, verbose = verbose)
    
    invisible(TRUE)
}


#' Summarize mutual_ci from written orthopair TSVs (direct scores only).
#'
#' @param orthopair_dir Directory containing `<genomeA>_<genomeB>.tsv` files.
#' @param verbose Logical.
#' @return Invisibly `NULL` if nothing to summarize; otherwise invisibly `NULL` after writing files.
#' @keywords internal
.summarizeOrthopairMutualCi <- function(orthopair_dir, verbose = TRUE) {
    tsv_files <- list.files(orthopair_dir, pattern = "\\.tsv$", full.names = TRUE)
    skip_names <- c(
        "orthopair_pairwise_mutual_ci_stats.tsv",
        "orthopair_genome_mean_mutual_ci_matrix.tsv"
    )
    tsv_files <- tsv_files[!basename(tsv_files) %in% skip_names]
    if (!length(tsv_files)) {
        if (verbose) message("[ortho] mutual_ci summary: no pair TSV files found")
        return(invisible(NULL))
    }
    
    long_parts <- list()
    for (fn in tsv_files) {
        pair_id <- sub("\\.tsv$", "", basename(fn))
        parts <- strsplit(pair_id, "_", fixed = TRUE)[[1L]]
        if (length(parts) != 2L || parts[1L] == parts[2L]) next
        dt <- tryCatch(
            data.table::fread(fn, sep = "\t", header = TRUE, select = "mutual_ci"),
            error = function(e) NULL
        )
        if (is.null(dt) || !("mutual_ci" %in% names(dt))) next
        mc <- suppressWarnings(as.numeric(dt$mutual_ci))
        mc <- mc[is.finite(mc)]
        if (!length(mc)) next
        long_parts[[length(long_parts) + 1L]] <- data.table::data.table(
            genome_a = parts[1L],
            genome_b = parts[2L],
            mutual_ci = mc
        )
    }
    
    if (!length(long_parts)) {
        if (verbose) message("[ortho] mutual_ci summary: no finite mutual_ci values in pair TSVs")
        return(invisible(NULL))
    }
    
    long_dt <- data.table::rbindlist(long_parts, use.names = TRUE)
    g1i <- suppressWarnings(as.integer(long_dt$genome_a))
    g2i <- suppressWarnings(as.integer(long_dt$genome_b))
    use_int <- !anyNA(g1i) && !anyNA(g2i)
    if (use_int) {
        long_dt[, genome_min := pmin(g1i, g2i)]
        long_dt[, genome_max := pmax(g1i, g2i)]
    } else {
        long_dt[, genome_min := pmin(genome_a, genome_b)]
        long_dt[, genome_max := pmax(genome_a, genome_b)]
    }
    
    sum_dt <- long_dt[, {
        n <- .N
        list(
            n_pairs = n,
            mean = mean(mutual_ci),
            sd = if (n > 1L) stats::sd(mutual_ci) else NA_real_,
            median = stats::median(mutual_ci),
            q10 = stats::quantile(mutual_ci, 0.10, na.rm = TRUE, names = FALSE, type = 7),
            q25 = stats::quantile(mutual_ci, 0.25, na.rm = TRUE, names = FALSE, type = 7),
            q75 = stats::quantile(mutual_ci, 0.75, na.rm = TRUE, names = FALSE, type = 7),
            q90 = stats::quantile(mutual_ci, 0.90, na.rm = TRUE, names = FALSE, type = 7)
        )
    }, by = .(genome_min, genome_max)]
    
    all_gid <- sort(unique(c(sum_dt$genome_min, sum_dt$genome_max)))
    if (use_int) {
        w <- max(4L, max(nchar(as.character(all_gid)), na.rm = TRUE))
        fmt <- paste0("%0", w, "d")
        sum_dt[, genome1 := sprintf(fmt, genome_min)]
        sum_dt[, genome2 := sprintf(fmt, genome_max)]
    } else {
        sum_dt[, genome1 := as.character(genome_min)]
        sum_dt[, genome2 := as.character(genome_max)]
    }
    
    stats_fn <- file.path(orthopair_dir, "orthopair_pairwise_mutual_ci_stats.tsv")
    data.table::fwrite(
        sum_dt[, .(genome1, genome2, n_pairs, mean, sd, median, q10, q25, q75, q90)],
        stats_fn,
        sep = "\t",
        quote = FALSE,
        na = ""
    )
    
    lab <- sort(unique(c(sum_dt$genome1, sum_dt$genome2)))
    mat <- matrix(NA_real_, nrow = length(lab), ncol = length(lab),
                  dimnames = list(lab, lab))
    for (k in seq_len(nrow(sum_dt))) {
        i <- sum_dt$genome1[k]
        j <- sum_dt$genome2[k]
        v <- sum_dt$mean[k]
        mat[i, j] <- v
        mat[j, i] <- v
    }
    mat_df <- cbind(genome = rownames(mat), as.data.frame(mat, stringsAsFactors = FALSE))
    rownames(mat_df) <- NULL
    mat_fn <- file.path(orthopair_dir, "orthopair_genome_mean_mutual_ci_matrix.tsv")
    data.table::fwrite(mat_df, mat_fn, sep = "\t", quote = FALSE, na = "")
    
    if (verbose) {
        message("[ortho] Wrote mutual_ci summary: ", basename(stats_fn))
        message("[ortho] Wrote mean mutual_ci matrix: ", basename(mat_fn))
    }
    invisible(NULL)
}

#' Create GFF lookup: genome ID -> GFF file path
#' RDS can be list(gff_df, cds_df); we use the first element for lookup.
.create_gff_lookup_dataframe <- function(gff_df_paths, genome_width, verbose) {
    lookup <- list()
    
    for (path in gff_df_paths) {
        tryCatch({
            gff <- readRDS(path)
            if (is.list(gff) && length(gff) >= 1L) gff <- gff[[1L]]
            if (is.data.frame(gff) && nrow(gff) > 0L) {
                tx_idx <- as.character(gff$tx_index[1L])
                if (nchar(tx_idx) >= genome_width) {
                    genome_id <- as.integer(substr(tx_idx, 1L, genome_width))
                    lookup[[as.character(genome_id)]] <- path
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
#' RDS can be list(gff_df, cds_df); we return the full list per path so callers can use both.
.load_gff_data_dataframe <- function(gff_df_paths, genome_width) {
    gff_list <- lapply(gff_df_paths, function(path) {
        gff <- readRDS(path)
        if (is.list(gff) && length(gff) >= 1L && is.data.frame(gff[[1L]])) return(gff)
        if (is.data.frame(gff)) return(list(gff, NULL))
        return(NULL)
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
    rbh <- fread(file = rbh_fn,
                 header = FALSE,
                 sep = "\t",
                 select = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L),
                 col.names = c("genome1_tx", "genome2_tx", "pident", "q2s_qcovs", 
                               "s2q_qcovs", "q2s_ci", "s2q_ci", "mutual_ci"),
                 na.strings = "",
                 colClasses = c("integer", "integer", "numeric", "integer", 
                                "integer", "numeric", "numeric", "numeric"),
                 verbose = FALSE)
    setDT(rbh)
    
    if (verbose) message("[ortho] Pair ", pair_id, ": Read ", nrow(rbh), " RBH records")
    
    # Filter valid pairs
    rbh <- rbh[!is.na(genome1_tx) & !is.na(genome2_tx) & genome1_tx != genome2_tx]
    if (nrow(rbh) == 0L) {
        if (verbose) message("[ortho] Pair ", pair_id, ": No valid RBH pairs")
        return(NULL)
    }
    
    # Extract genome IDs
    genome1_genome_pre <- as.integer(substr(as.character(rbh$genome1_tx), 1L, genome_width))
    genome2_genome_pre <- as.integer(substr(as.character(rbh$genome2_tx), 1L, genome_width))
    
    # Filter by mutual_ci if specified
    if (!is.null(mutual_ci_min)) {
        rbh <- rbh[mutual_ci >= mutual_ci_min]
        if (nrow(rbh) == 0L) return(NULL)
    }
    
    # Load GFF data for this pair
    if (!load_all_gff || is.null(gff_data_all)) {
        pair_genomes <- unique(c(genome1_genome_pre, genome2_genome_pre))
        gff_data_pair <- lapply(as.character(pair_genomes), function(gid) {
            path <- gff_lookup[[gid]]
            if (is.null(path) || !file.exists(path)) return(NULL)
            return(readRDS(path))
        })
        gff_data_pair <- gff_data_pair[!sapply(gff_data_pair, is.null)]
    } else {
        # Extract GFF data for this pair from pre-loaded data (each element can be list(gff_df, cds_df))
        pair_genomes <- unique(c(genome1_genome_pre, genome2_genome_pre))
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
    
    if (length(gff_data_pair) == 0L) {
        warning("Could not load GFF data for pair: ", pair_id)
        return(NULL)
    }
    
    # Split GFF data into genome1 (min genome id) and genome2 (max genome id)
    genome1_genome_id <- as.character(min(pair_genomes))
    genome2_genome_id <- as.character(max(pair_genomes))
    
    genome1_gff_full <- NULL
    genome2_gff_full <- NULL
    genome1_cds_full <- NULL
    genome2_cds_full <- NULL
    
    for (el in gff_data_pair) {
        gff_df <- if (is.list(el) && length(el) >= 1L) el[[1L]] else el
        cds_df_el <- if (is.list(el) && length(el) >= 2L) el[[2L]] else NULL
        if (is.null(gff_df) || !is.data.frame(gff_df) || nrow(gff_df) == 0L) next
        tx_idx <- as.character(gff_df$tx_index[1L])
        if (nchar(tx_idx) >= genome_width) {
            gid <- as.character(as.integer(substr(tx_idx, 1L, genome_width)))
            if (gid == genome1_genome_id) {
                genome1_gff_full <- gff_df
                genome1_cds_full <- cds_df_el
            } else if (gid == genome2_genome_id) {
                genome2_gff_full <- gff_df
                genome2_cds_full <- cds_df_el
            }
        }
    }
    
    if (is.null(genome1_gff_full) || is.null(genome2_gff_full)) {
        warning("Could not identify genome1/genome2 GFF for pair: ", pair_id)
        return(NULL)
    }
    
    # Order GFF data (as expected by linkGene2Genome)
    genome1_gff_full <- genome1_gff_full[order(genome1_gff_full$seqnames, genome1_gff_full$start), ]
    genome2_gff_full <- genome2_gff_full[order(genome2_gff_full$seqnames, genome2_gff_full$start), ]
    genome1_gff_full$gene_index <- as.integer(genome1_gff_full$gene_index)
    genome1_gff_full$tx_index <- as.integer(genome1_gff_full$tx_index)
    genome2_gff_full$gene_index <- as.integer(genome2_gff_full$gene_index)
    genome2_gff_full$tx_index <- as.integer(genome2_gff_full$tx_index)
    
    # g2g_graph list names genome1_df/genome2_df match R/05_functions_orthology.R (.linkGene2Genome).
    g2g_graph <- .linkGene2Genome_dataframe(genome1_gff_df = genome1_gff_full,
                                            genome2_gff_df = genome2_gff_full,
                                            genome1_cds_df = genome1_cds_full,
                                            genome2_cds_df = genome2_cds_full)
    
    # Convert RBH to format expected by R/05_functions_orthology.R
    genome1_tx_match <- match(rbh$genome1_tx, genome1_gff_full$tx_index)
    genome2_tx_match <- match(rbh$genome2_tx, genome2_gff_full$tx_index)
    
    rbh_df <- data.frame(genome1_tx = as.integer(rbh$genome1_tx),
                         genome2_tx = as.integer(rbh$genome2_tx),
                         genome1_gene = as.integer(genome1_gff_full$gene_index[genome1_tx_match]),
                         genome2_gene = as.integer(genome2_gff_full$gene_index[genome2_tx_match]),
                         ci_q2s = rbh$q2s_ci,
                         ci_s2q = rbh$s2q_ci,
                         pident = rbh$pident,
                         mutual_ci = rbh$mutual_ci,
                         stringsAsFactors = FALSE)
    
    rbh_df <- rbh_df[!is.na(rbh_df$genome1_gene) & !is.na(rbh_df$genome2_gene), ]
    rbh_df$pair_id <- paste(rbh_df$genome1_gene, rbh_df$genome2_gene, sep = "_")
    rbh_df <- subset(rbh_df, subset = !duplicated(rbh_df$pair_id))
    
    if (nrow(rbh_df) == 0L) {
        if (verbose) message("[ortho] Pair ", pair_id, ": No RBH after GFF merge")
        return(NULL)
    }
    
    # 
    # # Apply R/05_functions_orthology.R algorithm from linkGene2Genome onwards
    # # Note: .findAnchors checks length(rbh) == 1, so ensure rbh_df is a data.frame
    anchor <- .findAnchors(rbh = rbh_df, g2g_graph = g2g_graph)
    # 
    # if (all(is.na(anchor)) || (is.data.frame(anchor) && nrow(anchor) == 0L)) {
    #     if (verbose) message("[ortho] Pair ", pair_id, ": No anchors found")
    #     return(NULL)
    # }
    # 
    # t2a_graph <- .link2Anchor(g2g_graph = g2g_graph, anchor = anchor)
    # # .findSyntenicOrtho expects rbh as data.frame and will rename first 2 columns
    # orthopair <- .findSyntenicOrtho(rbh = rbh_df,
    #                                 anchor = anchor,
    #                                 g2g_graph = g2g_graph,
    #                                 t2a_graph = t2a_graph)
    # 
    # if (is.null(orthopair) || nrow(orthopair) == 0L) {
    #     if (verbose) message("[ortho] Pair ", pair_id, ": No syntenic orthologs")
    #     return(NULL)
    # }
    # 
    # orthopair <- .findSyntenyBlocks(orthopair = orthopair)
    
    orthopair <- .findShortestPath(rbh = rbh_df, anchor = anchor, g2g_graph = g2g_graph)
    orthopair <- .pickBestPair(orthopair = orthopair)
    orthopair <- .classifyOrthoPair(orthopair = orthopair)
    orthopair <- .filterOrthopair(orthopair = orthopair, g2g_graph = g2g_graph)
    orthopair <- .classifyOrthoPair(orthopair = orthopair)
    # .reformatOrthoPair() was modified!!!
    orthopair <- .reformatOrthoPair(orthopair = orthopair, g2g_graph = g2g_graph) 
    
    # Convert to data.table for consistency
    setDT(orthopair)
    
    if (verbose) {
        message("[ortho] Pair ", pair_id, ": ", nrow(orthopair), " ortholog pairs")
    }
    
    return(orthopair)
}

.findShortestPath <- function(rbh, anchor, g2g_graph){
    rbh$genome1_chr <- g2g_graph$genome1_df$seqnames[match(rbh$genome1_gene, g2g_graph$genome1_df$gene_index)]
    rbh$genome2_chr <- g2g_graph$genome2_df$seqnames[match(rbh$genome2_gene, g2g_graph$genome2_df$gene_index)]
    anchor$genome1_chr <- g2g_graph$genome1_df$seqnames[match(anchor$genome1_gene, g2g_graph$genome1_df$gene_index)]
    anchor$genome2_chr <- g2g_graph$genome2_df$seqnames[match(anchor$genome2_gene, g2g_graph$genome2_df$gene_index)]
    anchor$index <- seq_along(anchor$genome1_tx)
    rbh$is_anchor_pair <- rbh$pair_id %in% anchor$pair_id
    anchor_out <- rbh_shortest_path(qgeneid = anchor$genome1_gene, 
                                    sgeneid = anchor$genome2_gene,
                                    q_chr = anchor$genome1_chr, 
                                    s_chr = anchor$genome2_chr,
                                    anchor_qgeneid = anchor$genome1_gene,
                                    anchor_sgeneid = anchor$genome2_gene,
                                    anchor_q_chr = anchor$genome1_chr,
                                    anchor_s_chr = anchor$genome2_chr,
                                    omit_zero = TRUE)
    q <- quantile(anchor_out$best_dist, na.rm = T, c(0.25, 0.75))
    dist_threshold <- q[2] + diff(q) * 1.5
    out <- rbh_shortest_path(qgeneid = rbh$genome1_gene, 
                             sgeneid = rbh$genome2_gene,
                             q_chr = rbh$genome1_chr, 
                             s_chr = rbh$genome2_chr,
                             anchor_qgeneid = anchor$genome1_gene,
                             anchor_sgeneid = anchor$genome2_gene,
                             anchor_q_chr = anchor$genome1_chr,
                             anchor_s_chr = anchor$genome2_chr,
                             omit_zero = FALSE)
    rbh$best_dist <- out$best_dist
    rbh$rbh_shortest_path <- as.integer(factor(x = out$shortest_path,
                                               levels = unique(out$shortest_path)))
    rbh$index <- seq_along(rbh$genome1_tx)
    out <- subset(rbh, subset = best_dist < dist_threshold | is_anchor_pair, 
                  select = -c(best_dist, rbh_shortest_path, index))
    names(out)[1:4] <- c("genome1_tx", "genome2_tx", "genome1_gene", "genome2_gene")
    return(out)
}

#' Create g2g_graph structure from GFF data (equivalent to linkGene2Genome).
#' genome1_* is the query side, genome2_* the subject side; list element names
#' List names match R/05_functions_orthology.R `.linkGene2Genome()`.
.linkGene2Genome_dataframe <- function(genome1_gff_df, genome2_gff_df,
                                        genome1_cds_df = NULL, genome2_cds_df = NULL) {
    has_type <- "type" %in% names(genome1_gff_df)
    
    if (!is.null(genome1_cds_df) && !is.null(genome2_cds_df) &&
        is.data.frame(genome1_cds_df) && is.data.frame(genome2_cds_df)) {
        genome1_df <- as.data.frame(genome1_gff_df)
        genome2_df <- as.data.frame(genome2_gff_df)
        if ("transcript_id" %in% names(genome1_df) && !("ID" %in% names(genome1_df))) genome1_df$ID <- genome1_df$transcript_id
        if ("transcript_id" %in% names(genome2_df) && !("ID" %in% names(genome2_df))) genome2_df$ID <- genome2_df$transcript_id
        genome1_gff_cds <- as.data.frame(genome1_cds_df)
        genome2_gff_cds <- as.data.frame(genome2_cds_df)
        if (nrow(genome1_gff_cds) > 0L && "gene_id" %in% names(genome1_gff_cds)) {
            hit <- match(genome1_gff_cds$gene_id, genome1_gff_df$gene_id)
            if ("tx_index" %in% names(genome1_gff_df)) genome1_gff_cds$tx_index <- genome1_gff_df$tx_index[hit]
            if ("gene_index" %in% names(genome1_gff_df)) genome1_gff_cds$gene_index <- genome1_gff_df$gene_index[hit]
            if ("strand" %in% names(genome1_gff_df)) genome1_gff_cds$strand <- genome1_gff_df$strand[hit]
            if (!("Parent" %in% names(genome1_gff_cds)) && "transcript_id" %in% names(genome1_gff_df)) {
                genome1_gff_cds$Parent <- genome1_gff_df$transcript_id[hit]
            }
        }
        if (nrow(genome2_gff_cds) > 0L && "gene_id" %in% names(genome2_gff_cds)) {
            hit <- match(genome2_gff_cds$gene_id, genome2_gff_df$gene_id)
            if ("tx_index" %in% names(genome2_gff_df)) genome2_gff_cds$tx_index <- genome2_gff_df$tx_index[hit]
            if ("gene_index" %in% names(genome2_gff_df)) genome2_gff_cds$gene_index <- genome2_gff_df$gene_index[hit]
            if ("strand" %in% names(genome2_gff_df)) genome2_gff_cds$strand <- genome2_gff_df$strand[hit]
            if (!("Parent" %in% names(genome2_gff_cds)) && "transcript_id" %in% names(genome2_gff_df)) {
                genome2_gff_cds$Parent <- genome2_gff_df$transcript_id[hit]
            }
        }
    } else if (has_type) {
        q_tx_i <- genome1_gff_df$type %in% c("transcript", "mRNA")
        genome1_df <- as.data.frame(genome1_gff_df[q_tx_i, ])
        q_cds_i <- genome1_gff_df$type %in% "CDS"
        genome1_gff_cds <- genome1_gff_df[q_cds_i, ]
        if (nrow(genome1_gff_cds) > 0L && "Parent" %in% names(genome1_gff_cds)) {
            hit <- match(unlist(genome1_gff_cds$Parent), genome1_gff_df$ID)
            genome1_gff_cds$tx_index <- genome1_gff_df$tx_index[hit]
        }
        genome1_gff_cds <- as.data.frame(genome1_gff_cds)
        s_tx_i <- genome2_gff_df$type %in% c("transcript", "mRNA")
        genome2_df <- as.data.frame(genome2_gff_df[s_tx_i, ])
        s_cds_i <- genome2_gff_df$type %in% "CDS"
        genome2_gff_cds <- genome2_gff_df[s_cds_i, ]
        if (nrow(genome2_gff_cds) > 0L && "Parent" %in% names(genome2_gff_cds)) {
            hit <- match(unlist(genome2_gff_cds$Parent), genome2_gff_df$ID)
            genome2_gff_cds$tx_index <- genome2_gff_df$tx_index[hit]
        }
        genome2_gff_cds <- as.data.frame(genome2_gff_cds)
    } else {
        genome1_df <- as.data.frame(genome1_gff_df)
        genome2_df <- as.data.frame(genome2_gff_df)
        if ("transcript_id" %in% names(genome1_df) && !("ID" %in% names(genome1_df))) genome1_df$ID <- genome1_df$transcript_id
        if ("transcript_id" %in% names(genome2_df) && !("ID" %in% names(genome2_df))) genome2_df$ID <- genome2_df$transcript_id
        genome1_gff_cds <- data.frame()
        genome2_gff_cds <- data.frame()
    }
    
    list(genome1_df = genome1_df,
         genome2_df = genome2_df,
         genome1_gff = genome1_gff_cds,
         genome2_gff = genome2_gff_cds)
}
