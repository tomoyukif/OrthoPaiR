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
    if(length(rbh_files) > 0L) {
        pair_ids <- sub("\\.rbh$", "", basename(rbh_files))
        parts <- strsplit(pair_ids, "_", fixed = TRUE)
        same_genome <- lengths(parts) == 2L & vapply(parts, function(p) p[1L] == p[2L], logical(1L))
        rbh_files <- rbh_files[!same_genome]
    }
    if(length(rbh_files) == 0L) {
        warning("No RBH files found in: ", rbh_dir, " (or only same-genome comparisons)")
        return(list())
    }
    
    if(verbose) message("[ortho] Found ", length(rbh_files), " RBH pair files")
    
    # Get GFF data paths - search in working_dir/input/<index>_<genome_name>/ folders
    input_dir <- file.path(working_dir, "input")
    if(!dir.exists(input_dir)) {
        stop("Input directory not found: ", input_dir)
    }
    
    # Find all folders matching pattern <index>_<genome_name>
    input_folders <- list.dirs(input_dir, full.names = TRUE, recursive = FALSE)
    # Filter folders that match pattern (have underscore and contain gff_df.rds)
    gff_df_paths <- character(0)
    for(folder in input_folders) {
        gff_fn <- file.path(folder, "gff_df.rds")
        if(file.exists(gff_fn)) {
            gff_df_paths <- c(gff_df_paths, gff_fn)
        }
    }
    
    if(length(gff_df_paths) == 0L) {
        stop("No gff_df.rds files found in ", input_dir, 
             " subdirectories (expected pattern: <index>_<genome_name>/gff_df.rds)")
    }
    
    if(verbose) message("[ortho] Found ", length(gff_df_paths), " GFF files")
    
    # Create output directory for ortholog pairs
    orthopair_dir <- file.path(working_dir, "orthopair")
    dir.create(orthopair_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Create genome ID to GFF path mapping for efficient lookup
    gff_lookup <- .create_gff_lookup_dataframe(gff_df_paths, genome_width, verbose)
    if(length(gff_lookup) == 0L) {
        stop("Failed to create GFF lookup. Check that GFF files contain valid tx_index values.")
    }
    
    # Load GFF data once if requested (faster but uses more memory)
    gff_data_all <- NULL
    if(load_all_gff) {
        if(verbose) message("[ortho] Loading all GFF data (memory-intensive mode)")
        gff_data_all <- .load_gff_data_dataframe(gff_df_paths, genome_width)
        if(is.null(gff_data_all) || length(gff_data_all) == 0L) {
            stop("Failed to load GFF data. Check GFF file format.")
        }
        if(verbose) message("[ortho] Loaded GFF data from ", length(gff_df_paths), " genomes")
    }
    
    # Process each RBH file in parallel
    if(verbose) message("[ortho] Processing ", length(rbh_files), " genome pairs")
    
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
            if(!is.null(result) && nrow(result) > 0L) {
                fwrite(result, file = out_tsv, sep = "\t",
                                   quote = FALSE, row.names = FALSE)
            } else {
                # Write empty file with header
                if(!is.null(result)) {
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
    if(n_threads > 1L && .Platform$OS.type != "windows") {
        if(verbose) message("[ortho] Using parallel processing with ", n_threads, " cores")
        results_list <- mclapply(rbh_files, process_one_rbh_file, mc.cores = n_threads)
        
    } else {
        if(verbose) message("[ortho] Using sequential processing")
        results_list <- lapply(rbh_files, process_one_rbh_file)
    }
    
    # Extract results and errors
    errors <- lapply(results_list, function(x) x$error)
    pair_ids <- sapply(results_list, function(x) x$pair_id)
    
    # Report errors
    error_idx <- !sapply(errors, is.null)
    if(any(error_idx)) {
        error_pairs <- pair_ids[error_idx]
        error_msgs <- errors[error_idx]
        warning("[ortho] Errors occurred in ", sum(error_idx), " pairs:\n",
                paste(sprintf("  %s: %s", error_pairs, error_msgs), collapse = "\n"))
    }
    
    # Summary statistics
    n_total <- length(rbh_files)
    n_success <- sum(!error_idx)
    n_error <- sum(error_idx)
    
    if(verbose) {
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
    if(!length(tsv_files)) {
        if(verbose) message("[ortho] mutual_ci summary: no pair TSV files found")
        return(invisible(NULL))
    }
    
    long_parts <- list()
    for(fn in tsv_files) {
        pair_id <- sub("\\.tsv$", "", basename(fn))
        parts <- strsplit(pair_id, "_", fixed = TRUE)[[1L]]
        if(length(parts) != 2L || parts[1L] == parts[2L]) next
        dt <- tryCatch(
            fread(fn, sep = "\t", header = TRUE, select = "mutual_ci"),
            error = function(e) NULL
        )
        if(is.null(dt) || !("mutual_ci" %in% names(dt))) next
        mc <- suppressWarnings(as.numeric(dt$mutual_ci))
        mc <- mc[is.finite(mc)]
        if(!length(mc)) next
        long_parts[[length(long_parts) + 1L]] <- data.table(
            genome_a = parts[1L],
            genome_b = parts[2L],
            mutual_ci = mc
        )
    }
    
    if(!length(long_parts)) {
        if(verbose) message("[ortho] mutual_ci summary: no finite mutual_ci values in pair TSVs")
        return(invisible(NULL))
    }
    
    long_dt <- rbindlist(long_parts, use.names = TRUE)
    g1i <- suppressWarnings(as.integer(long_dt$genome_a))
    g2i <- suppressWarnings(as.integer(long_dt$genome_b))
    use_int <- !anyNA(g1i) && !anyNA(g2i)
    if(use_int) {
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
            sd = if(n > 1L) sd(mutual_ci) else NA_real_,
            median = median(mutual_ci),
            q10 = quantile(mutual_ci, 0.10, na.rm = TRUE, names = FALSE, type = 7),
            q25 = quantile(mutual_ci, 0.25, na.rm = TRUE, names = FALSE, type = 7),
            q75 = quantile(mutual_ci, 0.75, na.rm = TRUE, names = FALSE, type = 7),
            q90 = quantile(mutual_ci, 0.90, na.rm = TRUE, names = FALSE, type = 7)
        )
    }, by = .(genome_min, genome_max)]
    
    all_gid <- sort(unique(c(sum_dt$genome_min, sum_dt$genome_max)))
    if(use_int) {
        w <- max(4L, max(nchar(as.character(all_gid)), na.rm = TRUE))
        fmt <- paste0("%0", w, "d")
        sum_dt[, genome1 := sprintf(fmt, genome_min)]
        sum_dt[, genome2 := sprintf(fmt, genome_max)]
    } else {
        sum_dt[, genome1 := as.character(genome_min)]
        sum_dt[, genome2 := as.character(genome_max)]
    }
    
    stats_fn <- file.path(orthopair_dir, "orthopair_pairwise_mutual_ci_stats.tsv")
    fwrite(
        sum_dt[, .(genome1, genome2, n_pairs, mean, sd, median, q10, q25, q75, q90)],
        stats_fn,
        sep = "\t",
        quote = FALSE,
        na = ""
    )
    
    lab <- sort(unique(c(sum_dt$genome1, sum_dt$genome2)))
    mat <- matrix(NA_real_, nrow = length(lab), ncol = length(lab),
                  dimnames = list(lab, lab))
    for(k in seq_len(nrow(sum_dt))) {
        i <- sum_dt$genome1[k]
        j <- sum_dt$genome2[k]
        v <- sum_dt$mean[k]
        mat[i, j] <- v
        mat[j, i] <- v
    }
    mat_df <- cbind(genome = rownames(mat), as.data.frame(mat, stringsAsFactors = FALSE))
    rownames(mat_df) <- NULL
    mat_fn <- file.path(orthopair_dir, "orthopair_genome_mean_mutual_ci_matrix.tsv")
    fwrite(mat_df, mat_fn, sep = "\t", quote = FALSE, na = "")
    
    if(verbose) {
        message("[ortho] Wrote mutual_ci summary: ", basename(stats_fn))
        message("[ortho] Wrote mean mutual_ci matrix: ", basename(mat_fn))
    }
    invisible(NULL)
}

#' Create GFF lookup: genome ID -> GFF file path
#' RDS can be list(gff_df, cds_df); we use the first element for lookup.
.create_gff_lookup_dataframe <- function(gff_df_paths, genome_width, verbose) {
    lookup <- list()
    
    for(path in gff_df_paths) {
        tryCatch({
            gff <- readRDS(path)
            if(is.list(gff) && length(gff) >= 1L) gff <- gff[[1L]]
            if(is.data.frame(gff) && nrow(gff) > 0L) {
                tx_idx <- as.character(gff$tx_index[1L])
                if(nchar(tx_idx) >= genome_width) {
                    genome_id <- as.integer(substr(tx_idx, 1L, genome_width))
                    lookup[[as.character(genome_id)]] <- path
                }
            }
        }, error = function(e) NULL)
    }
    
    if(verbose && length(lookup) > 0L) {
        message("[ortho] Created GFF lookup for ", length(lookup), " genomes")
    }
    
    return(lookup)
}

#' Load GFF data for all genomes (for linkGene2Genome compatibility)
#' RDS can be list(gff_df, cds_df); we return the full list per path so callers can use both.
.load_gff_data_dataframe <- function(gff_df_paths, genome_width) {
    gff_list <- lapply(gff_df_paths, function(path) {
        gff <- readRDS(path)
        if(is.list(gff) && length(gff) >= 1L && is.data.frame(gff[[1L]])) return(gff)
        if(is.data.frame(gff)) return(list(gff, NULL))
        return(NULL)
    })
    
    # Filter NULLs
    gff_list <- gff_list[!sapply(gff_list, is.null)]
    if(length(gff_list) == 0L) {
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
    if(verbose) message("[ortho] Reading RBH file: ", basename(rbh_fn))
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
    
    if(verbose) message("[ortho] Pair ", pair_id, ": Read ", nrow(rbh), " RBH records")
    
    # Filter valid pairs
    rbh <- rbh[!is.na(genome1_tx) & !is.na(genome2_tx) & genome1_tx != genome2_tx]
    if(nrow(rbh) == 0L) {
        if(verbose) message("[ortho] Pair ", pair_id, ": No valid RBH pairs")
        return(NULL)
    }
    
    # Extract genome IDs
    genome1_genome_pre <- as.integer(substr(as.character(rbh$genome1_tx), 1L, genome_width))
    genome2_genome_pre <- as.integer(substr(as.character(rbh$genome2_tx), 1L, genome_width))
    
    # Filter by mutual_ci if specified
    if(!is.null(mutual_ci_min)) {
        rbh <- rbh[mutual_ci >= mutual_ci_min]
        if(nrow(rbh) == 0L) return(NULL)
    }
    
    # Load GFF data for this pair
    if(!load_all_gff || is.null(gff_data_all)) {
        pair_genomes <- unique(c(genome1_genome_pre, genome2_genome_pre))
        gff_data_pair <- lapply(as.character(pair_genomes), function(gid) {
            path <- gff_lookup[[gid]]
            if(is.null(path) || !file.exists(path)) return(NULL)
            return(readRDS(path))
        })
        gff_data_pair <- gff_data_pair[!sapply(gff_data_pair, is.null)]
    } else {
        # Extract GFF data for this pair from pre-loaded data (each element can be list(gff_df, cds_df))
        pair_genomes <- unique(c(genome1_genome_pre, genome2_genome_pre))
        gff_data_pair <- lapply(as.character(pair_genomes), function(gid) {
            path <- gff_lookup[[gid]]
            if(is.null(path)) return(NULL)
            idx <- which(sapply(gff_data_all, function(g) {
                if(is.null(g)) return(FALSE)
                el <- if(is.list(g) && length(g) >= 1L) g[[1L]] else g
                if(!is.data.frame(el) || nrow(el) == 0L) return(FALSE)
                tx_idx <- as.character(el$tx_index[1L])
                if(nchar(tx_idx) >= genome_width) {
                    gid_check <- as.integer(substr(tx_idx, 1L, genome_width))
                    return(gid_check == as.integer(gid))
                }
                return(FALSE)
            }))
            if(length(idx) > 0L) return(gff_data_all[[idx[1L]]])
            return(NULL)
        })
        gff_data_pair <- gff_data_pair[!sapply(gff_data_pair, is.null)]
    }
    
    if(length(gff_data_pair) == 0L) {
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
    
    for(el in gff_data_pair) {
        gff_df <- if(is.list(el) && length(el) >= 1L) el[[1L]] else el
        cds_df_el <- if(is.list(el) && length(el) >= 2L) el[[2L]] else NULL
        if(is.null(gff_df) || !is.data.frame(gff_df) || nrow(gff_df) == 0L) next
        tx_idx <- as.character(gff_df$tx_index[1L])
        if(nchar(tx_idx) >= genome_width) {
            gid <- as.character(as.integer(substr(tx_idx, 1L, genome_width)))
            if(gid == genome1_genome_id) {
                genome1_gff_full <- gff_df
                genome1_cds_full <- cds_df_el
            } else if(gid == genome2_genome_id) {
                genome2_gff_full <- gff_df
                genome2_cds_full <- cds_df_el
            }
        }
    }
    
    if(is.null(genome1_gff_full) || is.null(genome2_gff_full)) {
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
    
    if(nrow(rbh_df) == 0L) {
        if(verbose) message("[ortho] Pair ", pair_id, ": No RBH after GFF merge")
        return(NULL)
    }
    
    # 
    # # Apply R/05_functions_orthology.R algorithm from linkGene2Genome onwards
    # # Note: .findAnchors checks length(rbh) == 1, so ensure rbh_df is a data.frame
    anchor <- .findAnchors(rbh = rbh_df, g2g_graph = g2g_graph)
    # 
    # if(all(is.na(anchor)) || (is.data.frame(anchor) && nrow(anchor) == 0L)) {
    #     if(verbose) message("[ortho] Pair ", pair_id, ": No anchors found")
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
    # if(is.null(orthopair) || nrow(orthopair) == 0L) {
    #     if(verbose) message("[ortho] Pair ", pair_id, ": No syntenic orthologs")
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
    
    if(verbose) {
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
    
    if(!is.null(genome1_cds_df) && !is.null(genome2_cds_df) &&
        is.data.frame(genome1_cds_df) && is.data.frame(genome2_cds_df)) {
        genome1_df <- as.data.frame(genome1_gff_df)
        genome2_df <- as.data.frame(genome2_gff_df)
        if("transcript_id" %in% names(genome1_df) && !("ID" %in% names(genome1_df))) genome1_df$ID <- genome1_df$transcript_id
        if("transcript_id" %in% names(genome2_df) && !("ID" %in% names(genome2_df))) genome2_df$ID <- genome2_df$transcript_id
        genome1_gff_cds <- as.data.frame(genome1_cds_df)
        genome2_gff_cds <- as.data.frame(genome2_cds_df)
        if(nrow(genome1_gff_cds) > 0L && "gene_id" %in% names(genome1_gff_cds)) {
            hit <- match(genome1_gff_cds$gene_id, genome1_gff_df$gene_id)
            if("tx_index" %in% names(genome1_gff_df)) genome1_gff_cds$tx_index <- genome1_gff_df$tx_index[hit]
            if("gene_index" %in% names(genome1_gff_df)) genome1_gff_cds$gene_index <- genome1_gff_df$gene_index[hit]
            if("strand" %in% names(genome1_gff_df)) genome1_gff_cds$strand <- genome1_gff_df$strand[hit]
            if(!("Parent" %in% names(genome1_gff_cds)) && "transcript_id" %in% names(genome1_gff_df)) {
                genome1_gff_cds$Parent <- genome1_gff_df$transcript_id[hit]
            }
        }
        if(nrow(genome2_gff_cds) > 0L && "gene_id" %in% names(genome2_gff_cds)) {
            hit <- match(genome2_gff_cds$gene_id, genome2_gff_df$gene_id)
            if("tx_index" %in% names(genome2_gff_df)) genome2_gff_cds$tx_index <- genome2_gff_df$tx_index[hit]
            if("gene_index" %in% names(genome2_gff_df)) genome2_gff_cds$gene_index <- genome2_gff_df$gene_index[hit]
            if("strand" %in% names(genome2_gff_df)) genome2_gff_cds$strand <- genome2_gff_df$strand[hit]
            if(!("Parent" %in% names(genome2_gff_cds)) && "transcript_id" %in% names(genome2_gff_df)) {
                genome2_gff_cds$Parent <- genome2_gff_df$transcript_id[hit]
            }
        }
    } else if(has_type) {
        q_tx_i <- genome1_gff_df$type %in% c("transcript", "mRNA")
        genome1_df <- as.data.frame(genome1_gff_df[q_tx_i, ])
        q_cds_i <- genome1_gff_df$type %in% "CDS"
        genome1_gff_cds <- genome1_gff_df[q_cds_i, ]
        if(nrow(genome1_gff_cds) > 0L && "Parent" %in% names(genome1_gff_cds)) {
            hit <- match(unlist(genome1_gff_cds$Parent), genome1_gff_df$ID)
            genome1_gff_cds$tx_index <- genome1_gff_df$tx_index[hit]
        }
        genome1_gff_cds <- as.data.frame(genome1_gff_cds)
        s_tx_i <- genome2_gff_df$type %in% c("transcript", "mRNA")
        genome2_df <- as.data.frame(genome2_gff_df[s_tx_i, ])
        s_cds_i <- genome2_gff_df$type %in% "CDS"
        genome2_gff_cds <- genome2_gff_df[s_cds_i, ]
        if(nrow(genome2_gff_cds) > 0L && "Parent" %in% names(genome2_gff_cds)) {
            hit <- match(unlist(genome2_gff_cds$Parent), genome2_gff_df$ID)
            genome2_gff_cds$tx_index <- genome2_gff_df$tx_index[hit]
        }
        genome2_gff_cds <- as.data.frame(genome2_gff_cds)
    } else {
        genome1_df <- as.data.frame(genome1_gff_df)
        genome2_df <- as.data.frame(genome2_gff_df)
        if("transcript_id" %in% names(genome1_df) && !("ID" %in% names(genome1_df))) genome1_df$ID <- genome1_df$transcript_id
        if("transcript_id" %in% names(genome2_df) && !("ID" %in% names(genome2_df))) genome2_df$ID <- genome2_df$transcript_id
        genome1_gff_cds <- data.frame()
        genome2_gff_cds <- data.frame()
    }
    
    list(genome1_df = genome1_df,
         genome2_df = genome2_df,
         genome1_gff = genome1_gff_cds,
         genome2_gff = genome2_gff_cds)
}

.findAnchors <- function(rbh, g2g_graph){
    if(length(rbh) == 1){
        return(NA)
    }
    
    rbbh <- .getRBBH(rbh = rbh)
    
    root_hit <- match(rbbh$genome1_gene, g2g_graph$genome1_df$gene_index)
    leaf_hit <- match(rbbh$genome2_gene, g2g_graph$genome2_df$gene_index)
    anchor <- data.frame(root = g2g_graph$genome1_df$gene_index[root_hit],
                         root_anchor = g2g_graph$genome1_df$gene_index[root_hit],
                         root_anchor_chr = g2g_graph$genome1_df$seqnames[root_hit],
                         leaf_anchor = g2g_graph$genome2_df$gene_index[leaf_hit],
                         leaf_anchor_chr = g2g_graph$genome2_df$seqnames[leaf_hit],
                         subset(rbbh, select = c(genome1_tx:genome2_gene, ci_q2s:pair_id)))
    anchor <- unique(anchor)
    
    anchor <- .checkHighCopyGenes(anchor = anchor, rbh = rbh)
    
    q <- quantile(anchor$pident, c(0.25, 0.75))
    anchor_threshold <- q[1] - 1.5 * diff(q)
    anchor <- anchor[anchor$pident >= anchor_threshold, ]
    
    return(anchor)
}

.pickBestPair <- function(orthopair){
    orthopair <- orthopair[order(orthopair$mutual_ci, decreasing = TRUE), ]
    orthopair <- orthopair[!duplicated(orthopair$pair_id), ]
    return(orthopair)
}

.classifyOrthoPair <- function(orthopair){
    d <- subset(orthopair, select = c(genome1_gene, genome2_gene))
    d$genome1_gene <- paste0("g1_", d$genome1_gene)
    d$genome2_gene <- paste0("g2_", d$genome2_gene)
    genome1_gene_list <- unique(d$genome1_gene)
    genome2_gene_list <- unique(d$genome2_gene)
    g <- graph_from_data_frame(d = d, directed = FALSE)
    grp <- split(V(g)$name, components(g)$membership)
    names(grp) <- paste0(names(grp), "_")
    grp <- unlist(grp)
    grp <- data.frame(grp = names(grp), gene_id = grp)
    grp$grp <- as.numeric(sub("_.*", "", grp$grp))
    grp$genome1 <- grp$gene_id %in% genome1_gene_list
    grp$genome2 <- grp$gene_id %in% genome2_gene_list
    n_genome1 <- tapply(grp$genome1, grp$grp, sum)
    hit <- match(grp$grp, as.numeric(names(n_genome1)))
    grp$n_genome1 <- n_genome1[hit]
    n_genome2 <- tapply(grp$genome2, grp$grp, sum)
    hit <- match(grp$grp, as.numeric(names(n_genome2)))
    grp$n_genome2 <- n_genome2[hit]
    
    sog_1to1 <- which(grp$n_genome1 == 1 & grp$n_genome2 == 1)
    sog_1toM <- which(grp$n_genome1 == 1 & grp$n_genome2 != 1)
    sog_Mto1 <- which(grp$n_genome1 != 1 & grp$n_genome2 == 1)
    sog_MtoM <- which(grp$n_genome1 != 1 & grp$n_genome2 != 1)
    orthopair$class <- NA
    hit <- d$genome1_gene %in% grp$gene_id[sog_1to1]
    orthopair$class[hit] <- "1to1"
    hit <- d$genome1_gene %in% grp$gene_id[sog_1toM]
    orthopair$class[hit] <- "1toM"
    hit <- d$genome1_gene %in% grp$gene_id[sog_Mto1]
    orthopair$class[hit] <- "Mto1"
    hit <- d$genome1_gene %in% grp$gene_id[sog_MtoM]
    orthopair$class[hit] <- "MtoM"
    
    hit <- match(d$genome1_gene, grp$gene_id)
    orthopair$SOG <- grp$grp[hit]
    
    return(orthopair)
}

.filterOrthopair <- function(orthopair, g2g_graph){
    split_gene <- .splitGene(orthopair = orthopair, g2g_graph = g2g_graph)
    hit <- orthopair$genome1_tx %in% split_gene$genome1
    orthopair$genome1_gene[hit] <- -orthopair$genome1_tx[hit]
    hit <- orthopair$genome2_tx %in% split_gene$genome2
    orthopair$genome2_gene[hit] <- -orthopair$genome2_tx[hit]
    orthopair <- .classifyOrthoPair(orthopair = orthopair)
    
    orthopair <- .untangleOrthoPair(orthopair = orthopair)
    orthopair <- .classifyOrthoPair(orthopair = orthopair)
    
    sog_best_mci <- tapply(orthopair$ci_q2s, orthopair$SOG, max)
    q <- quantile(sog_best_mci, c(0.25, 0.75))
    ci_q2s_threshold <- q[1] - 1.5 * diff(q)
    ci_q2s_filter <- orthopair$ci_q2s >= ci_q2s_threshold
    
    sog_best_mci <- tapply(orthopair$ci_s2q, orthopair$SOG, max)
    q <- quantile(sog_best_mci, c(0.25, 0.75))
    ci_s2q_threshold <- q[1] - 1.5 * diff(q)
    ci_s2q_filter <- orthopair$ci_s2q >= ci_s2q_threshold
    
    sog_best_pidt <- tapply(orthopair$pident, orthopair$SOG, max)
    q <- quantile(sog_best_pidt, c(0.25, 0.75))
    pidt_threshold <- q[1] - 1.5 * diff(q)
    pidt_filter <- orthopair$pident >= pidt_threshold
    
    filter <- ci_q2s_filter | ci_s2q_filter | pidt_filter
    orthopair <- orthopair[filter, ]
    
    non_anchor <- !orthopair$is_anchor_pair
    sog_best_pidt <- tapply(orthopair$pident[non_anchor], 
                            orthopair$SOG[non_anchor],
                            max)
    q <- quantile(sog_best_pidt, c(0.25, 0.75))
    random_mut_threshold <- q[1] - 1.5 * diff(q)
    orthopair <- orthopair[!non_anchor | orthopair$pident >= random_mut_threshold, ]
    return(orthopair)
}

.reformatOrthoPair <- function(orthopair, g2g_graph){
    q_tx_hit <- match(orthopair$genome1_tx, g2g_graph$genome1_df$tx_index)
    orthopair$genome1_tx <- g2g_graph$genome1_df$ID[q_tx_hit]
    orthopair$genome1_start <- g2g_graph$genome1_df$start[q_tx_hit]
    orthopair$genome1_end <- g2g_graph$genome1_df$end[q_tx_hit]
    orthopair$genome1_strand <- g2g_graph$genome1_df$strand[q_tx_hit]
    orthopair$genome1_strand[orthopair$genome1_strand == "1"] <- "+"
    orthopair$genome1_strand[orthopair$genome1_strand == "2"] <- "-"
    q_split <- which(orthopair$genome1_gene < 0)
    orthopair$genome1_gene <- g2g_graph$genome1_df$gene_id[q_tx_hit]
    orthopair$original_genome1_gene <- orthopair$genome1_gene
    orthopair$genome1_gene[q_split] <- paste0(orthopair$genome1_tx[q_split], ":split")
    
    s_tx_hit <- match(orthopair$genome2_tx, g2g_graph$genome2_df$tx_index)
    orthopair$genome2_tx <- g2g_graph$genome2_df$ID[s_tx_hit]
    orthopair$genome2_start <- g2g_graph$genome2_df$start[s_tx_hit]
    orthopair$genome2_end <- g2g_graph$genome2_df$end[s_tx_hit]
    orthopair$genome2_strand <- g2g_graph$genome2_df$strand[s_tx_hit]
    orthopair$genome2_strand[orthopair$genome2_strand == "1"] <- "+"
    orthopair$genome2_strand[orthopair$genome2_strand == "2"] <- "-"
    s_split <- which(orthopair$genome2_gene < 0)
    orthopair$genome2_gene <- g2g_graph$genome2_df$gene_id[s_tx_hit]
    orthopair$original_genome2_gene <- orthopair$genome2_gene
    orthopair$genome2_gene[s_split] <- paste0(orthopair$genome2_tx[s_split], ":split")
    
    # orthopair <- subset(orthopair,
    #                     select = c(query_gene:subject_chr,
    #                                pident:mutual_ci,
    #                                genome1_is_anchor,
    #                                genome2_is_anchor, 
    #                                is_anchor_pair,
    #                                genome1_synteny_block:original_genome2_gene))
    # orthopair$genome1_is_anchor <- as.numeric(orthopair$genome1_is_anchor)
    # orthopair$genome2_is_anchor <- as.numeric(orthopair$genome2_is_anchor)
    # orthopair$is_anchor_pair <- as.numeric(orthopair$is_anchor_pair)
    # orthopair$genome1_synteny_block <- factor(orthopair$genome1_synteny_block)
    # orthopair$genome2_synteny_block <- factor(orthopair$genome2_synteny_block)
    # orthopair$SOG <- factor(orthopair$SOG)
    # orthopair$genome1_synteny_block <- as.numeric(orthopair$genome1_synteny_block)
    # orthopair$genome2_synteny_block <- as.numeric(orthopair$genome2_synteny_block)
    # orthopair$SOG <- as.numeric(orthopair$SOG)
    return(orthopair)
}

# .best1to1 <- function(orthopair){
#     sog_1to1 <- orthopair$class == "1to1"
#     best_pair <- orthopair[sog_1to1, ]
#     # if(nrow(tmp) == 0){
#     #     return(NULL)
#     # }
#     # best_pair <- tapply(seq_along(tmp$subject_gene),
#     #                     tmp$subject_gene, function(i){
#     #                         tmp_i <- tmp[i, ]
#     #                         return(tmp_i[which.max(tmp_i$mutual_ci), ])
#     #                     })
#     # best_pair <- do.call("rbind", best_pair)
#     best_pair <- best_pair[order(best_pair$SOG), ]
#     return(best_pair)
# }

#' @importFrom GenomicRanges findOverlaps

.link2Anchor <- function(g2g_graph, anchor){
    genome1_link2anchor <- .link2Anchor_core(g2g_graph = g2g_graph,
                                             anchor = anchor,
                                             dataset = "genome1")
    genome2_link2anchor <- .link2Anchor_core(g2g_graph = g2g_graph,
                                               anchor = anchor,
                                               dataset = "genome2")
    
    out <- list(genome1_link2anchor = genome1_link2anchor,
                genome2_link2anchor = genome2_link2anchor)
    return(out)
}

#' @importFrom GenomicRanges precede follow findOverlaps resize 
#' @importFrom S4Vectors queryHits subjectHits

.linkGene2Genome <- function(object){
    gff_ls <- .getGFFlist(object = object)
    
    q_tx_i <- gff_ls$genome1_gff$type %in% c("transcript", "mRNA")
    genome1_df <- gff_ls$genome1_gff[q_tx_i, ]
    q_cds_i <- gff_ls$genome1_gff$type %in% "CDS"
    genome1_gff <- gff_ls$genome1_gff[q_cds_i, ]
    hit <- match(unlist((genome1_gff$Parent)), gff_ls$genome1_gff$ID)
    genome1_gff$tx_index <- gff_ls$genome1_gff$tx_index[hit]
    
    s_tx_i <- gff_ls$genome2_gff$type %in% c("transcript", "mRNA")
    genome2_df <- gff_ls$genome2_gff[s_tx_i, ]
    s_cds_i <- gff_ls$genome2_gff$type %in% "CDS"
    genome2_gff <- gff_ls$genome2_gff[s_cds_i, ]
    hit <- match(unlist((genome2_gff$Parent)), gff_ls$genome2_gff$ID)
    genome2_gff$tx_index <- gff_ls$genome2_gff$tx_index[hit]
    out <- list(genome1_df = genome1_df,
                genome2_df = genome2_df,
                genome1_gff = genome1_gff,
                genome2_gff = genome2_gff)
    return(out)
}

#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#'

.findSyntenicOrtho <- function(rbh, 
                               anchor,
                               g2g_graph,
                               t2a_graph){
    colnames(rbh)[1:2] <- c("genome1_tx", "genome2_tx")
    root_hit <- match(rbh$genome1_tx, g2g_graph$genome1_df$tx_index)
    leaf_hit <- match(rbh$genome2_tx, g2g_graph$genome2_df$tx_index)
    orthopair <- data.frame(root_tx = g2g_graph$genome1_df$tx_index[root_hit],
                            root = g2g_graph$genome1_df$gene_index[root_hit],
                            expected_leaf = g2g_graph$genome2_df$gene_index[leaf_hit],
                            leaf_tx = g2g_graph$genome2_df$tx_index[leaf_hit])
    orthopair <- left_join(orthopair, 
                           subset(t2a_graph$genome1_link2anchor, 
                                  subset = !is.na(root_anchor)), 
                           "root",
                           relationship = "many-to-many")
    orthopair <- left_join(orthopair,
                           subset(anchor,
                                  select = c(root_anchor, leaf_anchor)),
                           "root_anchor",
                           relationship = "many-to-many")
    orthopair <- left_join(orthopair, 
                           subset(t2a_graph$genome2_link2anchor, 
                                  subset = !is.na(leaf_anchor)), 
                           "leaf_anchor",
                           relationship = "many-to-many")
    orthopair <- subset(orthopair, 
                        subset = leaf == expected_leaf & 
                            !is.na(root_anchor) & 
                            !is.na(leaf_anchor))
    root_tx_hit <- match(orthopair$root_tx, g2g_graph$genome1_df$tx_index)
    leaf_tx_hit <- match(orthopair$leaf_tx, g2g_graph$genome2_df$tx_index)
    orthopair <- data.frame(genome1_gene = orthopair$root,
                            genome1_tx = orthopair$root_tx,
                            genome1_chr = g2g_graph$genome1_df$seqnames[root_tx_hit],
                            genome2_gene = orthopair$leaf,
                            genome2_tx = orthopair$leaf_tx,
                            genome2_chr = g2g_graph$genome2_df$seqnames[leaf_tx_hit],
                            orthopair)
    orthopair <- orthopair[order(orthopair$root), ]
    orthopair <- left_join(orthopair, rbh, c("genome1_tx", "genome2_tx"))
    orthopair <- unique(orthopair)
    orthopair$is_anchor_pair <- orthopair$pair_id %in% anchor$pair_id
    return(orthopair)
}

.findSyntenyBlocks <- function(orthopair){
    orthopair <- orthopair[order(orthopair$leaf), ]
    genome2_synteny_block <- .findSyntenyBlocksCore(orthopair = orthopair,
                                                    dataset = "genome2")
    
    root_order <- order(orthopair$root)
    genome2_synteny_block <- genome2_synteny_block[root_order]
    orthopair <- orthopair[root_order, ]
    genome1_synteny_block <- .findSyntenyBlocksCore(orthopair = orthopair,
                                                  dataset = "genome1")
    orthopair <- cbind(orthopair,
                       genome1_synteny_block = genome1_synteny_block, 
                       genome2_synteny_block = genome2_synteny_block)
    
    return(orthopair)
}

.findSyntenyBlocksCore <- function(orthopair, dataset){
    if(dataset == "genome1"){
        is_anchor <- orthopair$genome1_is_anchor
        chr <- orthopair$genome1_chr
        anchor_id <- orthopair$root_anchor
        
    } else {
        is_anchor <- orthopair$genome2_is_anchor
        chr <- orthopair$genome2_chr
        anchor_id <- orthopair$leaf_anchor
    }
    
    synteny_block_id <- rep(NA, nrow(orthopair))
    chr_list <- unique(chr)
    id_offset <- 0
    for(i_chr in chr_list){
        i <- chr == i_chr
        anchor_index <- which(is_anchor[i])
        anchor_gap_pos <- which(diff(anchor_index) > 1)
        synteny_block <- data.frame(start = c(head(anchor_index, 1), 
                                              anchor_index[anchor_gap_pos + 1]),
                                    end = c(anchor_index[anchor_gap_pos],
                                            tail(anchor_index, 1)))
        synteny_block$width <- synteny_block$end - synteny_block$start + 1
        synteny_block$prev_start <- c(1, head(synteny_block$end, -1) + 1)
        synteny_block$prev_end <- synteny_block$start - 1
        synteny_block$gap_width <- synteny_block$prev_end - synteny_block$prev_start + 1
        
        block_id <- sapply(seq_along(synteny_block$start), function(j) {
            rep(j, synteny_block$width[j])
        })
        synteny_block_id[i][anchor_index] <- unlist(block_id) + id_offset
        
        gap_id <- lapply(seq_along(synteny_block$start), function(j) {
            rep(j, synteny_block$gap_width[j])
        })
        gap_id <- unlist(gap_id)
        if(length(gap_id) > 0){
            non_anchor_index <- which(!is_anchor[i])
            hit <- match(anchor_id[i][-anchor_index], 
                         anchor_id[i][anchor_index])
            gap_id <- c(gap_id, rep(max(gap_id) + 1, 
                                    length(non_anchor_index) - length(gap_id)))
            gap_block <- data.frame(gene_index = non_anchor_index, 
                                    anchor_index = anchor_index[hit],
                                    gap_id = gap_id)
            gap_block$block_id <- synteny_block_id[i][gap_block$anchor_index]
            
            gap_hit2anchor <- tapply(gap_block$block_id, gap_block$gap_id, unique)
            n_gap_hit2anchor <- sapply(gap_hit2anchor, length)
            briging_gap <- as.integer(names(gap_hit2anchor)[n_gap_hit2anchor > 1])
            briging_gap_block <- gap_block[gap_block$gap_id %in% briging_gap, ]
            briging_gap_block_min <- tapply(briging_gap_block$block_id,
                                            briging_gap_block$gap_id, 
                                            min)
            hit <- match(gap_block$gap_id, names(briging_gap_block_min))
            gap_block$block_id[!is.na(hit)] <- briging_gap_block_min[na.omit(hit)] + 0.5
            synteny_block_id[i][-anchor_index] <- gap_block$block_id
            
            gap_pos <- which(diff(synteny_block_id[i]) >= 1)
            block_len <- diff(c(0, gap_pos, length(synteny_block_id[i])))
            block_id <- unlist(sapply(seq_along(block_len),
                                      function(i) {
                                          rep(i, block_len[i])
                                      }))
            synteny_block_id[i] <- block_id + id_offset
        }
        id_offset <- max(synteny_block_id, na.rm = TRUE)
    }
    return(synteny_block_id)
}

#' @importFrom igraph graph_from_data_frame V components

.link2Anchor_core <- function(g2g_graph, anchor, dataset){
    if(dataset == "genome1"){
        gene_df <- g2g_graph$genome1_df
        anchor_id <- anchor$root_anchor
        
    } else {
        gene_df <- g2g_graph$genome2_df
        anchor_id <- anchor$leaf_anchor
    }
    
    anchor_df <- subset(gene_df,
                        subset = gene_index %in% anchor_id,
                        select = c(seqnames, gene_index))
    anchor_df <- unique(anchor_df)
    names(anchor_df) <- c("chr", "gene")
    anchor_df$anchor <- anchor_df$gene
    non_anchor_df <- subset(gene_df,
                            subset = !gene_index %in% anchor_id,
                            select = c(seqnames, gene_index))
    non_anchor_df <- unique(non_anchor_df)
    
    link2anchor_upper <- lapply(non_anchor_df$gene_index, "-", anchor_df$gene)
    link2anchor_upper <- lapply(link2anchor_upper, ">", 0)
    link2anchor_upper <- lapply(link2anchor_upper, which)
    link2anchor_upper_length <- sapply(link2anchor_upper, length)
    link2anchor_upper[link2anchor_upper_length == 0] <- Inf
    link2anchor_upper <- sapply(link2anchor_upper, max)
    link2anchor_upper[is.infinite(link2anchor_upper)] <- NA
    
    link2anchor_lower <- lapply(non_anchor_df$gene_index, "-", anchor_df$gene)
    link2anchor_lower <- lapply(link2anchor_lower, "<", 0)
    link2anchor_lower <- lapply(link2anchor_lower, which)
    link2anchor_lower_length <- sapply(link2anchor_lower, length)
    link2anchor_lower[link2anchor_lower_length == 0] <- Inf
    link2anchor_lower <- sapply(link2anchor_lower, min)
    link2anchor_lower[is.infinite(link2anchor_lower)] <- NA
    
    non_anchor_df <- data.frame(gene = c(non_anchor_df$gene_index,
                                         non_anchor_df$gene_index),
                                chr = c(non_anchor_df$seqnames,
                                        non_anchor_df$seqnames),
                                anchor = c(anchor_df$gene[link2anchor_upper],
                                           anchor_df$gene[link2anchor_lower]),
                                anchor_chr = c(anchor_df$chr[link2anchor_upper],
                                               anchor_df$chr[link2anchor_lower]))
    non_anchor_df <- subset(non_anchor_df, 
                            subset = chr == anchor_chr, 
                            select = -c(chr, anchor_chr))
    anchor_df$is_anchor <- TRUE
    non_anchor_df$is_anchor <- FALSE
    
    link2anchor <- rbind(subset(anchor_df, select = c(gene, anchor, is_anchor)),
                         non_anchor_df)
    link2anchor <- unique(link2anchor[order(link2anchor$gene), ])
    
    if(dataset == "genome1"){
        colnames(link2anchor) <- c("root", "root_anchor", "genome1_is_anchor")
        
    } else {
        colnames(link2anchor) <- c("leaf", "leaf_anchor", "genome2_is_anchor")
    }
    
    return(link2anchor)
}

#' @importFrom dplyr left_join

.splitGene <- function(orthopair, g2g_graph){
    split_1toM <- .split1toM(orthopair, gff = g2g_graph$genome1_df)
    split_Mto1 <- .splitMto1(orthopair, gff = g2g_graph$genome2_df)
    split_MtoM <- .splitMtoM(orthopair, g2g_graph = g2g_graph)
    
    out <- list(genome1 = c(split_1toM, split_MtoM$genome1),
                genome2 = c(split_Mto1, split_MtoM$genome2))
    return(out)
}

.split1toM <- function(orthopair, gff){
    gff <- GRanges(seqnames = gff$seqnames, 
                   ranges = IRanges(start = gff$start,
                                    end = gff$end), 
                   gene_id = gff$gene_id,
                   tx_index = gff$tx_index, 
                   gene_index = gff$gene_index)
    genome1_ol <- findOverlaps(gff, gff)
    genome1_ol <- as.data.frame(genome1_ol)
    genome1_ol$genome1_tx <- gff$tx_index[genome1_ol$queryHits]
    genome1_ol$genome1_ol_tx <- gff$tx_index[genome1_ol$subjectHits]
    valid <- gff$gene_index[genome1_ol$queryHits] == gff$gene_index[genome1_ol$subjectHits]
    genome1_ol <- subset(genome1_ol, subset = genome1_tx != genome1_ol_tx & valid)
    genome1_ol <- unique(subset(genome1_ol, select = genome1_tx:genome1_ol_tx))
    
    sog_1toM <- orthopair$class == "1toM"
    if(sum(sog_1toM) == 0){
        return(NULL)
    }
    orthopair_subset <- subset(orthopair, 
                               subset = sog_1toM,
                               select = c(genome1_gene:genome1_tx,
                                          genome2_gene:genome2_tx, 
                                          SOG))
    orthopair_subset$sog_tx <- paste(orthopair_subset$SOG, 
                                     orthopair_subset$genome2_tx,
                                     sep = "_")
    
    sog_tx_par_q_tx <- tapply(orthopair_subset$sog_tx, 
                              orthopair_subset$genome1_tx,
                              unique)
    n_sog_tx_par_q_tx <- sapply(sog_tx_par_q_tx, length)
    splitable <- n_sog_tx_par_q_tx == 1
    splitable <- sog_tx_par_q_tx[splitable]
    splitable <- data.frame(genome1_tx = as.numeric(names(splitable)), 
                            genome2_tx = as.numeric(sub("[0-9]+_", "", unlist(splitable))),
                            SOG = as.numeric(sub("_.+", "", unlist(splitable))))
    genome1_ol <- subset(genome1_ol, 
                         subset = genome1_ol_tx %in% orthopair_subset$genome1_tx)
    splitable <- left_join(splitable, genome1_ol, "genome1_tx")
    hit <- match(splitable$genome1_tx, gff$tx_index)
    splitable$genome1_gene <-  gff$gene_index[hit]
    hit <- match(splitable$genome1_ol_tx, gff$tx_index)
    splitable$genome1_ol_gene <-  gff$gene_index[hit]
    splitable <- subset(splitable,
                        subset = genome1_gene != genome1_ol_gene | is.na(genome1_ol_gene))
    out <- splitable$genome1_tx
    return(out)
}

.splitMto1 <- function(orthopair, gff){
    gff <- GRanges(seqnames = gff$seqnames, 
                   ranges = IRanges(start = gff$start,
                                    end = gff$end), 
                   gene_id = gff$gene_id,
                   tx_index = gff$tx_index, 
                   gene_index = gff$gene_index)
    genome2_ol <- findOverlaps(gff, gff)
    genome2_ol <- as.data.frame(genome2_ol)
    genome2_ol$genome2_tx <- gff$tx_index[genome2_ol$subjectHits]
    genome2_ol$genome2_ol_tx <- gff$tx_index[genome2_ol$subjectHits]
    valid <- gff$gene_index[genome2_ol$queryHits] == gff$gene_index[genome2_ol$subjectHits]
    genome2_ol <- subset(genome2_ol, subset = genome2_tx != genome2_ol_tx & valid)
    genome2_ol <- unique(subset(genome2_ol, select = genome2_tx:genome2_ol_tx))
    
    sog_Mto1 <- orthopair$class == "Mto1"
    if(sum(sog_Mto1) == 0){
        return(NULL)
    }
    orthopair_subset <- subset(orthopair, 
                               subset = sog_Mto1,
                               select = c(genome1_gene:genome1_tx,
                                          genome2_gene:genome2_tx, 
                                          SOG))
    orthopair_subset$sog_tx <- paste(orthopair_subset$SOG, 
                                     orthopair_subset$genome1_tx,
                                     sep = "_")
    
    sog_tx_par_s_tx <- tapply(orthopair_subset$sog_tx, 
                              orthopair_subset$genome2_tx,
                              unique)
    n_sog_tx_par_s_tx <- sapply(sog_tx_par_s_tx, length)
    splitable <- n_sog_tx_par_s_tx == 1
    splitable <- sog_tx_par_s_tx[splitable]
    splitable <- data.frame(genome2_tx = as.numeric(names(splitable)), 
                            genome1_tx = as.numeric(sub("[0-9]+_", "", unlist(splitable))),
                            SOG = as.numeric(sub("_.+", "", unlist(splitable))))
    genome2_ol <- subset(genome2_ol, 
                         subset = genome2_ol_tx %in% orthopair_subset$genome2_tx)
    splitable <- left_join(splitable, genome2_ol, "genome2_tx")
    hit <- match(splitable$genome2_tx, gff$tx_index)
    splitable$genome2_gene <-  gff$gene_index[hit]
    hit <- match(splitable$genome2_ol_tx, gff$tx_index)
    splitable$genome2_ol_gene <-  gff$gene_index[hit]
    splitable <- subset(splitable,
                        subset = genome2_gene != genome2_ol_gene | is.na(genome2_ol_gene))
    out <- splitable$genome2_tx
    return(out)
}

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
# MtoM だけsplitの定義が違う？
# 1toM と Mto1は、tx-variants のsplit
# MtoMは、MtoMを1toMやMto1splitしてる？？

.splitMtoM <- function(orthopair, g2g_graph){
    sog_MtoM <- orthopair$class == "MtoM"
    if(sum(sog_MtoM) == 0){
        return(NULL)
    }
    orthopair_subset <- subset(orthopair, 
                               subset = sog_MtoM,
                               select = c(genome1_gene:genome1_tx,
                                          genome2_gene:genome2_tx, 
                                          SOG))
    orthopair_subset$class <- "1toM"
    orthopair_subset$SOG <- as.numeric(factor(orthopair_subset$genome1_gene))
    n_member <- table(orthopair_subset$SOG)
    split_1toM <- .split1toM(subset(orthopair_subset,
                                    subset = SOG %in% names(n_member[n_member > 1])), 
                             gff = g2g_graph$genome1_gff)
    orthopair_subset$class <- "Mto1"
    orthopair_subset$SOG <- as.numeric(factor(orthopair_subset$genome2_gene))
    n_member <- table(orthopair_subset$SOG)
    split_Mto1 <- .splitMto1(subset(orthopair_subset,
                                    subset = SOG %in% names(n_member[n_member > 1])),
                             gff = g2g_graph$genome2_gff)
    out <- list(genome1 = split_1toM, genome2 = split_Mto1)
    return(out)
}

.checkHighCopyGenes <- function(anchor, rbh){
    subset_rbh <- subset(rbh, select = c(genome1_gene, genome2_gene))
    subset_rbh <- unique(subset_rbh)
    q_rbh <- rbh[rbh$genome1_gene %in% anchor$genome1_gene, ]
    n_q_rbh <- table(q_rbh$genome1_gene)
    q <- quantile(n_q_rbh, c(0.25, 0.75), na.rm = TRUE)
    whisker <- 1.5 * diff(q)
    q_threshold <- q[2] + whisker
    q_omit <- names(n_q_rbh)[n_q_rbh > q_threshold]
    
    s_rbh <- rbh[rbh$genome2_gene %in% anchor$genome2_gene, ]
    n_s_rbh <- table(s_rbh$genome2_gene)
    q <- quantile(n_s_rbh, c(0.25, 0.75), na.rm = TRUE)
    whisker <- 1.5 * diff(q)
    q_threshold <- q[2] + whisker
    s_omit <- names(n_s_rbh)[n_s_rbh > q_threshold]
    
    anchor <- subset(anchor, subset = !genome1_gene %in% q_omit | !genome2_gene %in% s_omit)
    return(anchor)
}

.untangleOrthoPair <- function(orthopair){
    orthopair$index <- seq_along(orthopair$genome1_gene)
    orthopair <- orthopair[order(orthopair$ci_q2s, decreasing = TRUE), ]
    q_best <- tapply(orthopair$index,
                     orthopair$genome1_gene, 
                     "[", 1)
    orthopair <- orthopair[order(orthopair$ci_s2q, decreasing = TRUE), ]
    s_best <- tapply(orthopair$index,
                     orthopair$genome2_gene, 
                     "[", 1)
    best <- unique(c(q_best, s_best))
    best_orthopair <- orthopair[orthopair$index %in% best, ]
    local_anchor <- intersect(q_best, s_best)
    rest_pair <- best[!best %in% local_anchor]
    anchor_orthopair <- orthopair[orthopair$index %in% local_anchor, ]
    rest_orthopair <- orthopair[orthopair$index %in% rest_pair, ]
    
    query_to_anchor <- rest_orthopair$genome1_gene %in% anchor_orthopair$genome1_gene
    subject_to_anchor <- rest_orthopair$genome2_gene %in% anchor_orthopair$genome2_gene
    anchored_orthopair <- rest_orthopair[query_to_anchor | subject_to_anchor, ]
    rest_pair <- rest_pair[!rest_pair %in% anchored_orthopair$index]
    rest_orthopair <- rest_orthopair[rest_orthopair$index %in% rest_pair, ]
    
    query_to_anchored <- rest_orthopair$genome1_gene %in% anchored_orthopair$genome1_gene
    subject_to_anchored <- rest_orthopair$genome2_gene %in% anchored_orthopair$genome2_gene
    
    nonanchored_orthopair <- rest_orthopair[!(query_to_anchored | subject_to_anchored), ]
    out <- rbind(best_orthopair, anchored_orthopair, nonanchored_orthopair)
    out <- unique(out[order(out$genome1_tx), ])
    return(out)
}

.getRBBH <- function(rbh){
    rbh$index <- seq_len(nrow(rbh))
    rbh <- rbh[order(rbh$ci_q2s, decreasing = TRUE), ]
    q_best <- rbh$index[!duplicated(rbh$genome1_gene)]
    rbh <- rbh[order(rbh$ci_s2q, decreasing = TRUE), ]
    s_best <- rbh$index[!duplicated(rbh$genome2_gene)]
    rbbh <- intersect(q_best, s_best)
    rbbh <- rbh[rbh$index %in% rbbh, ]
    return(rbbh)
}

.getGFFlist <- function(object = NULL, gff_fn = NULL){
    # Import GFF files for query and subject genomes
    if(is.null(gff_fn)){
        files <- h5read(object$h5, "files")
        query_gff <- readRDS(files$query_gff_df)
        subject_gff <- readRDS(files$subject_gff_df)
        
        # Order the GFF data
        query_gff <- query_gff[order(query_gff$seqnames, query_gff$start), ]
        subject_gff <- subject_gff[order(subject_gff$seqnames, subject_gff$start), ]
        
        # Return the ordered GFF data as a list
        out <- list(genome1_gff = query_gff, genome2_gff = subject_gff)
        
    } else {
        if(grepl("\\.rds$", gff_fn)){
            out <- readRDS(gff_fn)
        } else {
            out <- import.gff(gff_fn)
        }
        out <- .orderGFF(gff = out)
        out <- .mRNA2transcript(gff = out)
    }
    
    return(out)
}

.mRNA2transcript <- function(gff){
    type <- as.character(gff$type)
    type[type == "mRNA"] <- "transcript"
    gff$type <- factor(type)
    return(gff)
}

## Map legacy BLAST-style RBH column names to genome1_/genome2_ (idempotent).

.orderGFF <- function(gff){
    # Order GFF data by chromosome and start position
    gff_order <- order(as.character(seqnames(gff)), start(gff))
    gff <- gff[gff_order]
    return(gff)
}
