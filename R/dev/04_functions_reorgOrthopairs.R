#' Reorganise orthopair TSV results and create graph/table outputs
#'
#' This is a high-level wrapper that:
#' \itemize{
#'   \item reorganises per-pair OrthoPaiR TSV results in \code{working_dir/orthopair},
#'   \item optionally applies split-gene renaming based on \code{*:split} gene IDs,
#'   \item creates a genome-level orthology graph (GraphML),
#'   \item and summarises the graph into tabular outputs in the same directory.
#' }
#'
#' It is intended for use after pairwise OrthoPaiR runs have completed and
#' per-pair TSV files have been written by \code{orthopair()} into
#' \code{working_dir/orthopair}.
#'
#' Compared to \code{reorgOrthopairs()} in \code{R/08_functions_postsynog.R},
#' this implementation:
#' \itemize{
#'   \item does not depend on HDF5 results; it reads TSV files directly,
#'   \item derives genome IDs from TSV file names (\code{<genomeA>_<genomeB>.tsv}),
#'   \item and uses data.table throughout for speed.
#' }
#'
#' @param working_dir Base working directory. Must contain:
#'   \itemize{
#'     \item \code{orthopair/} – per-pair TSV files named \code{<genomeA>_<genomeB>.tsv}
#'       produced by \code{orthopair()},
#'     \item \code{input/} – per-genome folders with GFF/CDS (not strictly required
#'       for graph construction but expected by the overall pipeline).
#'   }
#' @param rename Logical; whether to apply split-gene renaming when reorganising
#'   orthopair information. Split genes are identified by \code{\"split\"}
#'   appearing in \code{query_gene} / \code{subject_gene}.
#' @param n_threads Integer; number of cores to use for internal parallel steps
#'   (Linux/macOS only; on Windows this falls back to serial execution).
#' @param overwrite Logical; reserved for future use (currently ignored).
#' @param verbose Logical; whether to print progress messages.
#'
#' @return The path to the output directory,
#'   \code{file.path(working_dir, \"reorg_out\")}.
#'
#' @importFrom data.table data.table rbindlist fwrite setcolorder dcast
#' @importFrom parallel mclapply
#' @importFrom igraph write_graph graph_from_data_frame components vertex_attr
#'   edge_attr ends E
#' @importFrom dplyr full_join
#' @export
reorgOrthopairs <- function(working_dir,
                            rename = TRUE,
                            n_threads = 1L,
                            overwrite = FALSE,
                            verbose = TRUE) {
    stopifnot(length(working_dir) == 1L)
    working_dir <- normalizePath(working_dir, mustWork = TRUE)
    
    orthopair_dir <- file.path(working_dir, "orthopair")
    if (!dir.exists(orthopair_dir)) {
        stop("orthopair directory not found under working_dir: ", orthopair_dir)
    }
    
    orthopair_files <- list.files(orthopair_dir, pattern = "\\.tsv$", full.names = TRUE)
    if (length(orthopair_files) == 0L) {
        stop("No orthopair TSV files found in ", orthopair_dir,
             " (expected files named <genomeA>_<genomeB>.tsv).")
    }
    
    out_dir <- file.path(working_dir, "reorg_out")
    dir.create(path = out_dir, showWarnings = FALSE, recursive = TRUE)
    
    if (verbose) {
        message("[reorg] Found ", length(orthopair_files), " orthopair TSV files")
    }
    
    reorg_list <- .getRenameList(orthopair_files = orthopair_files,
                                 rename = rename,
                                 n_threads = n_threads,
                                 verbose = verbose)
    
    orthopair_list <- .reorgOrthoPair(orthopair_files = orthopair_files,
                                      reorg_list = reorg_list,
                                      n_threads = n_threads,
                                      rename = rename,
                                      verbose = verbose)
    
    graph <- .createGraph(orthopair_list = orthopair_list,
                          reorg_list = reorg_list)
    graph_fn <- file.path(out_dir, "orthopair.graphml")
    write_graph(graph = graph, file = graph_fn, format = "graphml")
    
    .graph2df(graph = graph, out_dir = out_dir)
    
    if (verbose) {
        message("[reorg] Wrote graph and summary tables to: ", out_dir)
    }
    
    out_dir
}


## Internal helpers for TSV-based reorganisation -----------------------------

#' Build split-gene rename list and per-genome gene sets from orthopair TSVs
.getRenameList <- function(orthopair_files,
                           rename = FALSE,
                           n_threads = 1L,
                           verbose = TRUE) {
    idx <- seq_along(orthopair_files)
    use_parallel <- (n_threads > 1L && .Platform$OS.type != "windows")
    
    worker_fun <- function(i) {
        fn <- orthopair_files[i]
        pair_id <- sub("\\.tsv$", "", basename(fn))
        parts <- strsplit(pair_id, "_", fixed = TRUE)[[1L]]
        if (length(parts) != 2L) {
            stop("Invalid orthopair TSV name (expected <genomeA>_<genomeB>.tsv): ", fn)
        }
        genomes_i <- parts
        
        dt <- data.table::fread(fn, sep = "\t", header = TRUE)
        
        # Gene sets per genome (after any renaming)
        dt_query <- unique(dt[, .(original = original_query_gene, gene = query_gene)])
        dt_query[, genome := genomes_i[1]]

        dt_subject <- unique(dt[, .(original = original_subject_gene, gene = subject_gene)])
        dt_subject[, genome := genomes_i[2]]

        genes_dt <- rbind(dt_query, dt_subject)
        
        split_dt <- NULL
        if (rename) {
            # Identify split genes by pattern "split" in gene IDs
            query_split <- dt[grepl("split", query_gene),
                              .(original = original_query_gene, 
                                gene = query_gene, 
                                tx = query_tx)]
            subject_split <- dt[grepl("split", subject_gene),
                                .(original = original_subject_gene, 
                                  gene = subject_gene, 
                                  tx = subject_tx)]
            
            if (nrow(query_split) > 0L) {
                split_dt <- rbind(
                    split_dt,
                    data.table::data.table(
                        genome = genomes_i[1],
                        original = query_split$original, 
                        gene = query_split$gene,
                        tx = query_split$tx
                    )
                )
            }
            if (nrow(subject_split) > 0L) {
                split_dt <- rbind(
                    split_dt,
                    data.table::data.table(
                        genome = genomes_i[2],
                        original = subject_split$original, 
                        gene = subject_split$gene,
                        tx = subject_split$tx
                    )
                )
            }
        }
        
        list(genes = genes_dt, split = split_dt)
    }
    
    if (use_parallel) {
        res <- parallel::mclapply(
            X = idx,
            FUN = worker_fun,
            mc.cores = n_threads
        )
    } else {
        res <- lapply(idx, worker_fun)
    }
    
    # Combine per-genome gene sets
    genes_list <- lapply(res, function(x) x$genes)
    genes_list <- rbindlist(genes_list)
    genes_list <- unique(genes_list)
    
    # Combine split_list if rename = TRUE
    split_list <- NULL
    if (rename) {
        split_pieces <- lapply(res, function(x) x$split)
        split_pieces <- Filter(Negate(is.null), split_pieces)
        if (length(split_pieces)) {
            split_list <- unique(data.table::rbindlist(split_pieces))
            n_split <- table(split_list$original)
            valid_split <- names(n_split)[n_split > 1]
            split_list <- subset(split_list,
                                 subset = original %in% valid_split)
        }
    }
    
    list(
        genes_list = genes_list,
        split_list = split_list
    )
}


#' Reorganise orthopair TSVs into per-pair tables suitable for graph building
.reorgOrthoPair <- function(orthopair_files,
                            reorg_list,
                            n_threads = 1L,
                            rename = FALSE,
                            verbose = TRUE) {
    idx <- seq_along(orthopair_files)
    use_parallel <- (n_threads > 1L && .Platform$OS.type != "windows")
    
    worker_fun <- function(i) {
        fn <- orthopair_files[i]
        pair_id <- sub("\\.tsv$", "", basename(fn))
        parts <- strsplit(pair_id, "_", fixed = TRUE)[[1L]]
        if (length(parts) != 2L) {
            stop("Invalid orthopair TSV name (expected <genomeA>_<genomeB>.tsv): ", fn)
        }
        genomes_i <- parts
        
        if (verbose) {
            message("[reorg] Reading orthopair TSV: ", basename(fn))
        }
        opr <- data.table::fread(fn, sep = "\t", header = TRUE)
        
        # Ensure required columns exist
        required_cols <- c("query_gene", "subject_gene",
                           "query_tx", "subject_tx",
                           "mutual_ci", "class")
        missing_cols <- setdiff(required_cols, names(opr))
        if (length(missing_cols)) {
            stop("Missing required columns in ", fn, ": ",
                 paste(missing_cols, collapse = ", "))
        }
        
        if (!rename) {
            # Use original_* genes if available (pre-split), otherwise current genes
            if (all(c("original_query_gene", "original_subject_gene") %in% names(opr))) {
                target_col <- c("original_query_gene", "original_subject_gene",
                                "query_tx", "subject_tx",
                                "mutual_ci", "class")
                opr <- as.data.frame(opr)[, target_col]
                names(opr) <- c("query_gene", "subject_gene",
                                "query_tx", "subject_tx",
                                "mutual_ci", "class")
            } else {
                opr <- as.data.frame(opr)[, required_cols]
            }
        } else {
            # Use current query_gene/subject_gene and apply split renaming
            opr <- as.data.frame(opr)[, required_cols]
            meta <- list(genomes = genomes_i)
            opr <- .renameOrthoPair(opr = opr,
                                    meta = meta,
                                    reorg_list = reorg_list)
        }
        
        opr$query_id <- paste(genomes_i[1], opr$query_gene, sep = ":")
        opr$subject_id <- paste(genomes_i[2], opr$subject_gene, sep = ":")
        opr
    }
    
    res_list <- if (use_parallel) {
        parallel::mclapply(
            X = idx,
            FUN = worker_fun,
            mc.cores = n_threads
        )
    } else {
        lapply(idx, worker_fun)
    }
    
    res_list
}


.renameOrthoPair <- function(opr, meta, reorg_list){
    split_list <- reorg_list$split_list
    if (is.null(split_list) || !nrow(split_list)) return(opr)
    
    genome_1 <- split_list$genome == meta$genomes[1]
    if (any(genome_1)) {
        hit <- match(opr$query_tx, split_list$tx[genome_1])
        not_na <- !is.na(hit)
        if (any(not_na)) {
            opr$query_gene[not_na] <- split_list$gene[genome_1][hit[not_na]]
        }
    }
    genome_2 <- split_list$genome == meta$genomes[2]
    if (any(genome_2)) {
        hit <- match(opr$subject_tx, split_list$tx[genome_2])
        not_na <- !is.na(hit)
        if (any(not_na)) {
            opr$subject_gene[not_na] <- split_list$gene[genome_2][hit[not_na]]
        }
    }
    opr
}


## Graph creation & summarisation (copied from 08_functions_postsynog.R) ----

#' @import data.table
.createGraph <- function(orthopair_list, reorg_list){
    edges_list <- data.table::rbindlist(orthopair_list)
    vertex_list <- reorg_list$genes_list
    vertex_list$id <- paste(vertex_list$genome, vertex_list$gene, sep = ":")
    
    ## Build edge data.frame
    edges_df <- data.frame(from = edges_list$query_id,
                           to = edges_list$subject_id,
                           mutual_ci = edges_list$mutual_ci,
                           class = edges_list$class,
                           stringsAsFactors = FALSE)
    
    ## Single genome lookup vector
    gene_to_genome <- setNames(vertex_list$genome, vertex_list$id)
    edges_df$genome1 <- gene_to_genome[edges_df$from]
    edges_df$genome2 <- gene_to_genome[edges_df$to]
    
    ## Vertex table
    vertex_df <- data.frame(name = vertex_list$id,
                            genome = vertex_list$genome,
                            original = vertex_list$original,
                            stringsAsFactors = FALSE)
    
    ## Build graph in a single call
    graph_out <- igraph::graph_from_data_frame(d = edges_df,
                                               vertices = vertex_df,
                                               directed = FALSE)
    graph_out
}


#' @importFrom igraph edge_attr ends E vertex_attr
#' @importFrom dplyr full_join
.graph2df <- function(graph, out_dir) {
    comp <- igraph::components(graph)
    va <- igraph::vertex_attr(graph)
    ea <- igraph::edge_attr(graph)
    em <- igraph::ends(graph, igraph::E(graph), names = FALSE)
    rm(graph); gc(); gc()
    
    orig_attr <- if ("original" %in% names(va)) va$original else va$name
    dtv <- data.table::data.table(membership = comp$membership,
                                  genome = va$genome,
                                  gene = va$name,
                                  original = orig_attr)
    data.table::fwrite(dtv,
                       file.path(out_dir, "orthopair_list_long.tsv"),
                       sep = "\t",
                       quote = FALSE,
                       na = "")
    
    genome_levels <- unique(dtv$genome)
    data.table::fwrite(
        data.table::setcolorder(
            data.table::dcast(
                dtv[, .(gene = paste(gene, collapse=",")),
                    by = .(membership, genome)],
                membership ~ genome,
                value.var = "gene",
                fill = NA_character_
            ),
            c("membership", genome_levels)
        ),
        file.path(out_dir, "orthopair_list.tsv"),
        sep = "\t",
        quote = FALSE,
        na = ""
    )
    
    ## -------- Pairwise ortholog lists (un-collapsed) per genome pair --------
    pairwise_dir <- file.path(out_dir, "pairwise")
    dir.create(pairwise_dir, showWarnings = FALSE, recursive = TRUE)

    if (nrow(dtv) > 0L) {
        mem_ids <- unique(dtv$membership)
        pairwise_list <- vector("list", length(mem_ids))
        pw_idx <- 1L

        for (m in mem_ids) {
            sub <- dtv[membership == m]
            gvec <- unique(sub$genome)
            if (length(gvec) < 2L) next

            if (length(gvec) == 2L) {
                gpairs <- matrix(sort(gvec), ncol = 2L, byrow = TRUE)
            } else {
                gpairs <- t(utils::combn(sort(gvec), 2L))
            }

            for (k in seq_len(nrow(gpairs))) {
                gA <- gpairs[k, 1L]
                gB <- gpairs[k, 2L]
                genesA <- sub[genome == gA, .(gene, original)]
                genesB <- sub[genome == gB, .(gene, original)]
                if (!nrow(genesA) || !nrow(genesB)) next

                idxA <- seq_len(nrow(genesA))
                idxB <- seq_len(nrow(genesB))
                grid <- data.table::CJ(i = idxA, j = idxB, sorted = FALSE)

                pw_dt <- data.table::data.table(
                    membership     = m,
                    genome1        = gA,
                    gene1          = genesA$gene[grid$i],
                    original_gene1 = genesA$original[grid$i],
                    genome2        = gB,
                    gene2          = genesB$gene[grid$j],
                    original_gene2 = genesB$original[grid$j]
                )

                pairwise_list[[pw_idx]] <- pw_dt
                pw_idx <- pw_idx + 1L
            }
        }

        pairwise_list <- Filter(Negate(is.null), pairwise_list)
        if (length(pairwise_list)) {
            all_pairs <- data.table::rbindlist(pairwise_list, use.names = TRUE)
            gpairs2 <- unique(all_pairs[, .(genome1, genome2)])
            for (i in seq_len(nrow(gpairs2))) {
                gA <- gpairs2$genome1[i]
                gB <- gpairs2$genome2[i]
                sub_dt <- all_pairs[genome1 == gA & genome2 == gB]
                if (!nrow(sub_dt)) next
                out_fn <- file.path(pairwise_dir, paste0(gA, "_", gB, ".tsv"))
                data.table::fwrite(
                    sub_dt,
                    out_fn,
                    sep = "\t",
                    quote = FALSE,
                    na = ""
                )
            }
        }
    }

    eid <- dtv$membership[em[, 1]]
    g1 <- va$genome[em[, 1]]
    g2 <- va$genome[em[, 2]]
    rm(comp, va, em); gc(); gc()
    
    v_sum <- dtv[, .(nV = .N, nGenome = data.table::uniqueN(genome)), by = membership]
    
    dte <- data.table::data.table(membership = eid, mutual_ci = ea$mutual_ci)
    
    e_sum <- dte[, .(nE = .N,
                     max_mci = max(mutual_ci, na.rm = TRUE),
                     q1 = stats::quantile(mutual_ci, 0.25, na.rm = TRUE),
                     median_mci = stats::median(mutual_ci, na.rm = TRUE),
                     q3 = stats::quantile(mutual_ci, 0.75, na.rm = TRUE),
                     min_mci = min(mutual_ci, na.rm = TRUE),
                     mean_mci = mean(mutual_ci, na.rm = TRUE),
                     sd_mci = stats::sd(mutual_ci, na.rm = TRUE)),
                 by = membership]
    
    out <- dplyr::full_join(v_sum, e_sum, "membership")
    data.table::fwrite(out,
                       file.path(out_dir, "orthogroup_stats.tsv"),
                       sep = "\t",
                       quote = FALSE,
                       na = "")
    
    ## genomepair_edge_stats も、間接ペアを含むペア数を反映させる。
    ## all_pairs が存在すれば、そこから genome1/genome2 ごとの行数で nE を計算し、
    ## mutual_ci は NA（統計値は direct edges からのものを保持）とする。
    if (exists("all_pairs") && is.data.frame(all_pairs) && nrow(all_pairs) > 0L) {
        dt_pair <- data.table::data.table(
            genome1   = all_pairs$genome1,
            genome2   = all_pairs$genome2,
            mutual_ci = NA_real_
        )
    } else {
        genome1 <- pmin(g1, g2)
        genome2 <- pmax(g1, g2)
        dt_pair <- data.table::data.table(genome1 = genome1,
                                          genome2 = genome2,
                                          mutual_ci = ea$mutual_ci)
    }
    pair_sum <- dt_pair[, .(nE = .N,
                            max_mci = max(mutual_ci, na.rm = TRUE),
                            q1 = stats::quantile(mutual_ci, 0.75, na.rm = TRUE),
                            median_mci= stats::median(mutual_ci, na.rm = TRUE),
                            q3 = stats::quantile(mutual_ci, 0.25, na.rm = TRUE),
                            min_mci = min(mutual_ci, na.rm = TRUE),
                            mean_mci = mean(mutual_ci, na.rm = TRUE),
                            sd_mci = stats::sd(mutual_ci, na.rm = TRUE)),
                        by = .(genome1, genome2)][order(-nE)]
    
    data.table::fwrite(
        pair_sum,
        file.path(out_dir, "genomepair_edge_stats.tsv"),
        sep = "\t",
        quote = FALSE,
        na = ""
    )

    rm(dtv, dt_pair, pair_sum, g1, g2); gc(); gc()
}

