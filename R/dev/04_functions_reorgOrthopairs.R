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
#'   appearing in \code{genome1_gene} / \code{genome2_gene}.
#' @param n_threads Integer; number of cores to use for internal parallel steps
#'   (Linux/macOS only; on Windows this falls back to serial execution).
#' @param overwrite Logical; reserved for future use (currently ignored).
#' @param verbose Logical; whether to print progress messages.
#'
#' @return The path to the output directory,
#'   \code{file.path(working_dir, \"reorg_out\")}.
#'
#' @importFrom data.table data.table rbindlist fread fwrite setcolorder melt
#' @importFrom parallel mclapply
#' @importFrom igraph write_graph graph_from_data_frame vertex_attr decompose induced_subgraph is_connected
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
    
    .graph2df(graph = graph, out_dir = out_dir, working_dir = working_dir)
    
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
        dt_genome1 <- unique(dt[, .(original = original_genome1_gene, gene = genome1_gene)])
        dt_genome1[, genome := genomes_i[1]]
        
        dt_genome2 <- unique(dt[, .(original = original_genome2_gene, gene = genome2_gene)])
        dt_genome2[, genome := genomes_i[2]]
        
        genes_dt <- rbind(dt_genome1, dt_genome2)
        
        split_dt <- NULL
        if (rename) {
            # Identify split genes by pattern "split" in gene IDs
            genome1_split <- dt[grepl("split", genome1_gene),
                                .(original = original_genome1_gene, 
                                  gene = genome1_gene, 
                                  tx = genome1_tx)]
            genome2_split <- dt[grepl("split", genome2_gene),
                                .(original = original_genome2_gene, 
                                  gene = genome2_gene, 
                                  tx = genome2_tx)]
            
            if (nrow(genome1_split) > 0L) {
                split_dt <- rbind(
                    split_dt,
                    data.table::data.table(
                        genome = genomes_i[1],
                        original = genome1_split$original, 
                        gene = genome1_split$gene,
                        tx = genome1_split$tx
                    )
                )
            }
            if (nrow(genome2_split) > 0L) {
                split_dt <- rbind(
                    split_dt,
                    data.table::data.table(
                        genome = genomes_i[2],
                        original = genome2_split$original, 
                        gene = genome2_split$gene,
                        tx = genome2_split$tx
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
        required_cols <- c("genome1_gene", "genome2_gene",
                           "genome1_tx", "genome2_tx",
                           "mutual_ci", "class")
        missing_cols <- setdiff(required_cols, names(opr))
        if (length(missing_cols)) {
            stop("Missing required columns in ", fn, ": ",
                 paste(missing_cols, collapse = ", "))
        }
        
        if (!rename) {
            # Use original_* genes if available (pre-split), otherwise current genes
            if (all(c("original_genome1_gene", "original_genome2_gene") %in% names(opr))) {
                target_col <- c("original_genome1_gene", "original_genome2_gene",
                                "genome1_tx", "genome2_tx",
                                "mutual_ci", "class")
                opr <- as.data.frame(opr)[, target_col]
                names(opr) <- c("genome1_gene", "genome2_gene",
                                "genome1_tx", "genome2_tx",
                                "mutual_ci", "class")
            } else {
                opr <- as.data.frame(opr)[, required_cols]
            }
        } else {
            # Use current genome1_gene/genome2_gene and apply split renaming
            opr <- as.data.frame(opr)[, required_cols]
            meta <- list(genomes = genomes_i)
            opr <- .renameOrthoPair(opr = opr,
                                    meta = meta,
                                    reorg_list = reorg_list)
        }
        
        opr$genome1_id <- paste(genomes_i[1], opr$genome1_gene, sep = ":")
        opr$genome2_id <- paste(genomes_i[2], opr$genome2_gene, sep = ":")
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
        hit <- match(opr$genome1_tx, split_list$tx[genome_1])
        not_na <- !is.na(hit)
        if (any(not_na)) {
            opr$genome1_gene[not_na] <- split_list$gene[genome_1][hit[not_na]]
        }
    }
    genome_2 <- split_list$genome == meta$genomes[2]
    if (any(genome_2)) {
        hit <- match(opr$genome2_tx, split_list$tx[genome_2])
        not_na <- !is.na(hit)
        if (any(not_na)) {
            opr$genome2_gene[not_na] <- split_list$gene[genome_2][hit[not_na]]
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
    edges_df <- data.frame(from = edges_list$genome1_id,
                           to = edges_list$genome2_id,
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


## Map orthology genome labels (\code{index + 1000} from GFF) to species folder suffix
## from \code{working_dir/input/<index>_<species_name>/}.
.genomeNumericIdToSpeciesFromInput <- function(working_dir) {
    if (is.null(working_dir)) {
        return(character(0L))
    }
    input_dir <- file.path(working_dir, "input")
    if (!dir.exists(input_dir)) {
        return(character(0L))
    }
    subdirs <- list.dirs(input_dir, full.names = FALSE, recursive = FALSE)
    subdirs <- subdirs[nzchar(subdirs) & subdirs != "."]
    out <- character(0L)
    for (bn in subdirs) {
        parts <- strsplit(bn, "_", fixed = TRUE)[[1L]]
        if (length(parts) < 2L) next
        idx <- suppressWarnings(as.integer(parts[[1L]]))
        if (is.na(idx)) next
        sp <- paste(parts[-1L], collapse = "_")
        vid <- as.character(idx + 1000L)
        out[[vid]] <- sp
    }
    out
}


## Matrix rows = orthogroup entries; name columns by genome; align to full \code{genomes} order.
.groupMatrixToWideDT <- function(mat, genomes, nm2g) {
    if (!nrow(mat)) {
        empty <- stats::setNames(
            replicate(length(genomes), character(0L), simplify = FALSE),
            genomes
        )
        return(data.table::as.data.table(empty))
    }
    df <- as.data.frame(mat, stringsAsFactors = FALSE)
    for (j in seq_len(ncol(df))) {
        col_vals <- as.character(df[[j]])
        col_vals <- col_vals[!is.na(col_vals) & nzchar(col_vals)]
        gn <- if (length(col_vals)) {
            unique(unname(nm2g[col_vals]))[1L]
        } else {
            genomes[j]
        }
        colnames(df)[j] <- gn
    }
    dt <- data.table::as.data.table(df)
    miss <- setdiff(genomes, names(dt))
    for (m in miss) {
        dt[[m]] <- NA_character_
    }
    data.table::setcolorder(dt, genomes)
    dt
}


## Upper bound on \code{|V|} per connected component for exact subset enumeration (2^k).
.MAX_ENUM_VERTICES_PER_COMPONENT <- 100L


.dedupeVertexSets <- function(sets) {
    if (!length(sets)) {
        return(sets)
    }
    keys <- vapply(sets, function(s) paste(sort(s), collapse = "\x1f"), character(1L))
    sets[!duplicated(keys)]
}


.dropStrictSubsets <- function(sets) {
    sets <- .dedupeVertexSets(sets)
    n <- length(sets)
    if (n < 2L) {
        return(sets)
    }
    drop <- logical(n)
    for (i in seq_len(n)) {
        for (j in seq_len(n)) {
            if (i == j) next
            if (length(sets[[i]]) >= length(sets[[j]])) next
            if (all(sets[[i]] %in% sets[[j]])) {
                drop[i] <- TRUE
                break
            }
        }
    }
    sets[!drop]
}


.enumerateValidSetsInComponent <- function(subg,
                                           vnames,
                                           species_per_vertex,
                                           max_k = .MAX_ENUM_VERTICES_PER_COMPONENT) {
    k <- length(vnames)
    if (k < 1L) {
        return(list())
    }
    if (k > max_k) {
        warning(
            "Orthology component with ", k,
            " vertices skipped (enumeration cap is ", max_k, "). ",
            "Increase .MAX_ENUM_VERTICES_PER_COMPONENT if needed."
        )
        return(list())
    }
    res <- list()
    nout <- 0L
    for (mask in seq_len(2^k - 1L)) {
        idx <- which(as.logical(intToBits(mask))[seq_len(k)])
        sp <- species_per_vertex[idx]
        if (length(unique(sp)) != length(sp)) {
            next
        }
        sub_ind <- igraph::induced_subgraph(subg, vids = idx)
        if (!igraph::is_connected(sub_ind, mode = "weak")) {
            next
        }
        nout <- nout + 1L
        res[[nout]] <- sort(vnames[idx])
    }
    res
}


.vertexSetsToWideDT <- function(sets, genomes, nm2g) {
    if (!length(sets)) {
        empty <- stats::setNames(
            replicate(length(genomes), character(0L), simplify = FALSE),
            genomes
        )
        return(data.table::as.data.table(empty))
    }
    rows <- lapply(sets, function(S) {
        row <- stats::setNames(rep(NA_character_, length(genomes)), genomes)
        for (v in S) {
            g <- unname(nm2g[v])[1L]
            if (!is.na(g) && nzchar(g) && g %in% genomes) {
                row[[g]] <- sub("[0-9]*:", "", v)
            }
        }
        row
    })
    dt <- data.table::as.data.table(do.call(rbind, rows))
    data.table::setcolorder(dt, genomes)
    dt
}


#' @importFrom igraph vertex_attr decompose induced_subgraph is_connected
#'
#' Enumerates connected vertex sets with at most one vertex per species
#' (\code{lab_genome} via \code{nm2g}), drops sets strictly contained in another,
#' and writes \code{orthopair_list.tsv} (columns \code{sort(unique(lab_genome))},
#' cells = vertex ids). Pairwise tables use original gene symbols.
#'
#' @param working_dir Project root (parent of \code{reorg_out} and \code{input}). Used
#'   to map vertex genome IDs (\code{folder_index + 1000}) to \code{input} folder species names.
.graph2df <- function(graph, out_dir, working_dir = NULL) {
    va <- igraph::vertex_attr(graph)
    
    orig_attr <- if ("original" %in% names(va)) va$original else va$name
    name_to_orig <- stats::setNames(orig_attr, va$name)
    
    gid2sp <- .genomeNumericIdToSpeciesFromInput(working_dir)
    lab_genome <- as.character(va$genome)
    if (length(gid2sp)) {
        hit <- lab_genome %in% names(gid2sp)
        if (any(hit)) {
            lab_genome[hit] <- unname(gid2sp[lab_genome[hit]])
        }
    }
    nm2g <- stats::setNames(lab_genome, va$name)
    
    genomes <- sort(unique(lab_genome))
    n_genomes <- length(genomes)
    fn_ortho_list <- file.path(out_dir, "orthopair_list.tsv")
    
    comps <- igraph::decompose(graph, mode = "weak")
    n_comps <- sapply(comps, length)
    origins_in_subg <- lapply(comps, function(x) {
        unname(nm2g[V(x)$name])
    })
    n_origins <- sapply(origins_in_subg, function(x)length(unique(x)))
    less_than_n_genomes <- n_comps == n_origins & n_origins <= n_genomes
    comps_less_than_n_genomes <- comps[less_than_n_genomes]
    rest_comps <- comps[!less_than_n_genomes]
    comps_less_than_n_genomes <- lapply(comps_less_than_n_genomes, function(x){V(x)$name})
    wide_out <- .vertexSetsToWideDT(comps_less_than_n_genomes, gid2sp, nm2g)
    data.table::fwrite(wide_out,
                       fn_ortho_list,
                       sep = "\t",
                       quote = FALSE,
                       na = "",
                       append = FALSE)
    rm(less_than_n_genomes, comps_less_than_n_genomes, wide_out)
    
    for (subg in rest_comps) {
        x_ego <- lapply(seq_len(n_genomes), ego, graph = subg)
        x_ego <- do.call(c, x_ego)
        valid_x_ego <- vapply(x_ego, FUN.VALUE = logical(1L),
                              FUN = function(y){
                                  origins <- nm2g[names(y)]
                                  n_comps <- length(origins)
                                  n_origins <- length(unique(origins))
                                  return(n_comps == n_origins)
                              })
        valid_subgraphs <- lapply(x_ego[valid_x_ego], names)
        valid_subgraphs <- .dropStrictSubsets(valid_subgraphs)
        wide_out <- .vertexSetsToWideDT(valid_subgraphs, gid2sp, nm2g)
        data.table::fwrite(wide_out,
                           fn_ortho_list,
                           sep = "\t",
                           quote = FALSE,
                           na = "",
                           append = TRUE,
                           col.names = FALSE)
    }
    rm(graph, rest_comps, wide_out, x_ego, valid_x_ego, valid_subgraphs)
    gc()
    
    ortho_df <- data.table::fread(fn_ortho_list,
                                  sep = "\t",
                                  header = TRUE,
                                  quote = "",
                                  na.strings = "")
    
    pairwise_dir <- file.path(out_dir, "pairwise")
    dir.create(pairwise_dir, showWarnings = FALSE, recursive = TRUE)
    
    if (nrow(ortho_df) > 0L && length(genome_levels) >= 2L) {
        pairwise_genomes <- combn(seq_along(gid2sp), 2)
        for(i in seq_len(ncol(pairwise_genomes))){
            pgi <- pairwise_genomes[, i]
            pairwise_df <- subset(ortho_df, select = gid2sp[pgi])
            names(pairwise_df) <- c("genome1_gene", "genome2_gene")
            pairwise_df <- subset(pairwise_df, 
                                  subset = !is.na(genome1_gene) & !is.na(genome2_gene))
            pairwise_df <- unique(pairwise_df)
            out_fn <- file.path(pairwise_dir, 
                                paste0(paste(names(gid2sp[pgi]), collapse = "_"),
                                       ".tsv"))
            data.table::fwrite(pairwise_df,
                               out_fn,
                               sep = "\t",
                               quote = FALSE,
                               na = "")
        }
    }
    rm(ortho_df); gc(); gc()
}

