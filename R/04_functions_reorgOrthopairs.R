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
#' @param max_components Numeric; if finite, connected components with more vertices
#'   than this skip maximal colorful subgraph enumeration (fast; those components
#'   contribute no rows to \code{orthopair_list.tsv}). Default \code{Inf} keeps
#'   exact behaviour.
#'
#' @return The path to the output directory,
#'   \code{file.path(working_dir, \"reorg_out\")}.
#'
#' @importFrom data.table data.table rbindlist fread fwrite setcolorder melt
#' @importFrom parallel mclapply
#' @importFrom igraph write_graph graph_from_data_frame vertex_attr decompose induced_subgraph is_connected as_edgelist components
#' @export
reorgOrthopairs <- function(working_dir,
                            rename = TRUE,
                            n_threads = 1L,
                            overwrite = FALSE,
                            verbose = TRUE,
                            max_components = Inf) {
    stopifnot(length(working_dir) == 1L)
    stopifnot(
        length(max_components) == 1L,
        is.numeric(max_components),
        max_components >= 1 || is.infinite(max_components)
    )
    working_dir <- normalizePath(working_dir, mustWork = TRUE)
    
    orthopair_dir <- file.path(working_dir, "orthopair")
    if(!dir.exists(orthopair_dir)) {
        stop("orthopair directory not found under working_dir: ", orthopair_dir)
    }
    
    orthopair_files <- list.files(orthopair_dir, pattern = "[0-9]*_[0-9]*\\.tsv$", full.names = TRUE)
    if(length(orthopair_files) == 0L) {
        stop("No orthopair TSV files found in ", orthopair_dir,
             " (expected files named <genomeA>_<genomeB>.tsv).")
    }
    
    out_dir <- file.path(working_dir, "reorg_out")
    dir.create(path = out_dir, showWarnings = FALSE, recursive = TRUE)
    
    if(verbose) {
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
    rm(orthopair_list, reorg_list)
    gc()
    graph_fn <- file.path(out_dir, "orthopair.graphml")
    write_graph(graph = graph, file = graph_fn, format = "graphml")
    
    .graph2df(graph = graph,
              out_dir = out_dir,
              working_dir = working_dir,
              n_threads = n_threads,
              max_components = max_components)
    .writeReorgPairwiseOrphan(working_dir = working_dir, out_dir = out_dir)
    
    if(verbose) {
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
        if(length(parts) != 2L) {
            stop("Invalid orthopair TSV name (expected <genomeA>_<genomeB>.tsv): ", fn)
        }
        genomes_i <- parts
        
        dt <- fread(fn, sep = "\t", header = TRUE)
        
        # Gene sets per genome (after any renaming)
        dt_genome1 <- unique(dt[, .(original = original_genome1_gene, gene = genome1_gene)])
        dt_genome1[, genome := genomes_i[1]]
        
        dt_genome2 <- unique(dt[, .(original = original_genome2_gene, gene = genome2_gene)])
        dt_genome2[, genome := genomes_i[2]]
        
        genes_dt <- rbind(dt_genome1, dt_genome2)
        
        split_dt <- NULL
        if(rename) {
            # Identify split genes by pattern "split" in gene IDs
            genome1_split <- dt[grepl("split", genome1_gene),
                                .(original = original_genome1_gene, 
                                  gene = genome1_gene, 
                                  tx = genome1_tx)]
            genome2_split <- dt[grepl("split", genome2_gene),
                                .(original = original_genome2_gene, 
                                  gene = genome2_gene, 
                                  tx = genome2_tx)]
            
            if(nrow(genome1_split) > 0L) {
                split_dt <- rbind(
                    split_dt,
                    data.table(
                        genome = genomes_i[1],
                        original = genome1_split$original, 
                        gene = genome1_split$gene,
                        tx = genome1_split$tx
                    )
                )
            }
            if(nrow(genome2_split) > 0L) {
                split_dt <- rbind(
                    split_dt,
                    data.table(
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
    
    if(use_parallel) {
        res <- mclapply(
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
    if(rename) {
        split_pieces <- lapply(res, function(x) x$split)
        split_pieces <- Filter(Negate(is.null), split_pieces)
        if(length(split_pieces)) {
            split_list <- unique(rbindlist(split_pieces))
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
        if(length(parts) != 2L) {
            stop("Invalid orthopair TSV name (expected <genomeA>_<genomeB>.tsv): ", fn)
        }
        genomes_i <- parts
        
        if(verbose) {
            message("[reorg] Reading orthopair TSV: ", basename(fn))
        }
        opr <- fread(fn, sep = "\t", header = TRUE)
        
        # Ensure required columns exist
        required_cols <- c("genome1_gene", "genome2_gene",
                           "genome1_tx", "genome2_tx",
                           "mutual_ci", "class")
        missing_cols <- setdiff(required_cols, names(opr))
        if(length(missing_cols)) {
            stop("Missing required columns in ", fn, ": ",
                 paste(missing_cols, collapse = ", "))
        }
        
        if(!rename) {
            # Use original_* genes if available (pre-split), otherwise current genes
            if(all(c("original_genome1_gene", "original_genome2_gene") %in% names(opr))) {
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
    
    res_list <- if(use_parallel) {
        mclapply(
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
    if(is.null(split_list) || !nrow(split_list)) return(opr)
    
    genome_1 <- split_list$genome == meta$genomes[1]
    if(any(genome_1)) {
        hit <- match(opr$genome1_tx, split_list$tx[genome_1])
        not_na <- !is.na(hit)
        if(any(not_na)) {
            opr$genome1_gene[not_na] <- split_list$gene[genome_1][hit[not_na]]
        }
    }
    genome_2 <- split_list$genome == meta$genomes[2]
    if(any(genome_2)) {
        hit <- match(opr$genome2_tx, split_list$tx[genome_2])
        not_na <- !is.na(hit)
        if(any(not_na)) {
            opr$genome2_gene[not_na] <- split_list$gene[genome_2][hit[not_na]]
        }
    }
    opr
}


## Graph creation & summarisation (copied from 08_functions_postsynog.R) ----

#' @import data.table
.createGraph <- function(orthopair_list, reorg_list){
    edges_list <- rbindlist(orthopair_list)
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
    graph_out <- graph_from_data_frame(d = edges_df,
                                               vertices = vertex_df,
                                               directed = FALSE)
    graph_out
}


## Map orthology genome labels (\code{index + 1000} from GFF) to species folder suffix
## from \code{working_dir/input/<index>_<species_name>/}.
.genomeNumericIdToSpeciesFromInput <- function(working_dir) {
    if(is.null(working_dir)) {
        return(character(0L))
    }
    input_dir <- file.path(working_dir, "input")
    if(!dir.exists(input_dir)) {
        return(character(0L))
    }
    subdirs <- list.dirs(input_dir, full.names = FALSE, recursive = FALSE)
    subdirs <- subdirs[nzchar(subdirs) & subdirs != "."]
    out <- character(0L)
    for(bn in subdirs) {
        parts <- strsplit(bn, "_", fixed = TRUE)[[1L]]
        if(length(parts) < 2L) next
        idx <- suppressWarnings(as.integer(parts[[1L]]))
        if(is.na(idx)) next
        sp <- paste(parts[-1L], collapse = "_")
        vid <- as.character(idx + 1000L)
        out[[vid]] <- sp
    }
    out
}


.vertexSetsToWideDT <- function(sets, genomes, nm2g, nm2o) {
    original_cols <- paste0(genomes, "_original")
    if(!length(sets)) {
        empty_names <- c(genomes, original_cols)
        empty <- setNames(replicate(length(empty_names), character(0L), simplify = FALSE), empty_names)
        return(as.data.table(empty))
    }
    rows <- lapply(sets, function(S) {
        row_gene <- setNames(rep(NA_character_, length(genomes)), genomes)
        row_orig <- setNames(rep(NA_character_, length(original_cols)), original_cols)
        for(v in S) {
            g <- unname(nm2g[v])[1L]
            if(!is.na(g) && nzchar(g) && g %in% genomes) {
                gene_id <- sub("[0-9]*:", "", v)
                row_gene[[g]] <- gene_id
                orig_id <- unname(nm2o[v])[1L]
                if(is.na(orig_id) || !nzchar(orig_id)) orig_id <- gene_id
                row_orig[[paste0(g, "_original")]] <- orig_id
            }
        }
        c(row_gene, row_orig)
    })
    dt <- as.data.table(do.call(rbind, rows))
    setcolorder(dt, c(genomes, original_cols))
    dt
}


## Left-to-right lexicographic row order for wide orthogroup tables (gene symbols per column).
.sort_ortho_wide_dt <- function(dt, cols) {
    if(!nrow(dt) || !length(cols)) {
        return(dt)
    }
    cols <- cols[cols %in% names(dt)]
    if(!length(cols)) {
        return(dt)
    }
    setorderv(dt, cols = cols, na.last = TRUE)
    dt
}


## Scratch dir for parallel chunk TSVs: \code{working_dir/orthopair/graph2df_tmp}.
.graph2df_tmp_root <- function(working_dir, out_dir) {
    if(!is.null(working_dir)) {
        wd <- normalizePath(working_dir, mustWork = TRUE)
        root <- file.path(wd, "reorg_out", "graph2df_tmp")
    } else {
        root <- file.path(normalizePath(out_dir, mustWork = TRUE), "graph2df_tmp")
    }
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    root
}


.graph2df_vertex_maps <- function(graph, working_dir) {
    va <- vertex_attr(graph)
    genomes <- .genomeNumericIdToSpeciesFromInput(working_dir)
    lab_genome <- as.character(va$genome)
    if(length(genomes)) {
        hit <- lab_genome %in% names(genomes)
        if(any(hit)) {
            lab_genome[hit] <- unname(genomes[lab_genome[hit]])
        }
    }
    nm2g <- setNames(lab_genome, va$name)
    nm2o <- setNames(as.character(va$original), va$name)
    list(
        genomes = genomes,
        nm2g = nm2g,
        nm2o = nm2o,
        n_genomes = length(genomes)
    )
}


.graph2df_partition_components <- function(graph, nm2g, n_genomes) {
    comps <- components(graph)
    comps$origins <- nm2g[names(comps$membership)]
    origins <- tapply(comps$origins, comps$membership, unique)
    n_origins <- sapply(origins, length)
    less_than_n_genomes <- comps$csize == n_origins & n_origins <= n_genomes
    less_than_n_genomes_groups <- comps$membership %in% which(less_than_n_genomes)
    comps_small <- comps$membership[less_than_n_genomes_groups]
    rest_comps <- list(
        membership = comps$membership[!less_than_n_genomes_groups],
        origins = comps$origins[!less_than_n_genomes_groups],
        csize = comps$csize[!less_than_n_genomes]
    )
    small_groups <- tapply(names(comps_small), comps_small, c)
    list(
        small_groups = small_groups,
        rest_comps = rest_comps,
        comp_ids = unique(rest_comps$membership)
    )
}


.graph2df_write_small_block <- function(small_groups, genomes, nm2g, nm2o, fn_ortho_list) {
    out <- .vertexSetsToWideDT(small_groups, genomes, nm2g, nm2o)
    gene_cols <- genomes[genomes %in% names(out)]
    wide_colnames <- names(out)
    out <- .sort_ortho_wide_dt(out, gene_cols)
    fwrite(
        out,
        fn_ortho_list,
        sep = "\t",
        quote = FALSE,
        na = "",
        append = FALSE
    )
    wide_colnames
}


## Rows of \code{el} with both endpoints in \code{vids} (topology only; cheap vs. \code{induced_subgraph}).
.graph2df_edges_in_component <- function(el, vids) {
    keep <- el[, 1L] %in% vids & el[, 2L] %in% vids
    el[keep, , drop = FALSE]
}


## One list element per heavy component: small enough to fork without the full \code{igraph} object.
.graph2df_build_jobs <- function(comp_ids, rest_comps, full_el) {
    lapply(seq_along(comp_ids), function(cid) {
        target <- rest_comps$membership == comp_ids[cid]
        vids <- names(rest_comps$membership)[target]
        el <- .graph2df_edges_in_component(full_el, vids)
        genome_of <- rest_comps$origins[target]
        list(vids = vids, el = el, genome_of = genome_of)
    })
}


## Inputs for \code{cpp_em_max_cc_subg()} (0-based edge indices, 0-based colour ids).
.graph2df_job_to_cpp <- function(job) {
    vids <- job$vids
    el <- job$el
    genome_of <- job$genome_of
    if(!length(vids)) {
        return(NULL)
    }
    vids_sorted <- sort(vids)
    gn <- unname(genome_of[vids_sorted])
    ug <- unique(gn)
    color <- match(gn, ug) - 1L
    if(!nrow(el)) {
        el_idx <- matrix(integer(0), nrow = 0L, ncol = 2L)
    } else {
        u <- match(el[, 1L], vids_sorted) - 1L
        w <- match(el[, 2L], vids_sorted) - 1L
        el_idx <- cbind(as.integer(u), as.integer(w))
    }
    list(
        vids = vids_sorted,
        el = el_idx,
        color = as.integer(color)
    )
}


.graph2df_wide_from_job <- function(job, genomes, nm2g, nm2o, max_components = Inf) {
    nv <- length(job$vids)
    if(is.finite(max_components) && nv > max_components) {
        warning(
            "graph2df: skipped one connected component (|V| = ", nv,
            ", max_components = ", max_components,
            "). No orthogroup rows from this component.",
            call. = FALSE
        )
        return(.vertexSetsToWideDT(list(), genomes, nm2g, nm2o))
    }
    prep <- .graph2df_job_to_cpp(job)
    if(is.null(prep)) {
        return(.vertexSetsToWideDT(list(), genomes, nm2g, nm2o))
    }
    sets <- cpp_em_max_cc_subg(prep$vids, prep$el, prep$color)
    .vertexSetsToWideDT(sets, genomes, nm2g, nm2o)
}


.graph2df_merge_part_files <- function(paths_exist, wide_colnames, fn_ortho_list) {
    combined <- rbindlist(
        lapply(paths_exist, function(p) {
            t <- fread(
                p,
                sep = "\t",
                header = FALSE,
                quote = "",
                na.strings = ""
            )
            setnames(t, wide_colnames)
            t
        }),
        use.names = TRUE,
        fill = TRUE
    )
    combined <- .sort_ortho_wide_dt(combined, wide_colnames)
    fwrite(
        combined,
        fn_ortho_list,
        sep = "\t",
        quote = FALSE,
        na = "",
        append = TRUE,
        col.names = FALSE
    )
    unlink(paths_exist)
}


.graph2df_write_pairwise_dir <- function(fn_ortho_list, out_dir, genomes) {
    ortho_df <- fread(
        fn_ortho_list,
        sep = "\t",
        header = TRUE,
        quote = "",
        na.strings = ""
    )
    pairwise_dir <- file.path(out_dir, "pairwise")
    dir.create(pairwise_dir, showWarnings = FALSE, recursive = TRUE)
    if(nrow(ortho_df) > 0L && ncol(ortho_df) >= 2L) {
        pairwise_genomes <- combn(seq_along(genomes), 2)
        for(i in seq_len(ncol(pairwise_genomes))) {
            pgi <- pairwise_genomes[, i]
            pairwise_df <- subset(ortho_df, select = genomes[pgi])
            names(pairwise_df) <- c("genome1_gene", "genome2_gene")
            ocols <- paste0(genomes[pgi], "_original")
            if(all(ocols %in% names(ortho_df))) {
                pairwise_df$genome1_original_gene <- ortho_df[[ocols[1]]]
                pairwise_df$genome2_original_gene <- ortho_df[[ocols[2]]]
            } else {
                pairwise_df$genome1_original_gene <- pairwise_df$genome1_gene
                pairwise_df$genome2_original_gene <- pairwise_df$genome2_gene
            }
            pairwise_df <- subset(
                pairwise_df,
                subset = !is.na(genome1_gene) & !is.na(genome2_gene)
            )
            pairwise_df <- unique(pairwise_df)
            out_fn <- file.path(
                pairwise_dir,
                paste0(paste(names(genomes[pgi]), collapse = "_"), ".tsv")
            )
            fwrite(
                pairwise_df,
                out_fn,
                sep = "\t",
                quote = FALSE,
                na = ""
            )
        }
    }
    rm(ortho_df)
}


#' Enumerates connected vertex sets with at most one vertex per species
#' (\code{lab_genome} via \code{nm2g}), drops sets strictly contained in another,
#' and writes \code{orthopair_list.tsv} (columns \code{sort(unique(lab_genome))},
#' cells = vertex ids). Pairwise tables use original gene symbols.
#'
#' @param working_dir Project root (parent of \code{reorg_out} and \code{input}). Used
#'   to map vertex genome IDs (\code{folder_index + 1000}) to \code{input} folder species names.
#' @param n_threads Number of cores for per-component enumeration (Linux/macOS; serial on Windows or when \code{1}).
#' @param max_components See \code{max_components} on \code{reorgOrthopairs()}.
.graph2df <- function(graph,
                      out_dir,
                      working_dir = NULL,
                      n_threads = 1L,
                      max_components = Inf) {
    maps <- .graph2df_vertex_maps(graph, working_dir)
    genomes <- maps$genomes
    nm2g <- maps$nm2g
    nm2o <- maps$nm2o
    n_genomes <- maps$n_genomes
    fn_ortho_list <- file.path(out_dir, "orthopair_list.tsv")
    
    part <- .graph2df_partition_components(graph, nm2g, n_genomes)
    wide_colnames <- .graph2df_write_small_block(
        part$small_groups,
        genomes,
        nm2g,
        nm2o,
        fn_ortho_list
    )
    
    comp_ids <- part$comp_ids
    rest_comps <- part$rest_comps
    if(length(comp_ids)) {
        full_el <- as_edgelist(graph, names = TRUE)
        jobs <- .graph2df_build_jobs(comp_ids, rest_comps, full_el)
        rm(full_el, part)
        ## Drop the large graph before fork so children are not CoW-linked to it (~GB/thread).
        rm(graph)
        gc()
        
        use_parallel <- (n_threads > 1L && .Platform$OS.type != "windows")
        idx <- seq_along(comp_ids)
        
        if(use_parallel) {
            tmp_root <- .graph2df_tmp_root(working_dir, out_dir)
            on.exit(unlink(tmp_root, recursive = TRUE), add = TRUE)
            part_paths <- vapply(idx,
                                 function(i) {
                                     tempfile(pattern = sprintf("g2df_%04d_", i),
                                              tmpdir = tmp_root,
                                              fileext = ".tsv"
                                     )
                                 },
                                 character(1L))
            mclapply(idx,
                               function(cid) {
                                   out_dt <- .graph2df_wide_from_job(
                                       jobs[[cid]],
                                       genomes,
                                       nm2g,
                                       nm2o,
                                       max_components = max_components)
                                   pp <- part_paths[cid]
                                   if(nrow(out_dt) > 0L) {
                                       fwrite(
                                           out_dt,
                                           pp,
                                           sep = "\t",
                                           quote = FALSE,
                                           na = "",
                                           col.names = FALSE
                                       )
                                   }
                                   invisible(NULL)
                               },
                               mc.cores = n_threads,
                               mc.preschedule = FALSE
            )
            paths_exist <- part_paths[file.exists(part_paths)]
            if(length(paths_exist)) {
                .graph2df_merge_part_files(paths_exist, wide_colnames, fn_ortho_list)
            }
            rm(jobs, part_paths)
        } else {
            parts <- lapply(
                jobs,
                .graph2df_wide_from_job,
                genomes = genomes,
                nm2g = nm2g,
                nm2o = nm2o,
                max_components = max_components
            )
            combined <- rbindlist(parts, use.names = TRUE, fill = TRUE)
            rm(jobs, parts)
            if(nrow(combined) > 0L) {
                gene_cols <- genomes[genomes %in% names(combined)]
                combined <- .sort_ortho_wide_dt(combined, gene_cols)
                fwrite(
                    combined,
                    fn_ortho_list,
                    sep = "\t",
                    quote = FALSE,
                    na = "",
                    append = TRUE,
                    col.names = FALSE
                )
            }
            rm(combined)
        }
        rm(rest_comps, comp_ids, idx, use_parallel)
    } else {
        rm(graph, part)
    }
    gc()
    
    .graph2df_write_pairwise_dir(fn_ortho_list, out_dir, genomes)
    gc()
}

.writeReorgPairwiseOrphan <- function(working_dir, out_dir) {
    pairwise_dir <- file.path(out_dir, "pairwise")
    pairwise_files <- list.files(pairwise_dir, pattern = "\\.tsv$", full.names = TRUE)
    if(!length(pairwise_files)) {
        return(invisible(NULL))
    }
    input_dir <- file.path(working_dir, "input")
    input_folders <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)
    gid_map <- lapply(input_folders, function(d){
        bn <- basename(d)
        idx <- suppressWarnings(as.integer(sub("_.*", "", bn)))
        if(is.na(idx)) return(NULL)
        list(gid = as.character(idx + 1000L), gff = file.path(d, "gff_df.rds"))
    })
    gid_map <- Filter(Negate(is.null), gid_map)
    gff_lookup <- setNames(vapply(gid_map, `[[`, character(1), "gff"),
                           vapply(gid_map, `[[`, character(1), "gid"))
    get_all_gene <- function(gid){
        fn <- gff_lookup[[gid]]
        if(is.null(fn) || !file.exists(fn)) return(character(0))
        x <- readRDS(fn)
        if(is.list(x) && length(x) >= 1L) x <- x[[1L]]
        if(!is.data.frame(x) || !("gene_id" %in% names(x))) return(character(0))
        unique(as.character(x$gene_id[!is.na(x$gene_id)]))
    }
    
    out_dir_orphan <- file.path(working_dir, "reorg_orphan")
    dir.create(out_dir_orphan, showWarnings = FALSE, recursive = TRUE)
    for(fn in pairwise_files){
        pair_id <- sub("\\.tsv$", "", basename(fn))
        parts <- strsplit(pair_id, "_", fixed = TRUE)[[1L]]
        if(length(parts) != 2L) next
        dt <- fread(fn, sep = "\t", header = TRUE)
        if(!all(c("genome1_gene", "genome2_gene") %in% names(dt))) next
        g1 <- parts[1L]
        g2 <- parts[2L]
        all1 <- get_all_gene(g1)
        all2 <- get_all_gene(g2)
        found1 <- unique(as.character(dt$genome1_gene[!is.na(dt$genome1_gene)]))
        found2 <- unique(as.character(dt$genome2_gene[!is.na(dt$genome2_gene)]))
        out <- rbind(
            data.frame(genome = rep(g1, length(setdiff(all1, found1))),
                       gene = setdiff(all1, found1),
                       stringsAsFactors = FALSE),
            data.frame(genome = rep(g2, length(setdiff(all2, found2))),
                       gene = setdiff(all2, found2),
                       stringsAsFactors = FALSE)
        )
        if(!nrow(out)) {
            out <- data.frame(genome = character(0), gene = character(0), stringsAsFactors = FALSE)
        }
        fwrite(out, file.path(out_dir_orphan, paste0(pair_id, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    }
    invisible(TRUE)
}
