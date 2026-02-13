#' Reorganise orthopair results and create graph/table outputs
#'
#' This is a high-level wrapper that:
#' \itemize{
#'   \item reorganises per-pair OrthoPaiR HDF5 results,
#'   \item rewrites GFF information when \code{rename = TRUE} to merge split genes,
#'   \item creates a genome-level orthology graph (GraphML),
#'   \item and summarises the graph into tabular outputs in the same directory.
#' }
#'
#' It is intended for use after pairwise OrthoPaiR runs have completed and
#' their HDF5 paths are available.
#'
#' @param hdf5_fn Character vector of paths to pairwise OrthoPaiR HDF5 result
#'   files.
#' @param rename Logical; whether to apply split-gene renaming when
#'   reorganising GFF and orthopair information.
#' @param n_threads Integer; number of cores to use for the internal parallel
#'   steps (Linux/macOS only; on Windows this falls back to serial execution).
#' @param overwrite Logical; if \code{TRUE}, recompute the reorganised outputs
#'   and overwrite existing files in the output directory. If \code{FALSE},
#'   an existing \code{orthopair.graphml} will be reused if present.
#' @param verbose Logical; whether to print progress messages.
#'
#' @return The path to the output directory (typically
#'   \code{file.path(<base>, \"reorg_out\")}).
#'
#' @export
reorgOrthopairs <- function(hdf5_fn,
                            rename = TRUE,
                            n_threads = 1L,
                            overwrite = FALSE,
                            verbose = TRUE) {
    out_dir <- file.path(sub("hdf5_out.*", "", hdf5_fn[1]), "reorg_out")
    dir.create(path = file.path(out_dir, "gff"),
               showWarnings = FALSE,
               recursive = TRUE)
    
    ## Faster rename-list construction
    rename_list <- .getRenameList(hdf5_fn = hdf5_fn,
                                  rename = rename,
                                  n_threads = n_threads)
    
    reorg_list <- .reorgGFF(out_dir = out_dir,
                            rename_list = rename_list,
                            n_threads = n_threads)
    
    orthopair_list <- .reorgOrthoPair(hdf5_fn = hdf5_fn,
                                      reorg_list = reorg_list,
                                      n_threads = n_threads)
    
    graph <- .createGraph(orthopair_list = orthopair_list, 
                          reorg_list = reorg_list)
    graph_fn <- file.path(out_dir, "orthopair.graphml")
    write_graph(graph = graph, file = graph_fn, format = "graphml")
    
    .graph2df(graph = graph, out_dir = out_dir)
    
    return(out_dir)
}

## Internal helper: faster version of .getRenameList()
##
## - Uses lapply / mclapply to gather metadata and (optionally) split genes
##   in a single pass over hdf5_fn.
## - Fixes the dependency on a loop-local 'meta' variable in the original
##   implementation by keeping per-file meta information explicitly.
.getRenameList <- function(hdf5_fn, rename = FALSE, n_threads = 1L) {
    idx <- seq_along(hdf5_fn)
    
    use_parallel <- (n_threads > 1L && .Platform$OS.type != "windows")
    worker_fun <- function(i) {
        meta_i <- getMeta(hdf5_fn = hdf5_fn[i])
        
        ## Access components defensively to avoid '$' on atomic vectors
        genomes_i <- meta_i[["genomes"]]
        files_i <- meta_i[["files"]]
        
        if (!is.null(files_i)) {
            query_gff_i <- files_i[["query_gff"]]
            subject_gff_i <- files_i[["subject_gff"]]
        } else {
            query_gff_i <- NA_character_
            subject_gff_i <- NA_character_
        }
        
        gff_df <- data.frame(
            genomes = genomes_i,
            gff = c(query_gff_i, subject_gff_i),
            stringsAsFactors = FALSE
        )
        
        split_df <- NULL
        if (rename) {
            opr <- getOrthoPair(hdf5_fn = hdf5_fn[i])[[1]]
            
            query_split <- subset(
                opr,
                subset = grepl("split", query_gene),
                select = c(query_gene, query_tx)
            )
            subject_split <- subset(
                opr,
                subset = grepl("split", subject_gene),
                select = c(subject_gene, subject_tx)
            )
            
            if (nrow(query_split) > 0L) {
                split_df <- rbind(
                    split_df,
                    data.frame(
                        genome = genomes_i[1],
                        gene = query_split$query_gene,
                        tx = query_split$query_tx,
                        stringsAsFactors = FALSE
                    )
                )
            }
            if (nrow(subject_split) > 0L) {
                split_df <- rbind(
                    split_df,
                    data.frame(
                        genome = genomes_i[2],
                        gene = subject_split$subject_gene,
                        tx = subject_split$subject_tx,
                        stringsAsFactors = FALSE
                    )
                )
            }
        }
        
        list(gff = gff_df, split = split_df)
    }
    
    if (use_parallel) {
        res <- mclapply(
            X = idx,
            FUN = worker_fun,
            mc.cores = n_threads
        )
    } else {
        res <- lapply(idx, worker_fun)
    }
    
    ## Combine gff_list
    gff_list <- do.call(
        rbind,
        lapply(res, function(x) x$gff)
    )
    gff_list <- unique(gff_list)
    
    ## Combine split_list if rename = TRUE
    split_list <- NULL
    if (rename) {
        split_pieces <- lapply(res, function(x) x$split)
        split_pieces <- Filter(Negate(is.null), split_pieces)
        if (length(split_pieces)) {
            split_list <- unique(do.call(rbind, split_pieces))
        }
    }
    
    list(gff_list = gff_list, split_list = split_list)
}

## Faster version of .reorgGFF()
##
## Parallelises per-genome GFF reorganisation when possible; each genome’s GFF
## is handled independently and written to its own output file, so there are
## no write-contention issues.
.reorgGFF <- function(out_dir, rename_list, n_threads = 1L) {
    gff_paths <- rename_list$gff_list$gff
    genomes <- rename_list$gff_list$genomes
    split_list <- rename_list$split_list
    
    idx <- seq_along(gff_paths)
    use_parallel <- (n_threads > 1L && .Platform$OS.type != "windows")
    
    worker_fun <- function(i) {
        i_genome <- genomes[i]
        out_gff <- .setSplitGFF(gff_fn = gff_paths[i], 
                                i_genome = i_genome, 
                                split_list = split_list,
                                out_dir = out_dir)
        
        return(data.table(gene = unique(out_gff$gene_id),
                          genome = i_genome))
    }
    
    reorg_list <- list(split_list = split_list)
    reorg_list$genes_list <- if (use_parallel) {
        mclapply(
            X = idx,
            FUN = worker_fun,
            mc.cores = n_threads
        )
    } else {
        lapply(idx, worker_fun)
    }
    return(reorg_list)
}

.setSplitGFF <- function(gff_fn, i_genome, split_list, out_dir){
    gff <- readRDS(gff_fn)
    gff$Target <- gff$tx_index <- gff$gene_index <- NULL
    out_gff_fn <- file.path(out_dir, "gff", paste0(i_genome, ".gff.rds"))
    hit_genome <- split_list$genome == i_genome
    if (sum(hit_genome) == 0) {
        file.copy(from = gff_fn, to = out_gff_fn, overwrite = TRUE)
        out_gff <- gff
        
    } else {
        tx_hit <- match(gff$ID, split_list$tx[hit_genome])
        split_gene_gff <- split_gene_tx <- gff[!is.na(tx_hit)]
        split_gene_gff$type <- "gene"
        split_gene_gff$ID <- paste0(split_gene_gff$ID, ":split")
        split_gene_gff$Parent <- ""
        
        element_hit <- match(gff$Parent, split_gene_tx$ID)
        split_gene_tx$Parent <- paste0(split_gene_tx$ID, ":split")
        split_gene_tx$ID <- paste0(split_gene_gff$ID, ":split:tx")
        
        split_gene_element <- gff[!is.na(element_hit)]
        split_gene_element$ID <- paste0(split_gene_element$Parent,
                                        ":split:",
                                        split_gene_element$type)
        split_gene_element$Parent <- paste0(split_gene_element$Parent, 
                                            ":split:tx")
        out_gff <- c(gff, split_gene_gff, split_gene_tx, split_gene_element)
        out_gff <- .fixGFF(gff = out_gff)
        ## Save reorganised GFF as RDS instead of exporting to GFF3 text
        saveRDS(object = out_gff, file = out_gff_fn)
    }
    return(out_gff)
}

## Faster version of .reorgOrthoPair()
##
## This variant parallelises the expensive per-pair work (reading orthopairs,
## meta data, and applying rename mappings) while keeping all HDF5 writes and
## orphan calculations in a single process to avoid concurrency issues with
## the HDF5 backend.
.reorgOrthoPair <- function(hdf5_fn,
                            reorg_list,
                            n_threads = 1L) {
    idx <- seq_along(hdf5_fn)
    use_parallel <- (n_threads > 1L && .Platform$OS.type != "windows")
    
    worker_fun <- function(i) {
        opr <- getOrthoPair(hdf5_fn = hdf5_fn[i], score = TRUE, loc = FALSE)[[1]]
        meta <- getMeta(hdf5_fn = hdf5_fn[i])
        if(is.null(reorg_list$split_list)){
            target_col <- c("original_query_gene", "original_subject_gene",
                            "query_tx", "subject_tx",
                            "mutual_ci", "class")
            opr <- opr[, target_col]
            names(opr) <- c("query_gene", "subject_gene",
                            "query_tx", "subject_tx",
                            "mutual_ci", "class")
            
        } else {
            target_col <- c("query_gene", "subject_gene",
                            "query_tx", "subject_tx",
                            "mutual_ci", "class")
            opr <- opr[, target_col]
            opr <- .renameOrthoPair(opr = opr,
                                    meta = meta,
                                    reorg_list = reorg_list)
        }
        opr$query_id <- paste(meta$genomes[1], opr$query_gene, sep = ":")
        opr$subject_id <- paste(meta$genomes[2], opr$subject_gene, sep = ":")
        return(opr)
    }
    
    res_list <- if (use_parallel) {
        mclapply(
            X = idx,
            FUN = worker_fun,
            mc.cores = n_threads
        )
    } else {
        lapply(idx, worker_fun)
    }
    
    return(res_list)
}

.renameOrthoPair <- function(opr, meta, reorg_list){
    genome_1 <- reorg_list$split_list$genome == meta$genomes[1]
    if (any(genome_1)) {
        hit <- match(opr$query_tx, reorg_list$split_list$tx[genome_1])
        not_na <- !is.na(hit)
        if (any(not_na)) {
            opr$query_gene[not_na] <- reorg_list$split_list$gene[genome_1][hit[not_na]]
        }
    }
    genome_2 <- reorg_list$split_list$genome == meta$genomes[2]
    if (any(genome_2)) {
        hit <- match(opr$subject_tx, reorg_list$split_list$tx[genome_2])
        not_na <- !is.na(hit)
        if (any(not_na)) {
            opr$subject_gene[not_na] <- reorg_list$split_list$gene[genome_2][hit[not_na]]
        }
    }
    return(opr)
}


#' @import data.table
.createGraph <- function(orthopair_list, reorg_list){
    edges_list <- rbindlist(orthopair_list)
    vertex_list <- rbindlist(reorg_list$genes_list)
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
                            stringsAsFactors = FALSE)
    
    ## Build graph in a single call
    graph_out <- graph_from_data_frame(d = edges_df,
                                       vertices = vertex_df,
                                       directed = FALSE
    )
    return(graph_out)
}

#' @importFrom igraph edge_attr ends E vertex_attr
#' @importFrom dplyr full_join
#'
.graph2df <- function(graph, out_dir) {
    comp <- components(graph)
    va <- vertex_attr(graph)
    ea <- edge_attr(graph)
    em <- ends(graph, E(graph), names = FALSE)
    rm(graph); gc(); gc()
    
    dtv <- data.table(membership = comp$membership,
                      genome = va$genome,
                      gene = va$name)
    fwrite(dtv,
           file.path(out_dir, "orthopair_list_long.tsv"),
           sep = "\t",
           quote = FALSE,
           na = "")
    
    genome_levels <- unique(dtv$genome)
    fwrite(setcolorder(dcast(dtv[, 
                                 .(gene = paste(gene, collapse=",")), 
                                 by = .(membership, genome)],
                             membership ~ genome,
                             value.var = "gene",
                             fill = NA_character_),
                       c("membership", genome_levels)),
           file.path(out_dir, "orthopair_list.tsv"),
           sep = "\t",
           quote = FALSE,
           na = "")
    eid <- dtv$membership[em[, 1]]
    g1 <- va$genome[em[, 1]]
    g2 <- va$genome[em[, 2]]
    rm(comp, va, em); gc(); gc()
    
    v_sum <- dtv[, .(nV = .N, nGenome = uniqueN(genome)), by = membership]
    rm(dtv); gc(); gc()
    
    dte <- data.table(membership = eid, mutual_ci = ea$mutual_ci)
    
    e_sum <- dte[, .(nE = .N,
                     max_mci = max(mutual_ci, na.rm = TRUE),
                     q1 = quantile(mutual_ci, 0.25, na.rm = TRUE),
                     median_mci = median(mutual_ci, na.rm = TRUE),
                     q3 = quantile(mutual_ci, 0.75, na.rm = TRUE),
                     min_mci = min(mutual_ci, na.rm = TRUE),
                     mean_mci = mean(mutual_ci, na.rm = TRUE),
                     sd_mci = sd(mutual_ci, na.rm = TRUE)),
                 by = membership]
    
    out <- full_join(v_sum, e_sum, "membership")
    fwrite(out,
           file.path(out_dir, "orthogroup_stats.tsv"),
           sep = "\t",
           quote = FALSE,
           na = "")
    
    genome1 <- pmin(g1, g2)
    genome2 <- pmax(g1, g2)
    dt_pair <- data.table(genome1 = genome1,
                          genome2 = genome2,
                          mutual_ci = ea$mutual_ci)
    pair_sum <- dt_pair[, .(nE = .N,
                            max_mci = max(mutual_ci, na.rm = TRUE),
                            q1 = quantile(mutual_ci, 0.75, na.rm = TRUE),
                            median_mci= median(mutual_ci, na.rm = TRUE),
                            q3 = quantile(mutual_ci, 0.25, na.rm = TRUE),
                            min_mci = min(mutual_ci, na.rm = TRUE),
                            mean_mci = mean(mutual_ci, na.rm = TRUE),
                            sd_mci = sd(mutual_ci, na.rm = TRUE)),
                        by = .(genome1, genome2)][order(-nE)]
    
    fwrite(
        pair_sum,
        file.path(out_dir, "genomepair_edge_stats.tsv"),
        sep = "\t",
        quote = FALSE,
        na = ""
    )
    
    rm(dt_pair, pair_sum, g1, g2, genome1, genome2); gc(); gc()
}
