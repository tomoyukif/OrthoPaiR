#' @importFrom rtracklayer export.gff3
#' @importFrom Biostrings writeXStringSet
#' @export
reorgOrthopiars <- function(hdf5_fn,
                            out_dir = "./",
                            out_fn = "reorg_orthopair.h5",
                            rename = FALSE,
                            overwrite = TRUE,
                            verbose = TRUE){
    dir.create(path = out_dir, showWarnings = FALSE, recursive = TRUE)
    hdf5_out_fn <- .makeHDF5(hdf5_path = file.path(out_dir, out_fn), 
                             overwrite = overwrite)
    if(overwrite){
        .h5creategroup(hdf5_out_fn, "orthopair")
        .h5creategroup(hdf5_out_fn, "gff")
    }
    
    rename_list <- .getRenameList(hdf5_fn, rename)
    
    if(rename){
        .reorgGFF(out_dir = out_dir, rename_list = rename_list)
    }
    
    .reorgOrthoPair(hdf5_out_fn = hdf5_out_fn,
                    hdf5_fn = hdf5_fn, 
                    out_dir = out_dir,
                    rename_list = rename_list,
                    rename = rename)
    
    .h5overwrite(obj = "reorg_orthopair",
                 file = hdf5_out_fn,
                 name = "data_type")
    
    return(hdf5_out_fn)
}

.getRenameList <- function(hdf5_fn, rename){
    gff_list <- NULL
    for(i in seq_along(hdf5_fn)){
        meta <- getMeta(hdf5_fn = hdf5_fn[i])
        gff_list <- rbind(gff_list,
                          data.frame(genomes = meta$genomes, 
                                     gff = c(meta$files$query_gff, 
                                             meta$files$subject_gff)))
    }
    gff_list <- unique(gff_list)
    
    split_list <- NULL
    if(rename){
        for(i in seq_along(hdf5_fn)){
            opr <- getOrthoPair(hdf5_fn = hdf5_fn[i])
            query_split <- subset(opr[[1]], 
                                  subset = grepl("split", query_gene),
                                  select = c(query_gene, query_tx))
            subject_split <- subset(opr[[1]], 
                                    subset = grepl("split", subject_gene),
                                    select = c(subject_gene, subject_tx))
            if(nrow(query_split) > 0){
                split_list <- rbind(split_list, 
                                    data.frame(genome = meta$genomes[1], 
                                               gene = query_split$query_gene, 
                                               tx = query_split$query_tx))
            }
            if(nrow(subject_split) > 0){
                split_list <- rbind(split_list, 
                                    data.frame(genome = meta$genomes[2], 
                                               gene = subject_split$subject_gene, 
                                               tx = subject_split$subject_tx))
            }
            split_list <- unique(split_list)
        }
    }
    
    return(list(gff_list = gff_list, split_list = split_list))
}

.reorgGFF <- function(out_dir, rename_list){
    for(i in seq_along(rename_list$gff_list$gff)){
        gff <- readRDS(rename_list$gff_list$gff[i])
        gff$Target <- gff$tx_index <- gff$gene_index <- NULL
        i_genome <- rename_list$gff_list$genomes[i]
        hit_genome <- rename_list$split_list$genome == i_genome
        tx_hit <- match(gff$ID, rename_list$split_list$tx[hit_genome])
        gff_split_gene <- unique(gff$Parent[!is.na(tx_hit)])
        if(length(gff_split_gene) > 0){
            gff$Parent[!is.na(tx_hit)] <- rename_list$split_list$gene[hit_genome][na.omit(tx_hit)]
            split_gene_gff <- gff[!is.na(tx_hit)]
            split_gene_gff$type <- "gene"
            split_gene_gff$ID <- split_gene_gff$Parent
            split_gene_gff$Parent <- ""
            out <- split_gene_gff
            split_gene_tx <- gff[gff$Parent %in% gff_split_gene]
            gff <- gff[!gff$ID %in% gff_split_gene]
            gff <- gff[!gff$Parent %in% gff_split_gene]
            out <- c(out, gff)
            ol <- findOverlaps(split_gene_tx,
                               gff[gff$type %in% c("mRNA", "transcript")])
            non_ol_tx <- split_gene_tx[-queryHits(ol)]
            if(length(non_ol_tx) > 0){
                non_ol_gene_gff <- non_ol_tx
                non_ol_gene_gff$type <- "gene"
                non_ol_gene_gff$ID <- non_ol_gene_gff$Parent
                non_ol_gene_gff$Parent <- ""
                out <- c(out, non_ol_gene_gff, non_ol_tx)
            }
            retain_entry <- out$Parent %in% c(out$ID, "")
            out <- out[retain_entry]
            out <- .fixGFF(gff = out)
            
        } else {
            out <- gff
        }
        
        out_gff_fn <- file.path(out_dir, paste0(i_genome, ".gff"))
        rename_list$gff_list$gff[i] <- out_gff_fn
        export.gff3(object = out, out_gff_fn)
    }
    return(rename_list)
}

.reorgOrthoPair <- function(hdf5_out_fn, hdf5_fn, out_dir, rename_list, rename){
    gene_list <- NULL
    for(i in seq_along(rename_list$gff_list$genomes)){
        i_genome <- rename_list$gff_list$genomes[i]
        i_gff_fn <- rename_list$gff_list$gff[i]
        if(grepl("\\.rds$", i_gff_fn)){
            gff <- readRDS(i_gff_fn)
            
        } else {
            gff <- import.gff3(con = i_gff_fn)
        }
        gene_list <- rbind(gene_list, 
                           data.frame(genome = i_genome,
                                      gene = unique(gff$gene_id)))
        .h5overwrite(obj = i_gff_fn,
                     file = hdf5_out_fn,
                     name = paste0("gff/", i_genome))
    }
    .h5overwrite(obj = gene_list,
                 file = hdf5_out_fn,
                 name = "gene_list")
    
    orphan_list <- gene_list
    for(i in seq_along(hdf5_fn)){
        opr <- getOrthoPair(hdf5_fn = hdf5_fn[i], score = TRUE, loc = TRUE)[[1]]
        meta <- getMeta(hdf5_fn = hdf5_fn[i])
        
        if(rename){
            genome_1 <- rename_list$split_list$genome == meta$genomes[1]
            if(sum(genome_1) > 0){
                hit <- match(opr$query_tx, rename_list$split_list$tx[genome_1])
                opr$query_gene[!is.na(hit)] <- rename_list$split_list$gene[genome_1][na.omit(hit)]
            }
            genome_2 <- rename_list$split_list$genome == meta$genomes[2]
            if(sum(genome_2) > 0){
                hit <- match(opr$subject_tx, rename_list$split_list$tx[genome_2])
                opr$subject_gene[!is.na(hit)] <- rename_list$split_list$gene[genome_2][na.omit(hit)]
            }
            
        } else {
            opr$query_gene <- opr$original_query_gene
            opr$subject_gene <- opr$original_subject_gene
        }
        .h5overwrite(obj = opr,
                     file = hdf5_out_fn,
                     name = paste0("orthopair/", meta$pair_id))
        
        genome_1 <- orphan_list$genome == meta$genomes[1]
        if(sum(genome_1) > 0){
            hit <- orphan_list$gene[genome_1] %in% opr$query_gene
            orphan_list$gene[genome_1][hit] <- NA
        }
        genome_2 <- orphan_list$genome == meta$genomes[2]
        if(sum(genome_2) > 0){
            hit <- orphan_list$gene[genome_2] %in% opr$subject_gene
            orphan_list$gene[genome_2][hit] <- NA        
        }
    }
    orphan_list <- subset(orphan_list, subset = !is.na(gene))
    .h5overwrite(obj = orphan_list,
                 file = hdf5_out_fn,
                 name = "orphan_list")
}

################################################################################
#' @importFrom igraph make_empty_graph add_edges add_vertices V<-
#' @export
#' 
makeOrthoGraph <- function(hdf5_fn){
    h5 <- H5Fopen(hdf5_fn)
    on.exit(H5Fclose(h5))
    
    check <- h5$data_type != "reorg_orthopair"
    if(check){
        stop("This function onyl accepts an output hdf5 file created",
             " by the reorgOrthopiars() function.")
    }
    n_files <- length(h5$orthopair)
    orphan_gene <- h5$orphan_list
    gene_list <- h5$gene_list
    graph <- make_empty_graph(n = 0, directed = FALSE)
    orthopair_gene <- h5$orthopair
    check_orthopair <- sapply(orthopair_gene, function(x){x[1, 1] == "NA"})
    orthopair_gene <- orthopair_gene[!check_orthopair]
    orthopair_gene <- lapply(orthopair_gene, subset,
                             select = c(query_gene, subject_gene, 
                                        mutual_ci, class))
    orthopair_gene <- do.call("rbind", orthopair_gene)
    edges <- subset(orthopair_gene, select = c(query_gene, subject_gene))
    genome1<- gene_list$genome[match(edges$query_gene, gene_list$gene)]
    genome2<- gene_list$genome[match(edges$subject_gene, gene_list$gene)]
    edges <- as.vector(t(edges))
    vertices <- unique(edges)
    genome <- gene_list$genome[match(vertices, gene_list$gene)]
    graph <- add_vertices(graph = graph, 
                          nv = length(vertices),
                          name = vertices,
                          genome = genome)
    
    graph <- add_edges(graph = graph,
                       edges = edges, 
                       attr = list(mutual_ci = orthopair_gene$mutual_ci,
                                   class = orthopair_gene$class,
                                   genome1 = genome1,
                                   genome2 = genome2))
    graph <- add_vertices(graph = graph, 
                          nv = length(orphan_gene$gene),
                          name = orphan_gene$gene,
                          genome = orphan_gene$genome)
    return(graph)
}

#' @importFrom igraph components ego V subgraph vcount induced_subgraph vertex_attr vertex_attr<-
#' @importFrom tidyr pivot_wider
#' @export
graph2df <- function(hdf5_fn, graph, orthopair_fn, n_core = 1, n_batch = 50){
    h5 <- H5Fopen(hdf5_fn)
    on.exit(H5Fclose(h5))
    gene_list <- h5$gene_list
    gene_list$assign <- FALSE
    genomes <- sort(unique(gene_list$genome))
    n_genomes <- length(genomes)
    v <- vertex_attr(graph)
    v$name <- paste(v$genome, v$name, sep = "&")
    vertex_attr(graph) <- v
    ego_list <- ego(graph = graph, order = n_genomes)
    ego_list <- lapply(ego_list, names)
    ego_list_genomes <- lapply(ego_list, sub, pattern = "&.+", replace = "")
    ego_list_genomes_len <- sapply(ego_list_genomes, length)
    ego_list_genomes_unique <- lapply(ego_list_genomes, unique)
    ego_list_genomes_unique_len <- sapply(ego_list_genomes_unique, length)
    all_genomes <- ego_list_genomes_unique_len == n_genomes
    solo_entries <- ego_list_genomes_len == n_genomes
    is_full_pair <- all_genomes & solo_entries
    valid_group1 <- ego_list[is_full_pair]
    valid_group1 <- lapply(valid_group1, sort)
    valid_group1 <- do.call(rbind, valid_group1)
    valid_group1_dup <- duplicated(valid_group1)
    valid_group1 <- valid_group1[!valid_group1_dup, ]
    
    rest_list <- ego_list[!is_full_pair]
    rest_list_genomes <- lapply(rest_list, sub, pattern = "&.+", replace = "")
    rest_list_genomes <- lapply(rest_list_genomes, factor, levels = genomes)
    rest_list <- mapply(split, x = rest_list, f = rest_list_genomes, SIMPLIFY = FALSE)
    rest_list <- lapply(rest_list, function(x){
        x[sapply(x, length) == 0] <- NA
        return(x)
    })
    valid_group2 <- lapply(rest_list, expand.grid)
    n_valid_group2 <- lapply(valid_group2, nrow)
    valid_group2_sog <- mapply(rep, seq_along(valid_group2), n_valid_group2)
    valid_group2 <- lapply(valid_group2, t)
    valid_group2 <- unlist(valid_group2)
    valid_group2 <- matrix(valid_group2, nrow = n_genomes)
    valid_group2 <- t(valid_group2)
    valid_group2_dup <- duplicated(valid_group2)
    valid_group2 <- valid_group2[!valid_group2_dup, ]
    valid_group2_sog <- unlist(valid_group2_sog)[!valid_group2_dup] + nrow(valid_group1)
    out <- rbind(data.frame(valid_group1, SOG = seq_len(nrow(valid_group1))),
                 data.frame(valid_group2, SOG = valid_group2_sog))
    out_label <- apply(out[, -ncol(out)], 2, sub, pattern = "&.+", replace = "")
    out_label <- apply(out_label, 2, unique)
    out_label <- apply(out_label, 2, na.omit)
    out_values <- apply(out[, -ncol(out)], 2, sub, pattern = ".+&", replace = "")
    out[, -ncol(out)] <- out_values
    colnames(out)[-ncol(out)] <- out_label
    write.csv(out, file = orthopair_fn, row.names = FALSE)
    invisible(orthopair_fn)
}

.writeEntries <- function(x, file, gene_list){
    x <- .setTxID(x = x, gene_list = gene_list)
    if(file.exists(file)){
        write.table(x = x, 
                    file = file,
                    append = TRUE, 
                    quote = FALSE,
                    sep = ",", 
                    col.names = FALSE,
                    row.names = FALSE)
        
    } else {
        write.table(x = x, 
                    file = file,
                    append = FALSE, 
                    quote = FALSE,
                    sep = ",", 
                    col.names = TRUE,
                    row.names = FALSE)
    }
}

.is_mult_mp_in_grp <- function(x){
    mp <- grep("_MP[0-9]+", x, value = TRUE)
    out <- FALSE
    if(length(mp) > 1){
        check <- length(unique(sub("_.+", "", mp))) > 1
        if(check){
            out <- TRUE
        }
    }
    return(out)
}

.groupGraph <- function(graph){
    grp <- split(V(graph)$name, components(graph)$membership)
    mult_mp_in_grp <- sapply(grp, .is_mult_mp_in_grp)
    
    n_gene_in_grp <- lapply(grp, length)
    out <- sapply(seq_along(n_gene_in_grp), function(i){
        return(rep(i, n_gene_in_grp[[i]]))
    })
    out <- data.frame(gene_id = unlist(grp), SOG = unlist(out))
    hit <- match(out$SOG, as.numeric(names(mult_mp_in_grp)))
    out$mult_mp_in_grp <- mult_mp_in_grp[hit]
    return(out)
}

#' @importFrom igraph max_cliques
.getGraphCliques <- function(graph){
    graph_cliques <- max_cliques(graph = graph)
    graph_cliques <- lapply(graph_cliques, attributes)
    graph_cliques <- lapply(graph_cliques, "[[", "names")
    n_cliques <- sapply(graph_cliques, length)
    out <- list(cliques = graph_cliques, n = n_cliques)
    return(out)
}

.orgFullConnected <- function(graph, genomes, gene_list, n_len, full){
    graph_cliques <- .getGraphCliques(graph = graph)
    is_full_con <- graph_cliques$n == n_len
    if(sum(is_full_con) == 0){
        out <- list(rest_genes = unlist(graph_cliques$cliques[!is_full_con]),
                    full_con = NULL)
        return(out)
    }
    full_con <- graph_cliques$cliques[is_full_con]
    full_con <- data.frame(gene_id = unlist(full_con),
                           grp_id = rep(seq_along(full_con),
                                        each = n_len))
    hit <- match(full_con$gene_id, gene_list$gene)
    full_con$genome <- gene_list$genome[hit]
    full_df <- data.frame(id = paste(rep(unique(full_con$grp_id),
                                         each = length(genomes)),
                                     genomes,
                                     sep = "_"))
    full_con$id <- paste(full_con$grp_id, full_con$genome, sep = "_")
    full_con <- left_join(x = full_df, y = full_con, by = "id")
    
    full_con <- matrix(data = full_con$gene_id, 
                       ncol = length(genomes),
                       byrow = TRUE)
    full_con <- as.data.frame(full_con)
    names(full_con) <- genomes
    
    if(full){
        full_con$class <- "full_connected"
        
    } else {
        full_con$class <- "partial_connected"
    }
    
    rest_genes <- unlist(graph_cliques$cliques[!is_full_con])
    if(length(rest_genes) == 0){
        rest_genes <- NULL
    }
    out <- list(rest_genes = rest_genes,
                full_con = full_con)
    return(out)
}

.orgFullMembered <- function(graph, rest_genes, genomes, gene_list, n_len, full){
    rest_ego <- ego(graph = graph, order = 1, nodes = rest_genes)
    rest_ego <- lapply(rest_ego, names)
    n_rest_ego <- sapply(rest_ego, length)
    ego_id <- lapply(seq_along(n_rest_ego), function(i){
        return(rep(i, n_rest_ego[[i]]))
    })
    rest_ego <- data.frame(gene_id = unlist(rest_ego), ego_id = unlist(ego_id))
    hit <- match(rest_ego$gene_id, gene_list$gene)
    rest_ego$genome <- gene_list$genome[hit]
    order_rest_ego <- order(rest_ego$ego_id, rest_ego$genome)
    rest_ego <- rest_ego[order_rest_ego, ]
    n_rest_ego_nodes <- table(rest_ego$ego_id)
    full_mem_candidate <- as.numeric(n_rest_ego_nodes) == n_len
    if(sum(full_mem_candidate) == 0){ return(list(rest_genes = rest_genes)) }
    full_mem_candidate <- rest_ego[rest_ego$ego_id %in% as.numeric(names(n_rest_ego_nodes[full_mem_candidate])), ]
    
    full_ego_id <- rep(unique(full_mem_candidate$ego_id), each = length(genomes))
    full_df <- data.frame(id = paste(full_ego_id, genomes, sep = "_"),
                          full_ego_id = full_ego_id)
    full_mem_candidate$id <- paste(full_mem_candidate$ego_id, 
                                   full_mem_candidate$genome, sep = "_")
    full_mem_candidate <- left_join(x = full_df, y = full_mem_candidate, by = "id")
    
    n_entry_full_mem_candidate <- table(full_mem_candidate$full_ego_id)
    non_multientries <- n_entry_full_mem_candidate == length(genomes)
    if(sum(non_multientries) == 0){ return(list(rest_genes = rest_genes)) }
    valid_ego_id <- names(n_entry_full_mem_candidate)[non_multientries]
    full_mem_candidate <- subset(full_mem_candidate, 
                                 subset = full_ego_id %in% valid_ego_id)
    
    is_na_full_mem_candidate <- tapply(full_mem_candidate$ego_id, 
                                       full_mem_candidate$full_ego_id,
                                       is.na)
    n_entry_full_mem_candidate <- sapply(is_na_full_mem_candidate, sum)
    n_entry_full_mem_candidate <- length(genomes) - n_entry_full_mem_candidate
    full_mem <- n_entry_full_mem_candidate == n_len
    if(sum(full_mem) == 0){ return(list(rest_genes = rest_genes)) }
    valid_ego_id <- names(n_entry_full_mem_candidate[full_mem])
    full_mem <- full_mem_candidate[full_mem_candidate$full_ego_id %in% as.numeric(valid_ego_id), ]
    full_mem <- matrix(data = full_mem$gene_id, 
                       ncol = length(genomes),
                       byrow = TRUE)
    full_mem <- unique(full_mem)
    full_mem <- as.data.frame(full_mem)
    names(full_mem) <- genomes
    
    if(full){
        full_mem$class <- "full_membered"
        
    } else {
        full_mem$class <- "partial_membered"
    }
    
    rest_genes <- rest_genes[!rest_genes %in% unlist(full_mem)]
    if(length(rest_genes) == 0){
        rest_genes <- NULL
    }
    out <- list(rest_genes = rest_genes, full_mem = full_mem)
    return(out)
}

.orgFullChained <- function(graph, rest_genes, genomes, gene_list, n_len, full){
    rest_ego <- ego(graph = graph, order = length(genomes), nodes = rest_genes)
    rest_ego <- lapply(rest_ego, names)
    n_rest_ego <- sapply(rest_ego, length)
    ego_id <- lapply(seq_along(n_rest_ego), function(i){
        return(rep(i, n_rest_ego[[i]]))
    })
    rest_ego <- data.frame(gene_id = unlist(rest_ego), ego_id = unlist(ego_id))
    hit <- match(rest_ego$gene_id, gene_list$gene)
    rest_ego$genome <- gene_list$genome[hit]
    order_rest_ego <- order(rest_ego$ego_id, rest_ego$genome)
    rest_ego <- rest_ego[order_rest_ego, ]
    n_rest_ego_nodes <- table(rest_ego$ego_id)
    full_chain_candidate <- as.numeric(n_rest_ego_nodes) == n_len
    if(sum(full_chain_candidate) == 0){ return(list(rest_genes = rest_genes)) }
    full_chain_candidate <- rest_ego[rest_ego$ego_id %in% as.numeric(names(n_rest_ego_nodes[full_chain_candidate])), ]
    
    full_ego_id <- rep(unique(full_chain_candidate$ego_id), each = length(genomes))
    full_df <- data.frame(id = paste(full_ego_id, genomes, sep = "_"),
                          full_ego_id = full_ego_id)
    full_chain_candidate$id <- paste(full_chain_candidate$ego_id, 
                                     full_chain_candidate$genome, sep = "_")
    full_chain_candidate <- left_join(x = full_df, y = full_chain_candidate, by = "id")
    
    n_entry_full_chain_candidate <- table(full_chain_candidate$full_ego_id)
    non_multientries <- n_entry_full_chain_candidate == length(genomes)
    if(sum(non_multientries) == 0){ return(list(rest_genes = rest_genes)) }
    valid_ego_id <- names(n_entry_full_chain_candidate)[non_multientries]
    full_chain_candidate <- subset(full_chain_candidate, 
                                   subset = full_ego_id %in% valid_ego_id)
    
    is_na_full_chain_candidate <- tapply(full_chain_candidate$ego_id, 
                                         full_chain_candidate$full_ego_id,
                                         is.na)
    n_entry_full_chain_candidate <- sapply(is_na_full_chain_candidate, sum)
    n_entry_full_chain_candidate <- length(genomes) - n_entry_full_chain_candidate
    full_chain <- n_entry_full_chain_candidate == n_len
    if(sum(full_chain) == 0){ return(list(rest_genes = rest_genes)) }
    valid_ego_id <- names(n_entry_full_chain_candidate[full_chain])
    
    full_chain <- full_chain_candidate[full_chain_candidate$full_ego_id %in% as.numeric(valid_ego_id), ]
    full_chain <- matrix(data = full_chain$gene_id, 
                         ncol = length(genomes),
                         byrow = TRUE)
    full_chain <- unique(full_chain)
    full_chain <- as.data.frame(full_chain)
    names(full_chain) <- genomes
    
    if(full){
        full_chain$class <- "full_chained"
        
    } else {
        full_chain$class <- "partial_chained"
    }
    
    rest_genes <- rest_genes[!rest_genes %in% unlist(full_chain)]
    if(length(rest_genes) == 0){
        rest_genes <- NULL
    }
    out <- list(rest_genes = rest_genes, full_chain = full_chain)
    return(out)
}

.setTxID <- function(x = x, gene_list = gene_list){
    tx_out <- subset(x, select = -class)
    names(tx_out) <- paste(names(tx_out), "TX", sep = "_")
    for(i in seq_len(ncol(tx_out))){
        hit <- match(tx_out[, i], gene_list$gene)
        tx_out[, i] <- gene_list$tx_id[hit]
    }
    x <- cbind(x, tx_out)
}