#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#' @importFrom parallel mclapply
syntenicOrtho <- function(object){
    # Check if the input object is of class "OrthoPairDB"
    stopifnot(inherits(x = object, "OrthoPairDB"))
    
    # Open the HDF5 file
    h5 <- H5Fopen(object$h5)
    pair_id <- h5$pair_id
    # Ensure the HDF5 file is closed when the function exits
    on.exit(H5Fclose(h5))
    
    if(!H5Lexists(h5, "blast/rbh")){
        stop("Run rbh() to obtain RBH info.")
    }
    H5Fclose(h5)
    on.exit(NULL)
    g2g_graph <- .linkGene2Genome(object = object)
    rbh <- h5read(object$h5, "blast/rbh")
    rbh <- as.data.frame(rbh, stringsAsFactors = FALSE)
    rbh <- .normalize_rbh_id_cols(rbh)
    anchor <- .findAnchors(rbh = rbh, g2g_graph = g2g_graph)
    
    if(!all(is.na(anchor))){
        t2a_graph <- .link2Anchor(anchor = anchor, g2g_graph = g2g_graph)
        orthopair <- .findSyntenicOrtho(rbh = rbh,
                                        anchor = anchor,
                                        g2g_graph = g2g_graph, 
                                        t2a_graph = t2a_graph)
        orthopair <- .findSyntenyBlocks(orthopair = orthopair)
        orthopair <- .pickBestPair(orthopair = orthopair)
        orthopair <- .classifyOrthoPair(orthopair = orthopair)
        orthopair <- .filterOrthopair(orthopair = orthopair, 
                                      g2g_graph = g2g_graph)
        orthopair <- .classifyOrthoPair(orthopair = orthopair)
        orthopair <- .reformatOrthoPair(orthopair = orthopair, 
                                        g2g_graph = g2g_graph)
        
        collapsed_id <- .collapseOverlappingGene(orthopair = orthopair, 
                                                 g2g_graph = g2g_graph)
        orthopair$genome1_collapse <- collapsed_id$genome1_collapse
        orthopair$genome2_collapse <- collapsed_id$genome2_collapse
        
    } else {
        orthopair <- NA
    }
    
    orphan <- .getOrphan(orthopair = orthopair, g2g_graph = g2g_graph)
    .h5creategroup(file = object$h5, name = "orthopair")
    .h5creategroup(file = object$h5, name = "orphan")
    .h5overwrite(obj = orthopair, file = object$h5, paste0("orthopair/", pair_id))
    genome1_id <- sub("_.*", "", pair_id)
    genome2_id <- sub(".*_", "", pair_id)
    .h5overwrite(obj = orphan$genome1, file = object$h5, paste0("orphan/", genome1_id))
    .h5overwrite(obj = orphan$genome2, file = object$h5, paste0("orphan/", genome2_id))
    
    .h5overwrite(obj = "orthopair",
                 file = object$h5,
                 name = "data_type")
    .h5overwrite(obj = as.character(Sys.time()), file = object$h5, "timestamp/pairing")
}

.collapseOverlappingGene <- function(orthopair, g2g_graph){
    id_map_Mto1 <- .collapseMto1(orthopair = orthopair, g2g_graph = g2g_graph)
    id_map_1toM <- .collapse1toM(orthopair = orthopair, g2g_graph = g2g_graph)
    id_map_MtoM <- .collapseMtoM(orthopair = orthopair, g2g_graph = g2g_graph)
    id_map <- rbind(id_map_Mto1, id_map_1toM, id_map_MtoM)
    hit <- match(orthopair$genome1_gene, id_map$old)
    genome1_collapse <- id_map$new[hit]
    if(length(genome1_collapse) == 0){
        genome1_collapse <- 0
    }
    hit <- match(orthopair$genome2_gene, id_map$old)
    genome2_collapse <- id_map$new[hit]
    if(length(genome2_collapse) == 0){
        genome2_collapse <- 0
    }
    out <- list(genome1_collapse = genome1_collapse, genome2_collapse = genome2_collapse)
    return(out)
}

.collapse1toM <- function(orthopair, g2g_graph){
    id_map <- NULL
    orthopair_1toM <- orthopair$genome2_gene[orthopair$class == "1toM"]
    gene_1toM <- g2g_graph$genome2_gff$gene_id %in% orthopair_1toM
    gene_1toM <- g2g_graph$genome2_gff[gene_1toM, ]
    gene_1toM <- gene_1toM[gene_1toM$type == "CDS", ]
    gene_1toM <- GRanges(seqnames = gene_1toM$seqnames, 
                         ranges = IRanges(start = gene_1toM$start,
                                          end = gene_1toM$end),
                         gene_id = gene_1toM$gene_id)
    gene_1toM_ol <- as.data.frame(findOverlaps(gene_1toM, gene_1toM))
    gene_1toM_ol$queryHits <- gene_1toM$gene_id[gene_1toM_ol$queryHits]
    gene_1toM_ol$subjectHits <- gene_1toM$gene_id[gene_1toM_ol$subjectHits]
    gene_1toM_ol <- subset(gene_1toM_ol,
                           subset = queryHits != subjectHits)
    if(nrow(gene_1toM_ol) > 0){
        gene_1toM_ol <- unique(gene_1toM_ol)
        g <- graph_from_data_frame(d = gene_1toM_ol, directed = FALSE)
        grp <- split(V(g)$name, components(g)$membership)
        id_map <- lapply(grp, function(x){
            return(data.frame(old = x, 
                              new = paste(sort(x), 
                                          collapse = "/")))
        })
        id_map <- do.call("rbind", id_map)
    }
    return(id_map)
}

.collapseMto1 <- function(orthopair, g2g_graph){
    id_map <- NULL
    orthopair_Mto1 <- orthopair$genome1_gene[orthopair$class == "Mto1"]
    gene_Mto1 <- g2g_graph$genome1_gff$gene_id %in% orthopair_Mto1
    gene_Mto1 <- g2g_graph$genome1_gff[gene_Mto1, ]
    gene_Mto1 <- gene_Mto1[gene_Mto1$type == "CDS", ]
    gene_Mto1 <- GRanges(seqnames = gene_Mto1$seqnames, 
                         ranges = IRanges(start = gene_Mto1$start,
                                          end = gene_Mto1$end),
                         gene_id = gene_Mto1$gene_id)
    gene_Mto1_ol <- as.data.frame(findOverlaps(gene_Mto1, gene_Mto1))
    gene_Mto1_ol$queryHits <- gene_Mto1$gene_id[gene_Mto1_ol$queryHits]
    gene_Mto1_ol$subjectHits <- gene_Mto1$gene_id[gene_Mto1_ol$subjectHits]
    gene_Mto1_ol <- subset(gene_Mto1_ol,
                           subset = queryHits != subjectHits)
    if(nrow(gene_Mto1_ol) > 0){
        gene_Mto1_ol <- unique(gene_Mto1_ol)
        g <- graph_from_data_frame(d = gene_Mto1_ol, directed = FALSE)
        grp <- split(V(g)$name, components(g)$membership)
        id_map <- lapply(grp, function(x){
            return(data.frame(old = x, 
                              new = paste(sort(x), 
                                          collapse = "/")))
        })
        id_map <- do.call("rbind", id_map)
    }
    return(id_map)
}

.collapseMtoM <- function(orthopair, g2g_graph){
    id_map <- NULL
    orthopair_MtoM <- orthopair[orthopair$class == "MtoM", ]
    if(nrow(orthopair_MtoM) == 0){
        return(id_map)
    }
    orthopair_MtoM$class <- "Mto1"
    MtoM_id_map_Mto1 <- .collapseMto1(orthopair = orthopair_MtoM, 
                                      g2g_graph = g2g_graph)
    orthopair_MtoM$class <- "1toM"
    MtoM_id_map_1toM <- .collapse1toM(orthopair = orthopair_MtoM,
                                      g2g_graph = g2g_graph)
    id_map <- rbind(id_map, MtoM_id_map_Mto1, MtoM_id_map_1toM)
    return(id_map)
}

.getOrphan <- function(orthopair, g2g_graph){
    ## Transcript-level genome1_df / genome2_df; orthopair$genome1_tx is tx_index
    genome1_tx_hit <- g2g_graph$genome1_df$tx_index %in% orthopair$genome1_tx
    hit_gene_id <- g2g_graph$genome1_df$gene_id[genome1_tx_hit]
    genome1_gene_hit <- g2g_graph$genome1_df$gene_id %in% hit_gene_id
    genome1_orphan <- g2g_graph$genome1_df$gene_id[!genome1_gene_hit]
    
    genome2_tx_hit <- g2g_graph$genome2_df$tx_index %in% orthopair$genome2_tx
    hit_gene_id <- g2g_graph$genome2_df$gene_id[genome2_tx_hit]
    genome2_gene_hit <- g2g_graph$genome2_df$gene_id %in% hit_gene_id
    genome2_orphan <- g2g_graph$genome2_df$gene_id[!genome2_gene_hit]
    out <- list(genome1 = unique(genome1_orphan),
                genome2 = unique(genome2_orphan))
}

.prepGenomeGraph <- function(h5){
    # Prepare genome LCB nodes
    block_pairs <- list(lcb_1to1 = .order(df = h5$sibeliaz$lcb_pairs$lcb_1to1),
                        lcb_non_1to1 = .order(df = h5$sibeliaz$lcb_pairs$lcb_non_1to1))
    genome1_blocks <- rbind(subset(block_pairs$lcb_1to1,
                                   select = query_chr:query_end),
                            subset(block_pairs$lcb_non_1to1,
                                   select = query_chr:query_end))
    genome2_blocks <- rbind(subset(block_pairs$lcb_1to1,
                                   select = subject_chr:subject_end),
                            subset(block_pairs$lcb_non_1to1,
                                   select = subject_chr:subject_end))
    genome1_blocks$node <- apply(genome1_blocks, 1, paste, collapse = "_")
    genome1_blocks$node <- paste0("q_", genome1_blocks$node)
    genome2_blocks$node <- apply(genome2_blocks, 1, paste, collapse = "_")
    genome2_blocks$node <- paste0("s_", genome2_blocks$node)
    
    # Graph for LCBs
    genome_edge <- data.frame(genome1_genome = genome1_blocks$node,
                              genome2_genome = genome2_blocks$node)
    genome1_blocks <- unique(genome1_blocks)
    genome2_blocks <- unique(genome2_blocks)
    genome1_blocks$node_id <- seq_along(genome1_blocks$node)
    genome2_blocks$node_id <- seq_along(genome2_blocks$node)
    
    hit <- match(genome_edge$genome1_genome, genome1_blocks$node)
    genome_edge$genome1_genome <- genome1_blocks$node_id[hit]
    hit <- match(genome_edge$genome2_genome, genome2_blocks$node)
    genome_edge$genome2_genome <- genome2_blocks$node_id[hit]
    
    out <- list(genome1_blocks = genome1_blocks,
                genome2_blocks = genome2_blocks,
                genome_edge = genome_edge)
    return(out)
}

.order <- function(df, by = "genome1"){
    if(by == "genome1"){
        # Order by query chromosome and start position
        q_order <- order(df$query_chr, df$query_start)
        df <- df[q_order, ]
    } else {
        # Order by subject chromosome and start position
        s_order <- order(df$subject_chr, df$subject_start)
        df <- df[s_order, ]
    }
    return(df)
}

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
.makeGRanges <- function(df, genome){
    if(genome == "genome1"){
        # Handle cases where start position is greater than end position
        minus <- df$query_start > df$query_end
        tmp <- df$query_start[minus]
        df$query_start[minus] <- df$query_end[minus]
        df$query_end[minus] <- tmp
        
        # Create GRanges object for query genome
        gr <- GRanges(seqnames = df$query_chr,
                      ranges = IRanges(start = df$query_start,
                                       end = df$query_end),
                      strand = "*")
        
    } else {
        # Handle cases where start position is greater than end position
        minus <- df$subject_start > df$subject_end
        tmp <- df$subject_start[minus]
        df$subject_start[minus] <- df$subject_end[minus]
        df$subject_end[minus] <- tmp
        
        # Create GRanges object for subject genome
        gr <- GRanges(seqnames = df$subject_chr,
                      ranges = IRanges(start = df$subject_start,
                                       end = df$subject_end),
                      strand = "*")
    }
    return(gr)
}

#' @importFrom rtracklayer import.gff
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

.checkGFFentiry <- function(gff){
    entry_type <-  c("gene", "transcript", "mRNA", "five_prime_UTR", "exon",
                     "CDS", "three_prime_UTR")
    gff <- gff[gff$type %in% entry_type]
    if(!is.null(gff$gene_id)){
        if(all(!is.na(gff$gene_id))){
            return(gff)
        }
    }
    gff <- .setGeneID(gff = gff)
    return(gff)
}

.orderGFF <- function(gff){
    # Order GFF data by chromosome and start position
    gff_order <- order(as.character(seqnames(gff)), start(gff))
    gff <- gff[gff_order]
    return(gff)
}

.mRNA2transcript <- function(gff){
    type <- as.character(gff$type)
    type[type == "mRNA"] <- "transcript"
    gff$type <- factor(type)
    return(gff)
}

## Map legacy BLAST-style RBH column names to genome1_/genome2_ (idempotent).
.normalize_rbh_id_cols <- function(rbh) {
    if (is.null(rbh) || !length(rbh)) return(rbh)
    if (!is.data.frame(rbh)) rbh <- as.data.frame(rbh, stringsAsFactors = FALSE)
    nm <- names(rbh)
    map <- c(qseqid = "genome1_tx", sseqid = "genome2_tx",
             qgeneid = "genome1_gene", sgeneid = "genome2_gene")
    for (i in seq_along(map)) {
        old <- names(map)[i]
        new <- unname(map[i])
        j <- match(old, nm)
        if (!is.na(j) && !(new %in% nm)) nm[j] <- new
    }
    names(rbh) <- nm
    rbh
}

#' @importFrom dplyr left_join
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

.getSyntenicAnchors <- function(anchor){
    anchor <- anchor[order(anchor$leaf_anchor), ]
    leaf_anchor_neighbors <- data.frame(leaf = anchor$leaf_anchor,
                                        chr = anchor$leaf_anchor_chr,
                                        leaf_neighbor = c(NA,
                                                          head(anchor$leaf_anchor, -1),
                                                          c(tail(anchor$leaf_anchor, -1),
                                                            NA)),
                                        leaf_neighbor_chr = c(NA,
                                                              head(anchor$leaf_anchor_chr, -1),
                                                              c(tail(anchor$leaf_anchor_chr, -1),
                                                                NA)))
    leaf_anchor_neighbors <- subset(leaf_anchor_neighbors,
                                    subset = chr == leaf_neighbor_chr & !is.na(leaf_neighbor),
                                    select = -c(chr, leaf_neighbor_chr))
    
    anchor <- anchor[order(anchor$root_anchor), ]
    root_anchor_neighbors <- data.frame(root = anchor$root_anchor,
                                        chr = anchor$root_anchor_chr,
                                        root_neighbor = c(NA,
                                                          head(anchor$root_anchor, -1),
                                                          c(tail(anchor$root_anchor, -1),
                                                            NA)),
                                        root_neighbor_chr = c(NA,
                                                              head(anchor$root_anchor_chr, -1),
                                                              c(tail(anchor$root_anchor_chr, -1),
                                                                NA)))
    root_anchor_neighbors <- subset(root_anchor_neighbors,
                                    subset = chr == root_neighbor_chr & !is.na(root_neighbor),
                                    select = -c(chr, root_neighbor_chr))
    
    anchor <- left_join(anchor, root_anchor_neighbors, "root")
    anchor <- left_join(anchor, 
                        data.frame(root_neighbor = anchor$root_anchor, 
                                   leaf = anchor$leaf_anchor),
                        "root_neighbor")
    anchor <- left_join(anchor, leaf_anchor_neighbors, "leaf")
    anchor <- subset(anchor,
                     subset = leaf_anchor == leaf_neighbor,
                     select = root_anchor:pair_id)
    anchor <- unique(anchor)
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

.calcThreshold <- function(x, pair_id, top = FALSE){
    if(top){
        x <- unlist(tapply(x, pair_id, max))
        q <- quantile(x, c(0.25, 0.75))
        whisker <- 1.5 * diff(q)
        threshold <- q[2] + whisker
        
    } else {
        x <- unlist(tapply(x, pair_id, max))
        q <- quantile(x, c(0.25, 0.75))
        whisker <- 1.5 * diff(q)
        threshold <- q[1] - whisker
    }
    return(threshold)
}

#' @importFrom dplyr left_join
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

.evalSynteny <- function(orthopair, anchor){
    orthopair <- orthopair[order(orthopair$genome1_synteny_block), ]
    q_block <- subset(orthopair,
                      select = c(genome1_gene, genome1_synteny_block))
    q_block <- unique(q_block)
    n_q_block <- table(q_block$genome1_synteny_block)
    q_singleton_block <- as.numeric(names(n_q_block[n_q_block == 1]))
    q_singleton_anchor <- orthopair$genome1_synteny_block %in% q_singleton_block
    q_singleton_anchor <- orthopair$root_anchor[q_singleton_anchor]
    
    s_block <- subset(orthopair,
                      select = c(genome2_gene, genome2_synteny_block))
    s_block <- unique(s_block)
    n_s_block <- table(s_block$genome2_synteny_block)
    s_singleton_block <- as.numeric(names(n_s_block[n_s_block == 1]))
    s_singleton_anchor <- orthopair$genome2_synteny_block %in% s_singleton_block
    s_singleton_anchor <- orthopair$leaf_anchor[s_singleton_anchor]
    
    block_pair <- paste(orthopair$genome1_synteny_block,
                        orthopair$genome2_synteny_block, 
                        sep = "_")
    n_block_pair <- table(block_pair)
    singleton_block_pair <- names(n_block_pair[n_block_pair == 1])
    singleton_block_pair <- subset(orthopair, 
                                   subset = block_pair %in% singleton_block_pair,
                                   select = c(root, leaf))
    
    anchor <- anchor[!anchor$root_anchor %in% q_singleton_anchor, ]
    anchor <- anchor[!anchor$leaf_anchor %in% s_singleton_anchor, ]
    q_hit <- match(anchor$root_anchor, singleton_block_pair$root)
    s_hit <- match(anchor$leaf_anchor, singleton_block_pair$leaf)
    singleton_block_pair_hit <- which(q_hit == s_hit)
    anchor <- anchor[-singleton_block_pair_hit, ]
    
    return(anchor)
}

.pickNonSyntenicOrtho <- function(h5, orthopair, g2g_graph){
    rbbh <- as.data.frame(h5$blast$rbbh, stringsAsFactors = FALSE)
    rbbh <- .normalize_rbh_id_cols(rbbh)
    hit <- match(rbbh$genome1_tx, g2g_graph$genome1_df$tx_index)
    rbbh$genome1_gene <- g2g_graph$genome1_df$gene_id[hit]
    hit <- match(rbbh$genome2_tx, g2g_graph$genome2_df$tx_index)
    rbbh$genome2_gene <- g2g_graph$genome2_df$gene_id[hit]
    rbbh$pair_id <- paste(rbbh$genome1_gene, rbbh$genome2_gene, sep = "_")
    orthopair$rbbh <- FALSE
    orthopair$rbbh[orthopair$pair_id %in% rbbh$pair_id] <- TRUE
    orthopair$syntenic <- TRUE
    query_orphan <- !rbbh$genome1_gene %in% orthopair$genome1_gene
    subject_orphan <- !rbbh$genome2_gene %in% orthopair$genome2_gene
    rbbh <- subset(rbbh, subset = query_orphan & subject_orphan)
    rbbh_score <- .getRBHscore(rbh = rbbh)
    rbbh <- subset(rbbh, select = c(genome1_gene, genome1_tx, genome2_gene, genome2_tx, pair_id))
    names(rbbh) <- c("genome1_gene", "genome1_tx", "genome2_gene", "genome2_tx", "pair_id")
    rbbh$mutual_ci <- rbbh_score$mutual_ci
    rbbh$rbbh <- TRUE
    rbbh$syntenic <- FALSE
    orthopair <- rbind(orthopair, rbbh)
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

.splitGene <- function(orthopair, g2g_graph){
    split_1toM <- .split1toM(orthopair, gff = g2g_graph$genome1_df)
    split_Mto1 <- .splitMto1(orthopair, gff = g2g_graph$genome2_df)
    split_MtoM <- .splitMtoM(orthopair, g2g_graph = g2g_graph)
    
    out <- list(genome1 = c(split_1toM, split_MtoM$genome1),
                genome2 = c(split_Mto1, split_MtoM$genome2))
    return(out)
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

.examine1to1 <- function(orthopair, valid){
    sog_1to1 <- orthopair$class == "1to1"
    threshold <-  quantile(orthopair$mutual_ci[sog_1to1], 0.05)
    update_valid <- orthopair$mutual_ci[sog_1to1] >= threshold
    valid[sog_1to1] <- valid[sog_1to1] & update_valid
    return(valid)
}

.examine1toM <- function(orthopair, valid = NULL){
    if (is.null(valid)) valid <- rep(TRUE, nrow(orthopair))
    sog_1toM <- orthopair$class == "1toM"
    mci_diff <- tapply(seq_along(orthopair$SOG[sog_1toM]),
                       orthopair$SOG[sog_1toM], 
                       function(i){
                           i_og <- orthopair[sog_1toM, ][i, ]
                           i_max <- max(orthopair$mutual_ci[sog_1toM][i])
                           out_i <- i_max - orthopair$mutual_ci[sog_1toM][i]
                           out_i <- data.frame(genome2_gene = orthopair$genome2_gene[sog_1toM][i],
                                               mci_diff = out_i)
                           return(out_i)
                       })
    mci_diff <- do.call("rbind", mci_diff)
    threshold <-  quantile(mci_diff$mci_diff[mci_diff$mci_diff > 0], 0.95)
    update_valid <- mci_diff$mci_diff <= threshold
    update_valid <- orthopair$genome2_gene[sog_1toM] %in% mci_diff$genome2_gene[update_valid]
    valid[sog_1toM] <- valid[sog_1toM] & update_valid
    return(valid)
}

.examineMto1 <- function(orthopair, valid){
    sog_Mto1 <- orthopair$class == "Mto1"
    mci_diff <- tapply(seq_along(orthopair$SOG[sog_Mto1]),
                       orthopair$SOG[sog_Mto1], 
                       function(i){
                           i_max <- max(orthopair$mutual_ci[sog_Mto1][i])
                           out_i <- i_max - orthopair$mutual_ci[sog_Mto1][i]
                           out_i <- data.frame(genome1_gene = orthopair$genome1_gene[sog_Mto1][i],
                                               mci_diff = out_i)
                           return(out_i)
                       })
    mci_diff <- do.call("rbind", mci_diff)
    threshold <-  quantile(mci_diff$mci_diff[mci_diff$mci_diff > 0], 0.95)
    update_valid <- mci_diff$mci_diff <= threshold
    update_valid <- orthopair$genome1_gene[sog_Mto1] %in% mci_diff$genome1_gene[update_valid]
    valid[sog_Mto1] <- valid[sog_Mto1] & update_valid
    return(valid)
}

.examineMtoM <- function(orthopair){
    sog_MtoM <- orthopair$class == "MtoM"
    orthopair <- orthopair[order(orthopair$pident, decreasing = TRUE), ]
    orthopair$index <- seq_along(orthopair$genome1_gene)
    valid_pair_id <- tapply(seq_along(orthopair$SOG[sog_MtoM]),
                            orthopair$SOG[sog_MtoM], 
                            function(i){
                                i_og <- orthopair[sog_MtoM, ][i, ]
                                i_og <- i_og[order(i_og$ci_q2s, decreasing = TRUE), ]
                                q_best <- tapply(i_og$index,
                                                 i_og$genome1_gene, 
                                                 "[", 1)
                                i_og <- i_og[order(i_og$ci_s2q, decreasing = TRUE), ]
                                s_best <- tapply(i_og$index,
                                                 i_og$genome2_gene, 
                                                 "[", 1)
                                best <- unique(c(q_best, s_best))
                                i_og <- i_og[i_og$index %in% best, ]
                                i_anchor <- intersect(q_best, s_best)
                                rest_og <- best[!best %in% i_anchor]
                                i_anchor_og <- i_og[i_og$index %in% i_anchor, ]
                                i_og <- i_og[i_og$index %in% rest_og, ]
                                
                                q_to_anchor <- i_og$genome1_gene %in% i_anchor_og$genome1_gene
                                s_to_anchor <- i_og$genome2_gene %in% i_anchor_og$genome2_gene
                                i_to_anchor_og <- i_og[q_to_anchor | s_to_anchor, ]
                                rest_og <- rest_og[!rest_og %in% i_to_anchor_og$index]
                                i_og <- i_og[i_og$index %in% rest_og, ]
                                
                                q_to_anchored <- i_og$genome1_gene %in% i_to_anchor_og$genome1_gene
                                s_to_anchored <- i_og$genome2_gene %in% i_to_anchor_og$genome2_gene
                                
                                i_to_nonanchored <- i_og[!(q_to_anchored | s_to_anchored), ]
                                
                                q_id <- NULL
                                s_id <- NULL
                                out_i <- NULL
                                while(TRUE){
                                    q_dup <- i_og$genome1_gene %in% q_id
                                    s_dup <- i_og$genome2_gene %in% s_id
                                    target_pair <- !(q_dup & s_dup)
                                    if(all(!target_pair)){
                                        break
                                    }
                                    out_i <- c(out_i, i_og$pair_id[target_pair][1])
                                    q_id <- c(q_id, i_og$genome1_gene[target_pair][1])
                                    s_id <- c(s_id, i_og$genome2_gene[target_pair][1])
                                }
                                return(out_i)
                            })
    valid_pair_id <- do.call("c", valid_pair_id)
    pair_id_MtoM <- orthopair$pair_id[sog_MtoM]
    out <- pair_id_MtoM[!pair_id_MtoM %in% valid_pair_id]
    return(out)
}
# .examineMtoM <- function(orthopair){
#     sog_MtoM <- orthopair$class == "MtoM"
#     orthopair$order <- abs(orthopair$qgeneid - orthopair$sgeneid)
#     valid_pair_id <- tapply(seq_along(orthopair$SOG[sog_MtoM]),
#                             orthopair$SOG[sog_MtoM], 
#                             function(i){
#                                 i_og <- orthopair[sog_MtoM, ][i, ]
#                                 q_id <- NULL
#                                 s_id <- NULL
#                                 out_i <- NULL
#                                 while(TRUE){
#                                     q_dup <- i_og$qgeneid %in% q_id
#                                     s_dup <- i_og$sgeneid %in% s_id
#                                     target_pair <- !(q_dup | s_dup)
#                                     if(all(!target_pair)){
#                                         break
#                                     }
#                                     target_mci <- i_og$mutual_ci[target_pair][1]
#                                     target_pair[target_pair] <- i_og$mutual_ci[target_pair] == target_mci
#                                     if(sum(target_pair) > 1){
#                                         pair_order <- order(i_og$order[target_pair])
#                                         target_pair[target_pair][pair_order != 1] <- FALSE
#                                     }
#                                     out_i <- c(out_i, i_og$pair_id[target_pair])
#                                     q_id <- c(q_id, i_og$qgeneid[target_pair])
#                                     s_id <- c(s_id, i_og$sgeneid[target_pair])
#                                 }
#                                 return(out_i)
#                             })
#     valid_pair_id <- do.call("c", valid_pair_id)
#     pair_id_MtoM <- orthopair$pair_id[sog_MtoM]
#     out <- pair_id_MtoM[!pair_id_MtoM %in% valid_pair_id]
#     return(out)
# }
.pickBestPair <- function(orthopair){
    orthopair <- orthopair[order(orthopair$mutual_ci, decreasing = TRUE), ]
    orthopair <- orthopair[!duplicated(orthopair$pair_id), ]
    return(orthopair)
}
