#' Filter ortholog pairs based on Locally Collinear Blocks (LCBs)
#'
#' This function filters ortholog pairs based on their presence in Locally Collinear Blocks (LCBs).
#'
#' @param object A OrthoPairDB object.
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#' @importFrom parallel mclapply
syntenicOrtho <- function(object){
    # Check if the input object is of class "OrthoPairDB"
    stopifnot(inherits(x = object, "OrthoPairDB"))
    
    # Open the HDF5 file
    h5 <- H5Fopen(object$h5)
    # Ensure the HDF5 file is closed when the function exits
    on.exit(H5Fclose(h5))
    
    if(!H5Lexists(h5, "blast/rbh")){
        stop("Run rbh() to obtain RBH info.")
    }
    H5Fclose(h5)
    on.exit(NULL)
    g2g_graph <- .linkGene2Genome(object = object)
    rbh <- h5read(object$h5, "blast/rbh")
    anchor <- .findAnchors(rbh = rbh, g2g_graph = g2g_graph)
    
    if(!all(is.na(anchor))){
        t2a_graph <- .link2Anchor(anchor = anchor, g2g_graph = g2g_graph)
        orthopair <- .findSyntenicOrtho(rbh = rbh,
                                        anchor = anchor,
                                        g2g_graph = g2g_graph, 
                                        t2a_graph = t2a_graph)
        
        # orthopair <- .sortSyntenicOrtho(orthopair = orthopair, g2g_graph = g2g_graph)
        # .h5overwrite(obj = orthopair, file = object$h5, "orthopair_tx")
        orthopair <- .findSyntenyBlocks(orthopair = orthopair)
        # anchor2 <- .evalSynteny(orthopair = orthopair, anchor = anchor)
        # t2a_graph <- .link2Anchor(anchor = anchor, g2g_graph = g2g_graph)
        # orthopair <- .findSyntenicOrtho(rbh = rbh,
        #                                 anchor = anchor2,
        #                                 g2g_graph = g2g_graph,
        #                                 t2a_graph = t2a_graph)
        # orthopair <- .findSyntenyBlocks(orthopair = orthopair)
        
        orthopair <- .pickBestPair(orthopair = orthopair)
        # Check gene order
        # orthopair <- .filterOrthopair(orthopair = orthopair)
        # te <- .labelTE(orthopair = orthopair, 
        #                h5 = h5, 
        #                t2a_graph = t2a_graph)
        # orthopair$te <- orthopair$query_gene %in% te$query_te | 
        #     orthopair$query_gene %in% te$subject_te
        orthopair <- .classifyOrthoPair(orthopair = orthopair)
        # orthopair <- .filterOrthopair(orthopair = orthopair)
        # orthopair <- .examineOrthoPair(orthopair = orthopair, g2g_graph = g2g_graph)
        
        # split_orthopair <- .splitGene(orthopair = orthopair, g2g_graph = g2g_graph)
        # rest_orthopair <- split_orthopair$rest
        # orthopair <- split_orthopair$splited
        # 
        # if(!is.null(rest_orthopair)){
        #     while(TRUE){
        #         rest_orthopair <- .classifyOrthoPair(orthopair = rest_orthopair)
        #         rest_orthopair <- .splitGene(orthopair = rest_orthopair, g2g_graph = g2g_graph)
        #         if(is.null(rest_orthopair$rest)){
        #             orthopair <- rbind(orthopair,
        #                                rest_orthopair$splited)
        #             break
        #             
        #         } else {
        #             orthopair <- rbind(orthopair,
        #                                rest_orthopair$splited)
        #             rest_orthopair <- rest_orthopair$rest
        #         }
        #     }
        # }
        # orthopair <- .sortSyntenicOrtho(orthopair = orthopair, g2g_graph = g2g_graph)
        orthopair <- .filterOrthopair(orthopair = orthopair, 
                                      g2g_graph = g2g_graph)
        orthopair <- .classifyOrthoPair(orthopair = orthopair)
        orthopair <- .reformatOrthoPair(orthopair = orthopair, 
                                        g2g_graph = g2g_graph)
        
        collapsed_id <- .collapseOverlappingGene(orthopair = orthopair, 
                                                 g2g_graph = g2g_graph)
        orthopair$query_collapse <- collapsed_id$query_collapse
        orthopair$subject_collapse <- collapsed_id$subject_collapse
        
    } else {
        orthopair <- NA
    }
    
    orphan <- .getOrphan(orthopair = orthopair, g2g_graph = g2g_graph)
    .h5overwrite(obj = orthopair, file = object$h5, "orthopair_gene")
    .h5overwrite(obj = orphan$query, file = object$h5, "orphan_query")
    .h5overwrite(obj = orphan$subject, file = object$h5, "orphan_subject")
    
    .h5overwrite(obj = "orthopair",
                 file = object$h5,
                 name = "data_type")
    .h5overwrite(obj = as.character(Sys.time()), file = object$h5, "timestamp/pairing")
}

# .prepateTElist <- function(h5, g2g_graph){
#     anchor <- .findAnchors(h5 = h5,
#                            g2g_graph = g2g_graph)
#     t2a_graph <- .link2Anchor(anchor = anchor, g2g_graph = g2g_graph)
#     orthopair <- .findSyntenicOrtho(h5 = h5,
#                                     g2g_graph = g2g_graph,
#                                     t2a_graph = t2a_graph)
#     orthopair <- .pickBestPair(orthopair = orthopair)
#     orthopair <- .filterOrthopair(orthopair = orthopair)
#     out <- .labelTE(orthopair = orthopair,
#                     h5 = h5,
#                     g2g_graph = g2g_graph)
#     return(out)
# }

.collapseOverlappingGene <- function(orthopair, g2g_graph){
    id_map_Mto1 <- .collapseMto1(orthopair = orthopair, g2g_graph = g2g_graph)
    id_map_1toM <- .collapse1toM(orthopair = orthopair, g2g_graph = g2g_graph)
    id_map_MtoM <- .collapseMtoM(orthopair = orthopair, g2g_graph = g2g_graph)
    id_map <- rbind(id_map_Mto1, id_map_1toM, id_map_MtoM)
    hit <- match(orthopair$query_gene, id_map$old)
    query_collapse <- id_map$new[hit]
    if(length(query_collapse) == 0){
        query_collapse <- NA
    }
    hit <- match(orthopair$subject_gene, id_map$old)
    subject_collapse <- id_map$new[hit]
    if(length(subject_collapse) == 0){
        subject_collapse <- NA
    }
    out <- list(query_collapse = query_collapse, subject_collapse = subject_collapse)
    return(out)
}

.collapse1toM <- function(orthopair, g2g_graph){
    id_map <- NULL
    orthopair_1toM <- orthopair$subject_gene[orthopair$class == "1toM"]
    gene_1toM <- g2g_graph$subject_gff$gene_id %in% orthopair_1toM
    gene_1toM <- g2g_graph$subject_gff[gene_1toM, ]
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
    orthopair_Mto1 <- orthopair$query_gene[orthopair$class == "Mto1"]
    gene_Mto1 <- g2g_graph$query_gff$gene_id %in% orthopair_Mto1
    gene_Mto1 <- g2g_graph$query_gff[gene_Mto1, ]
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
    orthopair_MtoM$class <- "Mto1"
    MtoM_id_map_Mto1 <- .collapseMto1(orthopair = orthopair_MtoM, 
                                      g2g_graph = g2g_graph)
    orthopair_MtoM$class <- "1toM"
    MtoM_id_map_1toM <- .collapse1toM(orthopair = orthopair_MtoM,
                                      g2g_graph = g2g_graph)
    id_map <- rbind(id_map, MtoM_id_map_Mto1, MtoM_id_map_1toM)
    return(id_map)
}

.labelTE <- function(orthopair, h5, g2g_graph){
    rbh <- h5$blast$rbh
    hit <- match(rbh$qseqid, g2g_graph$query_gff$ID)
    rbh$qgeneid <- g2g_graph$query_gff$gene_id[hit]
    hit <- match(rbh$sseqid, g2g_graph$subject_gff$ID)
    rbh$sgeneid <- g2g_graph$subject_gff$gene_id[hit]
    rbh$pair_id <- paste(rbh$qgeneid, rbh$sgeneid, sep = "_")
    
    non_syntenic_hits <- {rbh$qgeneid %in% orthopair$query_gene |
            rbh$sgeneid %in% orthopair$subject_gene} &
        !rbh$pair_id %in% orthopair$pair_id
    rbh <- unique(subset(rbh[non_syntenic_hits, ],
                         select = c(qgeneid, sgeneid)))
    q_rbh <- table(rbh$qgeneid)
    s_rbh <- table(rbh$sgeneid)
    
    q_q <- quantile(q_rbh, c(0.25, 0.75))
    q_whisker <- 1.5 * diff(q_q)
    q_th <- q_q[2] + q_whisker
    s_q <- quantile(s_rbh, c(0.25, 0.75))
    s_whisker <- 1.5 * diff(s_q)
    s_th <- s_q[2] + s_whisker
    q_te <- names(q_rbh[q_rbh >= q_th])
    s_te <- names(s_rbh[s_rbh >= s_th])
    out <- list(query_te = q_te, subject_te = s_te)
    return(out)
}

.getOrphan <- function(orthopair, g2g_graph){
    query_tx_hit <- g2g_graph$query_gff$Parent %in% orthopair$query_tx
    hit_gene_id <- g2g_graph$query_gff$gene_id[query_tx_hit]
    query_gene_hit <- g2g_graph$query_gff$gene_id %in% hit_gene_id
    query_orphan <- g2g_graph$query_gff$gene_id[!query_gene_hit]
    
    subject_tx_hit <- g2g_graph$subject_gff$Parent %in% orthopair$subject_tx
    hit_gene_id <- g2g_graph$subject_gff$gene_id[subject_tx_hit]
    subject_gene_hit <- g2g_graph$subject_gff$gene_id %in% hit_gene_id
    subject_orphan <- g2g_graph$subject_gff$gene_id[!subject_gene_hit]
    out <- list(query = unique(query_orphan),
                subject = unique(subject_orphan))
}

.prepGenomeGraph <- function(h5){
    # Prepare genome LCB nodes
    block_pairs <- list(lcb_1to1 = .order(df = h5$sibeliaz$lcb_pairs$lcb_1to1),
                        lcb_non_1to1 = .order(df = h5$sibeliaz$lcb_pairs$lcb_non_1to1))
    query_blocks <- rbind(subset(block_pairs$lcb_1to1,
                                 select = query_chr:query_end),
                          subset(block_pairs$lcb_non_1to1,
                                 select = query_chr:query_end))
    subject_blocks <- rbind(subset(block_pairs$lcb_1to1,
                                   select = subject_chr:subject_end),
                            subset(block_pairs$lcb_non_1to1,
                                   select = subject_chr:subject_end))
    query_blocks$node <- apply(query_blocks, 1, paste, collapse = "_")
    query_blocks$node <- paste0("q_", query_blocks$node)
    subject_blocks$node <- apply(subject_blocks, 1, paste, collapse = "_")
    subject_blocks$node <- paste0("s_", subject_blocks$node)
    
    # Graph for LCBs
    genome_edge <- data.frame(query_genome = query_blocks$node,
                              subject_genome = subject_blocks$node)
    query_blocks <- unique(query_blocks)
    subject_blocks <- unique(subject_blocks)
    query_blocks$node_id <- seq_along(query_blocks$node)
    subject_blocks$node_id <- seq_along(subject_blocks$node)
    
    hit <- match(genome_edge$query_genome, query_blocks$node)
    genome_edge$query_genome <- query_blocks$node_id[hit]
    hit <- match(genome_edge$subject_genome, subject_blocks$node)
    genome_edge$subject_genome <- subject_blocks$node_id[hit]
    
    out <- list(query_blocks = query_blocks,
                subject_blocks = subject_blocks,
                genome_edge = genome_edge)
    return(out)
}

.order <- function(df, by = "query"){
    if(by == "query"){
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
    
    q_tx_i <- gff_ls$query_gff$type %in% c("transcript", "mRNA")
    query_df <- gff_ls$query_gff[q_tx_i, ]
    q_cds_i <- gff_ls$query_gff$type %in% "CDS"
    query_gff <- gff_ls$query_gff[q_cds_i, ]
    hit <- match(unlist((query_gff$Parent)), gff_ls$query_gff$ID)
    query_gff$tx_index <- gff_ls$query_gff$tx_index[hit]
    # query_df <- data.frame(chr = query_gff$seqnames,
    #                        tx = query_gff$ID, 
    #                        tx_id = seq_along(query_gff$ID),
    #                        index = query_gff$index,
    #                        gene = query_gff$gene_id)
    # query_df$gene_id <- as.numeric(factor(query_df$gene, 
    #                                       unique(query_df$gene)))
    
    s_tx_i <- gff_ls$subject_gff$type %in% c("transcript", "mRNA")
    subject_df <- gff_ls$subject_gff[s_tx_i, ]
    s_cds_i <- gff_ls$subject_gff$type %in% "CDS"
    subject_gff <- gff_ls$subject_gff[s_cds_i, ]
    hit <- match(unlist((subject_gff$Parent)), gff_ls$subject_gff$ID)
    subject_gff$tx_index <- gff_ls$subject_gff$tx_index[hit]
    # subject_df <- data.frame(chr = subject_gff$seqnames,
    #                          tx = subject_gff$ID, 
    #                          tx_id = seq_along(subject_gff$ID),
    #                          index = subject_gff$index,
    #                          gene = subject_gff$gene_id)
    # subject_df$gene_id <- as.numeric(factor(subject_df$gene, 
    #                                         unique(subject_df$gene)))
    out <- list(query_df = query_df,
                subject_df = subject_df,
                query_gff = query_gff,
                subject_gff = subject_gff)
    return(out)
}

#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#'
.makeGRanges <- function(df, genome){
    if(genome == "query"){
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


.getGFFlist <- function(object = NULL, gff_fn = NULL){
    # Import GFF files for query and subject genomes
    if(is.null(gff_fn)){
        files <- h5read(object$h5, "files")
        query_gff <- h5read(files$query_h5, "gff")
        subject_gff <- h5read(files$subject_h5, "gff")
        
        # Order the GFF data
        query_gff <- query_gff[order(query_gff$seqnames, query_gff$start), ]
        subject_gff <- subject_gff[order(subject_gff$seqnames, subject_gff$start), ]
        
        # Return the ordered GFF data as a list
        # out <- list(query_gff = .checkGFFentiry(query_gff),
        # subject_gff = .checkGFFentiry(subject_gff))
        out <- list(query_gff = query_gff, subject_gff = subject_gff)
        
    } else {
        out <- .importAllGFF(gff_fn)
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

#' @importFrom dplyr left_join
.findAnchors <- function(rbh, g2g_graph){
    if(length(rbh) == 1){
        return(NA)
    }
    
    rbbh <- .getRBBH(rbh = rbh)
    
    root_hit <- match(rbbh$qgeneid, g2g_graph$query_df$gene_index)
    leaf_hit <- match(rbbh$sgeneid, g2g_graph$subject_df$gene_index)
    anchor <- data.frame(root = g2g_graph$query_df$gene_index[root_hit],
                         root_anchor = g2g_graph$query_df$gene_index[root_hit],
                         root_anchor_chr = g2g_graph$query_df$seqnames[root_hit],
                         leaf_anchor = g2g_graph$subject_df$gene_index[leaf_hit],
                         leaf_anchor_chr = g2g_graph$subject_df$seqnames[leaf_hit],
                         subset(rbbh, select = c(pident, mutual_ci:pair_id)))
    anchor <- unique(anchor)
    # anchor2 <- .getSyntenicAnchors(anchor = anchor)
    
    anchor <- .checkHighCopyGenes(anchor = anchor, rbh = rbh)
    
    q <- quantile(anchor$pident, c(0.25, 0.75))
    anchor_threshold <- q[1] - 1.5 * diff(q)
    anchor <- anchor[anchor$pident >= anchor_threshold, ]
    
    return(anchor)
}

.getRBBH <- function(rbh){
    rbh$index <- seq_len(nrow(rbh))
    rbh <- rbh[order(rbh$ci_q2s, decreasing = TRUE), ]
    q_best <- rbh$index[!duplicated(rbh$qgeneid)]
    rbh <- rbh[order(rbh$ci_s2q, decreasing = TRUE), ]
    s_best <- rbh$index[!duplicated(rbh$sgeneid)]
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
    subset_rbh <- subset(rbh, select = c(qgeneid, sgeneid))
    subset_rbh <- unique(subset_rbh)
    q_rbh <- rbh[rbh$qgeneid %in% anchor$qgeneid, ]
    n_q_rbh <- table(q_rbh$qgeneid)
    q <- quantile(n_q_rbh, c(0.25, 0.75), na.rm = TRUE)
    whisker <- 1.5 * diff(q)
    q_threshold <- q[2] + whisker
    q_omit <- names(n_q_rbh)[n_q_rbh > q_threshold]
    
    s_rbh <- rbh[rbh$sgeneid %in% anchor$sgeneid, ]
    n_s_rbh <- table(s_rbh$sgeneid)
    q <- quantile(n_s_rbh, c(0.25, 0.75), na.rm = TRUE)
    whisker <- 1.5 * diff(q)
    q_threshold <- q[2] + whisker
    s_omit <- names(n_s_rbh)[n_s_rbh > q_threshold]
    
    anchor <- subset(anchor, subset = !qgeneid %in% q_omit | !sgeneid %in% s_omit)
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
    query_link2anchor <- .link2Anchor_core(g2g_graph = g2g_graph,
                                           anchor = anchor,
                                           dataset = "query")
    subject_link2anchor <- .link2Anchor_core(g2g_graph = g2g_graph,
                                             anchor = anchor,
                                             dataset = "subject")
    
    out <- list(query_link2anchor = query_link2anchor,
                subject_link2anchor = subject_link2anchor)
    return(out)
}

#' @importFrom GenomicRanges precede follow findOverlaps resize 
#' @importFrom S4Vectors queryHits subjectHits 
.link2Anchor_core <- function(g2g_graph, anchor, dataset){
    if(dataset == "query"){
        gene_df <- g2g_graph$query_df
        anchor_id <- anchor$root_anchor
        
    } else {
        gene_df <- g2g_graph$subject_df
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
    
    if(dataset == "query"){
        colnames(link2anchor) <- c("root", "root_anchor", "query_is_anchor")
        
    } else {
        colnames(link2anchor) <- c("leaf", "leaf_anchor", "subject_is_anchor")
    }
    
    return(link2anchor)
}

#' @importFrom dplyr left_join
.findSyntenicOrtho <- function(rbh, 
                               anchor,
                               g2g_graph,
                               t2a_graph,
                               rbh_threshold){
    colnames(rbh)[1:2] <- c("query_tx", "subject_tx")
    root_hit <- match(rbh$query_tx, g2g_graph$query_df$tx_index)
    leaf_hit <- match(rbh$subject_tx, g2g_graph$subject_df$tx_index)
    orthopair <- data.frame(root_tx = g2g_graph$query_df$tx_index[root_hit],
                            root = g2g_graph$query_df$gene_index[root_hit],
                            expected_leaf = g2g_graph$subject_df$gene_index[leaf_hit],
                            leaf_tx = g2g_graph$subject_df$tx_index[leaf_hit])
    orthopair <- left_join(orthopair, 
                           subset(t2a_graph$query_link2anchor, 
                                  subset = !is.na(root_anchor)), 
                           "root",
                           relationship = "many-to-many")
    orthopair <- left_join(orthopair,
                           subset(anchor,
                                  select = c(root_anchor, leaf_anchor)),
                           "root_anchor",
                           relationship = "many-to-many")
    orthopair <- left_join(orthopair, 
                           subset(t2a_graph$subject_link2anchor, 
                                  subset = !is.na(leaf_anchor)), 
                           "leaf_anchor",
                           relationship = "many-to-many")
    orthopair <- subset(orthopair, 
                        subset = leaf == expected_leaf & 
                            !is.na(root_anchor) & 
                            !is.na(leaf_anchor))
    root_tx_hit <- match(orthopair$root_tx, g2g_graph$query_df$tx_index)
    leaf_tx_hit <- match(orthopair$leaf_tx, g2g_graph$subject_df$tx_index)
    orthopair <- data.frame(query_gene = orthopair$root,
                            query_tx = orthopair$root_tx,
                            query_chr = g2g_graph$query_df$seqnames[root_tx_hit],
                            subject_gene = orthopair$leaf,
                            subject_tx = orthopair$leaf_tx,
                            subject_chr = g2g_graph$subject_df$seqnames[leaf_tx_hit],
                            orthopair)
    orthopair <- orthopair[order(orthopair$root), ]
    orthopair <- left_join(orthopair, rbh, c("query_tx", "subject_tx"))
    orthopair <- unique(orthopair)
    orthopair$is_anchor_pair <- orthopair$pair_id %in% anchor$pair_id
    return(orthopair)
}

.evalSynteny <- function(orthopair, anchor){
    orthopair <- orthopair[order(orthopair$query_synteny_block), ]
    q_block <- subset(orthopair,
                      select = c(query_gene, query_synteny_block))
    q_block <- unique(q_block)
    n_q_block <- table(q_block$query_synteny_block)
    q_singleton_block <- as.numeric(names(n_q_block[n_q_block == 1]))
    q_singleton_anchor <- orthopair$query_synteny_block %in% q_singleton_block
    q_singleton_anchor <- orthopair$root_anchor[q_singleton_anchor]
    
    s_block <- subset(orthopair,
                      select = c(subject_gene, subject_synteny_block))
    s_block <- unique(s_block)
    n_s_block <- table(s_block$subject_synteny_block)
    s_singleton_block <- as.numeric(names(n_s_block[n_s_block == 1]))
    s_singleton_anchor <- orthopair$subject_synteny_block %in% s_singleton_block
    s_singleton_anchor <- orthopair$leaf_anchor[s_singleton_anchor]
    
    block_pair <- paste(orthopair$query_synteny_block,
                        orthopair$subject_synteny_block, 
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
    rbbh <- h5$blast$rbbh
    hit <- match(rbbh$qseqid, g2g_graph$query_tx$tx)
    rbbh$qgeneid <- g2g_graph$query_gene$gene[hit]
    hit <- match(rbbh$sseqid, g2g_graph$subject_tx$tx)
    rbbh$sgeneid <- g2g_graph$subject_gene$gene[hit]
    rbbh$pair_id <- paste(rbbh$qgeneid, rbbh$sgeneid, sep = "_")
    orthopair$rbbh <- FALSE
    orthopair$rbbh[orthopair$pair_id %in% rbbh$pair_id] <- TRUE
    orthopair$syntenic <- TRUE
    query_orphan <- !rbbh$qgeneid %in% orthopair$query_gene
    subject_orphan <- !rbbh$sgeneid %in% orthopair$subject_gene
    rbbh <- subset(rbbh, subset = query_orphan & subject_orphan)
    rbbh_score <- .getRBHscore(rbh = rbbh)
    rbbh <- subset(rbbh, select = c(qgeneid, qseqid, sgeneid, sseqid, pair_id))
    names(rbbh) <- c("query_gene", "query_tx", "subject_gene", "subject_tx", "pair_id")
    rbbh$mutual_ci <- rbbh_score$mutual_ci
    rbbh$rbbh <- TRUE
    rbbh$syntenic <- FALSE
    orthopair <- rbind(orthopair, rbbh)
    return(orthopair)
}

.findSyntenyBlocks <- function(orthopair){
    orthopair <- orthopair[order(orthopair$leaf), ]
    subject_synteny_block <- .findSyntenyBlocksCore(orthopair = orthopair,
                                                    dataset = "subject")
    
    root_order <- order(orthopair$root)
    subject_synteny_block <- subject_synteny_block[root_order]
    orthopair <- orthopair[root_order, ]
    query_synteny_block <- .findSyntenyBlocksCore(orthopair = orthopair,
                                                  dataset = "query")
    orthopair <- cbind(orthopair,
                       query_synteny_block = query_synteny_block, 
                       subject_synteny_block = subject_synteny_block)
    
    return(orthopair)
}

.findSyntenyBlocksCore <- function(orthopair, dataset){
    if(dataset == "query"){
        is_anchor <- orthopair$query_is_anchor
        chr <- orthopair$query_chr
        anchor_id <- orthopair$root_anchor
        
    } else {
        is_anchor <- orthopair$subject_is_anchor
        chr <- orthopair$subject_chr
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
    d <- subset(orthopair, select = c(query_gene, subject_gene))
    d$query_gene <- paste0("q_", d$query_gene)
    d$subject_gene <- paste0("s_", d$subject_gene)
    query_gene_list <- unique(d$query_gene)
    subject_gene_list <- unique(d$subject_gene)
    g <- graph_from_data_frame(d = d, directed = FALSE)
    grp <- split(V(g)$name, components(g)$membership)
    names(grp) <- paste0(names(grp), "_")
    grp <- unlist(grp)
    grp <- data.frame(grp = names(grp), gene_id = grp)
    grp$grp <- as.numeric(sub("_.*", "", grp$grp))
    grp$query <- grp$gene_id %in% query_gene_list
    grp$subject <- grp$gene_id %in% subject_gene_list
    n_query <- tapply(grp$query, grp$grp, sum)
    hit <- match(grp$grp, as.numeric(names(n_query)))
    grp$n_query <- n_query[hit]
    n_subject <- tapply(grp$subject, grp$grp, sum)
    hit <- match(grp$grp, as.numeric(names(n_subject)))
    grp$n_subject <- n_subject[hit]
    
    sog_1to1 <- which(grp$n_query == 1 & grp$n_subject == 1)
    sog_1toM <- which(grp$n_query == 1 & grp$n_subject != 1)
    sog_Mto1 <- which(grp$n_query != 1 & grp$n_subject == 1)
    sog_MtoM <- which(grp$n_query != 1 & grp$n_subject != 1)
    orthopair$class <- NA
    hit <- d$query_gene %in% grp$gene_id[sog_1to1]
    orthopair$class[hit] <- "1to1"
    hit <- d$query_gene %in% grp$gene_id[sog_1toM]
    orthopair$class[hit] <- "1toM"
    hit <- d$query_gene %in% grp$gene_id[sog_Mto1]
    orthopair$class[hit] <- "Mto1"
    hit <- d$query_gene %in% grp$gene_id[sog_MtoM]
    orthopair$class[hit] <- "MtoM"
    
    hit <- match(d$query_gene, grp$gene_id)
    orthopair$SOG <- grp$grp[hit]
    
    return(orthopair)
}

.filterOrthopair <- function(orthopair, g2g_graph){
    split_gene <- .splitGene(orthopair = orthopair, g2g_graph = g2g_graph)
    hit <- orthopair$query_tx %in% split_gene$query
    orthopair$query_gene[hit] <- -orthopair$query_tx[hit]
    hit <- orthopair$subject_tx %in% split_gene$subject
    orthopair$subject_gene[hit] <- -orthopair$subject_tx[hit]
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
    
    # invalid <- .examineMtoM(orthopair = orthopair)
    # orthopair <- subset(orthopair, subset = !pair_id %in% invalid)
    # orthopair <- .classifyOrthoPair(orthopair = orthopair)
    # invalid <- .examine1toM(orthopair = orthopair)
    # invalid <- .examineMto1(orthopair = orthopair)
    # invalid <- .examine1to1(orthopair = orthopair)
    return(orthopair)
}

.splitGene <- function(orthopair, g2g_graph){
    split_1toM <- .split1toM(orthopair, gff = g2g_graph$query_gff)
    split_Mto1 <- .splitMto1(orthopair, gff = g2g_graph$subject_gff)
    split_MtoM <- .splitMtoM(orthopair, g2g_graph = g2g_graph)
    
    out <- list(query = c(split_1toM, split_MtoM$query),
                subject = c(split_Mto1, split_MtoM$subject))
    return(out)
}

.untangleOrthoPair <- function(orthopair){
    orthopair$index <- seq_along(orthopair$query_gene)
    orthopair <- orthopair[order(orthopair$ci_q2s, decreasing = TRUE), ]
    q_best <- tapply(orthopair$index,
                     orthopair$query_gene, 
                     "[", 1)
    orthopair <- orthopair[order(orthopair$ci_s2q, decreasing = TRUE), ]
    s_best <- tapply(orthopair$index,
                     orthopair$subject_gene, 
                     "[", 1)
    best <- unique(c(q_best, s_best))
    orthopair <- orthopair[orthopair$index %in% best, ]
    local_anchor <- intersect(q_best, s_best)
    rest_pair <- best[!best %in% local_anchor]
    anchor_orthopair <- orthopair[orthopair$index %in% local_anchor, ]
    rest_orthopair <- orthopair[orthopair$index %in% rest_pair, ]
    
    query_to_anchor <- rest_orthopair$query_gene %in% anchor_orthopair$query_gene
    subject_to_anchor <- rest_orthopair$subject_gene %in% anchor_orthopair$subject_gene
    anchored_orthopair <- rest_orthopair[query_to_anchor | subject_to_anchor, ]
    rest_pair <- rest_pair[!rest_pair %in% anchored_orthopair$index]
    rest_orthopair <- rest_orthopair[rest_orthopair$index %in% rest_pair, ]
    
    query_to_anchored <- rest_orthopair$query_gene %in% anchored_orthopair$query_gene
    subject_to_anchored <- rest_orthopair$subject_gene %in% anchored_orthopair$subject_gene
    
    nonanchored_orthopair <- rest_orthopair[!(query_to_anchored | subject_to_anchored), ]
    out <- rbind(anchor_orthopair, anchored_orthopair, nonanchored_orthopair)
}

.reformatOrthoPair <- function(orthopair, g2g_graph){
    q_tx_hit <- match(orthopair$query_tx, g2g_graph$query_gff$tx_index)
    orthopair$query_tx <- g2g_graph$query_gff$Parent[q_tx_hit]
    q_split <- which(orthopair$query_gene < 0)
    orthopair$query_gene <- g2g_graph$query_gff$gene_id[q_tx_hit]
    orthopair$original_query_gene <- orthopair$query_gene
    orthopair$query_gene[q_split] <- paste0(orthopair$query_tx[q_split], ":split")
    
    s_tx_hit <- match(orthopair$subject_tx, g2g_graph$subject_gff$tx_index)
    orthopair$subject_tx <- g2g_graph$subject_gff$Parent[s_tx_hit]
    s_split <- which(orthopair$subject_gene < 0)
    orthopair$subject_gene <- g2g_graph$subject_gff$gene_id[s_tx_hit]
    orthopair$original_subject_gene <- orthopair$subject_gene
    orthopair$subject_gene[s_split] <- paste0(orthopair$subject_tx[s_split], ":split")
    
    orthopair <- subset(orthopair, 
                        select = c(original_query_gene, 
                                   original_subject_gene,
                                   query_gene:subject_chr,
                                   pident:mutual_ci,
                                   is_anchor_pair:SOG))
    orthopair$query_synteny_block <- factor(orthopair$query_synteny_block)
    orthopair$subject_synteny_block <- factor(orthopair$subject_synteny_block)
    orthopair$SOG <- factor(orthopair$SOG)
    orthopair$query_synteny_block <- as.numeric(orthopair$query_synteny_block)
    orthopair$subject_synteny_block <- as.numeric(orthopair$subject_synteny_block)
    orthopair$SOG <- as.numeric(orthopair$SOG)
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
    gff$strand[gff$strand == 1] <- "+"
    gff$strand[gff$strand == 2] <- "-"
    gff <- GRanges(seqnames = gff$seqnames, 
                   ranges = IRanges(start = gff$start,
                                    end = gff$end), 
                   strand = gff$strand,
                   gene_id = gff$gene_id,
                   Parent = gff$Parent,
                   tx_index = gff$tx_index, 
                   gene_index = gff$gene_index)
    query_ol <- findOverlaps(gff, gff)
    query_ol <- as.data.frame(query_ol)
    query_ol$query_tx <- gff$tx_index[query_ol$queryHits]
    query_ol$query_ol_tx <- gff$tx_index[query_ol$subjectHits]
    valid <- gff$gene_index[query_ol$queryHits] == gff$gene_index[query_ol$subjectHits]
    query_ol <- subset(query_ol, subset = query_tx != query_ol_tx & valid)
    query_ol <- unique(subset(query_ol, select = query_tx:query_ol_tx))
    
    sog_1toM <- orthopair$class == "1toM"
    if(sum(sog_1toM) == 0){
        return(NULL)
    }
    orthopair_subset <- subset(orthopair, 
                               subset = sog_1toM,
                               select = c(query_gene:query_tx,
                                          subject_gene:subject_tx, 
                                          SOG))
    orthopair_subset$sog_tx <- paste(orthopair_subset$SOG, 
                                     orthopair_subset$subject_tx,
                                     sep = "_")
    
    sog_tx_par_q_tx <- tapply(orthopair_subset$sog_tx, 
                              orthopair_subset$query_tx,
                              unique)
    n_sog_tx_par_q_tx <- sapply(sog_tx_par_q_tx, length)
    splitable <- n_sog_tx_par_q_tx == 1
    splitable <- sog_tx_par_q_tx[splitable]
    splitable <- data.frame(query_tx = as.numeric(names(splitable)), 
                            subject_tx = as.numeric(sub("[0-9]+_", "", unlist(splitable))),
                            SOG = as.numeric(sub("_.+", "", unlist(splitable))))
    query_ol <- subset(query_ol, 
                       subset = query_ol_tx %in% orthopair_subset$query_tx)
    splitable <- left_join(splitable, query_ol, "query_tx")
    hit <- match(splitable$query_tx, gff$tx_index)
    splitable$query_gene <-  gff$gene_index[hit]
    hit <- match(splitable$query_ol_tx, gff$tx_index)
    splitable$query_ol_gene <-  gff$gene_index[hit]
    splitable <- subset(splitable,
                        subset = query_gene != query_ol_gene | is.na(query_ol_gene))
    out <- splitable$query_tx
    # tmp <- orthopair[sog_1toM]
    # best_pair <- tapply(seq_along(tmp$subject_gene),
    #                     tmp$subject_gene, function(i){
    #                         tmp_i <- tmp[i, ]
    #                         return(tmp_i[which.max(tmp_i$mutual_ci), ])
    #                     })
    # best_pair <- do.call("rbind", best_pair)
    # best_pair <- best_pair[order(best_pair$SOG), ]
    # 
    # best_pair_split <- tapply(seq_along(best_pair$SOG), best_pair$SOG, function(i){
    #     tmp_i <- best_pair[i, ]
    #     if(length(unique(tmp_i$query_tx)) == 1){
    #         return(rep(FALSE, length(i)))
    #     }
    #     best_pair_ol <- sapply(seq_along(tmp_i$query_tx), function(j){
    #         x <- tmp_i$query_tx[j]
    #         x_ol <- query_ol$query_ol_tx[query_ol$query_tx %in% x]
    #         return(any(x_ol %in% tmp_i$query_tx[-j]))
    #     })
    #     return(!best_pair_ol)
    # })
    # best_pair$split <- unlist(best_pair_split)
    # 
    # out <- best_pair
    # if(any(best_pair$split)){
    #     split_gene <- subset(best_pair, subset = split)
    #     check <- tapply(split_gene$query_tx, split_gene$query_gene, unique)
    #     check <- sapply(check, length)
    #     valid_split <- names(check[check > 1])
    #     split_gene <- subset(split_gene, subset = query_gene %in% valid_split)
    #     if(nrow(split_gene) > 1){
    #         split_gene$query_gene <- paste0(split_gene$query_tx, ":split_gene")
    #         non_split_gene <- subset(best_pair,
    #                                  subset = !subject_gene %in% split_gene$subject_gene)
    #         out <- rbind(split_gene, non_split_gene)
    #     }
    # }
    return(out)
}

.splitMto1 <- function(orthopair, gff){
    gff$strand[gff$strand == 1] <- "+"
    gff$strand[gff$strand == 2] <- "-"
    gff <- GRanges(seqnames = gff$seqnames, 
                   ranges = IRanges(start = gff$start,
                                    end = gff$end), 
                   strand = gff$strand,
                   gene_id = gff$gene_id,
                   Parent = gff$Parent,
                   tx_index = gff$tx_index, 
                   gene_index = gff$gene_index)
    subject_ol <- findOverlaps(gff, gff)
    subject_ol <- as.data.frame(subject_ol)
    subject_ol$subject_tx <- gff$tx_index[subject_ol$subjectHits]
    subject_ol$subject_ol_tx <- gff$tx_index[subject_ol$subjectHits]
    valid <- gff$gene_index[subject_ol$queryHits] == gff$gene_index[subject_ol$subjectHits]
    subject_ol <- subset(subject_ol, subset = subject_tx != subject_ol_tx & valid)
    subject_ol <- unique(subset(subject_ol, select = subject_tx:subject_ol_tx))
    
    sog_Mto1 <- orthopair$class == "Mto1"
    if(sum(sog_Mto1) == 0){
        return(NULL)
    }
    orthopair_subset <- subset(orthopair, 
                               subset = sog_Mto1,
                               select = c(query_gene:query_tx,
                                          subject_gene:subject_tx, 
                                          SOG))
    orthopair_subset$sog_tx <- paste(orthopair_subset$SOG, 
                                     orthopair_subset$query_tx,
                                     sep = "_")
    
    sog_tx_par_s_tx <- tapply(orthopair_subset$sog_tx, 
                              orthopair_subset$subject_tx,
                              unique)
    n_sog_tx_par_s_tx <- sapply(sog_tx_par_s_tx, length)
    splitable <- n_sog_tx_par_s_tx == 1
    splitable <- sog_tx_par_s_tx[splitable]
    splitable <- data.frame(subject_tx = as.numeric(names(splitable)), 
                            query_tx = as.numeric(sub("[0-9]+_", "", unlist(splitable))),
                            SOG = as.numeric(sub("_.+", "", unlist(splitable))))
    subject_ol <- subset(subject_ol, 
                         subset = subject_ol_tx %in% orthopair_subset$subject_tx)
    splitable <- left_join(splitable, subject_ol, "subject_tx")
    hit <- match(splitable$subject_tx, gff$tx_index)
    splitable$subject_gene <-  gff$gene_index[hit]
    hit <- match(splitable$subject_ol_tx, gff$tx_index)
    splitable$subject_ol_gene <-  gff$gene_index[hit]
    splitable <- subset(splitable,
                        subset = subject_gene != subject_ol_gene | is.na(subject_ol_gene))
    out <- splitable$subject_tx
    
    # best_pair <- tapply(seq_along(tmp$query_gene),
    #                     tmp$query_gene, function(i){
    #                         tmp_i <- tmp[i, ]
    #                         return(tmp_i[which.max(tmp_i$mutual_ci), ])
    #                     })
    # best_pair <- do.call("rbind", best_pair)
    # best_pair <- best_pair[order(best_pair$SOG), ]
    # 
    # best_pair_split <- tapply(seq_along(best_pair$SOG), best_pair$SOG, function(i){
    #     tmp_i <- best_pair[i, ]
    #     if(length(unique(tmp_i$subject_tx)) == 1){
    #         return(rep(FALSE, length(i)))
    #     }
    #     best_pair_ol <- sapply(seq_along(tmp_i$subject_tx), function(j){
    #         x <- tmp_i$subject_tx[j]
    #         x_ol <- subject_ol$subject_ol_tx[subject_ol$subject_tx %in% x]
    #         return(any(x_ol %in% tmp_i$subject_tx[-j]))
    #     })
    #     return(!best_pair_ol)
    # })
    # best_pair$split <- unlist(best_pair_split)
    # 
    # out <- best_pair
    # if(any(best_pair$split)){
    #     split_gene <- subset(best_pair, subset = split)
    #     check <- tapply(split_gene$subject_tx, split_gene$subject_gene, unique)
    #     check <- sapply(check, length)
    #     valid_split <- names(check[check > 1])
    #     split_gene <- subset(split_gene, subset = subject_gene %in% valid_split)
    #     if(nrow(split_gene) > 1){
    #         split_gene$subject_gene <- paste0(split_gene$subject_tx, ":split_gene")
    #         non_split_gene <- subset(best_pair,
    #                                  subset = !query_gene %in% split_gene$query_gene)
    #         out <- rbind(split_gene, non_split_gene)
    #     }
    # }
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
                               select = c(query_gene:query_tx,
                                          subject_gene:subject_tx, 
                                          SOG))
    orthopair_subset$class <- "1toM"
    orthopair_subset$SOG <- as.numeric(factor(orthopair_subset$query_gene))
    n_member <- table(orthopair_subset$SOG)
    split_1toM <- .split1toM(subset(orthopair_subset,
                                    subset = SOG %in% names(n_member[n_member > 1])), 
                             gff = g2g_graph$query_gff)
    orthopair_subset$class <- "Mto1"
    orthopair_subset$SOG <- as.numeric(factor(orthopair_subset$subject_gene))
    n_member <- table(orthopair_subset$SOG)
    split_Mto1 <- .splitMto1(subset(orthopair_subset,
                                    subset = SOG %in% names(n_member[n_member > 1])),
                             gff = g2g_graph$subject_gff)
    out <- list(query = split_1toM, subject = split_Mto1)
    # hit <- match(tmp$tx_pair_id, split_1toM$tx_pair_id)
    # tmp$split_1toM <- split_1toM$split[hit]
    # hit <- match(tmp$tx_pair_id, split_Mto1$tx_pair_id)
    # tmp$split_Mto1 <- split_Mto1$split[hit]
    # tmp <- tmp[order(tmp$SOG), ]
    # 
    # tmp$split <- tmp$split_1toM & tmp$split_Mto1
    # tmp$split[is.na(tmp$split)] <- FALSE
    # tmp$split <- tmp$split | {tmp$query_is_anchor & tmp$subject_is_anchor}
    # 
    # if(any(tmp$split)){
    #     split_gene <- subset(tmp, subset = split)
    #     push_back <- lapply(seq_along(split_gene$query_gene),
    #                         function(i){
    #                             tmp_i <- tmp[tmp$query_tx %in% split_gene$query_tx[i], ]
    #                             q_valid_i <- tmp_i[tmp_i$mutual_ci == split_gene$mutual_ci[i], ]
    #                             tmp_i <- tmp[tmp$subject_tx %in% split_gene$subject_tx[i], ]
    #                             s_valid_i <- tmp_i[tmp_i$mutual_ci == split_gene$mutual_ci[i], ]
    #                             return(unique(rbind(q_valid_i, s_valid_i)))
    #                         })
    #     push_back <- do.call("rbind", push_back)
    #     split_gene <- rbind(split_gene, push_back)
    #     split_gene <- subset(split_gene, select = -c(tx_pair_id:split))
    #     non_split_gene <- subset(tmp,
    #                              subset = !(query_gene %in% split_gene$query_gene &
    #                                             subject_gene %in% split_gene$subject_gene),
    #                              select = -c(tx_pair_id:split))
    #     out <- list(split_gene = unique(split_gene),
    #                 non_split_gene = unique(non_split_gene))
    #     
    # } else {
    #     tmp <- subset(tmp, select = -c(tx_pair_id:split))
    #     out <- list(split_gene = unique(tmp), non_split_gene = NULL)
    # }
    return(out)
}

.examine1to1 <- function(orthopair, valid){
    sog_1to1 <- orthopair$class == "1to1"
    threshold <-  quantile(orthopair$mutual_ci[sog_1to1], 0.05)
    update_valid <- orthopair$mutual_ci[sog_1to1] >= threshold
    valid[sog_1to1] <- valid[sog_1to1] & update_valid
    return(valid)
}

.examine1toM <- function(orthopair){
    sog_1toM <- orthopair$class == "1toM"
    mci_diff <- tapply(seq_along(orthopair$SOG[sog_1toM]),
                       orthopair$SOG[sog_1toM], 
                       function(i){
                           i_og <- orthopair[sog_1toM, ][i, ]
                           i_max <- max(orthopair$mutual_ci[sog_1toM][i])
                           out_i <- i_max - orthopair$mutual_ci[sog_1toM][i]
                           out_i <- data.frame(subject_gene = orthopair$subject_gene[sog_1toM][i],
                                               mci_diff = out_i)
                           return(out_i)
                       })
    mci_diff <- do.call("rbind", mci_diff)
    threshold <-  quantile(mci_diff$mci_diff[mci_diff$mci_diff > 0], 0.95)
    update_valid <- mci_diff$mci_diff <= threshold
    update_valid <- orthopair$subject_gene[sog_1toM] %in% mci_diff$subject_gene[update_valid]
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
                           out_i <- data.frame(query_gene = orthopair$query_gene[sog_Mto1][i],
                                               mci_diff = out_i)
                           return(out_i)
                       })
    mci_diff <- do.call("rbind", mci_diff)
    threshold <-  quantile(mci_diff$mci_diff[mci_diff$mci_diff > 0], 0.95)
    update_valid <- mci_diff$mci_diff <= threshold
    update_valid <- orthopair$query_gene[sog_Mto1] %in% mci_diff$query_gene[update_valid]
    valid[sog_Mto1] <- valid[sog_Mto1] & update_valid
    return(valid)
}

.examineMtoM <- function(orthopair){
    sog_MtoM <- orthopair$class == "MtoM"
    orthopair <- orthopair[order(orthopair$pident, decreasing = TRUE), ]
    orthopair$index <- seq_along(orthopair$query_gene)
    valid_pair_id <- tapply(seq_along(orthopair$SOG[sog_MtoM]),
                            orthopair$SOG[sog_MtoM], 
                            function(i){
                                i_og <- orthopair[sog_MtoM, ][i, ]
                                i_og <- i_og[order(i_og$ci_q2s, decreasing = TRUE), ]
                                q_best <- tapply(i_og$index,
                                                 i_og$query_gene, 
                                                 "[", 1)
                                i_og <- i_og[order(i_og$ci_s2q, decreasing = TRUE), ]
                                s_best <- tapply(i_og$index,
                                                 i_og$subject_gene, 
                                                 "[", 1)
                                best <- unique(c(q_best, s_best))
                                i_og <- i_og[i_og$index %in% best, ]
                                i_anchor <- intersect(q_best, s_best)
                                rest_og <- best[!best %in% i_anchor]
                                i_anchor_og <- i_og[i_og$index %in% i_anchor, ]
                                i_og <- i_og[i_og$index %in% rest_og, ]
                                
                                q_to_anchor <- i_og$query_gene %in% i_anchor_og$query_gene
                                s_to_anchor <- i_og$subject_gene %in% i_anchor_og$subject_gene
                                i_to_anchor_og <- i_og[q_to_anchor | s_to_anchor, ]
                                rest_og <- rest_og[!rest_og %in% i_to_anchor_og$index]
                                i_og <- i_og[i_og$index %in% rest_og, ]
                                
                                q_to_anchored <- i_og$query_gene %in% i_to_anchor_og$query_gene
                                s_to_anchored <- i_og$subject_gene %in% i_to_anchor_og$subject_gene
                                
                                i_to_nonanchored <- i_og[!(q_to_anchored | s_to_anchored), ]
                                
                                q_id <- NULL
                                s_id <- NULL
                                out_i <- NULL
                                while(TRUE){
                                    q_dup <- i_og$query_gene %in% q_id
                                    s_dup <- i_og$subject_gene %in% s_id
                                    target_pair <- !(q_dup & s_dup)
                                    if(all(!target_pair)){
                                        break
                                    }
                                    out_i <- c(out_i, i_og$pair_id[target_pair][1])
                                    q_id <- c(q_id, i_og$query_gene[target_pair][1])
                                    s_id <- c(s_id, i_og$subject_gene[target_pair][1])
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
