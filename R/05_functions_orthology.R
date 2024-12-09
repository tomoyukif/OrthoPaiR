#' Filter ortholog pairs based on Locally Collinear Blocks (LCBs)
#'
#' This function filters ortholog pairs based on their presence in Locally Collinear Blocks (LCBs).
#'
#' @param object A OrthoPairDB object.
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#' @importFrom parallel mclapply
syntenicOrtho <- function(object, rbbh_mci_threshold = 0.1, rbh_mci_threshold = 0.4){
    # Check if the input object is of class "OrthoPairDB"
    stopifnot(inherits(x = object, "OrthoPairDB"))
    
    # Open the HDF5 file
    h5 <- H5Fopen(object$h5)
    # Ensure the HDF5 file is closed when the function exits
    on.exit(H5Fclose(h5))
    
    # Check if the necessary groups exist in the HDF5 file
    # if(!H5Lexists(h5, "sibeliaz/lcb_pairs")){
    #     stop("Run getLCBpairs to obtain LCB pair info.")
    # }
    if(!H5Lexists(h5, "blast/rbbh")){
        stop("Run rbh with `best = TRUE` to obtain RBBH info.")
    }
    
    # genome_graph <- .prepGenomeGraph(h5 = h5)
    # genome_graph <- h5$sibeliaz$lcbgraph
    # g2g_graph <- .linkGene2Genome(h5 = h5, genome_graph = genome_graph)
    g2g_graph <- .linkGene2Genome(h5 = h5)
    anchor <- .findAnchors(h5 = h5,
                           # genome_graph = genome_graph,
                           g2g_graph = g2g_graph,
                           rbbh_mci_threshold = rbbh_mci_threshold)
    
    t2a_graph <- .linkTx2Anchor(anchor = anchor$anchor, g2g_graph = g2g_graph)
    orthopair <- .findSyntenicOrtho(h5 = h5, g2g_graph = g2g_graph, t2a_graph = t2a_graph, 
                                    rbbh_mci_threshold = anchor$rbbh_mci_threshold,
                                    rbh_mci_threshold = rbh_mci_threshold)
    # orthopair <- .pickNonSyntenicOrtho(h5 = h5, orthopair = orthopair, g2g_graph = g2g_graph)
    orthopair <- .sortSyntenicOrtho(orthopair = orthopair, g2g_graph = g2g_graph)
    
    .h5overwrite(obj = orthopair, file = object$h5, "orthopair_tx")
    
    orthopair <- .classifyOrthoPair(orthopair = orthopair)
    
    split_orthopair <- .splitGene(orthopair = orthopair, g2g_graph = g2g_graph)
    rest_orthopair <- split_orthopair$rest
    
    if(!is.null(rest_orthopair$rest)){
        while(TRUE){
            rest_orthopair <- .classifyOrthoPair(orthopair = rest_orthopair)
            rest_orthopair <- .splitGene(orthopair = rest_orthopair, g2g_graph = g2g_graph)
            if(is.null(rest_orthopair$rest)){
                split_orthopair$splited <- rbind(split_orthopair$splited,
                                                 rest_orthopair$rest)
                break
                
            } else {
                split_orthopair$splited <- rbind(split_orthopair$splited,
                                                 rest_orthopair$splited)
                rest_orthopair <- rest_orthopair$rest
            }
        }
    }
    orthopair_gene <- split_orthopair$splited
    orthopair_gene <- .pickBestPair(orthopair_gene = orthopair_gene)
    orthopair_gene <- .sortSyntenicOrtho(orthopair = orthopair_gene, g2g_graph = g2g_graph)
    orthopair_gene <- .classifyOrthoPair(orthopair = orthopair_gene)
    orthopair_gene <- subset(orthopair_gene, select = -gene_pair_id)
    orphan <- .getOrphan(orthopair_gene = orthopair_gene, g2g_graph = g2g_graph)
    .h5overwrite(obj = orthopair_gene, file = object$h5, "orthopair_gene")
    .h5overwrite(obj = orphan$query, file = object$h5, "orphan_query")
    .h5overwrite(obj = orphan$subject, file = object$h5, "orphan_subject")
    
    .h5overwrite(obj = "orthopair",
                 file = object$h5,
                 name = "data_type")
    .h5overwrite(obj = as.character(Sys.time()), file = object$h5, "timestamp/pairing")
}

.getOrphan <- function(orthopair_gene, g2g_graph){
    query_tx_hit <- g2g_graph$query_tx$tx %in% orthopair_gene$query_tx
    subject_tx_hit <- g2g_graph$subject_tx$tx %in% orthopair_gene$subject_tx
    hit_gene_id <- g2g_graph$query_gene$gene[query_tx_hit]
    query_gene_hit <- g2g_graph$query_gene$gene %in% hit_gene_id
    query_orphan <- g2g_graph$query_gene$gene[!query_gene_hit]
    hit_gene_id <- g2g_graph$subject_gene$gene[subject_tx_hit]
    subject_gene_hit <- g2g_graph$subject_gene$gene %in% hit_gene_id
    subject_orphan <- g2g_graph$subject_gene$gene[!subject_gene_hit]
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

# #' @importFrom GenomicRanges nearest
# .linkGene2Genome <- function(h5, genome_graph){
.linkGene2Genome <- function(h5){
    # query_gr <- unique(subset(genome_graph, 
    #                           select = c(query_genome, query_chr:query_end)))
    # query_gr_id <- query_gr$query_genome
    # query_gr <- .makeGRanges(df = query_gr, genome = "query")
    # query_gr$query_genome <- query_gr_id
    # 
    # subject_gr <- unique(subset(genome_graph, 
    #                             select = c(subject_genome, subject_chr:subject_end)))
    # subject_gr_id <- subject_gr$subject_genome
    # subject_gr <- .makeGRanges(df = subject_gr, genome = "subject")
    # subject_gr$subject_genome <- subject_gr_id
    # 
    gff_ls <- .getGFFlist(h5 = h5)
    q_tx_i <- gff_ls$query_gff$type %in% c("transcript", "mRNA")
    query_gff <- gff_ls$query_gff[q_tx_i]
    # q_g2b_precede <- precede(query_gff, query_gr, ignore.strand = TRUE)
    # q_g2b_follow <- follow(query_gff, query_gr, ignore.strand = TRUE)
    query_tx <- query_gff$ID
    query_tx <- data.frame(tx = query_tx, id = seq_along(query_tx))
    query_gene <- query_gff$gene_id
    query_gene <- data.frame(gene = query_gene,
                             id = as.numeric(factor(query_gene)))
    # query_g2b_edge <- data.frame(root = c(query_tx$id, query_tx$id),
    #                              query_genome = c(query_gr$query_genome[q_g2b_precede], 
    #                                               query_gr$query_genome[q_g2b_follow]))
    # query_g2b_edge <- unique(query_g2b_edge[order(query_g2b_edge$root), ])
    # query_g2b_edge <- subset(query_g2b_edge, subset = !is.na(query_g2b_edge))
    
    s_tx_i <- gff_ls$subject_gff$type %in% c("transcript", "mRNA")
    subject_gff <- gff_ls$subject_gff[s_tx_i]
    # s_g2b_precede <- precede(subject_gff, subject_gr, ignore.strand = TRUE)
    # s_g2b_follow <- follow(subject_gff, subject_gr, ignore.strand = TRUE)
    subject_tx <- subject_gff$ID
    subject_tx <- data.frame(tx = subject_tx, id = seq_along(subject_tx))
    subject_gene <- subject_gff$gene_id
    subject_gene <- data.frame(gene = subject_gene,
                               id = as.numeric(factor(subject_gene)))
    # subject_g2b_edge <- data.frame(subject_genome = c(subject_gr$subject_genome[s_g2b_precede], 
    #                                                   subject_gr$subject_genome[s_g2b_follow]),
    #                                leaf = c(subject_tx$id, subject_tx$id))
    # subject_g2b_edge <- unique(subject_g2b_edge[order(subject_g2b_edge$leaf), ])
    # subject_g2b_edge <- subset(subject_g2b_edge, subset = !is.na(subject_g2b_edge))
    
    out <- list(query_tx = query_tx,
                subject_tx = subject_tx,
                # query_g2b_edge = query_g2b_edge,
                # subject_g2b_edge = subject_g2b_edge,
                query_gene = query_gene,
                subject_gene = subject_gene,
                query_gff = gff_ls$query_gff,
                subject_gff = gff_ls$subject_gff)
    
    # gff_ls <- .getGFFlist(h5 = h5)
    # query_gr <- .makeGRanges(df = genome_graph, genome = "query")
    # query_gr$node_id <- genome_graph$query_blocks$node_id
    # subject_gr <- .makeGRanges(df = genome_graph$subject_blocks, genome = "subject")
    # subject_gr$node_id <- genome_graph$subject_blocks$node_id
    # 
    # q_tx_i <- gff_ls$query_gff$type %in% c("transcript", "mRNA")
    # query_gff <- gff_ls$query_gff[q_tx_i]
    # q_g2b <- nearest(query_gff, query_gr)
    # query_tx <- query_gff$ID
    # query_tx <- data.frame(tx = query_tx, id = seq_along(query_tx))
    # query_gene <- query_gff$gene_id
    # query_gene <- data.frame(gene = query_gene,
    #                          id = as.numeric(factor(query_gene)))
    # query_g2b_edge <- data.frame(root = query_tx$id,
    #                              query_genome = query_gr$node_id[q_g2b])
    # 
    # s_tx_i <- gff_ls$subject_gff$type %in% c("transcript", "mRNA")
    # subject_gff <- gff_ls$subject_gff[s_tx_i]
    # s_g2b <- nearest(subject_gff, subject_gr)
    # subject_tx <- subject_gff$ID
    # subject_tx <- data.frame(tx = subject_tx, id = seq_along(subject_tx))
    # subject_gene <- subject_gff$gene_id
    # subject_gene <- data.frame(gene = subject_gene,
    #                            id = as.numeric(factor(subject_gene)))
    # subject_g2b_edge <- data.frame(subject_genome = subject_gr$node_id[s_g2b],
    #                                leaf = subject_tx$id)
    # out <- list(query_tx = query_tx,
    #             subject_tx = subject_tx,
    #             query_g2b_edge = query_g2b_edge,
    #             subject_g2b_edge = subject_g2b_edge,
    #             query_gene = query_gene,
    #             subject_gene = subject_gene,
    #             query_gff = gff_ls$query_gff,
    #             subject_gff = gff_ls$subject_gff)
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


.getGFFlist <- function(h5, mp = FALSE){
    # Import GFF files for query and subject genomes
    if(mp){
        query_gff <- .importAllGFF(as.vector(h5$files$query_gff))
        mp_query_gff <- .importAllGFF(as.vector(h5$miniprot$s2q_gff))
        subject_gff <- .importAllGFF(as.vector(h5$files$subject_gff))
        mp_subject_gff <- .importAllGFF(as.vector(h5$miniprot$q2s_gff))
        
        # Order the GFF data
        query_gff <- .orderGFF(gff = query_gff)
        subject_gff <- .orderGFF(gff = subject_gff)
        mp_query_gff <- .orderGFF(gff = mp_query_gff)
        mp_subject_gff <- .orderGFF(gff = mp_subject_gff)
        
        query_gff <- .mRNA2transcript(gff = query_gff)
        subject_gff <- .mRNA2transcript(gff = subject_gff)
        mp_query_gff <- .mRNA2transcript(gff = mp_query_gff)
        mp_subject_gff <- .mRNA2transcript(gff = mp_subject_gff)
        
        mp_query_gff <- .setIDforElements(gff = mp_query_gff)
        mp_subject_gff <- .setIDforElements(gff = mp_subject_gff)
        
        # Return the ordered GFF data as a list
        out <- list(query_gff = .checkGFFentiry(query_gff), 
                    subject_gff = .checkGFFentiry(subject_gff),
                    mp_query_gff = mp_query_gff,
                    mp_subject_gff = mp_subject_gff)
        
    } else {
        query_gff <- .importAllGFF(as.vector(h5$files$query_gff))
        subject_gff <- .importAllGFF(as.vector(h5$files$subject_gff))
        
        # Order the GFF data
        query_gff <- .orderGFF(gff = query_gff)
        subject_gff <- .orderGFF(gff = subject_gff)
        
        # Return the ordered GFF data as a list
        out <- list(query_gff = .checkGFFentiry(query_gff),
                    subject_gff = .checkGFFentiry(subject_gff))
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

.setIDforElements <- function(gff){
    element_i <- !gff$type %in% c("gene", "transcript", "mRNA")
    gff$ID[element_i] <- paste(unlist(gff$Parent[element_i]),
                               as.character(gff$type[element_i]),
                               sep = ":")
    gff$Name <- gff$ID
    return(gff)
}

#' @importFrom dplyr left_join
.findAnchors <- function(h5, 
                         # genome_graph,
                         g2g_graph,
                         rbbh_mci_threshold){
    rbbh <- h5$blast$rbbh
    score <- .getRBHscore(rbh = rbbh)
    rbbh_mci_threshold <- quantile(x = score$mutual_ci, probs = rbbh_mci_threshold)
    rbbh <- subset(rbbh, subset = score$mutual_ci >= rbbh_mci_threshold)
    # genome_graph <- subset(genome_graph, 
    #                        select = c(query_genome, subject_genome))
    root_hit <- match(rbbh$qseqid, g2g_graph$query_tx$tx)
    leaf_hit <- match(rbbh$sseqid, g2g_graph$subject_tx$tx)
    anchor <- data.frame(root = g2g_graph$query_tx$id[root_hit],
                         leaf = g2g_graph$subject_tx$id[leaf_hit])
    # rbbh <- data.frame(root = g2g_graph$query_tx$id[root_hit],
    #                    expected_leaf = g2g_graph$subject_tx$id[leaf_hit])
    # rbbh <- left_join(x = rbbh, y = g2g_graph$query_g2b_edge, by = "root")
    # rbbh <- left_join(x = rbbh, y = genome_graph, by = "query_genome",
    #                   relationship = "many-to-many")
    # rbbh <- left_join(x = rbbh, y = g2g_graph$subject_g2b_edge, 
    #                   by = "subject_genome", relationship = "many-to-many")
    # anchor <- subset(rbbh, subset = leaf == expected_leaf, select = c(root, leaf))
    anchor <- unique(anchor)
    out <- list(anchor = anchor, rbbh_mci_threshold = rbbh_mci_threshold)
    return(out)
}

#' @importFrom dplyr left_join
#' @importFrom GenomicRanges precede follow
.linkTx2Anchor <- function(g2g_graph, anchor){
    tx2gene <- data.frame(root = g2g_graph$query_tx$id,
                          query_anchor = g2g_graph$query_gene$id)
    anchor <- left_join(anchor, tx2gene, "root")
    tx2gene <- data.frame(leaf = g2g_graph$subject_tx$id,
                          subject_anchor = g2g_graph$subject_gene$id)
    anchor <- left_join(anchor, tx2gene, "leaf")
    
    query_anchor_gene_i <- which(g2g_graph$query_gene$id %in% anchor$query_anchor)
    query_anchor_gene <- unique(g2g_graph$query_gene$gene[query_anchor_gene_i])
    q_tx_i <- g2g_graph$query_gff$type %in% c("transcript", "mRNA")
    query_gff <- g2g_graph$query_gff[q_tx_i]
    query_anchor_gff <- query_gff[query_gff$gene_id %in% query_anchor_gene]
    hit <- match(query_anchor_gff$gene_id, g2g_graph$query_gene$gene)
    query_anchor_gff$node_id <- g2g_graph$query_gene$id[hit]
    query_non_anchor_gff <- query_gff[!query_gff$gene_id %in% query_anchor_gene]
    hit <- match(query_non_anchor_gff$gene_id, g2g_graph$query_gene$gene)
    query_non_anchor_gff$node_id <- g2g_graph$query_gene$id[hit]
    query_tx2anchor_precede <- precede(query_non_anchor_gff, query_anchor_gff, ignore.strand = TRUE)
    query_tx2anchor_follow <- follow(query_non_anchor_gff, query_anchor_gff, ignore.strand = TRUE)
    nonanchor_hit <- match(query_non_anchor_gff$ID, g2g_graph$query_tx$tx)
    anchor_hit <- match(query_anchor_gff$ID, g2g_graph$query_tx$tx)
    query_tx2anchor <- data.frame(root = c(g2g_graph$query_tx$id[nonanchor_hit],
                                           g2g_graph$query_tx$id[nonanchor_hit],
                                           g2g_graph$query_tx$id[anchor_hit]),
                                  query_anchor = c(query_anchor_gff$node_id[query_tx2anchor_precede],
                                                   query_anchor_gff$node_id[query_tx2anchor_follow],
                                                   query_anchor_gff$node_id))
    query_tx2anchor <- unique(query_tx2anchor[order(query_tx2anchor$root), ])
    query_tx2anchor <- subset(query_tx2anchor, subset = !is.na(query_tx2anchor))
    
    subject_anchor_gene_i <- which(g2g_graph$subject_gene$id %in% anchor$subject_anchor)
    subject_anchor_gene <- unique(g2g_graph$subject_gene$gene[subject_anchor_gene_i])
    s_tx_i <- g2g_graph$subject_gff$type %in% c("transcript", "mRNA")
    subject_gff <- g2g_graph$subject_gff[s_tx_i]
    subject_anchor_gff <- subject_gff[subject_gff$gene_id %in% subject_anchor_gene]
    hit <- match(subject_anchor_gff$gene_id, g2g_graph$subject_gene$gene)
    subject_anchor_gff$node_id <- g2g_graph$subject_gene$id[hit]
    subject_non_anchor_gff <- subject_gff[!subject_gff$gene_id %in% subject_anchor_gene]
    hit <- match(subject_non_anchor_gff$gene_id, g2g_graph$subject_gene$gene)
    subject_non_anchor_gff$node_id <- g2g_graph$subject_gene$id[hit]
    subject_tx2anchor_precede <- precede(subject_non_anchor_gff, subject_anchor_gff, ignore.strand = TRUE)
    subject_tx2anchor_follow <- follow(subject_non_anchor_gff, subject_anchor_gff, ignore.strand = TRUE)
    nonanchor_hit <- match(subject_non_anchor_gff$ID, g2g_graph$subject_tx$tx)
    anchor_hit <- match(subject_anchor_gff$ID, g2g_graph$subject_tx$tx)
    subject_tx2anchor <- data.frame(leaf = c(g2g_graph$subject_tx$id[nonanchor_hit],
                                             g2g_graph$subject_tx$id[nonanchor_hit],
                                             g2g_graph$subject_tx$id[anchor_hit]),
                                    subject_anchor = c(subject_anchor_gff$node_id[subject_tx2anchor_precede],
                                                       subject_anchor_gff$node_id[subject_tx2anchor_follow],
                                                       subject_anchor_gff$node_id))
    subject_tx2anchor <- unique(subject_tx2anchor[order(subject_tx2anchor$leaf), ])
    subject_tx2anchor <- subset(subject_tx2anchor, subset = !is.na(subject_anchor))
    
    out <- list(anchor = anchor,
                query_tx2anchor = query_tx2anchor,
                subject_tx2anchor = subject_tx2anchor,
                query_tx = g2g_graph$query_tx,
                subject_tx = g2g_graph$subject_tx,
                query_gene = g2g_graph$query_gene,
                subject_gene = g2g_graph$subject_gene)
    return(out)
}

#' @importFrom dplyr left_join
.findSyntenicOrtho <- function(h5, 
                               g2g_graph,
                               t2a_graph,
                               rbbh_mci_threshold,
                               rbh_mci_threshold){
    rbh <- h5$blast$rbh
    rbh_score <- .getRBHscore(rbh = rbh)
    rbh_mci_threshold <- rbbh_mci_threshold * rbh_mci_threshold
    rbh <- subset(rbh, subset = rbh_score$mutual_ci >= rbh_mci_threshold)
    
    root_hit <- match(rbh$qseqid, t2a_graph$query_tx$tx)
    leaf_hit <- match(rbh$sseqid, t2a_graph$subject_tx$tx)
    rbh <- data.frame(root = g2g_graph$query_tx$id[root_hit],
                      expected_leaf = g2g_graph$subject_tx$id[leaf_hit])
    rbh <- left_join(rbh, t2a_graph$query_tx2anchor, "root",
                     relationship = "many-to-many")
    rbh <- left_join(rbh,
                     unique(subset(t2a_graph$anchor,
                                   select = query_anchor:subject_anchor)),
                     "query_anchor",
                     relationship = "many-to-many")
    rbh <- left_join(rbh, t2a_graph$subject_tx2anchor, "subject_anchor",
                     relationship = "many-to-many")
    orthopair <- subset(rbh, subset = leaf == expected_leaf)
    orthopair <- data.frame(query_gene = t2a_graph$query_gene$gene[orthopair$root],
                            query_tx = t2a_graph$query_tx$tx[orthopair$root],
                            subject_gene = t2a_graph$subject_gene$gene[orthopair$leaf],
                            subject_tx = t2a_graph$subject_tx$tx[orthopair$leaf])
    
    orthopair <- left_join(orthopair, rbh_score, c("query_tx", "subject_tx"))
    orthopair$gene_pair_id <- paste(orthopair$query_gene, orthopair$subject_gene, sep = "_")
    orthopair <- orthopair[order(orthopair$gene_pair_id), ]
    orthopair <- unique(orthopair)
    return(orthopair)
}

.getRBHscore <- function(rbh){
    # Calculate coverage identity for query to subject (q2s) and subject to query (s2q)
    q2s_covident <- rbh$q2s_pident * rbh$q2s_qcovs * 1e-4
    s2q_covident <- rbh$s2q_pident * rbh$s2q_qcovs * 1e-4
    
    # Calculate mutual coverage identity as the Euclidean distance of the individual coverage identities
    mutual_ci <- sqrt(q2s_covident^2 + s2q_covident^2) / sqrt(2)
    
    rbh_score <- data.frame(query_tx = rbh$qseqid,
                            subject_tx = rbh$sseqid,
                            mutual_ci = mutual_ci)
    
    return(rbh_score)
}

.pickNonSyntenicOrtho <- function(h5, orthopair, g2g_graph){
    rbbh <- h5$blast$rbbh
    hit <- match(rbbh$qseqid, g2g_graph$query_tx$tx)
    rbbh$qgeneid <- g2g_graph$query_gene$gene[hit]
    hit <- match(rbbh$sseqid, g2g_graph$subject_tx$tx)
    rbbh$sgeneid <- g2g_graph$subject_gene$gene[hit]
    rbbh$gene_pair_id <- paste(rbbh$qgeneid, rbbh$sgeneid, sep = "_")
    orthopair$rbbh <- FALSE
    orthopair$rbbh[orthopair$gene_pair_id %in% rbbh$gene_pair_id] <- TRUE
    orthopair$syntenic <- TRUE
    query_orphan <- !rbbh$qgeneid %in% orthopair$query_gene
    subject_orphan <- !rbbh$sgeneid %in% orthopair$subject_gene
    rbbh <- subset(rbbh, subset = query_orphan & subject_orphan)
    rbbh_score <- .getRBHscore(rbh = rbbh)
    rbbh <- subset(rbbh, select = c(qgeneid, qseqid, sgeneid, sseqid, gene_pair_id))
    names(rbbh) <- c("query_gene", "query_tx", "subject_gene", "subject_tx", "gene_pair_id")
    rbbh$mutual_ci <- rbbh_score$mutual_ci
    rbbh$rbbh <- TRUE
    rbbh$syntenic <- FALSE
    orthopair <- rbind(orthopair, rbbh)
    return(orthopair)
}

.sortSyntenicOrtho <- function(orthopair, g2g_graph){
    query_id_hit <- match(orthopair$query_tx, g2g_graph$query_tx$tx)
    orthopair <- orthopair[order(query_id_hit), ]
    return(orthopair)
}

#' @importFrom igraph graph_from_data_frame V components
.classifyOrthoPair <- function(orthopair){
    query_gene_list <- unique(orthopair$query_gene)
    subject_gene_list <- unique(orthopair$subject_gene)
    g <- graph_from_data_frame(d = subset(orthopair,
                                          select = c(query_gene,
                                                     subject_gene)),
                               directed = FALSE)
    grp <- split(V(g)$name, components(g)$membership)
    names(grp) <- paste0(names(grp), "_")
    grp <- unlist(grp)
    grp <- data.frame(grp = names(grp), gene_id = grp)
    grp$grp <- as.numeric(sub("_.?", "", grp$grp))
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
    hit <- orthopair$query_gene %in% grp$gene_id[sog_1to1]
    orthopair$class[hit] <- "1to1"
    hit <- orthopair$query_gene %in% grp$gene_id[sog_1toM]
    orthopair$class[hit] <- "1toM"
    hit <- orthopair$query_gene %in% grp$gene_id[sog_Mto1]
    orthopair$class[hit] <- "Mto1"
    hit <- orthopair$query_gene %in% grp$gene_id[sog_MtoM]
    orthopair$class[hit] <- "MtoM"
    
    hit <- match(orthopair$query_gene, grp$gene_id)
    orthopair$SOG <- grp$grp[hit]
    
    return(orthopair)
}

.splitGene <- function(orthopair, g2g_graph){
    orthopair$original_query_gene <- orthopair$query_gene
    orthopair$original_subject_gene <- orthopair$subject_gene
    best_1to1 <- .best1to1(orthopair = orthopair)
    split_1toM <- .split1toM(orthopair, gff = g2g_graph$query_gff)
    split_Mto1 <- .splitMto1(orthopair, gff = g2g_graph$subject_gff)
    split_MtoM <- .splitMtoM(orthopair, g2g_graph = g2g_graph)
    splited <- NULL
    if(!is.null(best_1to1)){
        splited <- rbind(splited, best_1to1)
    }
    if(!is.null(split_1toM)){
        splited <- rbind(splited, subset(split_1toM, select = -split))
    }
    if(!is.null(split_Mto1)){
        splited <- rbind(splited, subset(split_Mto1, select = -split))
    }
    if(!is.null(split_MtoM$split_gene)){
        splited <- rbind(splited, split_MtoM$split_gene)
    }
    if(!is.null(split_MtoM$non_split_gene)){
        rest <- split_MtoM$non_split_gene
        
    } else {
        rest <- NULL
    }
    out <- list(splited = splited, rest = rest)
    return(out)
}

.best1to1 <- function(orthopair){
    sog_1to1 <- orthopair$class == "1to1"
    tmp <- orthopair[sog_1to1, ]
    if(nrow(tmp) == 0){
        return(NULL)
    }
    best_pair <- tapply(seq_along(tmp$subject_gene),
                        tmp$subject_gene, function(i){
                            tmp_i <- tmp[i, ]
                            return(tmp_i[which.max(tmp_i$mutual_ci), ])
                        })
    best_pair <- do.call("rbind", best_pair)
    best_pair <- best_pair[order(best_pair$SOG), ]
    return(best_pair)
}

#' @importFrom GenomicRanges findOverlaps
.split1toM <- function(orthopair, gff){
    cds_i <- gff$type == "CDS"
    query_ol <- findOverlaps(gff[cds_i],
                             gff[cds_i])
    query_ol <- as.data.frame(query_ol)
    query_ol$query_tx <- unlist(gff$Parent[cds_i][query_ol$queryHits])
    query_ol$query_ol_tx <- unlist(gff$Parent[cds_i][query_ol$subjectHits])
    query_ol <- unique(subset(query_ol, select = query_tx:query_ol_tx))
    query_ol <- subset(query_ol, subset = query_tx != query_ol_tx)
    sog_1toM <- orthopair$class == "1toM"
    tmp <- orthopair[sog_1toM, ]
    if(nrow(tmp) == 0){
        return(NULL)
    }
    best_pair <- tapply(seq_along(tmp$subject_gene),
                        tmp$subject_gene, function(i){
                            tmp_i <- tmp[i, ]
                            return(tmp_i[which.max(tmp_i$mutual_ci), ])
                        })
    best_pair <- do.call("rbind", best_pair)
    best_pair <- best_pair[order(best_pair$SOG), ]
    
    best_pair_split <- tapply(seq_along(best_pair$SOG), best_pair$SOG, function(i){
        tmp_i <- best_pair[i, ]
        if(length(unique(tmp_i$query_tx)) == 1){
            return(rep(FALSE, length(i)))
        }
        best_pair_ol <- sapply(seq_along(tmp_i$query_tx), function(j){
            x <- tmp_i$query_tx[j]
            x_ol <- query_ol$query_ol_tx[query_ol$query_tx %in% x]
            return(any(x_ol %in% tmp_i$query_tx[-j]))
        })
        return(!best_pair_ol)
    })
    best_pair$split <- unlist(best_pair_split)
    
    out <- best_pair
    if(any(best_pair$split)){
        split_gene <- subset(best_pair, subset = split)
        check <- tapply(split_gene$query_tx, split_gene$query_gene, unique)
        check <- sapply(check, length)
        valid_split <- names(check[check > 1])
        split_gene <- subset(split_gene, subset = query_gene %in% valid_split)
        if(nrow(split_gene) > 1){
            split_gene$query_gene <- paste0(split_gene$query_tx, ":split_gene")
            non_split_gene <- subset(best_pair,
                                     subset = !subject_gene %in% split_gene$subject_gene)
            out <- rbind(split_gene, non_split_gene)
        }
    }
    return(out)
}

.splitMto1 <- function(orthopair, gff){
    cds_i <- gff$type == "CDS"
    subject_ol <- findOverlaps(gff[cds_i],
                               gff[cds_i])
    subject_ol <- as.data.frame(subject_ol)
    subject_ol$subject_tx <- unlist(gff$Parent[cds_i][subject_ol$queryHits])
    subject_ol$subject_ol_tx <- unlist(gff$Parent[cds_i][subject_ol$subjectHits])
    subject_ol <- unique(subset(subject_ol, select = subject_tx:subject_ol_tx))
    subject_ol <- subset(subject_ol, subset = subject_tx != subject_ol_tx)
    sog_Mto1 <- orthopair$class == "Mto1"
    tmp <- orthopair[sog_Mto1, ]
    if(nrow(tmp) == 0){
        return(NULL)
    }
    best_pair <- tapply(seq_along(tmp$query_gene),
                        tmp$query_gene, function(i){
                            tmp_i <- tmp[i, ]
                            return(tmp_i[which.max(tmp_i$mutual_ci), ])
                        })
    best_pair <- do.call("rbind", best_pair)
    best_pair <- best_pair[order(best_pair$SOG), ]
    
    best_pair_split <- tapply(seq_along(best_pair$SOG), best_pair$SOG, function(i){
        tmp_i <- best_pair[i, ]
        if(length(unique(tmp_i$subject_tx)) == 1){
            return(rep(FALSE, length(i)))
        }
        best_pair_ol <- sapply(seq_along(tmp_i$subject_tx), function(j){
            x <- tmp_i$subject_tx[j]
            x_ol <- subject_ol$subject_ol_tx[subject_ol$subject_tx %in% x]
            return(any(x_ol %in% tmp_i$subject_tx[-j]))
        })
        return(!best_pair_ol)
    })
    best_pair$split <- unlist(best_pair_split)
    
    out <- best_pair
    if(any(best_pair$split)){
        split_gene <- subset(best_pair, subset = split)
        check <- tapply(split_gene$subject_tx, split_gene$subject_gene, unique)
        check <- sapply(check, length)
        valid_split <- names(check[check > 1])
        split_gene <- subset(split_gene, subset = subject_gene %in% valid_split)
        if(nrow(split_gene) > 1){
            split_gene$subject_gene <- paste0(split_gene$subject_tx, ":split_gene")
            non_split_gene <- subset(best_pair,
                                     subset = !query_gene %in% split_gene$query_gene)
            out <- rbind(split_gene, non_split_gene)
        }
    }
    return(out)
}

.splitMtoM <- function(orthopair, g2g_graph){
    sog_MtoM <- orthopair$class == "MtoM"
    tmp <- orthopair[sog_MtoM, ]
    if(nrow(tmp) == 0){
        return(NULL)
    }
    tmp$tx_pair_id <- paste(tmp$query_tx, tmp$subject_tx, sep = "_")
    tmp_Mto1 <- tmp_1toM <- tmp
    tmp_1toM$class <- "1toM"
    tmp_Mto1$class <- "Mto1"
    split_1toM <- .split1toM(tmp_1toM, gff = g2g_graph$query_gff)
    split_Mto1 <- .splitMto1(tmp_Mto1, gff = g2g_graph$subject_gff)
    hit <- match(tmp$tx_pair_id, split_1toM$tx_pair_id)
    tmp$split_1toM <- split_1toM$split[hit]
    hit <- match(tmp$tx_pair_id, split_Mto1$tx_pair_id)
    tmp$split_Mto1 <- split_Mto1$split[hit]
    tmp <- tmp[order(tmp$SOG), ]
    
    tmp$split <- tmp$split_1toM & tmp$split_Mto1
    tmp$split[is.na(tmp$split)] <- FALSE
    
    if(any(tmp$split)){
        split_gene <- subset(tmp, subset = split)
        push_back <- lapply(seq_along(split_gene$query_gene),
                            function(i){
                                tmp_i <- tmp[tmp$query_tx %in% split_gene$query_tx[i], ]
                                q_valid_i <- tmp_i[tmp_i$mutual_ci == split_gene$mutual_ci[i], ]
                                tmp_i <- tmp[tmp$subject_tx %in% split_gene$subject_tx[i], ]
                                s_valid_i <- tmp_i[tmp_i$mutual_ci == split_gene$mutual_ci[i], ]
                                return(unique(rbind(q_valid_i, s_valid_i)))
                            })
        push_back <- do.call("rbind", push_back)
        split_gene <- rbind(split_gene, push_back)
        split_gene <- subset(split_gene, select = -c(tx_pair_id:split))
        non_split_gene <- subset(tmp,
                                 subset = !query_gene %in% split_gene$query_gene &
                                     !subject_gene %in% split_gene$subject_gene,
                                 select = -c(tx_pair_id:split))
        out <- list(split_gene = split_gene, non_split_gene = non_split_gene)
        
    } else {
        tmp <- subset(tmp, select = -c(tx_pair_id:split))
        out <- list(split_gene = tmp, non_split_gene = NULL)
    }
    return(out)
}

.pickBestPair <- function(orthopair_gene){
    out <- tapply(seq_along(orthopair_gene$gene_pair_id),
                  orthopair_gene$gene_pair_id,
                  function(i){
                      best_pair <- which.max(orthopair_gene$mutual_ci[i])
                      return(orthopair_gene[i[best_pair], ])
                  })
    out <- do.call("rbind", out)
    return(out)
}
