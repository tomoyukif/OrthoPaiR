#' Define a function to filter ortholog pairs based on LCBs
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#'
anchorOrtho <- function(object,
                        non1to1 = TRUE){
    stopifnot(inherits(x = object, "SynogDB"))
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    if(!H5Lexists(h5, "sibeliaz/lcb_pairs")){
        stop("Run getLCBpairs to obtain LCB pair info.")
    }
    if(!H5Lexists(h5, "blast/rbbh")){
        stop("Run rbh with `best = TRUE` to obtain RBBH info.")
    }

    obj <- .makeObject(object = object, h5 = h5)

    obj <- .orthoIn1to1lcb(obj = obj)
    obj <- .orthoInNon1to1lcb(obj = obj)
    obj <- .orthoIn1to1cbi(obj = obj)

    .h5overwrite(obj = obj$out, file = object$h5, "anchor")
}

.makeObject <- function(object, h5){
    pairs <- list(lcb_1to1 = .order(df = h5$sibeliaz$lcb_pairs$lcb_1to1),
                  lcb_non_1to1 = .order(df = h5$sibeliaz$lcb_pairs$lcb_non_1to1))

    gr <- .df2gr(pairs = pairs)

    gap_gr <- list(query = .gapGR(gr = gr$query_1to1,
                                  chrLen = object$genome$query),
                   subject = .gapGR(gr = gr$subject_1to1,
                                    chrLen = object$genome$subject))

    query_gff <- .importAllGFF(object$query_gff)
    subject_gff <- .importAllGFF(object$subject_gff)
    gff <- .subsetGFF(query_gff = query_gff, subject_gff = subject_gff)

    rbbh <- h5$blast$rbbh
    rbbh$index <- seq_along(rbbh$qseqid)

    out <- list(gr = gr,
                rbbh = rbbh,
                gap_gr = gap_gr,
                gff = gff)
    return(out)
}

.order <- function(df, by = "query"){
    if(by == "query"){
        q_order <- order(df$query_chr, df$query_start)
        df <- df[q_order, ]
    } else {
        s_order <- order(df$subject_chr,
                         df$subject_start)
        df <- df[s_order, ]
    }
    return(df)
}

.df2gr <- function(pairs){
    gr <- list(query_1to1 = .makeGRanges(df = pairs$lcb_1to1,
                                         genome = "query"),
               subject_1to1 = .makeGRanges(df = pairs$lcb_1to1,
                                           genome = "subject"),
               query_non1to1 = .makeGRanges(df = pairs$lcb_non_1to1,
                                            genome = "query"),
               subject_non1to1 = .makeGRanges(df = pairs$lcb_non_1to1,
                                              genome = "subject"))
    return(gr)
}

#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
.makeGRanges <- function(df, genome){
    if(genome == "query"){
        minus <- df$query_start > df$query_end
        tmp <- df$query_start[minus]
        df$query_start[minus] <- df$query_end[minus]
        df$query_end[minus] <- tmp
        gr <- GRanges(seqnames = df$query_chr,
                      ranges = IRanges(start = df$query_start,
                                       end = df$query_end),
                      strand = "*")

    } else {
        minus <- df$subject_start > df$subject_end
        tmp <- df$subject_start[minus]
        df$subject_start[minus] <- df$subject_end[minus]
        df$subject_end[minus] <- tmp
        gr <- GRanges(seqnames = df$subject_chr,
                      ranges = IRanges(start = df$subject_start,
                                       end = df$subject_end),
                      strand = "*")
    }
    return(gr)
}

.subsetGFF <- function(query_gff, subject_gff){
    query_gff <- query_gff[query_gff$type %in% c("transcript", "mRNA")]
    subject_gff <- subject_gff[subject_gff$type %in% c("transcript", "mRNA")]

    query_gff <- .orderGFF(gff = query_gff)
    subject_gff <- .orderGFF(gff = subject_gff)
    return(list(query = query_gff, subject = subject_gff))
}

.orderGFF <- function(gff){
    gff_order <- order(as.character(seqnames(gff)), start(gff))
    gff <- gff[gff_order]
    return(gff)
}

.gapGR <- function(gr, chrLen){
    n_gr <- length(gr)
    gap_gr <- data.frame(chr_start = as.character(seqnames(gr[-n_gr])),
                         chr_end = as.character(seqnames(gr[-1])),
                         start = end(gr[-n_gr]) + 1,
                         end = start(gr[-1]) - 1,
                         start_block = seq_along(gr)[-n_gr],
                         end_block = seq_along(gr)[-1])

    gap_gr <- rbind(subset(gap_gr, subset = chr_start == chr_end),
                    .gapGrStartEdge(gr = gr, gap_gr = gap_gr),
                    .gapGrEndEdge(gr = gr,
                                  gap_gr = gap_gr,
                                  chrLen = chrLen,
                                  n_gr = n_gr))

    gap_gr <- .flipMinusStrand(df = gap_gr)

    gap_gr <- GRanges(seqnames = gap_gr$chr_start,
                      ranges = IRanges(start = gap_gr$start,
                                       end = gap_gr$end),
                      start_block = gap_gr$start_block,
                      end_block = gap_gr$end_block)

    gap_gr <- gap_gr[width(gap_gr) > 1]
    return(gap_gr)
}

#' @importFrom GenomeInfoDb seqnames
.gapGrStartEdge <- function(gr, gap_gr){
    out <- subset(gap_gr, subset = chr_start != chr_end)
    out$chr_start <- out$chr_end
    out$start <- 1
    out$start_block <- NA
    chr <- as.character(seqnames(gr[1]))
    end <- start(gr[1]) - 1
    out <- rbind(data.frame(chr_start = chr,
                            start = 1,
                            chr_end = chr,
                            end = end,
                            start_block = NA,
                            end_block = 1),
                 out)
    return(out)
}

#' @importFrom GenomeInfoDb seqnames
.gapGrEndEdge <- function(gr, gap_gr, chrLen, n_gr){
    out <- subset(gap_gr, subset = chr_start != chr_end)
    out$chr_end <- out$chr_start
    out$end <- chrLen$length[match(out$chr_start, chrLen$names)]
    chr <- as.character(seqnames(gr[n_gr]))
    start <- end(gr[n_gr]) + 1
    out$end_block <- NA
    out <- rbind(out,
                 data.frame(chr_start = chr,
                            start = start,
                            chr_end = chr,
                            end = chrLen$length[match(chr, chrLen$names)],
                            start_block = n_gr,
                            end_block = NA))
    return(out)
}

.flipMinusStrand <- function(df){
    minus <- df$start > df$end
    tmp <- df$start[minus]
    df$start[minus] <- df$end[minus]
    df$end[minus] <- tmp
    return(df)
}

.orthoIn1to1lcb <- function(obj){
    out <- .orthoInBlock(gr_query = obj$gr$query_1to1,
                         gr_subject = obj$gr$subject_1to1,
                         rbbh = obj$rbbh,
                         gff_query = obj$gff$query,
                         gff_subject = obj$gff$subject)
    obj$rbbh <- subset(obj$rbbh, subset = !index %in% out$index)
    obj$out <- cbind(out, type = "1to1lcb")
    return(obj)
}

.orthoInNon1to1lcb <- function(obj){
    out <- .orthoInBlock(gr_query = obj$gr$query_non1to1,
                         gr_subject = obj$gr$subject_non1to1,
                         rbbh = obj$rbbh,
                         gff_query = obj$gff$query,
                         gff_subject = obj$gff$subject)
    obj$rbbh <- subset(obj$rbbh, subset = !index %in% out$index)
    obj$out <- rbind(obj$out, cbind(out, type = "non1to1lcb"))
    return(obj)
}

.orthoIn1to1cbi <- function(obj){
    out <- .orthoInBlock(gr_query = obj$gap_gr$query,
                         gr_subject = obj$gap_gr$subject,
                         rbbh = obj$rbbh,
                         gff_query = obj$gff$query,
                         gff_subject = obj$gff$subject,
                         interval = TRUE)
    obj$rbbh <- subset(obj$rbbh, subset = !index %in% out$index)
    obj$out <- rbind(obj$out, cbind(out, type = "1to1cbi"))
    return(obj)
}

.orthoInBlock <- function(gr_query,
                          gr_subject,
                          rbbh,
                          gff_query,
                          gff_subject,
                          interval = FALSE){
    q_g2b <-  .gene2block(gr = gr_query,
                          rbbh_id = rbbh$qseqid,
                          gff = gff_query)

    s_g2b <-  .gene2block(gr = gr_subject,
                          rbbh_id = rbbh$sseqid,
                          gff = gff_subject)
    if(interval){
        q_g2b <- .getInterval(g2b = q_g2b, gr = gr_query)
        s_g2b <- .getInterval(g2b = s_g2b, gr = gr_subject)
    }
    valid <- q_g2b == s_g2b
    out <- subset(rbbh, subset = valid)
    return(out)
}

#' @importFrom GenomicRanges findOverlaps
.gene2block <- function(gr, rbbh_id, gff){
    ol <- findOverlaps(gff, gr)
    out <- rep(NA, length(gff))
    out[queryHits(ol)] <- subjectHits(ol)
    gene_index <- match(rbbh_id, gff$ID)
    return(out[gene_index])
}

.getInterval <- function(g2b, gr){
    not_na <- !is.na(g2b)
    start_int <- end_int <- g2b[not_na]
    start_int <- gr$start_block[start_int]
    end_int <- gr$end_block[end_int]
    out <- rep(NA, length(g2b))
    out[not_na] <- paste(start_int, end_int, sep = "_")
    return(out)
}

################################################################################
#' Define a function to filter ortholog pairs based on gene synteny.
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
syntenyOrtho <- function(object,
                         omit_chr = "",
                         pident = 90,
                         evalue = 1e-50,
                         qcovs = 50){
    stopifnot(inherits(x = object, "SynogDB"))
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    if(!H5Lexists(h5, "anchor")){
        stop("Run anchorOrtho to obtain ortholog anchors info.")
    }
    if(!H5Lexists(h5, "blast/rbh")){
        stop("Run rbh with `best = FALSE` to obtain RBH info.")
    }

    gff_ls <- .getGFFlist(object = object)
    anchor_ls <- .makeAnchorList(gff_ls = gff_ls, h5 = h5, omit_chr = omit_chr)
    rbh_dist <- .getDist(gff_ls = gff_ls, h5 = h5, anchor_ls = anchor_ls)
    tp_list <- read.csv(object$positive_list)
    fp_list <- read.csv(object$negative_list)
    metrics <- .getMetrics(rbh_dist = rbh_dist,
                           tp_list = tp_list,
                           fp_list = fp_list)
    rbh_dist$so_valid <- rbh_dist$dist <= metrics$so_threshold
    syn_og <- .makeSynogDF(rbh = rbh_dist, h5 = h5)
    syn_og <- .getMutualScore(ortho = syn_og)
    syn_og$class <- .classifySO(ortho = syn_og)

    rbbh_og <- .addRBBH(rbbh = h5$blast$rbbh, syn_og = syn_og)
    rbbh_og <- .filterOG(target = rbbh_og, ortho = rbh_dist,
                         pident = 90, evalue = 1e-50, qcovs = 50)
    syn_og <- .mergeSynogDF(syn_og = syn_og, rbbh_og = rbbh_og, gff_ls = gff_ls)

    syn_og$class <- .classifySO(ortho = syn_og)
    orphan <- .getOrphan(gff_ls = gff_ls, syn_og = syn_og)
    out <- .makeOutput(syn_og = syn_og, orphan = orphan)
    .h5creategroup(object$h5,"synog_tx")
    .h5overwrite(obj = out$syn_og, file = object$h5, "synog_tx/orthopairs")
    .h5overwrite(obj = out$summary, file = object$h5, "synog_tx/summary")
    .h5overwrite(obj = out$orphan, file = object$h5, "synog_tx/orphan")
    .h5overwrite(obj = metrics$so_metrics, file = object$h5, "synog_tx/metrics")
    .h5overwrite(obj = metrics$so_threshold, file = object$h5, "synog_tx/so_threshold")
    invisible(metrics$metrics_plot)
}

.getGFFlist <- function(object){
    query_gff <- .importAllGFF(object$query_gff)
    subject_gff <- .importAllGFF(object$subject_gff)
    query_gff <- .orderGFF(gff = query_gff)
    subject_gff <- .orderGFF(gff = subject_gff)
    out <- list(query_gff =query_gff, subject_gff = subject_gff)
    return(out)
}

.makeAnchorList <- function(gff_ls, h5, omit_chr){
    anchor <- h5$anchor
    query_anchor <- .findNearestAnchor(gff = gff_ls$query_gff,
                                       anchor = anchor$qseqid,
                                       omit_chr = omit_chr)
    subject_anchor <- .findNearestAnchor(gff = gff_ls$subject_gff,
                                         anchor = anchor$sseqid,
                                         omit_chr = omit_chr)
    gene_id <- .setGeneID2Anchor(anchor = anchor, gff_ls = gff_ls)
    anchor <- cbind(anchor, gene_id)
    out <- .orgAnchor(anchor = anchor,
                      query_anchor = query_anchor,
                      subject_anchor = subject_anchor)
    return(out)
}

#' @importFrom GenomeInfoDb seqnames
#' @importFrom  BiocGenerics start
.findNearestAnchor <- function(gff, anchor, omit_chr){
    tx <- gff[gff$type %in% c("transcript", "mRNA")]
    anchor <- unlist(tx$Parent[match(anchor, tx$ID)])
    anchor <- na.omit(anchor)
    gff <- gff[gff$type == "gene"]
    gff <- gff[order(as.numeric(seqnames(gff)), start(gff))]
    gff <- data.frame(chr = as.character(seqnames(gff)),
                      start = start(gff),
                      ID = gff$ID)
    gff <- gff[!grepl(omit_chr, gff$chr), ]
    gff$anchor <- FALSE
    gff$anchor[gff$ID %in% anchor] <- TRUE
    gff$nearest_anchor <- vapply(X = seq_along(gff$anchor),
                                 FUN.VALUE = numeric(1),
                                 FUN = .getNearestAnchorIndex, gff = gff)
    gff$location <- seq_along(gff$chr)
    return(gff)
}

.getNearestAnchorIndex <- function(index, gff){
    target <- gff$chr %in% gff$chr[index] & gff$anchor
    near <- which.min(abs(gff$start[target] - gff$start[index]))
    near <- which(target)[near]
    return(near)
}

#' @importFrom GenomeInfoDb seqnames
.setGeneID2Anchor <- function(anchor, gff_ls){
    tx_i <- gff_ls$query_gff$type %in% c("transcript", "mRNA")
    tx <- gff_ls$query_gff[tx_i]
    hit <- match(anchor$qseqid, tx$ID)
    qgeneid <- unlist(tx$Parent[hit])
    qchr <- as.character(seqnames(tx[hit]))
    tx_i <- gff_ls$subject_gff$type %in% c("transcript", "mRNA")
    tx <- gff_ls$subject_gff[tx_i]
    hit <- match(anchor$sseqid, tx$ID)
    sgeneid <- unlist(tx$Parent[hit])
    schr <- as.character(seqnames(tx[hit]))
    out <- data.frame(qgeneid = qgeneid, sgeneid = sgeneid,
                      qchr = qchr, schr = schr)
    return(out)
}

.orgAnchor <- function(anchor, query_anchor, subject_anchor){
    anchor <- subset(anchor, select = qgeneid:schr)
    anchor <- unique(anchor)
    hit <- match(anchor$qgeneid, query_anchor$ID)
    anchor$query_anchor <- query_anchor$nearest_anchor[hit]
    tbl <- table(anchor$query_anchor)
    multiple <- anchor$query_anchor %in% names(tbl[tbl > 1])
    anchor$q2s_multiple <- FALSE
    anchor$q2s_multiple[multiple] <- TRUE
    hit <- match(anchor$sgeneid, subject_anchor$ID)
    anchor$subject_anchor <- subject_anchor$nearest_anchor[hit]
    tbl <- table(anchor$subject_anchor)
    multiple <- anchor$subject_anchor %in% names(tbl[tbl > 1])
    anchor$s2q_multiple <- FALSE
    anchor$s2q_multiple[multiple] <- TRUE
    out <- list(query = query_anchor, subject = subject_anchor, anchor = anchor)
    return(out)
}

.getDist <- function(gff_ls, h5, anchor_ls){
    rbh <- .extractRBH(h5 = h5, gff_ls = gff_ls)
    anchor2rbh <- .setAnchor2RBH(rbh = rbh, anchor_ls = anchor_ls)
    rbh <- cbind(rbh, anchor2rbh)
    q2s_dist <- .traceAnchorQ2S(rbh = rbh, anchor_ls = anchor_ls)
    s2q_dist <- .traceAnchorS2Q(rbh = rbh, anchor_ls = anchor_ls)
    rbh <- cbind(rbh, dist = rowMin(cbind(q2s_dist, s2q_dist)))
    return(rbh)
}

.extractRBH <- function(h5, gff_ls){
    out <- h5$blast$rbh
    q_hit <- match(out$qseqid, gff_ls$query_gff$ID)
    out$qgeneid <- gff_ls$query_gff$gene_id[q_hit]
    out$qchr <- as.character(seqnames(gff_ls$query_gff[q_hit]))
    s_hit <- match(out$sseqid, gff_ls$subject_gff$ID)
    out$sgeneid <- gff_ls$subject_gff$gene_id[s_hit]
    out$schr <- as.character(seqnames(gff_ls$subject_gff[s_hit]))
    out$pair_id <- paste(out$qgeneid, out$sgeneid, sep = "_")
    return(out)
}

.setAnchor2RBH <- function(rbh, anchor_ls){
    query_hit <- match(rbh$qgeneid, anchor_ls$query$ID)
    q2s_q_anchor <- anchor_ls$query$nearest_anchor[query_hit]
    s2q_q_location <- anchor_ls$query$location[query_hit]
    subject_hit <- match(rbh$sgeneid, anchor_ls$subject$ID)
    s2q_s_anchor <- anchor_ls$subject$nearest_anchor[subject_hit]
    q2s_s_location <- anchor_ls$subject$location[subject_hit]
    out <- data.frame(q2s_q_anchor = q2s_q_anchor,
                      s2q_q_location = s2q_q_location,
                      s2q_s_anchor = s2q_s_anchor,
                      q2s_s_location = q2s_s_location)
    return(out)
}

.traceAnchorQ2S <- function(rbh, anchor_ls){
    q2s_q_anchor_s <- match(rbh$q2s_q_anchor, anchor_ls$anchor$query_anchor)
    q2s_q_anchor_s[is.na(rbh$q2s_q_anchor)] <- NA
    multiple_anchor <- anchor_ls$anchor$q2s_multiple[q2s_q_anchor_s]
    multiple_anchor[is.na(multiple_anchor)] <- FALSE
    single_projection <- anchor_ls$anchor$subject_anchor[q2s_q_anchor_s]
    single_projection[multiple_anchor] <- NA
    q2s_dist <- abs(rbh$q2s_s_location - single_projection)
    single_projection_chr <- anchor_ls$anchor$schr[q2s_q_anchor_s]
    q2s_dist_chr_unmatch <- rbh$schr != single_projection_chr
    q2s_dist[q2s_dist_chr_unmatch] <- Inf
    multiple_projection_rbh <- rbh[multiple_anchor, ]
    q2s_dist_multi <- vapply(X = seq_along(multiple_projection_rbh$q2s_q_anchor),
                             FUN.VALUE = numeric(1), FUN = .solveMultiProjectionQ2S,
                             rbh = multiple_projection_rbh, anchor_ls = anchor_ls)
    q2s_dist[multiple_anchor] <- q2s_dist_multi
    q2s_dist[is.na(q2s_dist)] <- Inf
    return(q2s_dist)
}

.traceAnchorS2Q <- function(rbh, anchor_ls){
    s2q_s_anchor_q <- match(rbh$s2q_s_anchor, anchor_ls$anchor$subject_anchor)
    s2q_s_anchor_q[is.na(rbh$s2q_s_anchor)] <- NA
    multiple_anchor <- anchor_ls$anchor$s2q_multiple[s2q_s_anchor_q]
    multiple_anchor[is.na(multiple_anchor)] <- FALSE
    single_projection <- anchor_ls$anchor$query_anchor[s2q_s_anchor_q]
    single_projection[multiple_anchor] <- NA
    s2q_dist <- abs(rbh$s2q_q_location - single_projection)
    single_projection_chr <- anchor_ls$anchor$qchr[s2q_s_anchor_q]
    s2q_dist_chr_unmatch <- rbh$qchr != single_projection_chr
    s2q_dist[s2q_dist_chr_unmatch] <- Inf
    multiple_projection_rbh <- rbh[multiple_anchor, ]
    s2q_dist_multi <- vapply(X = seq_along(multiple_projection_rbh$s2q_s_anchor),
                             FUN.VALUE = numeric(1), FUN = .solveMultiProjectionS2Q,
                             rbh = multiple_projection_rbh, anchor_ls = anchor_ls)
    s2q_dist[multiple_anchor] <- s2q_dist_multi
    s2q_dist[is.na(s2q_dist)] <- Inf
    return(s2q_dist)
}

.solveMultiProjectionQ2S <- function(index, rbh, anchor_ls){
    q_anchor <- rbh$q2s_q_anchor[index]
    hit <- which(anchor_ls$anchor$query_anchor == q_anchor)
    projection <- anchor_ls$anchor$subject_anchor[hit]
    q2s_dist <- abs(rbh$q2s_s_location[index] - projection)
    projection_chr <- anchor_ls$anchor$schr[hit]
    q2s_dist_chr_unmatch <- rbh$schr[index] != projection_chr
    q2s_dist[q2s_dist_chr_unmatch] <- Inf
    out <- min(q2s_dist)
    return(out)
}

.solveMultiProjectionS2Q <- function(index, rbh, anchor_ls){
    s_anchor <- rbh$s2q_s_anchor[index]
    hit <- which(anchor_ls$anchor$subject_anchor == s_anchor)
    projection <- anchor_ls$anchor$query_anchor[hit]
    s2q_dist <- abs(rbh$s2q_q_location[index] - projection)
    projection_chr <- anchor_ls$anchor$qchr[hit]
    s2q_dist_chr_unmatch <- rbh$qchr[index] != projection_chr
    s2q_dist[s2q_dist_chr_unmatch] <- Inf
    out <- min(s2q_dist)
    return(out)
}

#' @import ggplot2
.getMetrics <- function(rbh_dist, tp_list, fp_list){
    n_list <- unique(c(fp_list$query, fp_list$subject))
    n_hit <- rbh_dist$qgeneid %in% n_list | rbh_dist$sgeneid %in% n_list
    n_hit_index_dist <- rbh_dist$dist[n_hit]
    n_hit_index_dist <- tapply(n_hit_index_dist,
                               rbh_dist$qgeneid[n_hit],
                               min)
    n_tbl <- table(n_hit_index_dist)
    n_tbl <- n_tbl[names(n_tbl) != "Inf"]

    p_id <- paste(tp_list$query, tp_list$subject, sep = "_")
    p_id <- unique(p_id)
    p_hit <- which(rbh_dist$pair_id %in% p_id)
    p_hit_index_dist <- rbh_dist$dist[p_hit]
    p_hit_index_dist <- tapply(p_hit_index_dist,
                               rbh_dist$pair_id[p_hit],
                               min)
    p_tbl <- table(p_hit_index_dist)
    p_tbl <- p_tbl[names(p_tbl) != "Inf"]

    non_inf_dist <- rbh_dist$dist[!is.infinite(rbh_dist$dist)]
    index <- seq(0, max(non_inf_dist, na.rm = TRUE))
    tp <- fp <- rep(0, length(index))
    fp[index %in% as.numeric(names(n_tbl))] <- n_tbl
    tp[index %in% as.numeric(names(p_tbl))] <- p_tbl
    fp <- cumsum(fp)
    tp <- cumsum(tp)
    fn <- nrow(tp_list) - tp
    tn <- nrow(fp_list) - fp
    precision <- tp / (tp + fp)
    recall <- tp / (tp + fn)
    specificity <- tn / (tn + fp)
    accuracy <- (tp + tn) / (tp + tn + fp + fn)
    f_measure <- (2 * recall * precision) / (recall + precision)

    df <- rbind(data.frame(Metrics = "precision",
                           Threshold = seq_along(precision),
                           Score = precision),
                data.frame(Metrics = "recall",
                           Threshold = seq_along(recall),
                           Score = recall),
                data.frame(Metrics = "specificity",
                           Threshold = seq_along(specificity),
                           Score = specificity),
                data.frame(Metrics = "accuracy",
                           Threshold = seq_along(accuracy),
                           Score = accuracy),
                data.frame(Metrics = "f_measure",
                           Threshold = seq_along(f_measure),
                           Score = f_measure))
    p <- ggplot(data = df) +
        geom_line(aes(x = Threshold, y = Score, group = Metrics, color = Metrics)) +
        xlim(0, 100) +
        scale_color_manual(values = c("blue",
                                      "darkgreen",
                                      "magenta",
                                      "skyblue",
                                      "darkorange"),
                           breaks = c("precision",
                                      "recall",
                                      "specificity",
                                      "accuracy",
                                      "f_measure"))
    out <- list()
    out$so_metrics <- data.frame(threshold = seq_along(precision),
                                 precision = precision,
                                 recall = recall,
                                 specificity = specificity,
                                 accuracy = accuracy,
                                 f_measure = f_measure)
    out$metrics_plot <- p
    out$so_threshold <- which.max(f_measure)
    return(out)
}

.makeSynogDF <- function(rbh, h5){
    out <- subset(rbh, subset = so_valid %in% TRUE,
                  select = -c(so_valid, qchr, schr, q2s_q_anchor,
                              s2q_q_location, s2q_s_anchor, q2s_s_location,
                              dist))
    out$syntenic <- TRUE
    out$rbbh <- FALSE
    out$pair_id <- paste(out$qseqid, out$sseqid, sep = "_")
    rbbh_id <- paste(h5$anchor$qseqid,
                     h5$anchor$sseqid, sep = "_")
    out$rbbh[out$pair_id %in% rbbh_id] <- TRUE
    return(out)
}

.mergeSynogDF <- function(syn_og, rbbh_og, gff_ls){
    if(nrow(rbbh_og) != 0){
        rbbh_og <- .setGeneIDsynog(rbbh_og = rbbh_og, gff_ls = gff_ls)
        rbbh_og$syntenic <- FALSE
        rbbh_og$rbbh <- TRUE
        rbbh_og$pair_id <- paste(rbbh_og$qseqid, rbbh_og$sseqid, sep = "_")
        rbbh_og <- .getMutualScore(ortho = rbbh_og)
        rbbh_og$class <- .classifySO(ortho = rbbh_og)
        syn_og <- rbind(syn_og, rbbh_og)
    }
    return(syn_og)
}

.makeOutput <- function(syn_og, orphan){
    out <- NULL
    out$syn_og <- syn_og
    out$syn_og$OG <- .numberingOG(ortho = syn_og)
    out$syn_og <- subset(out$syn_og, select = c(qseqid, sseqid, mutual_e,
                                                mutual_ci, syntenic, rbbh,
                                                class, OG, qgeneid, sgeneid))
    rownames(out$syn_og) <- NULL
    out$syn_og$syntenic <- as.character(out$syn_og$syntenic)
    out$syn_og$rbbh <- as.character(out$syn_og$rbbh)

    q_id <- tapply(out$syn_og$qseqid, out$syn_og$class, unique)
    n_q_id <- sapply(q_id, length)
    s_id <- tapply(out$syn_og$sseqid, out$syn_og$class, unique)
    n_s_id <- sapply(s_id, length)

    n_orphan <- lapply(orphan, nrow)

    q_tot <- sum(n_q_id) + n_orphan$query
    s_tot <- sum(n_s_id) + n_orphan$subject
    df <- data.frame(query = c(q_tot, sum(n_q_id), n_q_id, n_orphan$query),
                     subject = c(s_tot, sum(n_s_id), n_s_id, n_orphan$subject))
    out$summary <- df
    out$orphan <- orphan
    rownames(out$orphan$query) <- NULL
    rownames(out$orphan$subject) <- NULL
    return(out)
}

.classifySO <- function(ortho){
    tmp <- subset(ortho, select = c(qseqid, sseqid, pair_id))
    q_dup <- duplicated(tmp$qseqid)
    q_dup <- tmp$qseqid %in% tmp$qseqid[q_dup]
    s_dup <- duplicated(tmp$sseqid)
    s_dup <- tmp$sseqid %in% tmp$sseqid[s_dup]
    og_11 <- tmp[!(q_dup | s_dup), ]

    mult <- tmp[q_dup | s_dup, ]
    mult_q_dup <- mult$qseqid %in% mult$qseqid[duplicated(mult$qseqid)]
    mult_s_dup <- mult$sseqid %in% mult$sseqid[duplicated(mult$sseqid)]

    og_mm <- mult[mult_q_dup & mult_s_dup, ]
    og_mm_i <- mult$qseqid %in% og_mm$qseqid | mult$sseqid %in% og_mm$sseqid
    og_mm <- mult[og_mm_i, ]

    not_mm <- mult[!og_mm_i, ]
    not_mm_q_dup <- not_mm$qseqid %in% not_mm$qseqid[duplicated(not_mm$qseqid)]
    not_mm_s_dup <- not_mm$sseqid %in% not_mm$sseqid[duplicated(not_mm$sseqid)]
    og_1m <- not_mm[not_mm_q_dup, ]
    og_m1 <- not_mm[not_mm_s_dup, ]

    out <- og_11
    out$class <- "1to1"
    og_1m$class <-"1toM"
    out <- rbind(out, og_1m)
    og_m1$class <-"Mto1"
    out <- rbind(out, og_m1)
    og_mm$class <-"MtoM"
    out <- rbind(out, og_mm)

    hit <- match(ortho$pair_id, out$pair_id)
    return(out$class[hit])
}

.getOrphan <- function(gff_ls, syn_og){
    tx_i <- gff_ls$query_gff$type %in% c("transcript", "mRNA")
    tx <- gff_ls$query_gff[tx_i]
    hit <- tx$ID %in% syn_og$qseqid
    q_orphan <- data.frame(ID = tx$ID[!hit], gene_id = tx$gene_id[!hit])

    tx_i <- gff_ls$subject_gff$type %in% c("transcript", "mRNA")
    tx <- gff_ls$subject_gff[tx_i]
    hit <- tx$ID %in% syn_og$sseqid
    s_orphan <- data.frame(ID = tx$ID[!hit], gene_id = tx$gene_id[!hit])
    return(list(query = q_orphan, subject = s_orphan))
}

.getMutualScore <- function(ortho){
    q2s_covident <- ortho$q2s_pident * ortho$q2s_qcovs * 1e-2
    s2q_covident <- ortho$s2q_pident * ortho$s2q_qcovs * 1e-2
    mutual_e <- sqrt(ortho$q2s_evalue^2 + ortho$s2q_evalue^2)
    mutual_ci <- sqrt(q2s_covident^2 + s2q_covident^2)
    ortho$mutual_e <- mutual_e
    ortho$mutual_ci <- mutual_ci
    return(ortho)
}

.addRBBH <- function(rbbh, syn_og){
    candidate_rbbh <- !(rbbh$qseqid %in% syn_og$qseqid |
                            rbbh$sseqid %in% syn_og$sseqid)
    out <- rbbh[candidate_rbbh, ]
    return(out)
}

.filterOG <- function(target, ortho, pident, evalue, qcovs){
    ortho_id <- paste(ortho$qseqid, ortho$sseqid, sep = "_")
    target_id <- paste(target$qseqid, target$sseqid, sep = "_")
    target_metrics <- ortho[ortho_id %in% target_id, ]
    valid_og <- target_metrics$q2s_pident >= pident &
        target_metrics$s2q_pident >= pident &
        target_metrics$q2s_qcovs >= qcovs &
        target_metrics$s2q_qcovs >= qcovs &
        target_metrics$q2s_evalue <= evalue &
        target_metrics$s2q_evalue <= evalue
    valid_og <- target_metrics[valid_og, ]
    valid_id <- paste(valid_og$qseqid, valid_og$sseqid, sep = "_")
    valid_og <- target[target_id %in% valid_id, ]
    return(valid_og)
}

.setGeneIDsynog <- function(rbbh_og, gff_ls){
    q_hit <- match(rbbh_og$qseqid, gff_ls$query_gff$ID)
    rbbh_og$qgeneid <- gff_ls$query_gff$gene_id[q_hit]
    s_hit <- match(rbbh_og$sseqid, gff_ls$subject_gff$ID)
    rbbh_og$sgeneid <- gff_ls$subject_gff$gene_id[s_hit]
    return(rbbh_og)
}

.numberingOG <- function(ortho){
    tmp <- subset(ortho, select = c(qseqid, sseqid, pair_id))
    q_dup <- duplicated(tmp$qseqid)
    q_dup <- tmp$qseqid %in% tmp$qseqid[q_dup]
    s_dup <- duplicated(tmp$sseqid)
    s_dup <- tmp$sseqid %in% tmp$sseqid[s_dup]
    og_11 <- tmp[!(q_dup | s_dup), ]
    og_11 <- og_11[order(og_11$qseqid), ]
    og_11$OG <- seq_along(og_11$qseqid)

    mult <- tmp[q_dup | s_dup, ]
    mult_q_dup <- mult$qseqid %in% mult$qseqid[duplicated(mult$qseqid)]
    mult_s_dup <- mult$sseqid %in% mult$sseqid[duplicated(mult$sseqid)]

    og_mm <- mult[mult_q_dup & mult_s_dup, ]
    og_mm_i <- mult$qseqid %in% og_mm$qseqid | mult$sseqid %in% og_mm$sseqid
    og_mm <- mult[og_mm_i, ]

    not_mm <- mult[!og_mm_i, ]
    not_mm_q_dup <- not_mm$qseqid %in% not_mm$qseqid[duplicated(not_mm$qseqid)]
    not_mm_s_dup <- not_mm$sseqid %in% not_mm$sseqid[duplicated(not_mm$sseqid)]
    og_1m <- not_mm[not_mm_q_dup, ]
    og_m1 <- not_mm[not_mm_s_dup, ]

    og_1m <- og_1m[order(og_1m$qseqid), ]
    og_1m_index <- data.frame(id = unique(og_1m$qseqid),
                              OG = seq_along(unique(og_1m$qseqid)))
    og_1m$OG <- og_1m_index$OG[match(og_1m$qseqid, og_1m_index$id)]

    og_m1 <- og_m1[order(og_m1$qseqid), ]
    og_m1_index <- data.frame(id = unique(og_m1$sseqid),
                              OG = seq_along(unique(og_m1$sseqid)))
    og_m1$OG <- og_m1_index$OG[match(og_m1$sseqid, og_m1_index$id)]

    og_mm <- og_mm[order(og_mm$qseqid), ]
    og_mm$OG <- seq_along(og_mm$qseqid)
    og_mm_q_min <- tapply(og_mm$OG, og_mm$qseqid, min)
    og_mm$q_min <- og_mm_q_min[match(og_mm$qseqid, names(og_mm_q_min))]
    og_mm_s_min <- tapply(og_mm$q_min, og_mm$sseqid, min)
    og_mm$s_min <- og_mm_s_min[match(og_mm$sseqid, names(og_mm_s_min))]
    og_mm_q_min <- tapply(og_mm$s_min, og_mm$qseqid, min)
    og_mm$OG <- og_mm_q_min[match(og_mm$qseqid, names(og_mm_q_min))]
    og_mm <- subset(og_mm, select = -c(q_min, s_min))

    out <- og_11
    out$class <- "1to1"
    og_1m$OG <- og_1m$OG + max(out$OG)
    og_1m$class <-"1toM"
    out <- rbind(out, og_1m)
    og_m1$OG <- og_m1$OG + max(out$OG)
    og_m1$class <-"Mto1"
    out <- rbind(out, og_m1)
    og_mm$OG <- og_mm$OG + max(out$OG)
    og_mm$class <-"MtoM"
    out <- rbind(out, og_mm)

    hit <- match(ortho$pair_id, out$pair_id)
    return(out$OG[hit])
}

#'
#'
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
geneOrtho <- function(object){
    stopifnot(inherits(x = object, "SynogDB"))
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    if(!H5Lexists(h5, "synog_tx/orthopairs")){
        stop("Run syntenyOrtho to obtain ortholog anchors info.")
    }
    ortho <- h5$synog_tx$orthopairs
    ortho <- .findBestPair(object = object, ortho = ortho)
    ortho$pair_id <- paste(ortho$qgeneid, ortho$sgeneid, sep = "_")
    ortho$genewise_class <- .classifySOgene(ortho = ortho)
    omit_id <- .find1to1gene(target = ortho)
    ortho <- subset(ortho, subset = !ortho$pair_id %in% omit_id)
    class_og <- .numberingOGgene(ortho = ortho)
    ortho$genewise_class <- class_og$class
    ortho$genewise_OG <- class_og$OG
    ortho <- ortho[order(ortho$genewise_OG, ortho$qgeneid, ortho$sgeneid), ]
    rownames(ortho) <- NULL
    orphan <- .getOrphanGene(ortho = ortho, tx_orphan = h5$synog_tx$orphan)
    gene_summary <- .makeSummary(ortho = ortho, orphan = orphan)
    ortho$pair_id <- NULL
    .h5creategroup(object$h5,"synog_gene")
    .h5overwrite(obj = ortho, file = object$h5, "synog_gene/orthopairs")
    .h5overwrite(obj = gene_summary, file = object$h5, "synog_gene/summary")
    .h5overwrite(obj = orphan, file = object$h5, "synog_gene/orphan")
}

.classifySOgene <- function(ortho){
    tmp <- subset(ortho, select = c(qgeneid, sgeneid, pair_id))
    tmp <- unique(tmp)
    q_dup <- duplicated(tmp$qgeneid)
    q_dup <- tmp$qgeneid %in% tmp$qgeneid[q_dup]
    s_dup <- duplicated(tmp$sgeneid)
    s_dup <- tmp$sgeneid %in% tmp$sgeneid[s_dup]
    og_11 <- tmp[!(q_dup | s_dup), ]

    mult <- tmp[q_dup | s_dup, ]
    mult_q_dup <- mult$qgeneid %in% mult$qgeneid[duplicated(mult$qgeneid)]
    mult_s_dup <- mult$sgeneid %in% mult$sgeneid[duplicated(mult$sgeneid)]

    og_mm <- mult[mult_q_dup & mult_s_dup, ]
    og_mm_i <- mult$qgeneid %in% og_mm$qgeneid | mult$sgeneid %in% og_mm$sgeneid
    og_mm <- mult[og_mm_i, ]

    not_mm <- mult[!og_mm_i, ]
    not_mm_q_dup <- not_mm$qgeneid %in% not_mm$qgeneid[duplicated(not_mm$qgeneid)]
    not_mm_s_dup <- not_mm$sgeneid %in% not_mm$sgeneid[duplicated(not_mm$sgeneid)]
    og_1m <- not_mm[not_mm_q_dup, ]
    og_m1 <- not_mm[not_mm_s_dup, ]

    out <- og_11
    out$class <- "1to1"
    og_1m$class <-"1toM"
    out <- rbind(out, og_1m)
    og_m1$class <-"Mto1"
    out <- rbind(out, og_m1)
    og_mm$class <-"MtoM"
    out <- rbind(out, og_mm)

    hit <- match(ortho$pair_id, out$pair_id)

    return(out$class[hit])
}

.find1to1gene <- function(target){
    og_target <- target[target$genewise_class == "MtoM", ]
    score_index <- seq_along(og_target$qgeneid)
    q2s <- tapply(score_index, og_target$qgeneid, function(i){
        min_e <- min(og_target$mutual_e[i])
        min_e_i <- which(og_target$mutual_e[i] == min_e)
        max_ci <- max(og_target$mutual_ci[i][min_e_i])
        max_ci_i <- og_target$mutual_ci[i][min_e_i] == max_ci
        best_i <- min_e_i[max_ci_i]
        return(score_index[i][best_i])
    })
    q2s <- sapply(seq_along(q2s), function(i){
        out <- q2s[[i]]
        names(out) <- rep(names(q2s)[i], length(out))
        return(out)
    })
    q2s <- unlist(q2s)
    q2s <- data.frame(qgeneid = names(q2s), sgeneid = og_target$sgeneid[q2s])
    s2q <-  tapply(score_index, og_target$sgeneid, function(i){
        min_e <- min(og_target$mutual_e[i])
        min_e_i <- which(og_target$mutual_e[i] == min_e)
        max_ci <- max(og_target$mutual_ci[i][min_e_i])
        max_ci_i <- og_target$mutual_ci[i][min_e_i] == max_ci
        best_i <- min_e_i[max_ci_i]
        return(score_index[i][best_i])
    })
    s2q <- sapply(seq_along(s2q), function(i){
        out <- s2q[[i]]
        names(out) <- rep(names(s2q)[i], length(out))
        return(out)
    })
    s2q <- unlist(s2q)
    s2q <- data.frame(qgeneid = og_target$qgeneid[s2q], sgeneid = names(s2q))

    add_1to1 <- unique(rbind(q2s, s2q))
    add_1to1$pair_id <- paste(add_1to1$qgeneid, add_1to1$sgeneid, sep = "_")

    add_1to1_index <- which(og_target$pair_id %in% add_1to1$pair_id)
    q_hit <- og_target$qgeneid %in% og_target$qgeneid[add_1to1_index]
    s_hit <- og_target$sgeneid %in% og_target$sgeneid[add_1to1_index]
    out <- og_target$pair_id[q_hit | s_hit]
    out <- out[!out %in% add_1to1$pair_id]
    return(out)
}

.numberingOGgene <- function(ortho){
    tmp <- subset(ortho, select = c(qgeneid, sgeneid, pair_id))
    tmp <- unique(tmp)
    q_dup <- duplicated(tmp$qgeneid)
    q_dup <- tmp$qgeneid %in% tmp$qgeneid[q_dup]
    s_dup <- duplicated(tmp$sgeneid)
    s_dup <- tmp$sgeneid %in% tmp$sgeneid[s_dup]
    og_11 <- tmp[!(q_dup | s_dup), ]
    og_11$OG <- seq_along(og_11$qgeneid)

    mult <- tmp[q_dup | s_dup, ]
    mult_q_dup <- mult$qgeneid %in% mult$qgeneid[duplicated(mult$qgeneid)]
    mult_s_dup <- mult$sgeneid %in% mult$sgeneid[duplicated(mult$sgeneid)]

    og_mm <- mult[mult_q_dup & mult_s_dup, ]
    og_mm_i <- mult$qgeneid %in% og_mm$qgeneid | mult$sgeneid %in% og_mm$sgeneid
    og_mm <- mult[og_mm_i, ]

    not_mm <- mult[!og_mm_i, ]
    not_mm_q_dup <- not_mm$qgeneid %in% not_mm$qgeneid[duplicated(not_mm$qgeneid)]
    not_mm_s_dup <- not_mm$sgeneid %in% not_mm$sgeneid[duplicated(not_mm$sgeneid)]
    og_1m <- not_mm[not_mm_q_dup, ]
    og_m1 <- not_mm[not_mm_s_dup, ]

    og_1m <- og_1m[order(og_1m$qgeneid), ]
    og_1m_index <- data.frame(id = unique(og_1m$qgeneid),
                              OG = seq_along(unique(og_1m$qgeneid)))
    og_1m$OG <- og_1m_index$OG[match(og_1m$qgeneid, og_1m_index$id)]

    og_m1 <- og_m1[order(og_m1$qgeneid), ]
    og_m1_index <- data.frame(id = unique(og_m1$sgeneid),
                              OG = seq_along(unique(og_m1$sgeneid)))
    og_m1$OG <- og_m1_index$OG[match(og_m1$sgeneid, og_m1_index$id)]

    og_mm <- og_mm[order(og_mm$qgeneid), ]
    og_mm$OG <- seq_along(og_mm$qgeneid)
    og_mm_q_min <- tapply(og_mm$OG, og_mm$qgeneid, min)
    og_mm$q_min <- og_mm_q_min[match(og_mm$qgeneid, names(og_mm_q_min))]
    og_mm_s_min <- tapply(og_mm$q_min, og_mm$sgeneid, min)
    og_mm$s_min <- og_mm_s_min[match(og_mm$sgeneid, names(og_mm_s_min))]
    og_mm_q_min <- tapply(og_mm$s_min, og_mm$qgeneid, min)
    og_mm$OG <- og_mm_q_min[match(og_mm$qgeneid, names(og_mm_q_min))]
    og_mm <- subset(og_mm, select = -c(q_min, s_min))

    out <- og_11
    out$class <- "1to1"
    og_1m$class <-"1toM"
    og_1m$OG <- og_1m$OG + max(out$OG)
    out <- rbind(out, og_1m)
    og_m1$class <-"Mto1"
    og_m1$OG <- og_m1$OG + max(out$OG)
    out <- rbind(out, og_m1)
    og_mm$class <-"MtoM"
    og_mm$OG <- og_mm$OG + max(out$OG)
    out <- rbind(out, og_mm)

    hit <- match(ortho$pair_id, out$pair_id)

    return(list(class = out$class[hit], OG = out$OG[hit]))
}

.getOrphanGene <- function(ortho, tx_orphan){
    q_genes <- unique(tx_orphan$query$gene_id)
    q_orphan <- q_genes[!q_genes %in% ortho$qgeneid]
    q_orphan <- data.frame(qgeneid = q_orphan)
    s_genes <- unique(tx_orphan$subject$gene_id)
    s_orphan <- s_genes[!s_genes %in% ortho$sgeneid]
    s_orphan <- data.frame(sgeneid = s_orphan)
    return(list(query = q_orphan, subject = s_orphan))
}

.makeSummary <- function(ortho, orphan){
    q_id <- tapply(ortho$qgeneid, ortho$genewise_class, unique)
    n_q_id <- sapply(q_id, length)
    s_id <- tapply(ortho$sgeneid, ortho$genewise_class, unique)
    n_s_id <- sapply(s_id, length)

    n_orphan <- lapply(orphan, nrow)

    q_tot <- sum(n_q_id) + n_orphan$query
    s_tot <- sum(n_s_id) + n_orphan$subject
    out <- data.frame(query = c(q_tot, sum(n_q_id), n_q_id, n_orphan$query),
                      subject = c(s_tot, sum(n_s_id), n_s_id, n_orphan$subject))
    return(out)
}

#' @importFrom  BiocGenerics width
.findBestPair <- function(object, ortho){
    pair_id <- paste(ortho$qgeneid, ortho$sgeneid, sep = "_")

    q_gff <- .importAllGFF(object$query_gff)
    s_gff <- .importAllGFF(object$subject_gff)
    q_gff <- q_gff[q_gff$type %in% c("mRNA", "transcript")]
    s_gff <- s_gff[s_gff$type %in% c("mRNA", "transcript")]
    q_len <- width(q_gff[match(ortho$qseqid, q_gff$ID)])
    s_len <- width(s_gff[match(ortho$sseqid, s_gff$ID)])

    best_pair <- tapply(seq_along(pair_id), pair_id, function(index){
        check_rbbh <- ortho$rbbh[index] == "TRUE"
        if(sum(check_rbbh) == 1){
            return(index[check_rbbh])

        } else if(sum(check_rbbh) > 1){
            index <- index[check_rbbh]
        }

        check_e <- ortho$mutual_e[index] == min(ortho$mutual_e[index])
        if(sum(check_e) == 1){
            return(index[check_e])

        } else if(sum(check_e) > 1){
            index <- index[check_e]
        }

        check_ci <- ortho$mutual_ci[index] == max(ortho$mutual_ci[index])
        if(sum(check_ci) == 1){
            return(index[check_ci])

        } else if(sum(check_ci) > 1){
            index <- index[check_ci]
        }

        if(length(index) == 1){
            return(index)
        }

        len <- q_len[index] + s_len[index]
        check_len <- len == max(len)
        if(sum(check_len) == 1){
            return(index[check_len])

        } else {
            return(sample(index, 1))
        }
    })
    return(ortho[best_pair, ])
}


################################################################################
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
splitGenes <- function(object){
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    if(!H5Lexists(h5, "synog_gene/orthopairs")){
        stop("Run geneOrtho to obtain genewise ortholog info.")
    }

    ortho <- h5$synog_gene$orthopairs
    ortho$pair_id1 <- paste(ortho$qseqid, ortho$sgeneid, sep = "_")
    ortho$pair_id2 <- paste(ortho$qgeneid, ortho$sseqid, sep = "_")
    ortho$pair_id3 <- paste(ortho$qseqid, ortho$sseqid, sep = "_")

    gff <- .prepGFF(object)

    ortho$genewise_class <- .eval1toM(ortho = ortho, gff = gff)
    ortho$genewise_class <- .evalMto1(ortho = ortho, gff = gff)
    ortho$genewise_class <- .evalMtoM(ortho = ortho, gff = gff)

    new_id <- .replaceIDs(ortho = ortho)
    ortho$qgeneid <- new_id$qgeneid
    ortho$sgeneid <- new_id$sgeneid
    ortho$pair_id <- paste(ortho$qgeneid, ortho$sgeneid, sep = "_")

    class_og <- .numberingOGgene(ortho = ortho)
    ortho$genewise_class <- class_og$class
    ortho$genewise_OG <- class_og$OG
    ortho <- ortho[order(ortho$genewise_OG, ortho$qgeneid, ortho$sgeneid), ]
    rownames(ortho) <- NULL
    orphan <- h5$synog_gene$orphan
    gene_summary <- .makeSummary(ortho = ortho, orphan = orphan)
    ortho$pair_id <- ortho$pair_id1 <- ortho$pair_id2 <- ortho$pair_id3 <- NULL

    .h5creategroup(object$h5,"synog_gene_split")
    .h5overwrite(obj = ortho, file = object$h5, "synog_gene_split/orthopairs")
    .h5overwrite(obj = gene_summary, file = object$h5, "synog_gene_split/summary")
    .h5overwrite(obj = orphan, file = object$h5, "synog_gene_split/orphan")
}

.prepGFF <- function(object){
    query_gff <- .importAllGFF(object$query_gff)
    subject_gff <- .importAllGFF(object$subject_gff)
    query_gff <- query_gff[query_gff$type == "CDS"]
    query_gff$Parent <- unlist(query_gff$Parent)
    subject_gff <- subject_gff[subject_gff$type == "CDS"]
    subject_gff$Parent <- unlist(subject_gff$Parent)
    return(list(query_gff = query_gff, subject_gff = subject_gff))
}

.eval1toM <- function(ortho, gff = gff){
    out <- ortho$genewise_class
    tmp <- subset(ortho, select = c(qgeneid, qseqid, sgeneid, pair_id1),
                  subset = genewise_class == "1toM")

    og_11 <- .find1to1split(seqid = tmp$qseqid, geneid = tmp$sgeneid,
                            pair_id = tmp$pair_id1, gff = gff$query_gff)
    out[ortho$pair_id1 %in% og_11] <- "split_1toM"

    rest <- tmp[!tmp$pair_id1 %in% og_11, ]
    og_11 <- .pick1to1rest(rest = rest, pair_id = rest$pair_id1)
    out[ortho$pair_id1 %in% og_11] <- "split_1toM"

    return(out)
}

.evalMto1 <- function(ortho, gff = gff){
    out <- ortho$genewise_class
    tmp <- subset(ortho, select = c(qgeneid, sseqid, sgeneid, pair_id2),
                  subset = genewise_class == "Mto1")

    og_11 <- .find1to1split(seqid = tmp$sseqid, geneid = tmp$qgeneid,
                            pair_id = tmp$pair_id2, gff = gff$subject_gff)
    out[ortho$pair_id2 %in% og_11] <- "split_Mto1"

    rest <- tmp[!tmp$pair_id2 %in% og_11, ]
    og_11 <- .pick1to1rest(rest = rest, pair_id = rest$pair_id2)
    out[ortho$pair_id2 %in% og_11] <- "split_Mto1"

    return(out)
}

.evalMtoM <- function(ortho, gff = gff){
    out <- ortho$genewise_class
    tmp <- subset(ortho, select = c(qseqid, qgeneid, sseqid, sgeneid,
                                    pair_id1, pair_id2, pair_id3),
                  subset = genewise_class == "MtoM")

    og_11_1 <- .find1to1split(seqid = tmp$qseqid, geneid = tmp$sgeneid,
                              pair_id = tmp$pair_id3, gff = gff$query_gff)
    og_11_2 <- .find1to1split(seqid = tmp$sseqid, geneid = tmp$qgeneid,
                              pair_id = tmp$pair_id3, gff = gff$subject_gff)

    og_11 <- og_11_1[og_11_1 %in% og_11_2]
    out[ortho$pair_id3 %in% og_11] <- "split_MtoM"

    rest <- tmp[!tmp$pair_id3 %in% og_11, ]
    og_11 <- .pick1to1rest(rest = rest, pair_id = rest$pair_id3)
    out[ortho$pair_id3 %in% og_11] <- "split_MtoM"

    return(out)
}

#' @importFrom  GenomicRanges  findOverlaps
.find1to1split <- function(seqid, geneid, pair_id, gff){
    out <- pair_id
    tx_dup <- duplicated(seqid)
    tx_dup <- seqid %in% seqid[tx_dup]
    gene_dup <- duplicated(geneid)
    gene_dup <- geneid %in% geneid[gene_dup]
    seqid_1to1 <- seqid[!(tx_dup | gene_dup)]
    out <- out[!(tx_dup | gene_dup)]

    gff_1to1 <- gff[gff$Parent %in% seqid_1to1]
    ol <- findOverlaps(gff_1to1, gff_1to1)
    ol <- as.data.frame(ol)
    ol <- subset(ol, subset = queryHits != subjectHits)
    hit_id <- unique(gff_1to1$Parent[ol$queryHits])
    seqid_1to1_candidate <- seqid_1to1[!seqid_1to1 %in% hit_id]
    seqid_1to1_rest <- seqid[!seqid %in% seqid_1to1_candidate]
    out <- out[seqid_1to1 %in% seqid_1to1_candidate]

    gff_1to1_candidate <- gff[gff$Parent %in% seqid_1to1_candidate]
    gff_1to1_rest <- gff[gff$Parent %in% seqid_1to1_rest]
    ol <- findOverlaps(gff_1to1_candidate, gff_1to1_rest)
    invalid_seqid <- gff_1to1_candidate$Parent[queryHits(ol)]
    out <- out[!seqid_1to1_candidate %in% invalid_seqid]
    return(out)
}

.pick1to1rest <- function(rest, pair_id){
    out <- pair_id
    q_dup <- duplicated(rest$qgeneid)
    q_dup <- rest$qgeneid %in% rest$qgeneid[q_dup]
    s_dup <- duplicated(rest$sgeneid)
    s_dup <- rest$sgeneid %in% rest$sgeneid[s_dup]
    out <- out[!(q_dup | s_dup)]
    return(out)
}

.replaceIDs <- function(ortho){
    out <- list(qgeneid = ortho$qgeneid, sgeneid = ortho$sgeneid)
    hit <- ortho$genewise_class == "split_1toM"
    out$qgeneid[hit] <- ortho$qseqid[hit]

    hit <- ortho$genewise_class == "split_Mto1"
    out$sgeneid[hit] <- ortho$sseqid[hit]

    hit <- ortho$genewise_class == "split_MtoM"
    out$qgeneid[hit] <- ortho$qseqid[hit]
    out$sgeneid[hit] <- ortho$sseqid[hit]
    return(out)
}
