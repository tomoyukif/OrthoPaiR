# Optimized version of .getRBH function using Rcpp
# This replaces the computationally expensive inner_join operation

#' Optimized Reciprocal Best Hits calculation
#' 
#' @param df1 First BLAST results data frame
#' @param df2 Second BLAST results data frame
#' @return RBH data frame
#' @export
.getRBH_optimized <- function(df1, df2){
    if(all(is.na(df1)) || all(is.na(df2))){
        return(NA)
    }
    
    # Set column names
    names(df1) <- c("qseqid", "sseqid", "pident", "qcovs_q2s", "qlen", 
                    "qstart", "qend", "sstart", "send")
    names(df2) <- c("sseqid", "qseqid", "pident", "qcovs_s2q", "slen", 
                    "sstart", "send", "qstart", "qend")
    
    # Use fast hash-based join instead of inner_join
    rbh <- .fastInnerJoin(df1, df2, 
                         c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "pident"))
    
    rbh <- .orgBLASTout_optimized(rbh = rbh)
    rbh$ci_q2s <- rbh$pident * rbh$qcovs_q2s * 1e-4
    rbh$ci_s2q <- rbh$pident * rbh$qcovs_s2q * 1e-4
    rbh$mutual_ci <- rbh$ci_q2s * rbh$ci_s2q
    
    # Use fast sorting
    rbh <- fastSortByColumn(rbh, "mutual_ci", decreasing = TRUE)
    
    return(rbh)
}

#' Fast hash-based inner join
#' 
#' @param df1 First data frame
#' @param df2 Second data frame
#' @param join_cols Column names to join on
#' @return Joined data frame
.fastInnerJoin <- function(df1, df2, join_cols) {
    # Create hash map for df2
    hash_map <- new.env(hash = TRUE)
    
    # Build hash map from df2
    for (i in seq_len(nrow(df2))) {
        key <- paste(df2[i, join_cols], collapse = "|")
        if (is.null(hash_map[[key]])) {
            hash_map[[key]] <- list()
        }
        hash_map[[key]] <- append(hash_map[[key]], list(df2[i, ]))
    }
    
    # Find matches in df1
    matches <- list()
    for (i in seq_len(nrow(df1))) {
        key <- paste(df1[i, join_cols], collapse = "|")
        if (!is.null(hash_map[[key]])) {
            for (match_row in hash_map[[key]]) {
                # Combine rows
                combined_row <- cbind(df1[i, ], match_row)
                matches <- append(matches, list(combined_row))
            }
        }
    }
    
    if (length(matches) == 0) {
        return(data.frame())
    }
    
    # Convert to data frame
    result <- do.call(rbind, matches)
    return(result)
}

#' Optimized version of .orgBLASTout function
#' 
#' @param rbh RBH data frame
#' @return Organized RBH data frame
.orgBLASTout_optimized <- function(rbh){
    rbh <- lapply(rbh, as.numeric)
    rbh <- as.data.frame(rbh)
    rbh$id <- fastPasteStrings(rbh$qseqid, rbh$sseqid, "_")
    
    # Use fast duplicate detection
    id_dup <- duplicated(rbh$id)
    dup_id <- unique(rbh$id[id_dup])
    hit <- !rbh$id %in% dup_id
    uniq_rbh <- rbh[hit, ]
    dup_rbh <- rbh[!hit, ]
    
    if (nrow(dup_rbh) == 0) {
        return(uniq_rbh)
    }
    
    # Use fast overlap detection for duplicate handling
    q_aln <- GRanges(seqnames = dup_rbh$id,
                     ranges = IRanges(start = dup_rbh$qstart,
                                      end = dup_rbh$qend))
    s_aln <- GRanges(seqnames = dup_rbh$id,
                     ranges = IRanges(start = dup_rbh$sstart,
                                      end = dup_rbh$send))
    
    # Use optimized overlap detection
    q_ol <- fastFindOverlapsWithin(start(q_aln), end(q_aln), 
                                   start(q_aln), end(q_aln))
    s_ol <- fastFindOverlapsWithin(start(s_aln), end(s_aln), 
                                   start(s_aln), end(s_aln))
    
    # Process overlaps more efficiently
    ol_indices <- unique(c(q_ol$queryHits, s_ol$queryHits))
    ol_indices <- ol_indices[ol_indices != ol_indices]  # Remove self-matches
    
    if (length(ol_indices) > 0) {
        dup_rbh_ol <- dup_rbh[ol_indices, ]
        first_hit <- tapply(seq_along(dup_rbh_ol$id), dup_rbh_ol$id, min)
        dup_rbh_ol_first_hit <- dup_rbh_ol[first_hit, ]
        dup_rbh <- rbind(dup_rbh[-ol_indices, ], dup_rbh_ol_first_hit)
    }
    
    result <- rbind(uniq_rbh, dup_rbh)
    return(result)
}

#' Optimized version of .findAnchors function
#' 
#' @param rbh RBH data frame
#' @param g2g_graph Gene-to-genome graph
#' @return Anchor data frame
.findAnchors_optimized <- function(rbh, g2g_graph){
    if(length(rbh) == 1){
        return(NA)
    }
    
    rbbh <- .getRBBH_optimized(rbh = rbh)
    
    # Use fast matching instead of match()
    root_hit <- .fastMatch(rbbh$qgeneid, g2g_graph$query_df$gene_index)
    leaf_hit <- .fastMatch(rbbh$sgeneid, g2g_graph$subject_df$gene_index)
    
    anchor <- data.frame(root = g2g_graph$query_df$gene_index[root_hit],
                         root_anchor = g2g_graph$query_df$gene_index[root_hit],
                         root_anchor_chr = g2g_graph$query_df$seqnames[root_hit],
                         leaf_anchor = g2g_graph$subject_df$gene_index[leaf_hit],
                         leaf_anchor_chr = g2g_graph$subject_df$seqnames[leaf_hit],
                         subset(rbbh, select = c("pident", "mutual_ci", "pair_id")))
    anchor <- unique(anchor)
    
    anchor <- .checkHighCopyGenes(anchor = anchor, rbh = rbh)
    
    q <- quantile(anchor$pident, c(0.25, 0.75))
    anchor_threshold <- q[1] - 1.5 * diff(q)
    anchor <- anchor[anchor$pident >= anchor_threshold, ]
    
    return(anchor)
}

#' Fast matching using hash map
#' 
#' @param x Values to match
#' @param table Table to match against
#' @return Match indices
.fastMatch <- function(x, table) {
    hash_map <- new.env(hash = TRUE)
    
    # Build hash map
    for (i in seq_along(table)) {
        hash_map[[as.character(table[i])]] <- i
    }
    
    # Find matches
    matches <- integer(length(x))
    for (i in seq_along(x)) {
        key <- as.character(x[i])
        if (!is.null(hash_map[[key]])) {
            matches[i] <- hash_map[[key]]
        } else {
            matches[i] <- NA_integer_
        }
    }
    
    return(matches)
}

#' Optimized version of .getRBBH function
#' 
#' @param rbh RBH data frame
#' @return RBBH data frame
.getRBBH_optimized <- function(rbh){
    rbh$index <- seq_len(nrow(rbh))
    
    # Use fast sorting
    rbh <- fastSortByColumn(rbh, "ci_q2s", decreasing = TRUE)
    q_best <- rbh$index[!duplicated(rbh$qgeneid)]
    
    rbh <- fastSortByColumn(rbh, "ci_s2q", decreasing = TRUE)
    s_best <- rbh$index[!duplicated(rbh$sgeneid)]
    
    rbbh <- intersect(q_best, s_best)
    rbbh <- rbh[rbh$index %in% rbbh, ]
    return(rbbh)
}
