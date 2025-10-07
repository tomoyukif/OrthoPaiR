# Optimized overlap detection functions using Rcpp
# These functions replace the computationally expensive findOverlaps operations

#' Fast overlap detection for genomic ranges
#' 
#' @param start1 Start positions of query ranges
#' @param end1 End positions of query ranges  
#' @param start2 Start positions of subject ranges
#' @param end2 End positions of subject ranges
#' @return List with queryHits and subjectHits
#' @export
fastFindOverlaps <- function(start1, end1, start2, end2) {
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
        stop("Rcpp package is required for fast overlap detection")
    }
    
    # Use Rcpp implementation
    result <- .Call("fast_find_overlaps", start1, end1, start2, end2, PACKAGE = "OrthoPaiR")
    return(result)
}

#' Fast overlap detection with "within" type
#' 
#' @param start1 Start positions of query ranges
#' @param end1 End positions of query ranges
#' @param start2 Start positions of subject ranges  
#' @param end2 End positions of subject ranges
#' @return List with queryHits and subjectHits
#' @export
fastFindOverlapsWithin <- function(start1, end1, start2, end2) {
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
        stop("Rcpp package is required for fast overlap detection")
    }
    
    # Use Rcpp implementation
    result <- .Call("fast_find_overlaps_within", start1, end1, start2, end2, PACKAGE = "OrthoPaiR")
    return(result)
}

#' Fast duplicate removal for BLAST results
#' 
#' @param df Data frame to process
#' @param id_col Name of the ID column
#' @return Data frame with duplicates removed
#' @export
fastRemoveDuplicates <- function(df, id_col) {
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
        stop("Rcpp package is required for fast duplicate removal")
    }
    
    # Use Rcpp implementation
    result <- .Call("fast_remove_duplicates", df, id_col, PACKAGE = "OrthoPaiR")
    return(result)
}

#' Fast sorting by column
#' 
#' @param df Data frame to sort
#' @param sort_col Name of column to sort by
#' @param decreasing Logical, sort in decreasing order
#' @return Sorted data frame
#' @export
fastSortByColumn <- function(df, sort_col, decreasing = TRUE) {
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
        stop("Rcpp package is required for fast sorting")
    }
    
    # Use Rcpp implementation
    result <- .Call("fast_sort_by_column", df, sort_col, decreasing, PACKAGE = "OrthoPaiR")
    return(result)
}

#' Fast sequence similarity calculation
#' 
#' @param seq1 First sequence
#' @param seq2 Second sequence
#' @return Similarity score between 0 and 1
#' @export
fastSequenceSimilarity <- function(seq1, seq2) {
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
        stop("Rcpp package is required for fast sequence operations")
    }
    
    # Use Rcpp implementation
    result <- .Call("fast_sequence_similarity", seq1, seq2, PACKAGE = "OrthoPaiR")
    return(result)
}

#' Fast pairwise sequence scoring
#' 
#' @param sequences Character vector of sequences
#' @return Matrix of pairwise similarity scores
#' @export
fastPairwiseScores <- function(sequences) {
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
        stop("Rcpp package is required for fast sequence operations")
    }
    
    # Use Rcpp implementation
    result <- .Call("fast_pairwise_scores", sequences, PACKAGE = "OrthoPaiR")
    return(result)
}

#' Fast grouping of similar sequences
#' 
#' @param sequences Character vector of sequences
#' @param threshold Similarity threshold for grouping
#' @return List with group IDs and number of groups
#' @export
fastGroupSimilarSequences <- function(sequences, threshold = 0.8) {
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
        stop("Rcpp package is required for fast sequence operations")
    }
    
    # Use Rcpp implementation
    result <- .Call("fast_group_similar_sequences", sequences, threshold, PACKAGE = "OrthoPaiR")
    return(result)
}

#' Fast string concatenation
#' 
#' @param vec1 First character vector
#' @param vec2 Second character vector
#' @param separator Separator string
#' @return Character vector of concatenated strings
#' @export
fastPasteStrings <- function(vec1, vec2, separator = "_") {
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
        stop("Rcpp package is required for fast string operations")
    }
    
    # Use Rcpp implementation
    result <- .Call("fast_paste_strings", vec1, vec2, separator, PACKAGE = "OrthoPaiR")
    return(result)
}

#' Fast ID replacement using hash map
#' 
#' @param ids Character vector of IDs to replace
#' @param old_ids Character vector of old IDs
#' @param new_ids Character vector of new IDs
#' @return Character vector with replaced IDs
#' @export
fastReplaceIds <- function(ids, old_ids, new_ids) {
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
        stop("Rcpp package is required for fast string operations")
    }
    
    # Use Rcpp implementation
    result <- .Call("fast_replace_ids", ids, old_ids, new_ids, PACKAGE = "OrthoPaiR")
    return(result)
}

#' Fast counting of unique values
#' 
#' @param vec Character vector to count
#' @return List with unique values and their counts
#' @export
fastCountUnique <- function(vec) {
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
        stop("Rcpp package is required for fast counting")
    }
    
    # Use Rcpp implementation
    result <- .Call("fast_count_unique", vec, PACKAGE = "OrthoPaiR")
    return(result)
}
