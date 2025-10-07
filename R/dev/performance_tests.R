# Performance testing script for OrthoPaiR optimizations
# This script benchmarks the original vs optimized functions

library(OrthoPaiR)
library(microbenchmark)
library(ggplot2)

# Test data generation functions
.generateTestData <- function(n_ranges = 1000, n_sequences = 500) {
    # Generate test genomic ranges
    start1 <- sample(1:1000000, n_ranges, replace = TRUE)
    end1 <- start1 + sample(100:5000, n_ranges, replace = TRUE)
    start2 <- sample(1:1000000, n_ranges, replace = TRUE)
    end2 <- start2 + sample(100:5000, n_ranges, replace = TRUE)
    
    # Generate test sequences
    sequences <- replicate(n_sequences, {
        paste(sample(c("A", "T", "G", "C"), 100, replace = TRUE), collapse = "")
    })
    
    # Generate test data frames
    df1 <- data.frame(
        qseqid = paste0("Q", 1:n_ranges),
        sseqid = paste0("S", 1:n_ranges),
        pident = runif(n_ranges, 50, 100),
        qcovs_q2s = runif(n_ranges, 0, 100),
        qlen = sample(500:5000, n_ranges, replace = TRUE),
        qstart = start1,
        qend = end1,
        sstart = start2,
        send = end2
    )
    
    df2 <- data.frame(
        sseqid = paste0("S", 1:n_ranges),
        qseqid = paste0("Q", 1:n_ranges),
        pident = runif(n_ranges, 50, 100),
        qcovs_s2q = runif(n_ranges, 0, 100),
        slen = sample(500:5000, n_ranges, replace = TRUE),
        sstart = start2,
        send = end2,
        qstart = start1,
        qend = end1
    )
    
    return(list(
        ranges = list(start1 = start1, end1 = end1, start2 = start2, end2 = end2),
        sequences = sequences,
        df1 = df1,
        df2 = df2
    ))
}

# Benchmark overlap detection functions
.benchmarkOverlapDetection <- function(test_data) {
    cat("Benchmarking overlap detection functions...\n")
    
    # Test fast overlap detection
    fast_result <- microbenchmark(
        fastFindOverlaps(test_data$ranges$start1, test_data$ranges$end1,
                         test_data$ranges$start2, test_data$ranges$end2),
        times = 10
    )
    
    # Test fast overlap detection with "within" type
    fast_within_result <- microbenchmark(
        fastFindOverlapsWithin(test_data$ranges$start1, test_data$ranges$end1,
                              test_data$ranges$start2, test_data$ranges$end2),
        times = 10
    )
    
    return(list(
        fast_overlap = fast_result,
        fast_within = fast_within_result
    ))
}

# Benchmark sequence operations
.benchmarkSequenceOperations <- function(test_data) {
    cat("Benchmarking sequence operations...\n")
    
    # Test sequence similarity
    seq1 <- test_data$sequences[1]
    seq2 <- test_data$sequences[2]
    
    similarity_result <- microbenchmark(
        fastSequenceSimilarity(seq1, seq2),
        times = 100
    )
    
    # Test pairwise scoring
    pairwise_result <- microbenchmark(
        fastPairwiseScores(test_data$sequences[1:50]),
        times = 5
    )
    
    # Test grouping
    grouping_result <- microbenchmark(
        fastGroupSimilarSequences(test_data$sequences[1:100], 0.8),
        times = 5
    )
    
    return(list(
        similarity = similarity_result,
        pairwise = pairwise_result,
        grouping = grouping_result
    ))
}

# Benchmark data frame operations
.benchmarkDataFrameOperations <- function(test_data) {
    cat("Benchmarking data frame operations...\n")
    
    # Test duplicate removal
    duplicate_result <- microbenchmark(
        fastRemoveDuplicates(test_data$df1, "qseqid"),
        times = 10
    )
    
    # Test sorting
    sorting_result <- microbenchmark(
        fastSortByColumn(test_data$df1, "pident", decreasing = TRUE),
        times = 10
    )
    
    # Test string operations
    string_result <- microbenchmark(
        fastPasteStrings(test_data$df1$qseqid, test_data$df1$sseqid, "_"),
        times = 10
    )
    
    return(list(
        duplicate = duplicate_result,
        sorting = sorting_result,
        string = string_result
    ))
}

# Benchmark RBH operations
.benchmarkRBHOperations <- function(test_data) {
    cat("Benchmarking RBH operations...\n")
    
    # Test optimized RBH calculation
    rbh_result <- microbenchmark(
        .getRBH_optimized(test_data$df1, test_data$df2),
        times = 5
    )
    
    # Test optimized BLAST output organization
    org_result <- microbenchmark(
        .orgBLASTout_optimized(test_data$df1),
        times = 10
    )
    
    return(list(
        rbh = rbh_result,
        organization = org_result
    ))
}

# Main benchmarking function
runPerformanceTests <- function(n_ranges = 1000, n_sequences = 500) {
    cat("Starting OrthoPaiR performance tests...\n")
    cat("Test parameters: n_ranges =", n_ranges, ", n_sequences =", n_sequences, "\n\n")
    
    # Generate test data
    test_data <- .generateTestData(n_ranges, n_sequences)
    
    # Run benchmarks
    overlap_results <- .benchmarkOverlapDetection(test_data)
    sequence_results <- .benchmarkSequenceOperations(test_data)
    dataframe_results <- .benchmarkDataFrameOperations(test_data)
    rbh_results <- .benchmarkRBHOperations(test_data)
    
    # Compile results
    results <- list(
        overlap = overlap_results,
        sequence = sequence_results,
        dataframe = dataframe_results,
        rbh = rbh_results
    )
    
    # Print summary
    cat("\n=== PERFORMANCE TEST SUMMARY ===\n")
    cat("Overlap Detection (fast):\n")
    print(summary(overlap_results$fast_overlap))
    
    cat("\nSequence Similarity:\n")
    print(summary(sequence_results$similarity))
    
    cat("\nData Frame Operations:\n")
    print(summary(dataframe_results$duplicate))
    
    cat("\nRBH Operations:\n")
    print(summary(rbh_results$rbh))
    
    return(results)
}

# Memory usage testing
testMemoryUsage <- function() {
    cat("Testing memory usage...\n")
    
    # Test with different data sizes
    sizes <- c(100, 500, 1000, 2000)
    memory_results <- list()
    
    for (size in sizes) {
        cat("Testing with", size, "ranges...\n")
        test_data <- .generateTestData(size, size)
        
        # Measure memory before
        mem_before <- gc()
        
        # Run operations
        result <- fastFindOverlaps(test_data$ranges$start1, test_data$ranges$end1,
                                   test_data$ranges$start2, test_data$ranges$end2)
        
        # Measure memory after
        mem_after <- gc()
        
        memory_results[[as.character(size)]] <- list(
            before = mem_before,
            after = mem_after,
            result_size = length(result$queryHits)
        )
    }
    
    return(memory_results)
}

# Scalability testing
testScalability <- function() {
    cat("Testing scalability...\n")
    
    # Test with increasing data sizes
    sizes <- c(100, 500, 1000, 2000, 5000)
    scalability_results <- list()
    
    for (size in sizes) {
        cat("Testing with", size, "ranges...\n")
        test_data <- .generateTestData(size, size)
        
        # Time the operation
        start_time <- Sys.time()
        result <- fastFindOverlaps(test_data$ranges$start1, test_data$ranges$end1,
                                   test_data$ranges$start2, test_data$ranges$end2)
        end_time <- Sys.time()
        
        scalability_results[[as.character(size)]] <- list(
            size = size,
            time = as.numeric(end_time - start_time),
            overlaps = length(result$queryHits)
        )
    }
    
    return(scalability_results)
}

# Export functions
export("runPerformanceTests")
export("testMemoryUsage")
export("testScalability")
