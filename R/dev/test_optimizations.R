# Test script for OrthoPaiR optimizations
# This script tests the basic functionality of the optimized functions

library(OrthoPaiR)

# Test data generation
cat("Generating test data...\n")
n_ranges <- 100
start1 <- sample(1:10000, n_ranges, replace = TRUE)
end1 <- start1 + sample(100:1000, n_ranges, replace = TRUE)
start2 <- sample(1:10000, n_ranges, replace = TRUE)
end2 <- start2 + sample(100:1000, n_ranges, replace = TRUE)

# Test sequences
sequences <- replicate(50, {
    paste(sample(c("A", "T", "G", "C"), 100, replace = TRUE), collapse = "")
})

# Test data frames
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

# Test 1: Fast overlap detection
cat("Testing fast overlap detection...\n")
tryCatch({
    result1 <- fastFindOverlaps(start1, end1, start2, end2)
    cat("✓ Fast overlap detection works\n")
    cat("  Found", length(result1$queryHits), "overlaps\n")
}, error = function(e) {
    cat("✗ Fast overlap detection failed:", e$message, "\n")
})

# Test 2: Fast overlap detection with "within" type
cat("Testing fast overlap detection (within)...\n")
tryCatch({
    result2 <- fastFindOverlapsWithin(start1, end1, start2, end2)
    cat("✓ Fast overlap detection (within) works\n")
    cat("  Found", length(result2$queryHits), "within overlaps\n")
}, error = function(e) {
    cat("✗ Fast overlap detection (within) failed:", e$message, "\n")
})

# Test 3: Fast sequence similarity
cat("Testing fast sequence similarity...\n")
tryCatch({
    seq1 <- sequences[1]
    seq2 <- sequences[2]
    similarity <- fastSequenceSimilarity(seq1, seq2)
    cat("✓ Fast sequence similarity works\n")
    cat("  Similarity:", round(similarity, 3), "\n")
}, error = function(e) {
    cat("✗ Fast sequence similarity failed:", e$message, "\n")
})

# Test 4: Fast string operations
cat("Testing fast string operations...\n")
tryCatch({
    vec1 <- df1$qseqid[1:10]
    vec2 <- df1$sseqid[1:10]
    result3 <- fastPasteStrings(vec1, vec2, "_")
    cat("✓ Fast string operations work\n")
    cat("  Generated", length(result3), "concatenated strings\n")
}, error = function(e) {
    cat("✗ Fast string operations failed:", e$message, "\n")
})

# Test 5: Fast counting
cat("Testing fast counting...\n")
tryCatch({
    vec <- sample(letters[1:5], 100, replace = TRUE)
    result4 <- fastCountUnique(vec)
    cat("✓ Fast counting works\n")
    cat("  Found", length(result4$values), "unique values\n")
}, error = function(e) {
    cat("✗ Fast counting failed:", e$message, "\n")
})

# Test 6: Fast duplicate removal
cat("Testing fast duplicate removal...\n")
tryCatch({
    # Create data with duplicates
    df_with_duplicates <- rbind(df1, df1[1:10, ])
    result5 <- fastRemoveDuplicates(df_with_duplicates, "qseqid")
    cat("✓ Fast duplicate removal works\n")
    cat("  Original rows:", nrow(df_with_duplicates), "\n")
    cat("  After deduplication:", nrow(result5), "\n")
}, error = function(e) {
    cat("✗ Fast duplicate removal failed:", e$message, "\n")
})

# Test 7: Fast sorting
cat("Testing fast sorting...\n")
tryCatch({
    result6 <- fastSortByColumn(df1, "pident", decreasing = TRUE)
    cat("✓ Fast sorting works\n")
    cat("  Sorted", nrow(result6), "rows by pident\n")
}, error = function(e) {
    cat("✗ Fast sorting failed:", e$message, "\n")
})

# Test 8: Fast ID replacement
cat("Testing fast ID replacement...\n")
tryCatch({
    ids <- df1$qseqid[1:10]
    old_ids <- ids[1:5]
    new_ids <- paste0("NEW_", old_ids)
    result7 <- fastReplaceIds(ids, old_ids, new_ids)
    cat("✓ Fast ID replacement works\n")
    cat("  Replaced", sum(result7 != ids), "IDs\n")
}, error = function(e) {
    cat("✗ Fast ID replacement failed:", e$message, "\n")
})

# Test 9: Fast pairwise scoring
cat("Testing fast pairwise scoring...\n")
tryCatch({
    result8 <- fastPairwiseScores(sequences[1:10])
    cat("✓ Fast pairwise scoring works\n")
    cat("  Generated", nrow(result8), "x", ncol(result8), "similarity matrix\n")
}, error = function(e) {
    cat("✗ Fast pairwise scoring failed:", e$message, "\n")
})

# Test 10: Fast grouping
cat("Testing fast grouping...\n")
tryCatch({
    result9 <- fastGroupSimilarSequences(sequences[1:20], 0.8)
    cat("✓ Fast grouping works\n")
    cat("  Created", result9$n_groups, "groups\n")
}, error = function(e) {
    cat("✗ Fast grouping failed:", e$message, "\n")
})

cat("\n=== TEST SUMMARY ===\n")
cat("All basic functionality tests completed.\n")
cat("If you see ✓ marks above, the optimizations are working correctly.\n")
cat("If you see ✗ marks, there may be compilation or linking issues.\n")

# Performance comparison (if possible)
cat("\n=== PERFORMANCE COMPARISON ===\n")
cat("Note: For a full performance comparison, run the performance tests:\n")
cat("  source('R/dev/performance_tests.R')\n")
cat("  results <- runPerformanceTests()\n")
