#!/usr/bin/env Rscript
# Simple test script for OrthoPaiR core functions
# This script tests the main functions without external dependencies

library(OrthoPaiR)

# Set working directory
setwd("/home/ftom/01_wd/softDevel/OrthoPaiR")

# Define paths to subset data
subset_dir <- "inst/extdata/subset"
nb_genome <- file.path(subset_dir, "nb_subset_genome.fa")
nb_gff <- file.path(subset_dir, "nb_subset.gff")
nb_cds <- file.path(subset_dir, "nb_subset_cds.fa")
nb_prot <- file.path(subset_dir, "nb_subset_prot.fa")

wk21_genome <- file.path(subset_dir, "wk21_subset_genome.fa")
wk21_gff <- file.path(subset_dir, "wk21_subset.gff")
wk21_cds <- file.path(subset_dir, "wk21_subset_cds.fa")
wk21_prot <- file.path(subset_dir, "wk21_subset_prot.fa")

# Create test output directory
test_output_dir <- "test_output"
if (!dir.exists(test_output_dir)) {
    dir.create(test_output_dir, recursive = TRUE)
}

cat("=== OrthoPaiR Core Function Tests ===\n")
cat("Testing core functions with subset data\n\n")

# Test counter
test_count <- 0
passed_tests <- 0
failed_tests <- 0

# Function to run a test and track results
runTest <- function(test_name, test_function) {
    test_count <<- test_count + 1
    cat("Test", test_count, ":", test_name, "... ")
    
    tryCatch({
        result <- test_function()
        cat("PASSED\n")
        passed_tests <<- passed_tests + 1
        return(result)
    }, error = function(e) {
        cat("FAILED -", e$message, "\n")
        failed_tests <<- failed_tests + 1
        return(NULL)
    })
}

# Test 1: Check if subset files exist
runTest("Check subset files", function() {
    files_to_check <- c(nb_genome, nb_gff, nb_cds, nb_prot, 
                        wk21_genome, wk21_gff, wk21_cds, wk21_prot)
    
    for (file in files_to_check) {
        if (!file.exists(file)) {
            stop("File not found: ", file)
        }
    }
    
    return(TRUE)
})

# Test 2: fixInfiles
runTest("fixInfiles", function() {
    # Test with NB files
    nb_fixed <- fixInfiles(genome = nb_genome, 
                          gff = nb_gff, 
                          cds = nb_cds, 
                          prot = nb_prot,
                          autofix = TRUE)
    
    # Test with WK21 files
    wk21_fixed <- fixInfiles(genome = wk21_genome, 
                            gff = wk21_gff, 
                            cds = wk21_cds, 
                            prot = wk21_prot,
                            autofix = TRUE)
    
    return(list(nb = nb_fixed, wk21 = wk21_fixed))
})

# Test 3: formatGFF
runTest("formatGFF", function() {
    nb_formatted <- formatGFF(nb_gff, suffix = "_formatted.gff")
    wk21_formatted <- formatGFF(wk21_gff, suffix = "_formatted.gff")
    
    return(list(nb = nb_formatted, wk21 = wk21_formatted))
})

# Test 4: orgInputFiles
runTest("orgInputFiles", function() {
    # Create input file list
    in_list <- orgInputFiles(name = "NB",
                            genome = nb_genome,
                            gff = nb_gff,
                            cds = nb_cds,
                            prot = nb_prot)
    
    in_list <- orgInputFiles(object = in_list,
                            name = "WK21",
                            genome = wk21_genome,
                            gff = wk21_gff,
                            cds = wk21_cds,
                            prot = wk21_prot)
    
    return(in_list)
})

# Test 5: makeOrthoPairDB
runTest("makeOrthoPairDB", function() {
    hdf5_path <- file.path(test_output_dir, "test_orthopair.h5")
    
    object <- makeOrthoPairDB(query_genome = nb_genome,
                             subject_genome = wk21_genome,
                             query_gff = nb_gff,
                             subject_gff = wk21_gff,
                             query_cds = nb_cds,
                             subject_cds = wk21_cds,
                             query_prot = nb_prot,
                             subject_prot = wk21_prot,
                             hdf5_path = hdf5_path,
                             overwrite = TRUE)
    
    return(object)
})

# Test 6: rbh (Reciprocal Best Hits)
runTest("rbh", function() {
    hdf5_path <- file.path(test_output_dir, "test_orthopair.h5")
    
    # Create object first
    object <- makeOrthoPairDB(query_genome = nb_genome,
                             subject_genome = wk21_genome,
                             query_gff = nb_gff,
                             subject_gff = wk21_gff,
                             query_cds = nb_cds,
                             subject_cds = wk21_cds,
                             query_prot = nb_prot,
                             subject_prot = wk21_prot,
                             hdf5_path = hdf5_path,
                             overwrite = TRUE)
    
    # Run RBH analysis
    rbh_result <- rbh(object, n_threads = 1)
    
    return(rbh_result)
})

# Test 7: syntenicOrtho
runTest("syntenicOrtho", function() {
    hdf5_path <- file.path(test_output_dir, "test_orthopair.h5")
    
    # Create object and run RBH first
    object <- makeOrthoPairDB(query_genome = nb_genome,
                             subject_genome = wk21_genome,
                             query_gff = nb_gff,
                             subject_gff = wk21_gff,
                             query_cds = nb_cds,
                             subject_cds = wk21_cds,
                             query_prot = nb_prot,
                             subject_prot = wk21_prot,
                             hdf5_path = hdf5_path,
                             overwrite = TRUE)
    
    rbh(object, n_threads = 1)
    syntenic_result <- syntenicOrtho(object)
    
    return(syntenic_result)
})

# Test 8: summaryOrthoPair
runTest("summaryOrthoPair", function() {
    hdf5_path <- file.path(test_output_dir, "test_orthopair.h5")
    
    # Create object and run analysis first
    object <- makeOrthoPairDB(query_genome = nb_genome,
                             subject_genome = wk21_genome,
                             query_gff = nb_gff,
                             subject_gff = wk21_gff,
                             query_cds = nb_cds,
                             subject_cds = wk21_cds,
                             query_prot = nb_prot,
                             subject_prot = wk21_prot,
                             hdf5_path = hdf5_path,
                             overwrite = TRUE)
    
    rbh(object, n_threads = 1)
    syntenicOrtho(object)
    
    summary_result <- summaryOrthoPair(object)
    
    return(summary_result)
})

# Test 9: getOrthoPair
runTest("getOrthoPair", function() {
    hdf5_path <- file.path(test_output_dir, "test_orthopair.h5")
    
    # Create object and run analysis first
    object <- makeOrthoPairDB(query_genome = nb_genome,
                             subject_genome = wk21_genome,
                             query_gff = nb_gff,
                             subject_gff = wk21_gff,
                             query_cds = nb_cds,
                             subject_cds = wk21_cds,
                             query_prot = nb_prot,
                             subject_prot = wk21_prot,
                             hdf5_path = hdf5_path,
                             overwrite = TRUE)
    
    rbh(object, n_threads = 1)
    syntenicOrtho(object)
    
    orthopair_result <- getOrthoPair(object, gene = TRUE)
    
    return(orthopair_result)
})

# Test 10: getOrphan
runTest("getOrphan", function() {
    hdf5_path <- file.path(test_output_dir, "test_orthopair.h5")
    
    # Create object and run analysis first
    object <- makeOrthoPairDB(query_genome = nb_genome,
                             subject_genome = wk21_genome,
                             query_gff = nb_gff,
                             subject_gff = wk21_gff,
                             query_cds = nb_cds,
                             subject_cds = wk21_cds,
                             query_prot = nb_prot,
                             subject_prot = wk21_prot,
                             hdf5_path = hdf5_path,
                             overwrite = TRUE)
    
    rbh(object, n_threads = 1)
    syntenicOrtho(object)
    
    orphan_result <- getOrphan(object)
    
    return(orphan_result)
})

# Test 11: makeOrthoGraph
runTest("makeOrthoGraph", function() {
    hdf5_path <- file.path(test_output_dir, "test_orthopair.h5")
    
    # Create object and run analysis first
    object <- makeOrthoPairDB(query_genome = nb_genome,
                             subject_genome = wk21_genome,
                             query_gff = nb_gff,
                             subject_gff = wk21_gff,
                             query_cds = nb_cds,
                             subject_cds = wk21_cds,
                             query_prot = nb_prot,
                             subject_prot = wk21_prot,
                             hdf5_path = hdf5_path,
                             overwrite = TRUE)
    
    rbh(object, n_threads = 1)
    syntenicOrtho(object)
    
    # Create ortholog graph
    graph_result <- makeOrthoGraph(object)
    
    return(graph_result)
})

# Test 12: graph2df
runTest("graph2df", function() {
    hdf5_path <- file.path(test_output_dir, "test_orthopair.h5")
    
    # Create object and run analysis first
    object <- makeOrthoPairDB(query_genome = nb_genome,
                             subject_genome = wk21_genome,
                             query_gff = nb_gff,
                             subject_gff = wk21_gff,
                             query_cds = nb_cds,
                             subject_cds = wk21_cds,
                             query_prot = nb_prot,
                             subject_prot = wk21_prot,
                             hdf5_path = hdf5_path,
                             overwrite = TRUE)
    
    rbh(object, n_threads = 1)
    syntenicOrtho(object)
    
    # Create graph and convert to dataframe
    graph <- makeOrthoGraph(object)
    graph_df <- graph2df(object, graph)
    
    return(graph_df)
})

# Test 13: reorgOrthopiars
runTest("reorgOrthopiars", function() {
    hdf5_path <- file.path(test_output_dir, "test_orthopair.h5")
    
    # Create object and run analysis first
    object <- makeOrthoPairDB(query_genome = nb_genome,
                             subject_genome = wk21_genome,
                             query_gff = nb_gff,
                             subject_gff = wk21_gff,
                             query_cds = nb_cds,
                             subject_cds = wk21_cds,
                             query_prot = nb_prot,
                             subject_prot = wk21_prot,
                             hdf5_path = hdf5_path,
                             overwrite = TRUE)
    
    rbh(object, n_threads = 1)
    syntenicOrtho(object)
    
    # Reorganize orthopairs
    reorg_result <- reorgOrthopiars(object, out_dir = test_output_dir, overwrite = TRUE)
    
    return(reorg_result)
})

# Test 14: .getPromoterSeq (internal function)
runTest(".getPromoterSeq", function() {
    # Test promoter sequence extraction
    promoter_result <- .getPromoterSeq(gff = nb_gff, 
                                      genome = nb_genome, 
                                      out_fn = file.path(test_output_dir, "nb_promoter.fa"))
    
    return(promoter_result)
})

# Print final results
cat("\n=== Test Results Summary ===\n")
cat("Total tests:", test_count, "\n")
cat("Passed:", passed_tests, "\n")
cat("Failed:", failed_tests, "\n")
cat("Success rate:", round(passed_tests/test_count * 100, 1), "%\n")

if (failed_tests > 0) {
    cat("\nSome tests failed. This may be due to:\n")
    cat("- Insufficient data in subset files\n")
    cat("- System-specific issues\n")
    cat("- Missing dependencies\n")
} else {
    cat("\nAll core tests passed successfully!\n")
}

cat("\nTest output saved in:", test_output_dir, "\n")
cat("HDF5 file:", file.path(test_output_dir, "test_orthopair.h5"), "\n")
