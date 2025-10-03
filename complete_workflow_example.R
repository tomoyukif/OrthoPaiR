#!/usr/bin/env Rscript
# Complete OrthoPaiR workflow example
# This script demonstrates the full pipeline using subset data

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

# Create output directory
output_dir <- "orthopair_example_output"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

cat("=== OrthoPaiR Complete Workflow Example ===\n")
cat("Using subset data: NB (query) vs WK21 (subject)\n\n")

# Step 1: Fix input files
cat("Step 1: Fixing input files...\n")
nb_fixed <- fixInfiles(genome = nb_genome, 
                      gff = nb_gff, 
                      cds = nb_cds, 
                      prot = nb_prot,
                      autofix = TRUE)

wk21_fixed <- fixInfiles(genome = wk21_genome, 
                        gff = wk21_gff, 
                        cds = wk21_cds, 
                        prot = wk21_prot,
                        autofix = TRUE)

cat("✓ Input files fixed\n\n")

# Step 2: Organize input files
cat("Step 2: Organizing input files...\n")
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

cat("✓ Input files organized\n\n")

# Step 3: Create OrthoPairDB object
cat("Step 3: Creating OrthoPairDB object...\n")
hdf5_path <- file.path(output_dir, "nb_wk21_orthopair.h5")

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

cat("✓ OrthoPairDB object created\n")
cat("  HDF5 file:", hdf5_path, "\n\n")

# Step 4: Run Reciprocal Best Hits (RBH) analysis
cat("Step 4: Running RBH analysis...\n")
rbh_result <- rbh(object, n_threads = 1)
cat("✓ RBH analysis completed\n\n")

# Step 5: Run syntenic ortholog analysis
cat("Step 5: Running syntenic ortholog analysis...\n")
syntenic_result <- syntenicOrtho(object)
cat("✓ Syntenic ortholog analysis completed\n\n")

# Step 6: Get summary statistics
cat("Step 6: Getting summary statistics...\n")
summary_stats <- summaryOrthoPair(object)
print(summary_stats)
cat("\n")

# Step 7: Get ortholog pairs
cat("Step 7: Retrieving ortholog pairs...\n")
orthopairs <- getOrthoPair(object, gene = TRUE)
cat("✓ Retrieved", nrow(orthopairs), "ortholog pairs\n\n")

# Step 8: Get orphan genes
cat("Step 8: Retrieving orphan genes...\n")
orphans <- getOrphan(object)
cat("✓ Retrieved", length(orphans$query), "query orphans and", 
    length(orphans$subject), "subject orphans\n\n")

# Step 9: Create ortholog graph
cat("Step 9: Creating ortholog graph...\n")
graph <- makeOrthoGraph(object)
cat("✓ Ortholog graph created with", vcount(graph), "vertices\n\n")

# Step 10: Convert graph to dataframe
cat("Step 10: Converting graph to dataframe...\n")
graph_df <- graph2df(object, graph)
cat("✓ Graph converted to dataframe\n\n")

# Step 11: Reorganize orthopairs
cat("Step 11: Reorganizing orthopairs...\n")
reorg_result <- reorgOrthopiars(object, out_dir = output_dir, overwrite = TRUE)
cat("✓ Orthopairs reorganized\n\n")

# Step 12: Extract promoter sequences
cat("Step 12: Extracting promoter sequences...\n")
nb_promoter <- .getPromoterSeq(gff = nb_gff, 
                              genome = nb_genome, 
                              out_fn = file.path(output_dir, "nb_promoter.fa"))

wk21_promoter <- .getPromoterSeq(gff = wk21_gff, 
                                genome = wk21_genome, 
                                out_fn = file.path(output_dir, "wk21_promoter.fa"))

cat("✓ Promoter sequences extracted\n\n")

# Step 13: Run sequence comparison (if possible)
cat("Step 13: Running sequence comparison...\n")
tryCatch({
    compare_result <- compareOrthoSeq(object, n_threads = 1, verbose = FALSE)
    cat("✓ Sequence comparison completed\n\n")
}, error = function(e) {
    cat("⚠ Sequence comparison skipped:", e$message, "\n\n")
})

# Step 14: Save results
cat("Step 14: Saving results...\n")

# Save summary statistics
write.csv(summary_stats, file.path(output_dir, "summary_statistics.csv"), row.names = TRUE)

# Save ortholog pairs
write.csv(orthopairs, file.path(output_dir, "ortholog_pairs.csv"), row.names = FALSE)

# Save graph dataframe
write.csv(graph_df, file.path(output_dir, "ortholog_graph.csv"), row.names = FALSE)

# Save orphan genes
write.csv(orphans$query, file.path(output_dir, "query_orphans.csv"), row.names = FALSE)
write.csv(orphans$subject, file.path(output_dir, "subject_orphans.csv"), row.names = FALSE)

cat("✓ Results saved to", output_dir, "\n\n")

# Step 15: Display final summary
cat("=== Final Summary ===\n")
cat("Analysis completed successfully!\n")
cat("Input species: NB (query) vs WK21 (subject)\n")
cat("Total ortholog pairs found:", nrow(orthopairs), "\n")
cat("Query orphans:", length(orphans$query), "\n")
cat("Subject orphans:", length(orphans$subject), "\n")
cat("Graph vertices:", vcount(graph), "\n")
cat("Graph edges:", ecount(graph), "\n")

cat("\nOutput files:\n")
cat("- HDF5 database:", hdf5_path, "\n")
cat("- Summary statistics:", file.path(output_dir, "summary_statistics.csv"), "\n")
cat("- Ortholog pairs:", file.path(output_dir, "ortholog_pairs.csv"), "\n")
cat("- Ortholog graph:", file.path(output_dir, "ortholog_graph.csv"), "\n")
cat("- Query orphans:", file.path(output_dir, "query_orphans.csv"), "\n")
cat("- Subject orphans:", file.path(output_dir, "subject_orphans.csv"), "\n")
cat("- Promoter sequences:", file.path(output_dir, "nb_promoter.fa"), "\n")
cat("- Promoter sequences:", file.path(output_dir, "wk21_promoter.fa"), "\n")

cat("\n=== Workflow completed successfully! ===\n")
