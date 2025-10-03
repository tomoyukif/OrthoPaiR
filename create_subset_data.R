#!/usr/bin/env Rscript
# Script to create subset data files for OrthoPaiR testing
# This script creates small sample datasets from the original input files

library(Biostrings)
library(rtracklayer)
library(GenomicRanges)

# Set working directory
setwd("/home/ftom/01_wd/softDevel/OrthoPaiR")

# Create output directory for subset data
subset_dir <- "inst/extdata/subset"
if (!dir.exists(subset_dir)) {
    dir.create(subset_dir, recursive = TRUE)
}

cat("Creating subset data files for OrthoPaiR testing...\n")

# Function to create subset GFF file
createSubsetGFF <- function(input_gff, output_gff, n_genes = 50) {
    cat("Processing", basename(input_gff), "...\n")
    
    # Read GFF file
    gff <- import.gff3(input_gff)
    
    # Get unique genes
    gene_entries <- gff[gff$type == "gene"]
    if (length(gene_entries) == 0) {
        # If no gene entries, look for other feature types
        gene_entries <- gff[gff$type %in% c("transcript", "mRNA")]
    }
    
    # Select subset of genes
    n_select <- min(n_genes, length(gene_entries))
    selected_genes <- gene_entries[1:n_select]
    
    # Get all features related to selected genes
    selected_gene_ids <- selected_genes$ID
    if (is.null(selected_gene_ids)) {
        selected_gene_ids <- selected_genes$gene_id
    }
    
    # Find all related features
    related_features <- gff[gff$gene_id %in% selected_gene_ids | 
                           gff$ID %in% selected_gene_ids |
                           unlist(gff$Parent) %in% selected_gene_ids]
    
    # Export subset GFF
    export.gff3(related_features, output_gff)
    
    cat("  Created", basename(output_gff), "with", length(related_features), "features\n")
    return(related_features)
}

# Function to create subset FASTA file with matching sequences
createSubsetFASTA <- function(input_fasta, output_fasta, gff_features, n_sequences = 100) {
    cat("Processing", basename(input_fasta), "...\n")
    
    # Read FASTA file
    if (grepl("prot", input_fasta)) {
        sequences <- readAAStringSet(input_fasta)
    } else {
        sequences <- readDNAStringSet(input_fasta)
    }
    
    # Get sequence names from GFF features
    if (grepl("prot", input_fasta)) {
        # For protein files, look for transcript IDs
        gff_seq_names <- unique(unlist(gff_features$Parent))
    } else {
        # For CDS files, look for transcript IDs
        gff_seq_names <- unique(unlist(gff_features$Parent))
    }
    
    # Find matching sequences
    matching_sequences <- sequences[names(sequences) %in% gff_seq_names]
    
    if (length(matching_sequences) == 0) {
        # If no matches, take first n sequences
        n_select <- min(n_sequences, length(sequences))
        selected_sequences <- sequences[1:n_select]
    } else {
        # Use matching sequences
        n_select <- min(n_sequences, length(matching_sequences))
        selected_sequences <- matching_sequences[1:n_select]
    }
    
    # Write subset FASTA
    if (grepl("prot", input_fasta)) {
        writeXStringSet(selected_sequences, output_fasta)
    } else {
        writeXStringSet(selected_sequences, output_fasta)
    }
    
    cat("  Created", basename(output_fasta), "with", length(selected_sequences), "sequences\n")
    return(selected_sequences)
}

# Function to create subset genome file
createSubsetGenome <- function(input_genome, output_genome, max_length = 1000000) {
    cat("Processing", basename(input_genome), "...\n")
    
    # Read genome file
    genome <- readDNAStringSet(input_genome)
    
    # Select first few chromosomes/scaffolds
    n_select <- min(3, length(genome))
    selected_genome <- genome[1:n_select]
    
    # Truncate sequences if too long
    for (i in 1:length(selected_genome)) {
        if (width(selected_genome[i]) > max_length) {
            selected_genome[i] <- subseq(selected_genome[i], 1, max_length)
        }
    }
    
    # Write subset genome
    writeXStringSet(selected_genome, output_genome)
    
    cat("  Created", basename(output_genome), "with", length(selected_genome), "sequences\n")
    return(selected_genome)
}

# Create subset files for NB (query species)
cat("\n=== Creating NB subset files ===\n")
nb_gff_subset <- createSubsetGFF("inst/extdata/nb.gff", 
                                 file.path(subset_dir, "nb_subset.gff"), 
                                 n_genes = 30)
nb_cds_subset <- createSubsetFASTA("inst/extdata/nb_cds.fa", 
                                   file.path(subset_dir, "nb_subset_cds.fa"), 
                                   nb_gff_subset,
                                   n_sequences = 50)
nb_prot_subset <- createSubsetFASTA("inst/extdata/nb_prot.fa", 
                                    file.path(subset_dir, "nb_subset_prot.fa"), 
                                    nb_gff_subset,
                                    n_sequences = 50)
nb_genome_subset <- createSubsetGenome("inst/extdata/nb_genome.fa", 
                                       file.path(subset_dir, "nb_subset_genome.fa"), 
                                       max_length = 500000)

# Create subset files for WK21 (subject species)
cat("\n=== Creating WK21 subset files ===\n")
wk21_gff_subset <- createSubsetGFF("inst/extdata/wk21.gff", 
                                   file.path(subset_dir, "wk21_subset.gff"), 
                                   n_genes = 30)
wk21_cds_subset <- createSubsetFASTA("inst/extdata/wk21_cds.fa", 
                                     file.path(subset_dir, "wk21_subset_cds.fa"), 
                                     wk21_gff_subset,
                                     n_sequences = 50)
wk21_prot_subset <- createSubsetFASTA("inst/extdata/wk21_prot.fa", 
                                      file.path(subset_dir, "wk21_subset_prot.fa"), 
                                      wk21_gff_subset,
                                      n_sequences = 50)
wk21_genome_subset <- createSubsetGenome("inst/extdata/wk21_genome.fa", 
                                         file.path(subset_dir, "wk21_subset_genome.fa"), 
                                         max_length = 500000)

# Create a summary file
summary_file <- file.path(subset_dir, "subset_summary.txt")
cat("OrthoPaiR Subset Data Summary\n", file = summary_file)
cat("=============================\n", file = summary_file, append = TRUE)
cat("Created on:", as.character(Sys.time()), "\n", file = summary_file, append = TRUE)
cat("\n", file = summary_file, append = TRUE)

cat("NB (Query Species) Subset:\n", file = summary_file, append = TRUE)
cat("- GFF features:", length(nb_gff_subset), "\n", file = summary_file, append = TRUE)
cat("- CDS sequences:", length(nb_cds_subset), "\n", file = summary_file, append = TRUE)
cat("- Protein sequences:", length(nb_prot_subset), "\n", file = summary_file, append = TRUE)
cat("- Genome sequences:", length(nb_genome_subset), "\n", file = summary_file, append = TRUE)
cat("\n", file = summary_file, append = TRUE)

cat("WK21 (Subject Species) Subset:\n", file = summary_file, append = TRUE)
cat("- GFF features:", length(wk21_gff_subset), "\n", file = summary_file, append = TRUE)
cat("- CDS sequences:", length(wk21_cds_subset), "\n", file = summary_file, append = TRUE)
cat("- Protein sequences:", length(wk21_prot_subset), "\n", file = summary_file, append = TRUE)
cat("- Genome sequences:", length(wk21_genome_subset), "\n", file = summary_file, append = TRUE)

cat("\n=== Subset data creation completed ===\n")
cat("Files created in:", subset_dir, "\n")
cat("Summary saved to:", summary_file, "\n")
