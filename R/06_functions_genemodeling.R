#' Map Proteins Between Genomes
#'
#' This function maps proteins between the query and subject genomes using the specified miniprot binary and number of cores.
#'
#' @param object A SynogDB object.
#' @param out_prefix A character string specifying the prefix for the output files.
#' @param miniprot_bin A character string specifying the path to the miniprot binary.
#' @param n_threads An integer specifying the number of cores to use for the mapping.
#' @param len_diff A numeric value specifying the maximum allowable length difference for proteins. Default is 0.5.
#'
#' @return None. The function performs the mapping and writes the results to the output files specified by `out_prefix`.
#'
#' @export
mapProt <- function(object,
                    out_dir,
                    miniprot_bin = "miniprot",
                    conda = "conda",
                    condaenv = NULL,
                    n_threads,
                    len_diff = 0.5,
                    overlap = TRUE){
    stopifnot(inherits(x = object, "SynogDB"))

    # Call the mapping engine function with specified parameters
    .mapEngine(object = object,
               subject_prot = object$subject_prot,
               query_prot = object$query_prot,
               out_dir = out_dir,
               miniprot_bin = miniprot_bin,
               conda = conda,
               condaenv = condaenv,
               n_threads = n_threads,
               overlap = overlap,
               len_diff = len_diff)

    .createFASTA(object = object, out_dir = out_dir)
    object <- .updateFiles(object = object)

    .syntenyOrthoMiniprot(object = object)
}

#' Execute Protein Mapping Between Genomes
#'
#' This function maps proteins between the query and subject genomes using the specified miniprot binary and number of cores. It performs both directions of mapping and organizes the results.
#'
#' @importFrom rtracklayer export.gff3
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
.mapEngine <- function(object, subject_prot, query_prot,
                       out_dir, miniprot_bin, conda,
                       condaenv, n_threads,
                       overlap, len_diff){
    dir.create(path = out_dir, showWarnings = FALSE, recursive = TRUE)

    # Define output file paths for both directions of mapping
    s2q_out <- file.path(out_dir, "query_miniprot_out")
    q2s_out <- file.path(out_dir, "subject_miniprot_out")

    # Run miniprot mapping from subject to query genome
    .miniprot(query_fn = subject_prot, genome_fn = object$query_genome,
              out_prefix = s2q_out, miniprot_bin = miniprot_bin,
              conda = conda, condaenv = condaenv,
              n_threads = n_threads)

    # Run miniprot mapping from query to subject genome
    .miniprot(query_fn = query_prot, genome_fn = object$subject_genome,
              out_prefix = q2s_out, miniprot_bin = miniprot_bin,
              conda = conda, condaenv = condaenv,
              n_threads = n_threads)

    # Create HDF5 group for protein mapping results and save GFF files
    .h5creategroup(object$h5,"protmap")
    .h5overwrite(obj = paste0(q2s_out, ".gff"),
                 file = object$h5, "protmap/q2s_gff")
    .h5overwrite(obj = paste0(s2q_out, ".gff"),
                 file = object$h5, "protmap/s2q_gff")

    # Import existing GFF files for query and subject genomes
    query_gff <- .importAllGFF(object$query_gff)
    subject_gff <- .importAllGFF(object$subject_gff)

    # Organize GFF results from subject-to-query and query-to-subject mappings
    s2q_gff <- .orgGFF(gff1 = import.gff(paste0(s2q_out, ".gff")),
                       gff2 = query_gff,
                       gff3 = subject_gff,
                       prefix = "query_",
                       overlap = overlap,
                       len_diff = len_diff)
    q2s_gff <- .orgGFF(gff1 = import.gff(paste0(q2s_out, ".gff")),
                       gff2 = subject_gff,
                       gff3 = query_gff,
                       prefix = "subject_",
                       overlap = overlap,
                       len_diff = len_diff)

    # Merge organized GFF results with original query and subject GFFs
    s2q_gff <- .mergeGFF(gff1 = s2q_gff, gff2 = query_gff)
    q2s_gff <- .mergeGFF(gff1 = q2s_gff, gff2 = subject_gff)

    # Fix and filter GFF results, retaining only necessary columns
    s2q_gff <- .fixGFF(gff = s2q_gff)
    q2s_gff <- .fixGFF(gff = q2s_gff)
    m_s2q_gff <- mcols(s2q_gff)
    hit <- names(m_s2q_gff) %in% c("source", "type", "score", "phase",
                                   "ID", "Name", "gene_id", "Parent", "Target")
    mcols(s2q_gff) <- m_s2q_gff[, hit]
    m_q2s_gff <- mcols(q2s_gff)
    hit <- names(m_q2s_gff) %in% c("source", "type", "score", "phase",
                                   "ID", "Name", "gene_id", "Parent", "Target")
    mcols(q2s_gff) <- m_q2s_gff[, hit]
    s2q_gff$Name <- s2q_gff$ID
    q2s_gff$Name <- q2s_gff$ID

    # Export the final GFF results to output files
    export.gff3(s2q_gff, paste0(s2q_out, ".gff"))
    export.gff3(q2s_gff, paste0(q2s_out, ".gff"))
}

.miniprot <- function(query_fn, genome_fn, out_prefix,
                      miniprot_bin,
                      conda, condaenv,
                      n_threads = 1){

    if(!is.null(condaenv)){
        .condaExe(conda = conda, env = condaenv, command = miniprot_bin,
                  args = paste("-t", n_threads,
                               "-d", paste0(out_prefix, ".mpi"),
                               genome_fn))
        .condaExe(conda = conda, env = condaenv, command = miniprot_bin,
                  args = paste("-t", n_threads, "--gff",
                               paste0(out_prefix, ".mpi"),
                               query_fn, ">", paste0(out_prefix, ".gff")))

    } else {
        system2(command = miniprot_bin,
                args = paste("-t", n_threads,
                             "-d", paste0(out_prefix, ".mpi"),
                             genome_fn))
        system2(command = miniprot_bin,
                args = paste("-t", n_threads,
                             "--gff",
                             paste0(out_prefix, ".mpi"),
                             query_fn, ">", paste0(out_prefix, ".gff")))
    }
}


#' Write Protein Sequences to a FASTA File
#'
#' This function writes protein sequences to a FASTA file based on the specified IDs.
#'
#' @importFrom Biostrings readAAStringSet writeXStringSet
.writeProtFASTA <- function(id, prot_fn, gff_fn, fn){
    # Read protein sequences from the FASTA file
    aa <- readAAStringSet(prot_fn)

    # Import GFF data
    gff <- .importAllGFF(gff_fn)

    # Identify transcripts and match them to the provided IDs
    tx_i <- gff$type %in% c("mRNA", "transcript")
    tx_hit <- gff$ID[tx_i] %in% id
    hit_tx <- gff$ID[tx_i][tx_hit]

    # Filter the protein sequences to include only those matching the IDs
    hit_aa <- aa[names(aa) %in% hit_tx]

    # Write the filtered protein sequences to the output FASTA file
    writeXStringSet(hit_aa, fn)
}
#' Write Protein Sequences to a FASTA File
#'
#' This function writes protein sequences to a FASTA file based on the specified IDs.
#'
#' @importFrom Biostrings readAAStringSet writeXStringSet
.writeProtFASTA <- function(id, prot_fn, gff_fn, fn){
    # Read protein sequences from the FASTA file
    aa <- readAAStringSet(prot_fn)

    # Import GFF data
    gff <- .importAllGFF(gff_fn)

    # Identify transcripts and match them to the provided IDs
    tx_i <- gff$type %in% c("mRNA", "transcript")
    tx_hit <- gff$ID[tx_i] %in% id
    hit_tx <- gff$ID[tx_i][tx_hit]

    # Filter the protein sequences to include only those matching the IDs
    hit_aa <- aa[names(aa) %in% hit_tx]

    # Write the filtered protein sequences to the output FASTA file
    writeXStringSet(hit_aa, fn)
}
#' Organize GFF Data
#'
#' This function organizes and validates GFF data, filtering transcripts and organizing genes, transcripts, and CDS features.
#'
#' @importFrom BiocGenerics start
.orgGFF <- function(gff1, gff2, gff3, prefix, overlap = FALSE, len_diff){
    # Adjust levels for transcript types
    lv <- levels(gff1$type)
    lv[lv == "mRNA"] <- "transcript"
    levels(gff1$type) <- lv

    # Filter relevant types and sort GFF data
    gff1 <- gff1[gff1$type %in% c("gene", "transcript", "five_prime_UTR", "CDS", "three_prime_UTR")]
    gff1 <- gff1[order(as.numeric(seqnames(gff1)), start(gff1))]
    gff2 <- gff2[order(as.numeric(seqnames(gff2)), start(gff2))]

    # Validate GFF data
    gff_valid <- .validTX(gff1 = gff1, gff2 = gff2, gff3 = gff3, overlap = overlap, len_diff = len_diff)

    # Organize transcripts and CDS features
    gff_tx <- .orgTX(gff = gff_valid, prefix = prefix)
    gff_cds <- .orgCDS(gff1 = gff1, gff_tx = gff_tx)
    gff_gene <- .orgGene(gff_tx = gff_tx)

    # Remove old IDs from transcripts and genes
    gff_tx$old_id <- gff_gene$old_id <- NULL

    # Combine genes, transcripts, and CDS features
    out <- c(gff_gene, gff_tx, gff_cds)

    # Sort the final output
    out <- out[order(as.numeric(seqnames(out)), start(out), as.numeric(out$type))]

    return(out)
}
#' Validate Transcripts
#'
#' This function validates transcripts by checking overlaps and length differences with reference GFF data.
#'
#' @importFrom S4Vectors queryHits
#' @importFrom GenomicRanges findOverlaps
.validTX <- function(gff1, gff2, gff3, overlap, len_diff){
    if(overlap){
        # Get unique CDS blocks from gff1
        gff1_block <- .getCDSblock(gff = gff1)
        gff1_block_uniq <- gff1_block[!duplicated(gff1_block)]

        # Get CDS blocks from gff2 and remove overlaps from gff1_block_uniq
        gff2_block <- .getCDSblock(gff = gff2)
        gff1_block_uniq <- gff1_block_uniq[!gff1_block_uniq %in% gff2_block]

        # Get unique transcripts from gff1
        out <- .getUniqTx(gff1 = gff1, gff1_block_uniq = gff1_block_uniq)

    } else {
        # Find overlaps between gff1 and gff2 transcripts
        tx_i1 <- gff1$type == "transcript"
        tx_i2 <- gff2$type == "transcript"
        ol <- findOverlaps(gff1[tx_i1], gff2[tx_i2])
        hit <- unique(queryHits(ol))

        # Exclude overlapping transcripts
        out <- gff1[tx_i1][!seq_len(sum(tx_i1)) %in% hit]
    }

    # Check length differences with gff3
    out <- .checkLength(gff = out, gff3 = gff3, len_diff = len_diff)


    return(out)
}
#' Get CDS Blocks
#'
#' This function extracts the concatenated CDS blocks for each transcript in the GFF data.
#'
#' @importFrom BiocGenerics start end
.getCDSblock <- function(gff){
    # Identify CDS and transcript indices
    gff_cds_i <- gff$type == "CDS"
    gff_tx_i <- gff$type == "transcript"

    # Extract start and end positions of CDS
    gff_cds_start <- start(gff[gff_cds_i])
    gff_cds_end <- end(gff[gff_cds_i])

    # Create exon strings combining start and end positions
    gff_cds_exon <- paste(gff_cds_start, gff_cds_end, sep = "-")

    # Map CDS to their parent transcripts
    map_to_tx <- match(unlist(gff$Parent[gff_cds_i]), gff$ID[gff_tx_i])

    # Identify the first occurrence of each transcript
    first_i <- !duplicated(map_to_tx)

    # Get the chromosome names for CDS
    gff_cds_chr <- as.character(seqnames(gff[gff_cds_i]))

    # Add chromosome information to the first occurrence of each transcript's CDS block
    gff_cds_exon[first_i] <- paste(gff_cds_chr[first_i], gff_cds_exon[first_i], sep = ":")

    # Concatenate CDS blocks for each transcript
    out <- tapply(gff_cds_exon, map_to_tx, paste, collapse = ",")

    return(out)
}
#' Get Unique Transcripts
#'
#' This function extracts unique transcripts and their associated elements (e.g., CDS) from a GFF object.
#'
#' @importFrom BiocGenerics unlist
.getUniqTx <- function(gff1, gff1_block_uniq){
    # Identify transcript indices
    gff1_tx_i <- gff1$type == "transcript"

    # Extract unique transcripts using the unique CDS blocks
    out_tx <- gff1[gff1_tx_i][as.numeric(names(gff1_block_uniq))]

    # Identify and extract elements associated with the unique transcripts
    out_element <- unlist(gff1$Parent[!gff1_tx_i]) %in% out_tx$ID
    out_element <- gff1[!gff1_tx_i][out_element]

    # Combine unique transcripts and their associated elements
    out <- c(out_tx, out_element)

    return(out)
}
#' Check Length of Transcripts
#'
#' This function checks the length of transcripts and filters those that meet the length difference criteria.
#'
#' @importFrom BiocGenerics unlist
.checkLength <- function(gff, gff3, len_diff){
    # Calculate CDS lengths for the given GFF data
    gff_cds_len <- .getCDSlen(gff = gff)
    index <- as.numeric(names(gff_cds_len))
    gff_tx_i <- gff$type %in% "transcript"
    names(gff_cds_len) <- sub("\\s.+", "", gff$Target[gff_tx_i][index])

    # Calculate CDS lengths for the reference GFF data
    gff3_cds_len <- .getCDSlen(gff = gff3)
    index <- as.numeric(names(gff3_cds_len))
    gff3_tx_i <- gff3$type %in% c("transcript", "mRNA")
    names(gff3_cds_len) <- gff3$ID[gff3_tx_i][index]

    # Match CDS lengths by transcript IDs
    id_hit <- match(names(gff_cds_len), names(gff3_cds_len))
    gff3_cds_len <- gff3_cds_len[id_hit]

    # Determine which transcripts are longer
    longer <- gff_cds_len
    is_gff3_longer <- longer < gff3_cds_len
    longer[is_gff3_longer] <- gff3_cds_len[is_gff3_longer]

    # Calculate valid transcripts based on length difference criteria
    valid <- abs(gff_cds_len - gff3_cds_len) / longer <= len_diff
    valid_tx <- gff[gff_tx_i][as.vector(valid)]

    # Extract non-transcript elements associated with valid transcripts
    non_tx_gff <- gff[!gff_tx_i]
    out <- c(valid_tx, non_tx_gff[unlist(non_tx_gff$Parent) %in% valid_tx$ID])

    return(out)
}

#' Calculate CDS Lengths
#'
#' This function calculates the lengths of coding sequences (CDS) for each transcript in a GFF file.
#'
#' @importFrom BiocGenerics start end
.getCDSlen <- function(gff){
    # Identify indices for CDS and transcript elements
    gff_cds_i <- gff$type == "CDS"
    gff_tx_i <-  gff$type == "transcript"

    # Extract start and end positions for CDS elements
    gff_cds_start <- start(gff[gff_cds_i])
    gff_cds_end <- end(gff[gff_cds_i])

    # Calculate the length of each CDS element
    gff_cds_len <- gff_cds_end - gff_cds_start

    # Map CDS elements to their corresponding transcripts
    map_to_tx <- match(unlist(gff$Parent[gff_cds_i]), gff$ID[gff_tx_i])

    # Sum the lengths of CDS elements for each transcript
    out <- tapply(gff_cds_len, map_to_tx, sum)

    return(out)
}

#' Organize Transcripts
#'
#' This function organizes transcripts in a GFF file by assigning gene IDs and updating transcript IDs.
#'
#' @importFrom BiocGenerics start
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges reduce findOverlaps
#' @importFrom S4Vectors subjectHits
.orgTX <- function(gff, prefix){
    # Identify indices for transcript elements
    tx_i <- gff$type == "transcript"

    # Extract and order transcript elements by chromosome and start position
    gff_tx <- gff[tx_i]
    gff_tx <- gff_tx[order(as.numeric(seqnames(gff_tx)), start(gff_tx))]

    # Remove duplicate transcripts
    gff_tx <- unique(gff_tx)

    # Reduce transcripts to unique loci
    uni_loci <- reduce(gff_tx)

    # Map transcripts to unique loci
    map_to_loci <- findOverlaps(gff_tx, uni_loci)

    # Assign gene IDs based on unique loci
    gff_tx$gene_id <- paste0(prefix, "G", sprintf("%05d", subjectHits(map_to_loci)))

    # Store old transcript IDs
    gff_tx$old_id <- gff_tx$ID

    # Order transcripts by gene ID
    gff_tx <- gff_tx[order(gff_tx$gene_id)]

    # Update transcript IDs based on gene ID
    gff_tx$ID <- unlist(tapply(gff_tx$gene_id, gff_tx$gene_id, function(x){
        return(paste0(sub("G", "T", x[1]), ".", sprintf("%02d", seq_len(length(x)))))
    }))

    # Assign parent gene IDs to transcripts
    gff_tx$Parent <- lapply(gff_tx$gene_id, c)

    return(gff_tx)
}
#' Organize CDS Elements
#'
#' This function organizes CDS elements in a GFF file by assigning gene IDs and updating parent IDs.
#'
.orgCDS <- function(gff1, gff_tx){
    # Identify indices for non-transcript elements
    cds_i <- gff1$type != "transcript"

    # Extract CDS elements that have parents in the organized transcript data
    gff_cds <- gff1[cds_i][unlist(gff1$Parent[cds_i]) %in% gff_tx$old_id]

    # Assign gene IDs based on organized transcript data
    gff_cds$gene_id <- gff_tx$gene_id[match(unlist(gff_cds$Parent), gff_tx$old_id)]

    # Update parent IDs based on organized transcript data
    gff_cds$Parent <- lapply(gff_tx$ID[match(unlist(gff_cds$Parent), gff_tx$old_id)], c)

    # Update CDS IDs based on parent IDs
    gff_cds$ID <- paste0(unlist(gff_cds$Parent), ":CDS")

    return(gff_cds)
}

#' Organize Gene Elements
#'
#' This function organizes gene elements in a GFF file by creating gene entries from the organized transcript data.
#'
.orgGene <- function(gff_tx){
    # Extract unique gene entries from the transcript data
    gff_gene <- gff_tx[!duplicated(gff_tx$gene_id)]

    # Set the type of each entry to "gene"
    gff_gene$type <- "gene"

    # Assign gene IDs to the ID field
    gff_gene$ID <- gff_gene$gene_id

    # Remove parent information for gene entries
    gff_gene$Parent <- lapply(seq_along(gff_gene), function(i) character())

    return(gff_gene)
}
#' Merge GFF Annotations
#'
#' This function merges GFF annotations by combining two GFF files, removing duplicate gene and transcript entries, and creating a unified annotation.
#'
#' @importFrom BiocGenerics start
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges reduce findOverlaps
#' @importFrom S4Vectors subjectHits
.mergeGFF <- function(gff1, gff2 = NULL, ann_priority){
    # Combine the two GFF objects
    gff <- c(gff1, gff2)

    # Filter for relevant types
    gff <- gff[gff$type %in% c("gene", "transcript",
                               "five_prime_UTR", "CDS", "three_prime_UTR")]

    # Order and reduce gene entries
    gene_i <- gff$type == "gene"
    gff_gene <- gff[gene_i]
    gff_gene <- gff_gene[order(as.numeric(seqnames(gff_gene)), start(gff_gene))]
    uni_loci <- reduce(gff_gene)
    map_to_loci <- findOverlaps(gff_gene, uni_loci)

    # Identify overlapping entries and create ID map
    ol_list <- tapply(queryHits(map_to_loci), subjectHits(map_to_loci), c)
    ol_list <- ol_list[sapply(ol_list, length) > 1]
    id_map <- .getIDmap(gff_gene = gff_gene, ol_list = ol_list)

    # Map IDs in the GFF object based on the ID map
    gff <- .mapID(gff = gff, id_map = id_map)

    # Order the final merged GFF object
    gff <- gff[order(as.numeric(seqnames(gff)), start(gff), as.numeric(gff$type))]

    return(gff)
}
#' Create ID Map for Gene Overlaps
#'
#' This function creates an ID map for overlapping genes by selecting a representative gene ID for each set of overlapping genes.
#'
.getIDmap <- function(gff_gene, ol_list){
    # Create ID map for each set of overlapping genes
    out <- lapply(ol_list, function(i){
        gene_id <- gff_gene$ID[i]
        # Select the representative gene ID that does not start with 'query_' or 'subject_'
        set_gene_id <- gene_id[!grepl("^query_|^subject_", gene_id)][1]
        id_map <- cbind(gene_id, set_gene_id)
        return(id_map)
    })
    # Combine all ID maps into a single data frame
    out <- do.call("rbind", out)
    return(out)
}


.mapID <- function(gff, id_map){
    gff_parent <- sapply(gff$Parent, function(x){
        if(length(x) == 0){
            return(NA)
        } else {
            return(x)
        }
    })
    hit <- match(gff_parent, id_map[, 1])
    gff$Parent[!is.na(hit)] <- lapply(id_map[na.omit(hit), 2], c)
    gff$gene_id[!is.na(hit)] <- id_map[na.omit(hit), 2]

    rm_gene <- id_map[id_map[, 1] != id_map[, 2], ]
    gff <- gff[!gff$ID %in% rm_gene[, 1]]
    return(gff)
}


#' Fix and Standardize GFF Annotations
#'
#' This function applies a series of fixes and standardizations to GFF annotations.
#'
.fixGFF <- function(gff){
    gff <- .setGeneID(gff = gff)
    gff <- .fixGFFrange(gff = gff)
    gff <- .fixGFFexon(gff = gff)
    gff <- .fixGFFphase(gff = gff)
    return(gff)
}

#' Set Gene IDs for GFF Annotations
#'
#' This function sets the gene_id attribute for GFF annotations based on their parent-child relationships.
#'
.setGeneID <- function(gff){
    # For transcripts or mRNA, set gene_id to their ID
    tx_i <- gff$type %in% c("transcript", "mRNA")
    hit <- match(unlist(gff$Parent[tx_i]), gff$ID)
    gff$gene_id[tx_i] <- gff$ID[hit]

    # For other elements, set gene_id based on their parent transcript/mRNA
    element_i <- !gff$type %in% c("gene", "transcript", "mRNA")
    hit <- match(unlist(gff$Parent[element_i]), gff$ID[tx_i])
    gff$gene_id[element_i] <- gff$gene_id[tx_i][hit]
    return(gff)
}

#' Fix GFF Range
#'
#' This function adjusts the start and end positions of transcripts and genes to cover the complete range of their member elements.
#'
.fixGFFrange <- function(gff){
    # Fix the start and end positions of each transcript to
    # cover the whole range of member elements (CDS, exon, and UTRs)
    tx_i <- which(gff$type %in% c("transcript", "mRNA"))
    element_i <- !gff$type %in% c("gene", "transcript", "mRNA")

    # Determine minimum start and maximum end positions for member elements
    min_start <- tapply(start(gff[element_i]), unlist(gff$Parent[element_i]), min)
    max_end <- tapply(end(gff[element_i]), unlist(gff$Parent[element_i]), max)

    # Match and update transcript start and end positions
    hit <- match(gff$ID[tx_i], names(min_start))
    tx_start <- min_start[hit]
    tx_end <- max_end[hit]
    not_na_start <- !is.na(tx_start)
    not_na_end <- !is.na(tx_end)
    start(gff[tx_i[not_na_start]]) <- tx_start[not_na_start]
    end(gff[tx_i[not_na_end]]) <- tx_end[not_na_end]

    # Fix the start and end positions of each gene to
    # cover the whole range of member transcripts
    gene_i <- which(gff$type == "gene")
    tx_i <- which(gff$type %in% c("transcript", "mRNA"))

    # Determine minimum start and maximum end positions for transcripts
    min_start <- tapply(start(gff[tx_i]), gff$gene_id[tx_i], min)
    max_end <- tapply(end(gff[tx_i]), gff$gene_id[tx_i], max)

    # Match and update gene start and end positions
    hit <- match(gff$gene_id[gene_i], names(min_start))
    gene_start <- min_start[hit]
    gene_end <- max_end[hit]
    start(gff[gene_i]) <- gene_start
    end(gff[gene_i]) <- gene_end

    return(gff)
}


.fixGFFexon <- function(gff){
    gff <- gff[gff$type != "exon"]
    gff_exon <- gff[gff$type %in% c("CDS", "five_prime_UTR", "three_prime_UTR")]
    gff_exon$type <- "exon"
    gff_exon$Name <- gff_exon$ID <- paste0(unlist(gff_exon$Parent), ":exon")
    gff_exon$score <- gff_exon$phase <- NA

    tx_i <- gff$type %in% c("transcript", "mRNA")
    non_cds <- gff[tx_i][!gff$ID[tx_i] %in% unlist(gff_exon$Parent)]
    if(length(non_cds) > 0){
        non_cds$type <- "exon"
        non_cds$Parent <- lapply(non_cds$ID, c)
        non_cds$Name <- non_cds$ID <- paste0(non_cds$ID, ":exon")
        non_cds$score <- non_cds$phase <- NA
        gff_exon <- c(gff_exon, non_cds)
    }

    out <- c(gff, gff_exon)
    out$type <- factor(out$type, levels = c("gene", "transcript", "mRNA",
                                            "five_prime_UTR", "exon",
                                            "CDS", "three_prime_UTR"))
    out$type <- droplevels(out$type)
    out <- out[order(as.numeric(seqnames(out)), start(out), as.numeric(out$type))]
    return(out)
}


#' Fix GFF Phase
#'
#' This function adjusts the phase of CDS features in GFF annotations.
#'
#' @importFrom BiocGenerics start
.fixGFFphase <- function(gff){
    # Extract CDS features
    gff_cds <- gff[gff$type == "CDS"]

    # Adjust phase for plus strand
    gff_cds_plus <- .phasePlus(gff_cds = gff_cds)

    # Adjust phase for minus strand
    gff_cds_minus <- .phaseMinus(gff_cds = gff_cds)

    # Combine non-CDS features with adjusted CDS features
    out <- c(gff[gff$type != "CDS"], gff_cds_plus, gff_cds_minus)

    # Set factor levels for feature types
    out$type <- factor(out$type, levels = c("gene", "transcript", "mRNA",
                                            "five_prime_UTR", "exon",
                                            "CDS", "three_prime_UTR"))
    out$type <- droplevels(out$type)

    # Order the features by sequence name, start position, and type
    out <- out[order(as.numeric(seqnames(out)), start(out), as.numeric(out$type))]
    return(out)
}
#' Adjust Phase for CDS Features on Plus Strand
#'
#' This function adjusts the phase for CDS features on the plus strand in GFF annotations.
#'
#' @importFrom BiocGenerics start width
#' @importFrom GenomeInfoDb seqnames
.phasePlus <- function(gff_cds){
    # Filter CDS features on the plus strand
    gff_cds_plus <- gff_cds[as.character(strand(gff_cds)) == "+"]
    gff_cds_plus <- gff_cds_plus[order(as.numeric(seqnames(gff_cds_plus)), start(gff_cds_plus))]

    # Extract parent IDs for CDS features
    gff_cds_plus_parent <- unlist(gff_cds_plus$Parent)

    # Initialize phase for the first CDS in each transcript
    target_i <- which(!duplicated(gff_cds_plus_parent))
    gff_cds_plus$phase[target_i] <- 0

    # Calculate the phase for the next CDS feature
    next_phase <- (3 - (width(gff_cds_plus[target_i]) - gff_cds_plus$phase[target_i]) %% 3) %% 3
    names(next_phase) <- gff_cds_plus_parent[target_i]
    gff_cds_plus_parent[target_i] <- "NA"

    # Iterate over the remaining CDS features to set their phases
    while(TRUE){
        target_i <- which(!duplicated(gff_cds_plus_parent))[-1]
        if(length(target_i) == 0){
            break
        }
        next_phase <- next_phase[names(next_phase) %in% gff_cds_plus_parent[target_i]]
        gff_cds_plus$phase[target_i] <- next_phase
        next_phase <- (3 - (width(gff_cds_plus[target_i]) - gff_cds_plus$phase[target_i]) %% 3) %% 3
        names(next_phase) <- gff_cds_plus_parent[target_i]
        gff_cds_plus_parent[target_i] <- "NA"
    }

    return(gff_cds_plus)
}

#' Adjust Phase for CDS Features on Minus Strand
#'
#' This function adjusts the phase for CDS features on the minus strand in GFF annotations.
#'
#' @importFrom BiocGenerics end width
#' @importFrom GenomeInfoDb seqnames
.phaseMinus <- function(gff_cds){
    # Filter CDS features on the minus strand
    gff_cds_minus <- gff_cds[as.character(strand(gff_cds)) == "-"]
    gff_cds_minus <- gff_cds_minus[order(as.numeric(seqnames(gff_cds_minus)),
                                         end(gff_cds_minus), decreasing = TRUE)]

    # Extract parent IDs for CDS features
    gff_cds_minus_parent <- unlist(gff_cds_minus$Parent)

    # Initialize phase for the first CDS in each transcript
    target_i <- which(!duplicated(gff_cds_minus_parent))
    gff_cds_minus$phase[target_i] <- 0

    # Calculate the phase for the next CDS feature
    next_phase <- (3 - (width(gff_cds_minus[target_i]) - gff_cds_minus$phase[target_i]) %% 3) %% 3
    names(next_phase) <- gff_cds_minus_parent[target_i]
    gff_cds_minus_parent[target_i] <- "NA"

    # Iterate over the remaining CDS features to set their phases
    while(TRUE){
        target_i <- which(!duplicated(gff_cds_minus_parent))[-1]
        if(length(target_i) == 0){
            break
        }
        next_phase <- next_phase[names(next_phase) %in% gff_cds_minus_parent[target_i]]
        gff_cds_minus$phase[target_i] <- next_phase
        next_phase <- (3 - (width(gff_cds_minus[target_i]) - gff_cds_minus$phase[target_i]) %% 3) %% 3
        names(next_phase) <- gff_cds_minus_parent[target_i]
        gff_cds_minus_parent[target_i] <- "NA"
    }

    return(gff_cds_minus)
}

################################################################################

#' Create FASTA files
#'
#' This function generates FASTA files for query and subject CDS from the SynogDB object.
#'
#' @param object A SynogDB object.
#' @param out_dir Output directory for FASTA files.
#'
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#' @importFrom Biostrings writeXStringSet
.createFASTA <- function(object, out_dir){
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    # Create CDS sequences for query and subject genomes
    q_cds <- .makeCDS(gff = as.vector(h5$protmap$s2q_gff),
                      genome = object$query_genome)
    s_cds <- .makeCDS(gff = as.vector(h5$protmap$q2s_gff),
                      genome = object$subject_genome)

    # Define filenames for the CDS FASTA files
    q_cds_fn <- sub("\\.gff", "_cds.fa", as.vector(h5$protmap$s2q_gff))
    s_cds_fn <- sub("\\.gff", "_cds.fa", as.vector(h5$protmap$q2s_gff))

    # Write the CDS sequences to FASTA files
    writeXStringSet(q_cds, q_cds_fn)
    writeXStringSet(s_cds, s_cds_fn)

    # Overwrite the HDF5 file with the new CDS FASTA filenames
    .h5overwrite(obj = q_cds_fn, file = object$h5, "protmap/s2q_cds")
    .h5overwrite(obj = s_cds_fn, file = object$h5, "protmap/q2s_cds")
}

#' Create CDS sequences from GFF and genome files
#'
#' This function generates CDS sequences from provided GFF and genome files.
#'
#' @param gff Path to the GFF file.
#' @param genome Path to the genome file.
#'
#' @return A DNAStringSet object containing CDS sequences.
#' @importFrom Biostrings readDNAStringSet
#' @importFrom GenomicFeatures cdsBy extractTranscriptSeqs
#' @importFrom txdbmaker makeTxDbFromGFF
#' @import BSgenome
#'
.makeCDS <- function(gff, genome){
    # Create a TxDb object from the GFF file
    txdb <- makeTxDbFromGFF(file = gff)

    # Read the genome file as a DNAStringSet object
    genome <- readDNAStringSet(filepath = genome)

    # Extract CDS sequences from the TxDb object
    cds_db <- cdsBy(x = txdb, by = "tx", use.names = TRUE)
    cds <- extractTranscriptSeqs(x = genome, transcripts = cds_db)

    # Order CDS sequences by their names
    cds <- cds[order(names(cds))]
    return(cds)
}

#' Update SynogDB files
#'
#' This function updates the GFF and CDS files in a SynogDB object from the HDF5 file.
#'
#' @param object A SynogDB object.
#'
#' @return The updated SynogDB object.
#' @importFrom rhdf5 H5Fopen H5Fclose
.updateFiles <- function(object){
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))

    # Update GFF and CDS file paths from the HDF5 file
    object$query_gff <- as.vector(h5$protmap$s2q_gff)
    object$subject_gff <- as.vector(h5$protmap$q2s_gff)
    object$query_cds <- as.vector(h5$protmap$s2q_cds)
    object$subject_cds <- as.vector(h5$protmap$q2s_cds)
    return(object)
}

.syntenyOrthoMiniprot <- function(object){
    # Open the HDF5 file
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))

    # Check for the existence of syntenic ortholog pairs
    if(!H5Lexists(h5, "synog_tx/orthopairs")){
        stop("Run syntenyOrtho to obtain ortholog anchors info.")
    }

    # Get the GFF lists
    gff_ls <- .getGFFlist(object = object, h5 = h5)

    out <- .filterGFF(object = object, gff_ls = gff_ls, h5 = h5)
#
#     lcb_gr <- .makeLCBgr(object = object, h5 = h5)
#
#     .remakeAnchor(gff_ls = gff_ls, lcb_gr = lcb_gr, h5 = h5)
#
#     # Create anchor lists
#     anchor_ls <- .makeAnchorList(gff_ls = gff_ls,
#                                  h5 = h5,
#                                  omit_chr = omit_chr)
#
#     # Calculate distances for RBH
#     mp_dist <- .getDist(gff_ls = gff_ls, anchor_ls = anchor_ls, miniprot = TRUE)
#     dist_threshold <- as.vector(h5$synog_tx$dist_threshold)
#     mp_dist$so_valid <- mp_dist$dist <= dist_threshold
#
#     # Create a data frame for syntenic orthologs
#     syn_og <- .makeSynogDF(rbh = mp_dist, h5 = h5)
#     syn_og$rbbh <- TRUE
#
#     # Add Reciprocal Best BLAST Hits (RBBH) and filter
#     rbbh_og <- .addRBBH(rbbh = subset(mp_dist, !is.infinite(dist)),
#                         syn_og = syn_og)
#     rbbh_og <- subset(rbbh_og, select = c(qseqid,
#                                           sseqid,
#                                           qgeneid,
#                                           sgeneid))
#     rbbh_og$rbbh <- TRUE
#     rbbh_og$syntenic <- FALSE
#
#     # Combine the syntenic orthologs with the RBBH orthologs
#     syn_og <- rbind(subset(syn_og, select = c(qseqid,
#                                               sseqid,
#                                               syntenic,
#                                               rbbh,
#                                               qgeneid,
#                                               sgeneid)),
#                     rbbh_og)
#     syn_og$pair_id <- paste(syn_og$qseqid, syn_og$sseqid, sep = "_")
#     syn_og$class <- NA
#     syn_og$mutual_e <- NA
#     syn_og$mutual_ci <- NA
#     old_syn_og <- subset(h5$synog_tx$orthopairs, select = -OG)
#     old_syn_og$pair_id <- paste(old_syn_og$qseqid, old_syn_og$sseqid, sep = "_")
#     syn_og <- rbind(syn_og, old_syn_og)
#     syn_og$class <- .classifySO(ortho = syn_og)
#
#     # Identify orphan genes
#     orphan <- .getOrphan(gff_ls = gff_ls, syn_og = syn_og)
#
#     # Create the final output
#     out <- .makeOutput(syn_og = syn_og, orphan = orphan)
#
#     # Save results to the HDF5 file
#     .h5creategroup(object$h5,"synog_tx")
#     .h5overwrite(obj = out$syn_og, file = object$h5, "synog_tx/orthopairs")
#     .h5overwrite(obj = out$summary, file = object$h5, "synog_tx/summary")
#     .h5overwrite(obj = out$orphan, file = object$h5, "synog_tx/orphan")
#     .h5overwrite(obj = dist_threshold, file = object$h5, "synog_tx/dist_threshold")

    # Export the final GFF results to output files
    export.gff3(out$gff_ls$query_gff, as.vector(h5$protmap$s2q_gff))
    export.gff3(out$gff_ls$subject_gff, as.vector(h5$protmap$q2s_gff))

    hit <- names(out$query_cds) %in% out$gff_ls$query_gff$ID
    writeXStringSet(out$query_cds[hit], h5$protmap$s2q_cds)
    hit <- names(out$subject_cds) %in% out$gff_ls$subject_gff$ID
    writeXStringSet(out$subject_cds[hit], h5$protmap$q2s_cds)
}


# .makeLCBgr <- function(object, h5){
#     # Order the LCB pairs
#     pairs <- list(lcb_1to1 = .order(df = h5$sibeliaz$lcb_pairs$lcb_1to1),
#                   lcb_non_1to1 = .order(df = h5$sibeliaz$lcb_pairs$lcb_non_1to1))
#
#     # Convert data frames to GRanges objects
#     gr <- .df2gr(pairs = pairs)
#
#     # Create GRanges objects for gaps in the query and subject genomes
#     gap_gr <- list(query = .gapGR(gr = gr$query_1to1,
#                                   chrLen = object$genome$query),
#                    subject = .gapGR(gr = gr$subject_1to1,
#                                     chrLen = object$genome$subject))
#     out <- list(gr = gr,
#                 gap_gr = gap_gr)
#     return(out)
# }
#
# .remakeAnchor <- function(gff_ls, lcb_gr, h5){
#     df1 <- .makeMiniprotPairDF(gff1 = gff_ls$query_gff,
#                                gff2 = gff_ls$subject_gff)
#     df1 <- .filterMPtx(df = df1, gff_ls = gff_ls, query = TRUE)
#
#     df2 <- .makeMiniprotPairDF(gff1 = gff_ls$subject_gff,
#                                gff2 = gff_ls$query_gff)
#     df2 <- .filterMPtx(df = df2, gff_ls = gff_ls, query = FALSE)
#
#     new_anchor <- .getNewAnchors(df1 = df1,
#                                  df2 = df2,
#                                  h5 = h5,
#                                  lcb_gr = lcb_gr,
#                                  gff_ls = gff_ls)
#
#     out <- rbind(subset(h5$anchor, select = c(qseqid, sseqid)),
#                  subset(new_anchor, select = c(qseqid, sseqid)))
#
#     # Overwrite the "anchor" group in the HDF5 file with the filtered orthologs
#     .h5overwrite(obj = out, file = object$h5, "anchor")
# }
#
# .filterMPtx <- function(df, gff_ls, query = TRUE){
#     if(query){
#         df <- subset(df, subset = grepl("query_", gene))
#         df$tx_single <- .isSingleExon(tx = df$tx,
#                                       gff = gff_ls$query_gff)
#         df$target_tx_single <- .isSingleExon(tx = df$target_tx,
#                                              gff = gff_ls$subject_gff)
#         df <- subset(df, subset = !tx_single & !target_tx_single)
#
#     } else {
#         df <- subset(df, subset = grepl("subject_", gene))
#         df$tx_single <- .isSingleExon(tx = df$tx,
#                                       gff = gff_ls$subject_gff)
#         df$target_tx_single <- .isSingleExon(tx = df$target_tx,
#                                              gff = gff_ls$query_gff)
#         df <- subset(df, subset = !tx_single & !target_tx_single)
#     }
#     return(df)
# }
#
# .isSingleExon <- function(tx, gff){
#     cds_i <- gff$type == "CDS"
#     query_cds_parents <- unlist(gff$Parent[cds_i])
#     n_query_cds_parents <- table(query_cds_parents)
#     hit <- match(tx, names(n_query_cds_parents))
#     out <- n_query_cds_parents[hit] == 1
#     return(out)
# }
#
# .getNewAnchors <- function(df1, df2, h5, lcb_gr, gff_ls){
#     df1 <- .omitPairedGenes(df = df1,
#                             geneid = h5$synog_gene$orthopairs$sgeneid,
#                             syntenic = h5$synog_gene$orthopairs$syntenic)
#     df2 <- .omitPairedGenes(df = df2,
#                             geneid = h5$synog_gene$orthopairs$qgeneid,
#                             syntenic = h5$synog_gene$orthopairs$syntenic)
#
#     rbbh <- data.frame(qseqid = c(df1$tx, df2$target_tx),
#                        sseqid = c(df1$target_tx, df2$tx))
#     rbbh$index <- seq_len(nrow(rbbh))
#
#     obj <- list(gr = lcb_gr$gr,
#                 gap_gr = lcb_gr$gap_gr,
#                 rbbh = rbbh,
#                 gff = list(query = gff_ls$query_gff,
#                            subject = gff_ls$subject_gff))
#
#     # Filter orthologs in 1-to-1 LCBs
#     obj <- .orthoIn1to1lcb(obj = obj)
#
#     # Filter orthologs in non-1-to-1 LCBs if specified
#     obj <- .orthoInNon1to1lcb(obj = obj)
#
#     # Filter orthologs in 1-to-1 CBI
#     obj <- .orthoIn1to1cbi(obj = obj)
#
#     return(obj$out)
# }
#
# .omitPairedGenes <- function(df, geneid, syntenic){
#     hit <- match(df$target_gene, geneid)
#     df$paired[!is.na(hit)] <- TRUE
#     df$syntenic <- as.logical(syntenic[hit])
#     df$synog <- df$paired & df$syntenic
#     df$synog[is.na(df$synog)] <- FALSE
#     df <- subset(df, subset = !synog)
#     return(df)
# }

.filterGFF <- function(object, gff_ls, h5){
    query_cds <- readDNAStringSet(object$query_cds)
    subject_cds <- readDNAStringSet(object$subject_cds)

    query_df <- .makeQueryDF(query_gff = gff_ls$query_gff,
                             subject_gff = gff_ls$subject_gff,
                             query_cds = query_cds)
    subject_df <- .makeSubjectDF(subject_gff = gff_ls$subject_gff,
                                 query_gff = gff_ls$query_gff,
                                 subject_cds = subject_cds)
    query_valid <- .filterMPgenes(df = query_df, pattern = "^query_")
    subject_valid <- .filterMPgenes(df = subject_df, pattern = "^subject_")
    gff_ls$query_gff <- .removeInvalidMP(gff = gff_ls$query_gff,
                                         valid_id = query_valid,
                                         pattern = "^query_")
    gff_ls$subject_gff <- .removeInvalidMP(gff = gff_ls$subject_gff,
                                           valid_id = subject_valid,
                                           pattern = "^subject_")
    out <- list(gff_ls = gff_ls,
                query_cds = query_cds,
                subject_cds = subject_cds)
    return(out)
}


.makeQueryDF <- function(query_gff, subject_gff, query_cds){
    query_cds_init <- .checkInitCodon(cds = query_cds)
    query_cds_term <- .checkTermCodon(cds = query_cds)
    query_gff_tx_index <- query_gff$type %in% c("transcript", "mRNA")
    query_df <- data.frame(tx = query_gff$ID[query_gff_tx_index],
                           gene = query_gff$gene_id[query_gff_tx_index])
    hit <- match(query_df$tx, names(query_cds_init))
    query_df$init <- query_cds_init[hit]
    query_df$term <- query_cds_term[hit]
    query_df$valid <- query_df$init & query_df$term
    query_df$len <- width(query_cds)[hit]
    hit <- match(query_df$tx, query_gff$ID)
    query_df$target <- sub("\\s.+", "", query_gff$Target[hit])
    hit <- match(query_df$target, subject_gff$ID)
    query_df$target <- subject_gff$gene_id[hit]
    query_df$pair_id <- paste(query_df$gene, query_df$target, sep = "_")
    return(query_df)
}

.makeSubjectDF <- function(subject_gff, query_gff, subject_cds){
    subject_cds_init <- .checkInitCodon(cds = subject_cds)
    subject_cds_term <- .checkTermCodon(cds = subject_cds)
    subject_gff_tx_index <- subject_gff$type %in% c("transcript", "mRNA")
    subject_df <- data.frame(tx = subject_gff$ID[subject_gff_tx_index],
                             gene = subject_gff$gene_id[subject_gff_tx_index])
    hit <- match(subject_df$tx, names(subject_cds_init))
    subject_df$init <- subject_cds_init[hit]
    subject_df$term <- subject_cds_term[hit]
    subject_df$valid <- subject_df$init & subject_df$term
    subject_df$len <- width(subject_cds)[hit]
    hit <- match(subject_df$tx, subject_gff$ID)
    subject_df$target <- sub("\\s.+", "", subject_gff$Target[hit])
    hit <- match(subject_df$target, query_gff$ID)
    subject_df$target <- query_gff$gene_id[hit]
    subject_df$pair_id <- paste(subject_df$gene, subject_df$target, sep = "_")
    return(subject_df)
}

.checkInitCodon <- function(cds){
    return(substr(cds, 1, 3) == "ATG")
}

.checkTermCodon <- function(cds){
    len <- width(cds)
    return(substr(cds, len - 2, len) %in% c("TAA", "TGA", "TAG"))
}

.filterMPgenes <- function(df, pattern = "^query_"){
    df <- df[order(df$gene), ]
    mp_tx_index <- grepl(pattern, df$tx)
    eval <- data.frame(id = unique(df$gene))

    mp_valid_max_len <- tapply(df$len[mp_tx_index][df$valid[mp_tx_index]],
                               df$gene[mp_tx_index][df$valid[mp_tx_index]],
                               max)
    hit <- match(eval$id, names(mp_valid_max_len))
    eval$mp_valid_max_len <- mp_valid_max_len[hit]
    eval$mp_valid_max_len[is.na(eval$mp_valid_max_len)] <- 0

    non_mp_valid_max_len <- tapply(df$len[!mp_tx_index][df$valid[!mp_tx_index]],
                                   df$gene[!mp_tx_index][df$valid[!mp_tx_index]],
                                   max)
    hit <- match(eval$id, names(non_mp_valid_max_len))
    eval$non_mp_valid_max_len <- non_mp_valid_max_len[hit]
    eval$non_mp_valid_max_len[is.na(eval$non_mp_valid_max_len)] <- 0

    max_tx_index <- tapply(df$len[df$valid],
                           df$gene[df$valid],
                           which.max)
    max_tx_index <- sapply(max_tx_index, function(x){
        if(length(x) == 0){
            x <- NA
        }
        return(x)
    })
    hit <- match(eval$id, names(max_tx_index))
    eval$max_tx_index <- max_tx_index[hit]

    max_tx_is_na <- is.na(eval$max_tx_index)
    hit <- match(grep(pattern, eval$id[max_tx_is_na], value = TRUE),
                 df$gene)
    max_tx_index <- tapply(df$len[hit],
                           df$gene[hit],
                           which.max)
    hit <- match(eval$id[max_tx_is_na], names(max_tx_index))
    eval$max_tx_index[max_tx_is_na] <- max_tx_index[hit]

    start_index <- match(df$gene, df$gene)
    hit <- match(df$gene, eval$id)
    df$index <- eval$max_tx_index[hit] + start_index - 1
    out <- grep(pattern, df$tx[df$index], value = TRUE)
    return(out)
}

.removeInvalidMP <- function(gff, valid_id, pattern){
    parent <- sapply(gff$Parent, function(x){
        if(length(x) == 0){
            x <- NA
        }
        return(x)
    })
    valid_tx <- gff$ID %in% valid_id
    valid_elemetns <- parent %in% valid_id
    non_mp_entries <- !grepl(pattern, gff$ID)
    valid_entries <- valid_tx | valid_elemetns | non_mp_entries
    return(gff[valid_entries])
}

