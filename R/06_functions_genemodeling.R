#' Map Proteins Between Genomes
#'
#' This function maps proteins between the query and subject genomes using the specified miniprot binary and number of cores.
#'
#' @param object A SynogDB object.
#' @param out_prefix A character string specifying the prefix for the output files.
#' @param miniprot_bin A character string specifying the path to the miniprot binary.
#' @param n_core An integer specifying the number of cores to use for the mapping.
#' @param len_diff A numeric value specifying the maximum allowable length difference for proteins. Default is 0.5.
#'
#' @return None. The function performs the mapping and writes the results to the output files specified by `out_prefix`.
#'
#' @export
mapProt <- function(object, out_prefix, miniprot_bin = "miniprot",
                    conda = "conda",
                    condaenv = NULL, n_core, len_diff = 0.5){
    stopifnot(inherits(x = object, "SynogDB"))

    # Call the mapping engine function with specified parameters
    .mapEngine(object = object,
               subject_prot = object$subject_prot,
               query_prot = object$query_prot,
               out_prefix = out_prefix,
               miniprot_bin = miniprot_bin,
               conda = conda,
               condaenv = condaenv,
               n_core = n_core,
               overlap = TRUE,
               len_diff = len_diff)
}

#' Map Orphan Proteins Between Genomes
#'
#' This function maps orphan proteins between the query and subject genomes using the specified miniprot binary and number of cores.
#'
#' @param object A SynogDB object.
#' @param out_prefix A character string specifying the prefix for the output files.
#' @param miniprot_bin A character string specifying the path to the miniprot binary.
#' @param n_core An integer specifying the number of cores to use for the mapping.
#' @param len_diff A numeric value specifying the maximum allowable length difference for proteins. Default is 0.5.
#'
#' @return None. The function performs the mapping and writes the results to the output files specified by `out_prefix`.
#'
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
mapOrphan <- function(object, out_prefix, miniprot_bin, conda,
                      condaenv, n_core, len_diff = 0.5){
    stopifnot(inherits(x = object, "SynogDB"))

    # Open the HDF5 file and ensure it gets closed on exit
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))

    # Check if genewise ortholog information exists
    if(!H5Lexists(h5, "synog_gene/orthopairs")){
        stop("Run geneOrtho to obtain genewise ortholog info.")
    }

    # Define file paths for query and subject orphan protein sequences
    q_aa_fn <- paste0(out_prefix, "miniprot_query_prot.fa")
    s_aa_fn <- paste0(out_prefix, "miniprot_subject_prot.fa")

    # Write orphan protein sequences to FASTA files
    .writeProtFASTA(id = h5$synog_tx$orphan$query$qseqid,
                    prot_fn = object$query_prot,
                    gff_fn = object$query_gff,
                    fn = q_aa_fn)

    .writeProtFASTA(id = h5$synog_tx$orphan$subject$sseqid,
                    prot_fn = object$subject_prot,
                    gff_fn = object$subject_gff,
                    fn = s_aa_fn)

    # Call the mapping engine function with specified parameters
    .mapEngine(object = object,
               subject_prot = s_aa_fn,
               query_prot = q_aa_fn,
               out_prefix = out_prefix,
               miniprot_bin = miniprot_bin,
               conda = conda,
               condaenv = condaenv,
               n_core = n_core,
               overlap = TRUE,
               len_diff = len_diff)
}

#' Execute Protein Mapping Between Genomes
#'
#' This function maps proteins between the query and subject genomes using the specified miniprot binary and number of cores. It performs both directions of mapping and organizes the results.
#'
#' @importFrom rtracklayer export.gff3
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
.mapEngine <- function(object, subject_prot, query_prot,
                       out_prefix, miniprot_bin, conda,
                       condaenv, n_core,
                       overlap, len_diff){
    # Define output file paths for both directions of mapping
    s2q_out <- paste0(out_prefix, "query_miniprot_out")
    q2s_out <- paste0(out_prefix, "subject_miniprot_out")

    # Run miniprot mapping from subject to query genome
    .miniprot(query_fn = subject_prot, genome_fn = object$query_genome,
              out_prefix = s2q_out, miniprot_bin = miniprot_bin,
              conda = conda, condaenv = condaenv,
              n_core = n_core)

    # Run miniprot mapping from query to subject genome
    .miniprot(query_fn = query_prot, genome_fn = object$subject_genome,
              out_prefix = q2s_out, miniprot_bin = miniprot_bin,
              conda = conda, condaenv = condaenv,
              n_core = n_core)

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
                                   "ID", "Name", "gene_id", "Parent")
    mcols(s2q_gff) <- m_s2q_gff[, hit]
    m_q2s_gff <- mcols(q2s_gff)
    hit <- names(m_q2s_gff) %in% c("source", "type", "score", "phase",
                                   "ID", "Name", "gene_id", "Parent")
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
                      n_core = 1){

    if(!is.null(condaenv)){
        .condaExe(conda = conda, env = condaenv, command = miniprot_bin,
                  args = paste("-t", n_core,
                               "-d", paste0(out_prefix, ".mpi"),
                               genome_fn))
        .condaExe(conda = conda, env = condaenv, command = miniprot_bin,
                  args = paste("-t", n_core, "--gff",
                               paste0(out_prefix, ".mpi"),
                               query_fn, ">", paste0(out_prefix, ".gff")))

    } else {
        system2(command = miniprot_bin,
                args = paste("-t", n_core,
                             "-d", paste0(out_prefix, ".mpi"),
                             genome_fn))
        system2(command = miniprot_bin,
                args = paste("-t", n_core,
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
