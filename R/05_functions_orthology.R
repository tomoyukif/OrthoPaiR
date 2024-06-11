#' Filter ortholog pairs based on Locally Collinear Blocks (LCBs)
#'
#' This function filters ortholog pairs based on their presence in Locally Collinear Blocks (LCBs).
#'
#' @param object A SynogDB object.
#' @param non1to1 Logical indicating whether to include non-1-to-1 LCB pairs (default is TRUE).
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
anchorOrtho <- function(object,
                        non1to1 = TRUE){
    # Check if the input object is of class "SynogDB"
    stopifnot(inherits(x = object, "SynogDB"))

    # Open the HDF5 file
    h5 <- H5Fopen(object$h5)
    # Ensure the HDF5 file is closed when the function exits
    on.exit(H5Fclose(h5))

    # Check if the necessary groups exist in the HDF5 file
    if(!H5Lexists(h5, "sibeliaz/lcb_pairs")){
        stop("Run getLCBpairs to obtain LCB pair info.")
    }
    if(!H5Lexists(h5, "blast/rbbh")){
        stop("Run rbh with `best = TRUE` to obtain RBBH info.")
    }

    # Create a working object from the SynogDB object and HDF5 file
    obj <- .makeObject(object = object, h5 = h5)

    # Filter orthologs in 1-to-1 LCBs
    obj <- .orthoIn1to1lcb(obj = obj)
    # Filter orthologs in non-1-to-1 LCBs if specified
    if (non1to1) {
        obj <- .orthoInNon1to1lcb(obj = obj)
    }
    # Filter orthologs in 1-to-1 CBI
    obj <- .orthoIn1to1cbi(obj = obj)

    # Overwrite the "anchor" group in the HDF5 file with the filtered orthologs
    .h5overwrite(obj = obj$out, file = object$h5, "anchor")
}#' Filter ortholog pairs based on Locally Collinear Blocks (LCBs)
#'
#' This function filters ortholog pairs based on their presence in Locally Collinear Blocks (LCBs).
#'
#' @param object A SynogDB object.
#' @param non1to1 Logical indicating whether to include non-1-to-1 LCB pairs (default is TRUE).
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
anchorOrtho <- function(object,
                        non1to1 = TRUE){
    # Check if the input object is of class "SynogDB"
    stopifnot(inherits(x = object, "SynogDB"))

    # Open the HDF5 file
    h5 <- H5Fopen(object$h5)
    # Ensure the HDF5 file is closed when the function exits
    on.exit(H5Fclose(h5))

    # Check if the necessary groups exist in the HDF5 file
    if(!H5Lexists(h5, "sibeliaz/lcb_pairs")){
        stop("Run getLCBpairs to obtain LCB pair info.")
    }
    if(!H5Lexists(h5, "blast/rbbh")){
        stop("Run rbh with `best = TRUE` to obtain RBBH info.")
    }

    # Create a working object from the SynogDB object and HDF5 file
    obj <- .makeObject(object = object, h5 = h5)

    # Filter orthologs in 1-to-1 LCBs
    obj <- .orthoIn1to1lcb(obj = obj)
    # Filter orthologs in non-1-to-1 LCBs if specified
    if (non1to1) {
        obj <- .orthoInNon1to1lcb(obj = obj)
    }
    # Filter orthologs in 1-to-1 CBI
    obj <- .orthoIn1to1cbi(obj = obj)

    # Overwrite the "anchor" group in the HDF5 file with the filtered orthologs
    .h5overwrite(obj = obj$out, file = object$h5, "anchor")
}

#' Create a working object for ortholog filtering
#'
#' This function creates a working object containing necessary data for ortholog filtering based on LCBs.
#'
#' @param object A SynogDB object.
#' @param h5 An HDF5 file handle.
#' @return A list containing the working object data.
.makeObject <- function(object, h5){
    # Order the LCB pairs
    pairs <- list(lcb_1to1 = .order(df = h5$sibeliaz$lcb_pairs$lcb_1to1),
                  lcb_non_1to1 = .order(df = h5$sibeliaz$lcb_pairs$lcb_non_1to1))

    # Convert data frames to GRanges objects
    gr <- .df2gr(pairs = pairs)

    # Create GRanges objects for gaps in the query and subject genomes
    gap_gr <- list(query = .gapGR(gr = gr$query_1to1,
                                  chrLen = object$genome$query),
                   subject = .gapGR(gr = gr$subject_1to1,
                                    chrLen = object$genome$subject))

    # Import GFF files for query and subject genomes
    query_gff <- .importAllGFF(object$query_gff)
    subject_gff <- .importAllGFF(object$subject_gff)
    gff <- .subsetGFF(query_gff = query_gff, subject_gff = subject_gff)

    # Get the Reciprocal Best BLAST Hits (RBBH) from the HDF5 file
    rbbh <- h5$blast$rbbh
    rbbh$index <- seq_along(rbbh$qseqid)

    # Create and return the working object
    out <- list(gr = gr,
                rbbh = rbbh,
                gap_gr = gap_gr,
                gff = gff)
    return(out)
}

#' Order LCB pairs
#'
#' This function orders Locally Collinear Block (LCB) pairs based on query or subject genome coordinates.
#'
#' @param df A data.frame containing LCB pairs.
#' @param by A character string indicating whether to order by "query" or "subject" (default is "query").
#'
#' @return An ordered data.frame of LCB pairs.
.order <- function(df, by = "query"){
    if(by == "query"){
        # Order by query chromosome and start position
        q_order <- order(df$query_chr, df$query_start)
        df <- df[q_order, ]
    } else {
        # Order by subject chromosome and start position
        s_order <- order(df$subject_chr, df$subject_start)
        df <- df[s_order, ]
    }
    return(df)
}

#' Convert LCB pairs to GRanges objects
#'
#' This function converts Locally Collinear Block (LCB) pairs to GRanges objects for both query and subject genomes.
#'
#' @param pairs A list containing LCB pairs data frames.
#'
#' @return A list of GRanges objects for query and subject genomes.
#'
.df2gr <- function(pairs){
    gr <- list(
        # Create GRanges object for 1-to-1 LCB pairs in the query genome
        query_1to1 = .makeGRanges(df = pairs$lcb_1to1, genome = "query"),
        # Create GRanges object for 1-to-1 LCB pairs in the subject genome
        subject_1to1 = .makeGRanges(df = pairs$lcb_1to1, genome = "subject"),
        # Create GRanges object for non-1-to-1 LCB pairs in the query genome
        query_non1to1 = .makeGRanges(df = pairs$lcb_non_1to1, genome = "query"),
        # Create GRanges object for non-1-to-1 LCB pairs in the subject genome
        subject_non1to1 = .makeGRanges(df = pairs$lcb_non_1to1, genome = "subject")
    )
    return(gr)
}

#' Convert LCB pairs to GRanges object
#'
#' This function converts Locally Collinear Block (LCB) pairs to a GRanges object for either the query or subject genome.
#'
#' @param df A data.frame containing LCB pairs.
#' @param genome A character string indicating whether to create GRanges for "query" or "subject".
#' @return A GRanges object for the specified genome.
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#'
.makeGRanges <- function(df, genome){
    if(genome == "query"){
        # Handle cases where start position is greater than end position
        minus <- df$query_start > df$query_end
        tmp <- df$query_start[minus]
        df$query_start[minus] <- df$query_end[minus]
        df$query_end[minus] <- tmp

        # Create GRanges object for query genome
        gr <- GRanges(seqnames = df$query_chr,
                      ranges = IRanges(start = df$query_start,
                                       end = df$query_end),
                      strand = "*")

    } else {
        # Handle cases where start position is greater than end position
        minus <- df$subject_start > df$subject_end
        tmp <- df$subject_start[minus]
        df$subject_start[minus] <- df$subject_end[minus]
        df$subject_end[minus] <- tmp

        # Create GRanges object for subject genome
        gr <- GRanges(seqnames = df$subject_chr,
                      ranges = IRanges(start = df$subject_start,
                                       end = df$subject_end),
                      strand = "*")
    }
    return(gr)
}

#' Subset and order GFF files
#'
#' This function subsets GFF files to include only transcripts and mRNA, and orders them.
#'
#' @param query_gff A GRanges object containing the query GFF data.
#' @param subject_gff A GRanges object containing the subject GFF data.
#'
#' @return A list containing the subsetted and ordered query and subject GFF data.
#'
.subsetGFF <- function(query_gff, subject_gff){
    # Subset the GFF data to include only "transcript" and "mRNA" types
    query_gff <- query_gff[query_gff$type %in% c("transcript", "mRNA")]
    subject_gff <- subject_gff[subject_gff$type %in% c("transcript", "mRNA")]

    # Order the GFF data
    query_gff <- .orderGFF(gff = query_gff)
    subject_gff <- .orderGFF(gff = subject_gff)

    return(list(query = query_gff, subject = subject_gff))
}

#' Order GFF data
#'
#' This function orders GFF entries based on chromosome and start position.
#'
#' @param gff A GRanges object containing the GFF data.
#'
#' @return An ordered GRanges object.
#'
.orderGFF <- function(gff){
    # Order GFF data by chromosome and start position
    gff_order <- order(as.character(seqnames(gff)), start(gff))
    gff <- gff[gff_order]
    return(gff)
}

#' Create GRanges object for genomic gaps
#'
#' This function creates a GRanges object for genomic gaps between Locally Collinear Blocks (LCBs).
#'
#' @param gr A GRanges object containing LCBs.
#' @param chrLen A data.frame containing chromosome lengths.
#'
#' @return A GRanges object representing the gaps between LCBs.
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#'
.gapGR <- function(gr, chrLen){
    n_gr <- length(gr)
    # Create a data.frame for the gaps between LCBs
    gap_gr <- data.frame(chr_start = as.character(seqnames(gr[-n_gr])),
                         chr_end = as.character(seqnames(gr[-1])),
                         start = end(gr[-n_gr]) + 1,
                         end = start(gr[-1]) - 1,
                         start_block = seq_along(gr)[-n_gr],
                         end_block = seq_along(gr)[-1])

    # Include gaps at the start and end of chromosomes
    gap_gr <- rbind(
        subset(gap_gr, subset = chr_start == chr_end),
        .gapGrStartEdge(gr = gr, gap_gr = gap_gr),
        .gapGrEndEdge(gr = gr, gap_gr = gap_gr, chrLen = chrLen, n_gr = n_gr)
    )

    # Flip minus strand gaps
    gap_gr <- .flipMinusStrand(df = gap_gr)

    # Create a GRanges object from the gap data
    gap_gr <- GRanges(seqnames = gap_gr$chr_start,
                      ranges = IRanges(start = gap_gr$start, end = gap_gr$end),
                      start_block = gap_gr$start_block,
                      end_block = gap_gr$end_block)

    # Filter out gaps with width less than or equal to 1
    gap_gr <- gap_gr[width(gap_gr) > 1]
    return(gap_gr)
}

#' Create GRanges for start edges of gaps
#'
#' This function creates a GRanges object for the gaps at the start edges of chromosomes.
#'
#' @param gr A GRanges object containing LCBs.
#' @param gap_gr A data.frame containing gap information.
#'
#' @return A data.frame with updated gap information for start edges.
#'
#' @importFrom GenomeInfoDb seqnames
#'
.gapGrStartEdge <- function(gr, gap_gr){
    # Subset gaps where chromosome start and end differ
    out <- subset(gap_gr, subset = chr_start != chr_end)
    if(nrow(out) == 0){
        return(NULL)
    }
    out$chr_start <- out$chr_end
    out$start <- 1
    out$start_block <- NA

    # Add gap at the start of the first chromosome
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

#' Create GRanges for end edges of gaps
#'
#' This function creates a GRanges object for the gaps at the end edges of chromosomes.
#'
#' @param gr A GRanges object containing LCBs.
#' @param gap_gr A data.frame containing gap information.
#' @param chrLen A data.frame containing chromosome lengths.
#' @param n_gr The number of GRanges in the input.
#'
#' @return A data.frame with updated gap information for end edges.
#'
#' @importFrom GenomeInfoDb seqnames
#'
.gapGrEndEdge <- function(gr, gap_gr, chrLen, n_gr){
    # Subset gaps where chromosome start and end differ
    out <- subset(gap_gr, subset = chr_start != chr_end)
    if(nrow(out) == 0){
        return(NULL)
    }
    out$chr_end <- out$chr_start
    out$end <- chrLen$length[match(out$chr_start, chrLen$names)]

    # Add gap at the end of the last chromosome
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

#' Flip start and end positions for minus strand gaps
#'
#' This function flips the start and end positions for gaps on the minus strand.
#'
#' @param df A data.frame containing gap information.
#'
#' @return A data.frame with flipped start and end positions for minus strand gaps.
#'
.flipMinusStrand <- function(df){
    minus <- df$start > df$end
    tmp <- df$start[minus]
    df$start[minus] <- df$end[minus]
    df$end[minus] <- tmp
    return(df)
}

#' Filter orthologs in 1-to-1 LCBs
#'
#' This function filters ortholog pairs that are in 1-to-1 Locally Collinear Blocks (LCBs).
#'
#' @param obj A working object containing necessary data.
#'
#' @return A working object with filtered orthologs.
#'
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

#' Filter orthologs in non-1-to-1 LCBs
#'
#' This function filters ortholog pairs that are in non-1-to-1 Locally Collinear Blocks (LCBs).
#'
#' @param obj A working object containing necessary data.
#'
#' @return A working object with filtered orthologs.
#'
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

#' Filter orthologs in 1-to-1 CBI
#'
#' This function filters ortholog pairs that are in 1-to-1 Chromosome Break Intervals (CBI).
#'
#' @param obj A working object containing necessary data.
#'
#' @return A working object with filtered orthologs.
#'
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

#' Filter ortholog pairs in blocks
#'
#' This function filters ortholog pairs that are within specified blocks.
#'
#' @param gr_query A GRanges object for query genome.
#' @param gr_subject A GRanges object for subject genome.
#' @param rbbh A data.frame containing reciprocal best BLAST hits.
#' @param gff_query A GRanges object containing query GFF data.
#' @param gff_subject A GRanges object containing subject GFF data.
#' @param interval Logical indicating whether to use intervals (default is FALSE).
#'
#' @return A data.frame containing the filtered ortholog pairs.
#'
.orthoInBlock <- function(gr_query,
                          gr_subject,
                          rbbh,
                          gff_query,
                          gff_subject,
                          interval = FALSE){
    # Map genes to blocks for query and subject genomes
    q_g2b <-  .gene2block(gr = gr_query, rbbh_id = rbbh$qseqid, gff = gff_query)
    s_g2b <-  .gene2block(gr = gr_subject, rbbh_id = rbbh$sseqid, gff = gff_subject)

    # If intervals are specified, get the intervals for the blocks
    if(interval){
        q_g2b <- .getInterval(g2b = q_g2b, gr = gr_query)
        s_g2b <- .getInterval(g2b = s_g2b, gr = gr_subject)
    }

    # Find valid ortholog pairs
    valid <- q_g2b == s_g2b
    out <- subset(rbbh, subset = valid)
    return(out)
}

#' Map genes to blocks
#'
#' This function maps genes to their corresponding blocks (LCBs or CBIs) based on overlaps.
#'
#' @importFrom GenomicRanges findOverlaps
.gene2block <- function(gr, rbbh_id, gff){
    # Find overlaps between GFF entries and blocks
    ol <- findOverlaps(gff, gr)

    # Initialize output vector with NA values
    out <- rep(NA, length(gff))

    # Assign block indices to overlapping genes
    out[queryHits(ol)] <- subjectHits(ol)

    # Match RBBH IDs to GFF IDs and return the block assignments
    gene_index <- match(rbbh_id, gff$ID)
    return(out[gene_index])
}

#' Get intervals for blocks
#'
#' This function gets intervals for blocks, considering non-NA values and creating interval strings.
#'
.getInterval <- function(g2b, gr){
    # Identify non-NA block assignments
    not_na <- !is.na(g2b)

    # Get start and end block indices for non-NA values
    start_int <- end_int <- g2b[not_na]
    start_int <- gr$start_block[start_int]
    end_int <- gr$end_block[end_int]

    # Initialize output vector with NA values
    out <- rep(NA, length(g2b))

    # Create interval strings for non-NA values
    out[not_na] <- paste(start_int, end_int, sep = "_")
    return(out)
}

################################################################################
#' Filter ortholog pairs based on gene synteny
#'
#' This function filters ortholog pairs based on gene synteny, using specified thresholds for percentage identity, e-value, and query coverage.
#'
#' @param object A SynogDB object.
#' @param omit_chr A character string specifying chromosomes to omit (default is "").
#' @param pident Minimum percentage identity for BLAST hits (default is 90).
#' @param evalue Maximum e-value for BLAST hits (default is 1e-50).
#' @param qcovs Minimum query coverage for BLAST hits (default is 50).
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
syntenyOrtho <- function(object,
                         omit_chr = "",
                         pident = 90,
                         evalue = 1e-50,
                         qcovs = 50){
    stopifnot(inherits(x = object, "SynogDB"))

    # Open the HDF5 file
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))

    # Check for necessary data in the HDF5 file
    if(!H5Lexists(h5, "anchor")){
        stop("Run anchorOrtho to obtain ortholog anchors info.")
    }
    if(!H5Lexists(h5, "blast/rbh")){
        stop("Run rbh with `best = FALSE` to obtain RBH info.")
    }

    # Get the GFF lists
    gff_ls <- .getGFFlist(object = object)
    # Create anchor lists
    anchor_ls <- .makeAnchorList(gff_ls = gff_ls, h5 = h5, omit_chr = omit_chr)
    # Calculate distances for RBH
    rbh_dist <- .getDist(gff_ls = gff_ls, h5 = h5, anchor_ls = anchor_ls)

    rbh_dist <- .findSingleExonGene(rbh_dist = rbh_dist, gff_ls = gff_ls)
    dist_threshold <- .findDistThreshold(rbh_dist = rbh_dist)
    rbh_dist$so_valid <- rbh_dist$dist <= dist_threshold
    message("The threshold of syntenic gene distance was set to ",
            dist_threshold)

    # Create a data frame for syntenic orthologs
    syn_og <- .makeSynogDF(rbh = rbh_dist, h5 = h5)
    syn_og <- .getMutualScore(ortho = syn_og)
    syn_og$class <- .classifySO(ortho = syn_og)

    # Add Reciprocal Best BLAST Hits (RBBH) and filter
    rbbh_og <- .addRBBH(rbbh = h5$blast$rbbh, syn_og = syn_og)
    rbbh_og <- .filterOG(target = rbbh_og, ortho = rbh_dist, pident = pident, evalue = evalue, qcovs = qcovs)
    # Merge syntenic ortholog data frames
    syn_og <- .mergeSynogDF(syn_og = syn_og, rbbh_og = rbbh_og, gff_ls = gff_ls)

    syn_og$class <- .classifySO(ortho = syn_og)
    # Identify orphan genes
    orphan <- .getOrphan(gff_ls = gff_ls, syn_og = syn_og)
    # Create the final output
    out <- .makeOutput(syn_og = syn_og, orphan = orphan)
    # Save results to the HDF5 file
    .h5creategroup(object$h5,"synog_tx")
    .h5overwrite(obj = out$syn_og, file = object$h5, "synog_tx/orthopairs")
    .h5overwrite(obj = out$summary, file = object$h5, "synog_tx/summary")
    .h5overwrite(obj = out$orphan, file = object$h5, "synog_tx/orphan")
    .h5overwrite(obj = dist_threshold, file = object$h5, "synog_tx/dist_threshold")
}

#' Get GFF List
#'
#' This function imports and orders GFF files for query and subject genomes.
.getGFFlist <- function(object){
    # Import GFF files for query and subject genomes
    query_gff <- .importAllGFF(object$query_gff)
    subject_gff <- .importAllGFF(object$subject_gff)

    # Order the GFF data
    query_gff <- .orderGFF(gff = query_gff)
    subject_gff <- .orderGFF(gff = subject_gff)

    # Return the ordered GFF data as a list
    out <- list(query_gff = query_gff, subject_gff = subject_gff)
    return(out)
}

#' Make Anchor List
#'
#' This function creates a list of nearest anchors for query and subject genomes.
.makeAnchorList <- function(gff_ls, h5, omit_chr){
    anchor <- h5$anchor

    # Find nearest anchors for query and subject genomes
    query_anchor <- .findNearestAnchor(gff = gff_ls$query_gff, anchor = anchor$qseqid, omit_chr = omit_chr)
    subject_anchor <- .findNearestAnchor(gff = gff_ls$subject_gff, anchor = anchor$sseqid, omit_chr = omit_chr)

    # Set gene IDs to anchors
    gene_id <- .setGeneID2Anchor(anchor = anchor, gff_ls = gff_ls)

    # Combine anchor data with gene IDs
    anchor <- cbind(anchor, gene_id)

    # Organize anchor data
    out <- .orgAnchor(anchor = anchor, query_anchor = query_anchor, subject_anchor = subject_anchor)
    return(out)
}

.findSingleExonGene <- function(rbh_dist, gff_ls){
    query_tx_ids <- unique(rbh_dist$qseqid)
    cds_i <- gff_ls$query_gff$type == "CDS"
    query_cds_parents <- unlist(gff_ls$query_gff$Parent[cds_i])
    n_query_cds_parents <- table(query_cds_parents)
    hit <- match(rbh_dist$qseqid, names(n_query_cds_parents))
    rbh_dist$q_sigle_exon <- n_query_cds_parents[hit] == 1

    subject_tx_ids <- unique(rbh_dist$sseqid)
    cds_i <- gff_ls$subject_gff$type == "CDS"
    subject_cds_parents <- unlist(gff_ls$subject_gff$Parent[cds_i])
    n_subject_cds_parents <- table(subject_cds_parents)
    hit <- match(rbh_dist$sseqid, names(n_subject_cds_parents))
    rbh_dist$s_sigle_exon <- n_subject_cds_parents[hit] == 1

    return(rbh_dist)
}

.findDistThreshold <- function(rbh_dist){
    rbh_dist <- .getMutualScore(ortho = rbh_dist)
    rbh_dist <- subset(rbh_dist, subset = !is.infinite(dist))
    dup_q <- rbh_dist$qseqid[duplicated(rbh_dist$qseqid)]
    dup_s <- rbh_dist$sseqid[duplicated(rbh_dist$sseqid)]
    dist_zero <- rbh_dist$dist == 0
    dup_ids <- rbh_dist$qseqid %in% dup_q | rbh_dist$sseqid %in% dup_s
    single_exon <- rbh_dist$q_sigle_exon | rbh_dist$s_sigle_exon
    dist_eval_filt <- !(dist_zero | single_exon)
    dist_vs_ci <- sapply(seq_len(100),
                         function(i){
                             cond <- dist_eval_filt & rbh_dist$dist <= i
                             return(mean(rbh_dist$mutual_ci[cond]))
                         })
    out <- which.max(dist_vs_ci)
    return(out)
}

#' Find Nearest Anchor
#'
#' This function finds the nearest anchor for each gene in the GFF data.
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start
.findNearestAnchor <- function(gff, anchor, omit_chr){
    # Subset GFF data to transcripts and mRNAs
    tx <- gff[gff$type %in% c("transcript", "mRNA")]

    # Find the parent gene for each anchor
    anchor <- unlist(tx$Parent[match(anchor, tx$ID)])
    anchor <- na.omit(anchor)

    # Subset GFF data to genes and order by chromosome and start position
    gff <- gff[gff$type == "gene"]
    gff <- gff[order(as.numeric(seqnames(gff)), start(gff))]

    # Convert GFF data to a data frame and exclude specified chromosomes
    gff <- data.frame(chr = as.character(seqnames(gff)), start = start(gff), ID = gff$ID)
    gff <- gff[!grepl(omit_chr, gff$chr), ]

    # Initialize anchor and nearest anchor columns
    gff$anchor <- FALSE
    gff$anchor[gff$ID %in% anchor] <- TRUE
    gff$nearest_anchor <- vapply(X = seq_along(gff$anchor), FUN.VALUE = numeric(1), FUN = .getNearestAnchorIndex, gff = gff)

    # Add a location column
    gff$location <- seq_along(gff$chr)
    return(gff)
}

#' Get Nearest Anchor Index
#'
#' This function finds the nearest anchor index for a given gene.
#'
.getNearestAnchorIndex <- function(index, gff){
    # Find the nearest anchor on the same chromosome
    target <- gff$chr %in% gff$chr[index] & gff$anchor
    near <- which.min(abs(gff$start[target] - gff$start[index]))
    near <- which(target)[near]
    return(near)
}

#' Get Nearest Anchor Index
#'
#' This function finds the nearest anchor index for a given gene.
#'
.getNearestAnchorIndex <- function(index, gff){
    # Identify target genes on the same chromosome that are anchors
    target <- gff$chr %in% gff$chr[index] & gff$anchor

    # Find the nearest anchor based on the start position
    near <- which.min(abs(gff$start[target] - gff$start[index]))
    near <- which(target)[near]
    return(near)
}

#' Set Gene IDs to Anchors
#'
#' This function sets the gene IDs for query and subject anchors.
#'
#' @importFrom GenomeInfoDb seqnames
.setGeneID2Anchor <- function(anchor, gff_ls){
    # Get gene IDs for query anchors
    tx_i <- gff_ls$query_gff$type %in% c("transcript", "mRNA")
    tx <- gff_ls$query_gff[tx_i]
    hit <- match(anchor$qseqid, tx$ID)
    qgeneid <- unlist(tx$Parent[hit])
    qchr <- as.character(seqnames(tx[hit]))

    # Get gene IDs for subject anchors
    tx_i <- gff_ls$subject_gff$type %in% c("transcript", "mRNA")
    tx <- gff_ls$subject_gff[tx_i]
    hit <- match(anchor$sseqid, tx$ID)
    sgeneid <- unlist(tx$Parent[hit])
    schr <- as.character(seqnames(tx[hit]))

    # Return a data frame with gene IDs and chromosome information
    out <- data.frame(qgeneid = qgeneid, sgeneid = sgeneid, qchr = qchr, schr = schr)
    return(out)
}

#' Organize Anchor Data
#'
#' This function organizes anchor data for query and subject genomes.
#'
.orgAnchor <- function(anchor, query_anchor, subject_anchor){
    # Subset and unique anchor data
    anchor <- subset(anchor, select = qgeneid:schr)
    anchor <- unique(anchor)

    # Set query anchor information
    hit <- match(anchor$qgeneid, query_anchor$ID)
    anchor$query_anchor <- query_anchor$nearest_anchor[hit]
    tbl <- table(anchor$query_anchor)
    multiple <- anchor$query_anchor %in% names(tbl[tbl > 1])
    anchor$q2s_multiple <- FALSE
    anchor$q2s_multiple[multiple] <- TRUE

    # Set subject anchor information
    hit <- match(anchor$sgeneid, subject_anchor$ID)
    anchor$subject_anchor <- subject_anchor$nearest_anchor[hit]
    tbl <- table(anchor$subject_anchor)
    multiple <- anchor$subject_anchor %in% names(tbl[tbl > 1])
    anchor$s2q_multiple <- FALSE
    anchor$s2q_multiple[multiple] <- TRUE

    # Return a list with query, subject, and anchor data
    out <- list(query = query_anchor, subject = subject_anchor, anchor = anchor)
    return(out)
}

#' Get Distance between Anchors and RBH
#'
#' This function calculates the distance between anchors and reciprocal best hits (RBH).
#'
#' @importFrom Biobase rowMin
#'
.getDist <- function(gff_ls, h5, anchor_ls){
    # Extract RBH data
    rbh <- .extractRBH(h5 = h5, gff_ls = gff_ls)

    # Set anchor information for RBH
    anchor2rbh <- .setAnchor2RBH(rbh = rbh, anchor_ls = anchor_ls)
    rbh <- cbind(rbh, anchor2rbh)

    # Trace anchor distances for query to subject and subject to query
    q2s_dist <- .traceAnchorQ2S(rbh = rbh, anchor_ls = anchor_ls)
    s2q_dist <- .traceAnchorS2Q(rbh = rbh, anchor_ls = anchor_ls)

    # Calculate minimum distance
    rbh <- cbind(rbh, dist = rowMin(cbind(q2s_dist, s2q_dist)))
    return(rbh)
}

#' Extract Reciprocal Best Hits (RBH)
#'
#' This function extracts RBH information from the HDF5 object.
#'
.extractRBH <- function(h5, gff_ls){
    out <- h5$blast$rbh

    # Match query gene IDs and chromosome information
    q_hit <- match(out$qseqid, gff_ls$query_gff$ID)
    out$qgeneid <- gff_ls$query_gff$gene_id[q_hit]
    out$qchr <- as.character(seqnames(gff_ls$query_gff[q_hit]))

    # Match subject gene IDs and chromosome information
    s_hit <- match(out$sseqid, gff_ls$subject_gff$ID)
    out$sgeneid <- gff_ls$subject_gff$gene_id[s_hit]
    out$schr <- as.character(seqnames(gff_ls$subject_gff[s_hit]))

    # Create pair ID for each RBH
    out$pair_id <- paste(out$qgeneid, out$sgeneid, sep = "_")
    return(out)
}

#' Set Anchor Information to RBH
#'
#' This function sets the anchor information to the reciprocal best hits (RBH).
#'
.setAnchor2RBH <- function(rbh, anchor_ls){
    # Match gene IDs to nearest anchors for query and subject
    query_hit <- match(rbh$qgeneid, anchor_ls$query$ID)
    q2s_q_anchor <- anchor_ls$query$nearest_anchor[query_hit]
    s2q_q_location <- anchor_ls$query$location[query_hit]

    subject_hit <- match(rbh$sgeneid, anchor_ls$subject$ID)
    s2q_s_anchor <- anchor_ls$subject$nearest_anchor[subject_hit]
    q2s_s_location <- anchor_ls$subject$location[subject_hit]

    # Return data frame with anchor information
    out <- data.frame(q2s_q_anchor = q2s_q_anchor,
                      s2q_q_location = s2q_q_location,
                      s2q_s_anchor = s2q_s_anchor,
                      q2s_s_location = q2s_s_location)
    return(out)
}

#' Trace Anchor Distance from Query to Subject
#'
#' This function calculates the distance from query anchors to subject anchors.
#'
.traceAnchorQ2S <- function(rbh, anchor_ls){
    # Match query anchors to anchor data
    q2s_q_anchor_s <- match(rbh$q2s_q_anchor, anchor_ls$anchor$query_anchor)
    q2s_q_anchor_s[is.na(rbh$q2s_q_anchor)] <- NA

    # Identify multiple anchors
    multiple_anchor <- anchor_ls$anchor$q2s_multiple[q2s_q_anchor_s]
    multiple_anchor[is.na(multiple_anchor)] <- FALSE

    # Calculate single projection distance
    single_projection <- anchor_ls$anchor$subject_anchor[q2s_q_anchor_s]
    single_projection[multiple_anchor] <- NA
    q2s_dist <- abs(rbh$q2s_s_location - single_projection)

    # Handle chromosome mismatch
    single_projection_chr <- anchor_ls$anchor$schr[q2s_q_anchor_s]
    q2s_dist_chr_unmatch <- rbh$schr != single_projection_chr
    q2s_dist[q2s_dist_chr_unmatch] <- Inf

    # Calculate multiple projection distances
    multiple_projection_rbh <- rbh[multiple_anchor, ]
    q2s_dist_multi <- vapply(seq_along(multiple_projection_rbh$q2s_q_anchor),
                             .solveMultiProjectionQ2S,
                             FUN.VALUE = numeric(1),
                             rbh = multiple_projection_rbh, anchor_ls = anchor_ls)
    q2s_dist[multiple_anchor] <- q2s_dist_multi
    q2s_dist[is.na(q2s_dist)] <- Inf
    return(q2s_dist)
}

#' Trace Anchor Distance from Subject to Query
#'
#' This function calculates the distance from subject anchors to query anchors.
#'
.traceAnchorS2Q <- function(rbh, anchor_ls){
    # Match subject anchors to anchor data
    s2q_s_anchor_q <- match(rbh$s2q_s_anchor, anchor_ls$anchor$subject_anchor)
    s2q_s_anchor_q[is.na(rbh$s2q_s_anchor)] <- NA

    # Identify multiple anchors
    multiple_anchor <- anchor_ls$anchor$s2q_multiple[s2q_s_anchor_q]
    multiple_anchor[is.na(multiple_anchor)] <- FALSE

    # Calculate single projection distance
    single_projection <- anchor_ls$anchor$query_anchor[s2q_s_anchor_q]
    single_projection[multiple_anchor] <- NA
    s2q_dist <- abs(rbh$s2q_q_location - single_projection)

    # Handle chromosome mismatch
    single_projection_chr <- anchor_ls$anchor$qchr[s2q_s_anchor_q]
    s2q_dist_chr_unmatch <- rbh$qchr != single_projection_chr
    s2q_dist[s2q_dist_chr_unmatch] <- Inf

    # Calculate multiple projection distances
    multiple_projection_rbh <- rbh[multiple_anchor, ]
    s2q_dist_multi <- vapply(seq_along(multiple_projection_rbh$s2q_s_anchor),
                             .solveMultiProjectionS2Q,
                             FUN.VALUE = numeric(1),
                             rbh = multiple_projection_rbh, anchor_ls = anchor_ls)
    s2q_dist[multiple_anchor] <- s2q_dist_multi
    s2q_dist[is.na(s2q_dist)] <- Inf
    return(s2q_dist)
}

#' Solve Multiple Projections from Query to Subject
#'
#' This function solves multiple projections from query to subject anchors.
#'
.solveMultiProjectionQ2S <- function(index, rbh, anchor_ls){
    # Get query anchor and match to anchor data
    q_anchor <- rbh$q2s_q_anchor[index]
    hit <- which(anchor_ls$anchor$query_anchor == q_anchor)

    # Calculate distances for multiple projections
    projection <- anchor_ls$anchor$subject_anchor[hit]
    q2s_dist <- abs(rbh$q2s_s_location[index] - projection)
    projection_chr <- anchor_ls$anchor$schr[hit]
    q2s_dist_chr_unmatch <- rbh$schr[index] != projection_chr
    q2s_dist[q2s_dist_chr_unmatch] <- Inf

    # Return minimum distance
    out <- min(q2s_dist)
    return(out)
}

#' Solve Multiple Projections from Subject to Query
#'
#' This function solves multiple projections from subject to query anchors.
#'
.solveMultiProjectionS2Q <- function(index, rbh, anchor_ls){
    # Extract the subject anchor for the current index
    s_anchor <- rbh$s2q_s_anchor[index]

    # Identify the corresponding projections on the query genome
    hit <- which(anchor_ls$anchor$subject_anchor == s_anchor)
    projection <- anchor_ls$anchor$query_anchor[hit]

    # Calculate the distances between the subject anchor and the projections on the query genome
    s2q_dist <- abs(rbh$s2q_q_location[index] - projection)

    # Identify projection chromosome for comparison
    projection_chr <- anchor_ls$anchor$qchr[hit]

    # Set distances to Inf where the chromosomes do not match
    s2q_dist_chr_unmatch <- rbh$qchr[index] != projection_chr
    s2q_dist[s2q_dist_chr_unmatch] <- Inf

    # Find the minimum distance
    out <- min(s2q_dist)

    return(out)
}

#' Get Metrics for Ortholog Prediction
#'
#' This function calculates various metrics (precision, recall, specificity, accuracy, and F-measure)
#' for ortholog prediction based on the distances to anchors.
#'
#' @importFrom ggplot2 ggplot geom_line aes xlim scale_color_manual
.getMetrics <- function(rbh_dist, tp_list, fp_list){
    # Identify non-infinite distances and their metrics
    n_list <- unique(c(fp_list$query, fp_list$subject))
    n_hit <- rbh_dist$qgeneid %in% n_list | rbh_dist$sgeneid %in% n_list
    n_hit_index_dist <- rbh_dist$dist[n_hit]
    n_hit_index_dist <- tapply(n_hit_index_dist, rbh_dist$qgeneid[n_hit], min)
    n_tbl <- table(n_hit_index_dist)
    n_tbl <- n_tbl[names(n_tbl) != "Inf"]

    # Identify true positive hits and their metrics
    p_id <- paste(tp_list$query, tp_list$subject, sep = "_")
    p_id <- unique(p_id)
    p_hit <- which(rbh_dist$pair_id %in% p_id)
    p_hit_index_dist <- rbh_dist$dist[p_hit]
    p_hit_index_dist <- tapply(p_hit_index_dist, rbh_dist$pair_id[p_hit], min)
    p_tbl <- table(p_hit_index_dist)
    p_tbl <- p_tbl[names(p_tbl) != "Inf"]

    # Calculate metrics
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

    # Create data frame for plotting
    df <- rbind(data.frame(Metrics = "precision", Threshold = seq_along(precision), Score = precision),
                data.frame(Metrics = "recall", Threshold = seq_along(recall), Score = recall),
                data.frame(Metrics = "specificity", Threshold = seq_along(specificity), Score = specificity),
                data.frame(Metrics = "accuracy", Threshold = seq_along(accuracy), Score = accuracy),
                data.frame(Metrics = "f_measure", Threshold = seq_along(f_measure), Score = f_measure))

    # Plot metrics
    p <- ggplot(data = df) +
        geom_line(aes(x = Threshold, y = Score, group = Metrics, color = Metrics)) +
        xlim(0, 100) +
        scale_color_manual(values = c("blue", "darkgreen", "magenta", "skyblue", "darkorange"),
                           breaks = c("precision", "recall", "specificity", "accuracy", "f_measure"))

    # Return metrics and plot
    out <- list()
    out$so_metrics <- data.frame(threshold = seq_along(precision), precision = precision, recall = recall,
                                 specificity = specificity, accuracy = accuracy, f_measure = f_measure)
    out$metrics_plot <- p
    out$so_threshold <- which.max(f_measure)
    return(out)
}

#' Create Syntenic Ortholog Data Frame
#'
#' This function creates a data frame of syntenic orthologs based on RBH information.
#'
.makeSynogDF <- function(rbh, h5){
    # Filter RBH data to retain only valid syntenic orthologs
    out <- subset(rbh, subset = so_valid %in% TRUE,
                  select = -c(so_valid, qchr, schr, q2s_q_anchor, s2q_q_location, s2q_s_anchor, q2s_s_location, dist))

    # Mark these orthologs as syntenic
    out$syntenic <- TRUE

    # Initialize the RBBH column as FALSE
    out$rbbh <- FALSE

    # Create a unique pair identifier for each ortholog pair
    out$pair_id <- paste(out$qseqid, out$sseqid, sep = "_")

    # Create a unique pair identifier for RBBH from the anchor information in the HDF5 file
    rbbh_id <- paste(h5$anchor$qseqid, h5$anchor$sseqid, sep = "_")

    # Mark the pairs that are also RBBH
    out$rbbh[out$pair_id %in% rbbh_id] <- TRUE

    return(out)
}

#' Merge Syntenic Ortholog Data Frames
#'
#' This function merges syntenic ortholog data frames with RBBH orthologs.
#'
.mergeSynogDF <- function(syn_og, rbbh_og, gff_ls){
    syn_og <- subset(syn_og, select = -c(q_sigle_exon, s_sigle_exon))
    # Check if RBBH orthologs data frame is not empty
    if(nrow(rbbh_og) != 0){
        # Set gene IDs for the RBBH orthologs
        rbbh_og <- .setGeneIDsynog(rbbh_og = rbbh_og, gff_ls = gff_ls)

        # Mark these orthologs as not syntenic and set RBBH flag to TRUE
        rbbh_og$syntenic <- FALSE
        rbbh_og$rbbh <- TRUE

        # Create a unique pair identifier for each ortholog pair
        rbbh_og$pair_id <- paste(rbbh_og$qseqid, rbbh_og$sseqid, sep = "_")

        # Calculate mutual scores for the ortholog pairs
        rbbh_og <- .getMutualScore(ortho = rbbh_og)

        # Classify the ortholog pairs
        rbbh_og$class <- .classifySO(ortho = rbbh_og)

        # Combine the syntenic orthologs with the RBBH orthologs
        syn_og <- rbind(syn_og, rbbh_og)
    }
    return(syn_og)
}

#' Generate Output for Syntenic Ortholog Pairs
#'
#' This function generates the output for syntenic ortholog pairs, including a summary and orphan genes.
#'
.makeOutput <- function(syn_og, orphan){
    out <- list()  # Initialize the output list
    out$syn_og <- syn_og  # Assign syntenic orthologs to the output list
    out$syn_og$OG <- .numberingOG(ortho = syn_og)  # Number the ortholog groups
    out$syn_og <- subset(out$syn_og, select = c(qseqid, sseqid, mutual_e,
                                                mutual_ci, syntenic, rbbh,
                                                class, OG, qgeneid, sgeneid))  # Subset relevant columns
    rownames(out$syn_og) <- NULL  # Reset row names
    out$syn_og$syntenic <- as.character(out$syn_og$syntenic)  # Convert syntenic to character
    out$syn_og$rbbh <- as.character(out$syn_og$rbbh)  # Convert rbbh to character

    # Calculate the number of unique IDs for each class
    q_id <- tapply(out$syn_og$qseqid, out$syn_og$class, unique)
    n_q_id <- sapply(q_id, length)
    s_id <- tapply(out$syn_og$sseqid, out$syn_og$class, unique)
    n_s_id <- sapply(s_id, length)

    n_orphan <- lapply(orphan, nrow)  # Calculate the number of orphan genes

    # Calculate total numbers for query and subject
    q_tot <- sum(n_q_id) + n_orphan$query
    s_tot <- sum(n_s_id) + n_orphan$subject

    # Create a summary data frame
    df <- data.frame(query = c(q_tot, sum(n_q_id), n_q_id, n_orphan$query),
                     subject = c(s_tot, sum(n_s_id), n_s_id, n_orphan$subject))
    out$summary <- df  # Assign the summary to the output list
    out$orphan <- orphan  # Assign the orphan genes to the output list
    rownames(out$orphan$query) <- NULL  # Reset row names for orphan query genes
    rownames(out$orphan$subject) <- NULL  # Reset row names for orphan subject genes
    return(out)  # Return the output list
}

#' Classify Orthologs
#'
#' This function classifies ortholog pairs into categories: 1-to-1, 1-to-many, many-to-1, and many-to-many.
#'
.classifySO <- function(ortho){
    # Subset relevant columns for classification
    tmp <- subset(ortho, select = c(qseqid, sseqid, pair_id))

    # Identify duplicated query and subject sequences
    q_dup <- duplicated(tmp$qseqid)
    q_dup <- tmp$qseqid %in% tmp$qseqid[q_dup]
    s_dup <- duplicated(tmp$sseqid)
    s_dup <- tmp$sseqid %in% tmp$sseqid[s_dup]

    # Classify 1-to-1 orthologs
    og_11 <- tmp[!(q_dup | s_dup), ]

    # Identify and classify many-to-many orthologs
    mult <- tmp[q_dup | s_dup, ]
    mult_q_dup <- mult$qseqid %in% mult$qseqid[duplicated(mult$qseqid)]
    mult_s_dup <- mult$sseqid %in% mult$sseqid[duplicated(mult$sseqid)]
    og_mm <- mult[mult_q_dup & mult_s_dup, ]
    og_mm_i <- mult$qseqid %in% og_mm$qseqid | mult$sseqid %in% og_mm$sseqid
    og_mm <- mult[og_mm_i, ]

    # Identify and classify 1-to-many and many-to-1 orthologs
    not_mm <- mult[!og_mm_i, ]
    not_mm_q_dup <- not_mm$qseqid %in% not_mm$qseqid[duplicated(not_mm$qseqid)]
    not_mm_s_dup <- not_mm$sseqid %in% not_mm$sseqid[duplicated(not_mm$sseqid)]
    og_1m <- not_mm[not_mm_q_dup, ]
    og_m1 <- not_mm[not_mm_s_dup, ]

    # Combine all classifications into one data frame
    out <- NULL
    if(nrow(og_11) != 0){
        og_11$class <- "1to1"
        out <- rbind(out, og_11)
    }
    if(nrow(og_1m) != 0){
        og_1m$class <- "1toM"
        out <- rbind(out, og_1m)
    }
    if(nrow(og_m1) != 0){
        og_m1$class <- "Mto1"
        out <- rbind(out, og_m1)
    }
    if(nrow(og_mm) != 0){
        og_mm$class <- "MtoM"
        out <- rbind(out, og_mm)
    }

    # Match the classifications to the original ortholog pairs
    hit <- match(ortho$pair_id, out$pair_id)
    return(out$class[hit])
}

#' Get Orphan Transcripts
#'
#' This function identifies orphan transcripts that do not have orthologs.
#'
.getOrphan <- function(gff_ls, syn_og){
    # Identify orphan transcripts in the query genome
    tx_i <- gff_ls$query_gff$type %in% c("transcript", "mRNA")
    tx <- gff_ls$query_gff[tx_i]
    hit <- tx$ID %in% syn_og$qseqid
    q_orphan <- data.frame(ID = tx$ID[!hit], gene_id = tx$gene_id[!hit])

    # Identify orphan transcripts in the subject genome
    tx_i <- gff_ls$subject_gff$type %in% c("transcript", "mRNA")
    tx <- gff_ls$subject_gff[tx_i]
    hit <- tx$ID %in% syn_og$sseqid
    s_orphan <- data.frame(ID = tx$ID[!hit], gene_id = tx$gene_id[!hit])

    return(list(query = q_orphan, subject = s_orphan))
}

#' Calculate Mutual Scores for Orthologs
#'
#' This function calculates mutual coverage identity (CI) and mutual e-value for ortholog pairs.
#'
.getMutualScore <- function(ortho){
    # Calculate coverage identity for query to subject (q2s) and subject to query (s2q)
    q2s_covident <- ortho$q2s_pident * ortho$q2s_qcovs * 1e-2
    s2q_covident <- ortho$s2q_pident * ortho$s2q_qcovs * 1e-2

    # Calculate mutual e-value as the Euclidean distance of the individual e-values
    mutual_e <- sqrt(ortho$q2s_evalue^2 + ortho$s2q_evalue^2)

    # Calculate mutual coverage identity as the Euclidean distance of the individual coverage identities
    mutual_ci <- sqrt(q2s_covident^2 + s2q_covident^2)

    # Add mutual scores to the ortholog data frame
    ortho$mutual_e <- mutual_e
    ortho$mutual_ci <- mutual_ci

    return(ortho)
}

#' Add Reciprocal Best BLAST Hits (RBBH) to Syntenic Orthologs
#'
#' This function identifies candidate Reciprocal Best BLAST Hits (RBBH) that are not already present in the syntenic ortholog list.
#'
.addRBBH <- function(rbbh, syn_og){
    # Identify candidate RBBH that are not present in the syntenic ortholog list
    candidate_rbbh <- !(rbbh$qseqid %in% syn_og$qseqid | rbbh$sseqid %in% syn_og$sseqid)

    # Subset the RBBH data frame to include only the candidate RBBH
    out <- rbbh[candidate_rbbh, ]

    return(out)
}

#' Filter Ortholog Groups
#'
#' This function filters ortholog groups based on identity percentage, e-value, and query coverage.
#'
.filterOG <- function(target, ortho, pident, evalue, qcovs){
    # Create unique identifiers for ortholog pairs in both target and ortho data frames
    ortho_id <- paste(ortho$qseqid, ortho$sseqid, sep = "_")
    target_id <- paste(target$qseqid, target$sseqid, sep = "_")

    # Subset ortho data frame to include only the target pairs
    target_metrics <- ortho[ortho_id %in% target_id, ]

    # Apply filtering criteria based on percentage identity, query coverage, and e-value
    valid_og <- target_metrics$q2s_pident >= pident &
        target_metrics$s2q_pident >= pident &
        target_metrics$q2s_qcovs >= qcovs &
        target_metrics$s2q_qcovs >= qcovs &
        target_metrics$q2s_evalue <= evalue &
        target_metrics$s2q_evalue <= evalue

    # Subset the valid ortholog groups
    valid_og <- target_metrics[valid_og, ]

    # Create unique identifiers for the valid ortholog pairs
    valid_id <- paste(valid_og$qseqid, valid_og$sseqid, sep = "_")

    # Subset the target data frame to include only the valid ortholog pairs
    valid_og <- target[target_id %in% valid_id, ]

    return(valid_og)
}

#' Set Gene IDs for Syntenic Orthologs
#'
#' This function sets the gene IDs for query and subject sequences in the syntenic orthologs data frame.
#'
.setGeneIDsynog <- function(rbbh_og, gff_ls){
    # Match query sequence IDs with GFF to obtain gene IDs
    q_hit <- match(rbbh_og$qseqid, gff_ls$query_gff$ID)
    rbbh_og$qgeneid <- gff_ls$query_gff$gene_id[q_hit]

    # Match subject sequence IDs with GFF to obtain gene IDs
    s_hit <- match(rbbh_og$sseqid, gff_ls$subject_gff$ID)
    rbbh_og$sgeneid <- gff_ls$subject_gff$gene_id[s_hit]

    return(rbbh_og)
}

#' Number Ortholog Groups
#'
#' This function assigns unique identifiers to ortholog groups (OGs) based on the ortholog relationships between query and subject sequences.
#'
.numberingOG <- function(ortho){
    # Subset relevant columns
    tmp <- subset(ortho, select = c(qseqid, sseqid, pair_id))

    # Identify duplicated query and subject sequences
    q_dup <- duplicated(tmp$qseqid)
    q_dup <- tmp$qseqid %in% tmp$qseqid[q_dup]
    s_dup <- duplicated(tmp$sseqid)
    s_dup <- tmp$sseqid %in% tmp$sseqid[s_dup]

    # Identify 1-to-1 ortholog pairs
    og_11 <- tmp[!(q_dup | s_dup), ]
    og_11 <- og_11[order(og_11$qseqid), ]
    og_11$OG <- seq_along(og_11$qseqid)

    # Identify many-to-many (M-to-M) ortholog pairs
    mult <- tmp[q_dup | s_dup, ]
    mult_q_dup <- mult$qseqid %in% mult$qseqid[duplicated(mult$qseqid)]
    mult_s_dup <- mult$sseqid %in% mult$sseqid[duplicated(mult$sseqid)]

    og_mm <- mult[mult_q_dup & mult_s_dup, ]
    og_mm_i <- mult$qseqid %in% og_mm$qseqid | mult$sseqid %in% og_mm$sseqid
    og_mm <- mult[og_mm_i, ]

    # Identify 1-to-many (1-to-M) ortholog pairs
    not_mm <- mult[!og_mm_i, ]
    not_mm_q_dup <- not_mm$qseqid %in% not_mm$qseqid[duplicated(not_mm$qseqid)]
    not_mm_s_dup <- not_mm$sseqid %in% not_mm$sseqid[duplicated(not_mm$sseqid)]
    og_1m <- not_mm[not_mm_q_dup, ]
    og_m1 <- not_mm[not_mm_s_dup, ]

    # Assign unique identifiers to 1-to-many ortholog pairs
    og_1m <- og_1m[order(og_1m$qseqid), ]
    og_1m_index <- data.frame(id = unique(og_1m$qseqid), OG = seq_along(unique(og_1m$qseqid)))
    og_1m$OG <- og_1m_index$OG[match(og_1m$qseqid, og_1m_index$id)]

    # Assign unique identifiers to many-to-one ortholog pairs
    og_m1 <- og_m1[order(og_m1$qseqid), ]
    og_m1_index <- data.frame(id = unique(og_m1$sseqid), OG = seq_along(unique(og_m1$sseqid)))
    og_m1$OG <- og_m1_index$OG[match(og_m1$sseqid, og_m1_index$id)]

    # Assign unique identifiers to many-to-many ortholog pairs
    og_mm <- og_mm[order(og_mm$qseqid), ]
    og_mm$OG <- seq_along(og_mm$qseqid)
    og_mm_q_min <- tapply(og_mm$OG, og_mm$qseqid, min)
    og_mm$q_min <- og_mm_q_min[match(og_mm$qseqid, names(og_mm_q_min))]
    og_mm_s_min <- tapply(og_mm$q_min, og_mm$sseqid, min)
    og_mm$s_min <- og_mm_s_min[match(og_mm$sseqid, names(og_mm$s_min))]
    og_mm_q_min <- tapply(og_mm$s_min, og_mm$qseqid, min)
    og_mm$OG <- og_mm_q_min[match(og_mm$qseqid, names(og_mm_q_min))]
    og_mm <- subset(og_mm, select = -c(q_min, s_min))

    # Combine all ortholog pairs and assign classes
    out <- og_11
    out$class <- "1to1"
    og_1m$OG <- og_1m$OG + max(out$OG)
    og_1m$class <- "1toM"
    out <- rbind(out, og_1m)
    og_m1$OG <- og_m1$OG + max(out$OG)
    og_m1$class <- "Mto1"
    out <- rbind(out, og_m1)
    og_mm$OG <- og_mm$OG + max(out$OG)
    og_mm$class <- "MtoM"
    out <- rbind(out, og_mm)

    # Match ortholog pairs with assigned ortholog group identifiers
    hit <- match(ortho$pair_id, out$pair_id)
    return(out$OG[hit])
}

#' Identify Gene Orthologs
#'
#' This function identifies gene orthologs based on syntenic ortholog pairs and assigns ortholog groups.
#'
#' @param object A SynogDB object.
#'
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
geneOrtho <- function(object){
    # Check that object is of class SynogDB
    stopifnot(inherits(x = object, "SynogDB"))

    # Open HDF5 file
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))

    # Check for the existence of syntenic ortholog pairs
    if(!H5Lexists(h5, "synog_tx/orthopairs")){
        stop("Run syntenyOrtho to obtain ortholog anchors info.")
    }

    # Extract ortholog pairs
    ortho <- h5$synog_tx$orthopairs

    # Find the best gene pairs based on ortholog relationships
    ortho <- .findBestPair(object = object, ortho = ortho)

    # Create pair IDs and classify orthologs
    ortho$pair_id <- paste(ortho$qgeneid, ortho$sgeneid, sep = "_")
    ortho$genewise_class <- .classifySOgene(ortho = ortho)

    # Identify and omit 1-to-1 gene pairs
    omit_id <- .find1to1gene(target = ortho)
    ortho <- subset(ortho, subset = !ortho$pair_id %in% omit_id)

    # Number ortholog groups based on gene-wise classification
    class_og <- .numberingOGgene(ortho = ortho)
    ortho$genewise_class <- class_og$class
    ortho$genewise_OG <- class_og$OG

    # Order ortholog pairs by ortholog group and gene IDs
    ortho <- ortho[order(ortho$genewise_OG, ortho$qgeneid, ortho$sgeneid), ]
    rownames(ortho) <- NULL

    # Identify orphan genes
    orphan <- .getOrphanGene(ortho = ortho, tx_orphan = h5$synog_tx$orphan)

    # Create gene summary
    gene_summary <- .makeSummary(ortho = ortho, orphan = orphan)

    # Remove pair IDs
    ortho$pair_id <- NULL

    # Save results to HDF5 file
    .h5creategroup(object$h5, "synog_gene")
    .h5overwrite(obj = ortho, file = object$h5, "synog_gene/orthopairs")
    .h5overwrite(obj = gene_summary, file = object$h5, "synog_gene/summary")
    .h5overwrite(obj = orphan, file = object$h5, "synog_gene/orphan")
}

#' Classify Syntenic Ortholog Gene Pairs
#'
#' This function classifies syntenic ortholog gene pairs into categories: 1-to-1, 1-to-M, M-to-1, and M-to-M.
#'
.classifySOgene <- function(ortho){
    # Select relevant columns and ensure uniqueness
    tmp <- subset(ortho, select = c(qgeneid, sgeneid, pair_id))
    tmp <- unique(tmp)

    # Identify duplicated query and subject gene IDs
    q_dup <- duplicated(tmp$qgeneid)
    q_dup <- tmp$qgeneid %in% tmp$qgeneid[q_dup]
    s_dup <- duplicated(tmp$sgeneid)
    s_dup <- tmp$sgeneid %in% tmp$sgeneid[s_dup]

    # Identify 1-to-1 orthologs (no duplications in either gene ID)
    og_11 <- tmp[!(q_dup | s_dup), ]

    # Handle multiple gene duplications
    mult <- tmp[q_dup | s_dup, ]
    mult_q_dup <- mult$qgeneid %in% mult$qgeneid[duplicated(mult$qgeneid)]
    mult_s_dup <- mult$sgeneid %in% mult$sgeneid[duplicated(mult$sgeneid)]

    # Identify M-to-M orthologs (multiple duplications in both gene IDs)
    og_mm <- mult[mult_q_dup & mult_s_dup, ]
    og_mm_i <- mult$qgeneid %in% og_mm$qgeneid | mult$sgeneid %in% og_mm$sgeneid
    og_mm <- mult[og_mm_i, ]

    # Identify other duplications not in M-to-M
    not_mm <- mult[!og_mm_i, ]
    not_mm_q_dup <- not_mm$qgeneid %in% not_mm$qgeneid[duplicated(not_mm$qgeneid)]
    not_mm_s_dup <- not_mm$sgeneid %in% not_mm$sgeneid[duplicated(not_mm$sgeneid)]
    og_1m <- not_mm[not_mm_q_dup, ]
    og_m1 <- not_mm[not_mm_s_dup, ]

    # Combine classifications
    out <- og_11
    out$class <- "1to1"
    og_1m$class <- "1toM"
    out <- rbind(out, og_1m)
    og_m1$class <- "Mto1"
    out <- rbind(out, og_m1)
    og_mm$class <- "MtoM"
    out <- rbind(out, og_mm)

    # Match classifications to the original ortholog pairs
    hit <- match(ortho$pair_id, out$pair_id)
    return(out$class[hit])
}

#' Find 1-to-1 Gene Pairs from M-to-M Classes
#'
#' This function identifies the best 1-to-1 gene pairs from the M-to-M ortholog class by selecting pairs with the lowest mutual_e value and the highest mutual_ci value.
#'
.find1to1gene <- function(target){
    # Subset the target data frame to include only M-to-M gene pairs
    og_target <- target[target$genewise_class == "MtoM", ]
    score_index <- seq_along(og_target$qgeneid)

    # Identify the best 1-to-1 gene pairs for each query gene ID
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

    # Identify the best 1-to-1 gene pairs for each subject gene ID
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

    # Combine and deduplicate 1-to-1 gene pairs
    add_1to1 <- unique(rbind(q2s, s2q))
    add_1to1$pair_id <- paste(add_1to1$qgeneid, add_1to1$sgeneid, sep = "_")

    # Find indices of added 1-to-1 pairs and their query and subject gene hits
    add_1to1_index <- which(og_target$pair_id %in% add_1to1$pair_id)
    q_hit <- og_target$qgeneid %in% og_target$qgeneid[add_1to1_index]
    s_hit <- og_target$sgeneid %in% og_target$sgeneid[add_1to1_index]

    # Identify pairs to omit from the 1-to-1 classification
    out <- og_target$pair_id[q_hit | s_hit]
    out <- out[!out %in% add_1to1$pair_id]

    return(out)
}

#' Number Orthologous Groups (OGs) at the Gene Level
#'
#' This function assigns a unique Orthologous Group (OG) number to each gene pair based on their ortholog classification (1-to-1, 1-to-M, M-to-1, or M-to-M).
#'
.numberingOGgene <- function(ortho){
    # Subset and remove duplicate pairs
    tmp <- subset(ortho, select = c(qgeneid, sgeneid, pair_id))
    tmp <- unique(tmp)

    # Identify duplicated query and subject gene IDs
    q_dup <- duplicated(tmp$qgeneid)
    q_dup <- tmp$qgeneid %in% tmp$qgeneid[q_dup]
    s_dup <- duplicated(tmp$sgeneid)
    s_dup <- tmp$sgeneid %in% tmp$sgeneid[s_dup]

    # Create OGs for 1-to-1 gene pairs
    og_11 <- tmp[!(q_dup | s_dup), ]
    og_11$OG <- seq_along(og_11$qgeneid)

    # Subset pairs with duplicated query or subject gene IDs
    mult <- tmp[q_dup | s_dup, ]
    mult_q_dup <- mult$qgeneid %in% mult$qgeneid[duplicated(mult$qgeneid)]
    mult_s_dup <- mult$sgeneid %in% mult$sgeneid[duplicated(mult$sgeneid)]

    # Create OGs for M-to-M gene pairs
    og_mm <- mult[mult_q_dup & mult_s_dup, ]
    og_mm_i <- mult$qgeneid %in% og_mm$qgeneid | mult$sgeneid %in% og_mm$sgeneid
    og_mm <- mult[og_mm_i, ]

    # Create OGs for 1-to-M and M-to-1 gene pairs
    not_mm <- mult[!og_mm_i, ]
    not_mm_q_dup <- not_mm$qgeneid %in% not_mm$qgeneid[duplicated(not_mm$qgeneid)]
    not_mm_s_dup <- not_mm$sgeneid %in% not_mm$sgeneid[duplicated(not_mm$sgeneid)]
    og_1m <- not_mm[not_mm_q_dup, ]
    og_m1 <- not_mm[not_mm_s_dup, ]

    # Assign OG numbers to 1-to-M gene pairs
    og_1m <- og_1m[order(og_1m$qgeneid), ]
    og_1m_index <- data.frame(id = unique(og_1m$qgeneid),
                              OG = seq_along(unique(og_1m$qgeneid)))
    og_1m$OG <- og_1m_index$OG[match(og_1m$qgeneid, og_1m_index$id)]

    # Assign OG numbers to M-to-1 gene pairs
    og_m1 <- og_m1[order(og_m1$qgeneid), ]
    og_m1_index <- data.frame(id = unique(og_m1$sgeneid),
                              OG = seq_along(unique(og_m1$sgeneid)))
    og_m1$OG <- og_m1_index$OG[match(og_m1$sgeneid, og_m1_index$id)]

    # Assign OG numbers to M-to-M gene pairs and resolve conflicts
    og_mm <- og_mm[order(og_mm$qgeneid), ]
    og_mm$OG <- seq_along(og_mm$qgeneid)
    og_mm_q_min <- tapply(og_mm$OG, og_mm$qgeneid, min)
    og_mm$q_min <- og_mm_q_min[match(og_mm$qgeneid, names(og_mm_q_min))]
    og_mm_s_min <- tapply(og_mm$q_min, og_mm$sgeneid, min)
    og_mm$s_min <- og_mm_s_min[match(og_mm$sgeneid, names(og_mm$s_min))]
    og_mm_q_min <- tapply(og_mm$s_min, og_mm$qgeneid, min)
    og_mm$OG <- og_mm_q_min[match(og_mm$qgeneid, names(og_mm$q_min))]
    og_mm <- subset(og_mm, select = -c(q_min, s_min))

    # Combine all OGs and classes into the final output
    out <- og_11
    out$class <- "1to1"
    og_1m$class <- "1toM"
    og_1m$OG <- og_1m$OG + max(out$OG)
    out <- rbind(out, og_1m)
    og_m1$class <- "Mto1"
    og_m1$OG <- og_m1$OG + max(out$OG)
    out <- rbind(out, og_m1)
    og_mm$class <- "MtoM"
    og_mm$OG <- og_mm$OG + max(out$OG)
    out <- rbind(out, og_mm)

    # Match the final OGs and classes with the original ortholog pairs
    hit <- match(ortho$pair_id, out$pair_id)

    return(list(class = out$class[hit], OG = out$OG[hit]))
}

#' Identify Orphan Genes
#'
#' This function identifies genes that do not have orthologous pairs.
#'
.getOrphanGene <- function(ortho, tx_orphan){
    # Extract unique gene IDs for query and subject from tx_orphan
    q_genes <- unique(tx_orphan$query$gene_id)
    s_genes <- unique(tx_orphan$subject$gene_id)

    # Identify query genes that are not present in ortho$qgeneid
    q_orphan <- q_genes[!q_genes %in% ortho$qgeneid]
    q_orphan <- data.frame(qgeneid = q_orphan)

    # Identify subject genes that are not present in ortho$sgeneid
    s_orphan <- s_genes[!s_genes %in% ortho$sgeneid]
    s_orphan <- data.frame(sgeneid = s_orphan)

    # Return a list containing data frames for query and subject orphan genes
    return(list(query = q_orphan, subject = s_orphan))
}

#' Generate Summary of Orthologous and Orphan Genes
#'
#' This function creates a summary data frame that provides counts of orthologous gene pairs classified by type and counts of orphan genes.
#'
.makeSummary <- function(ortho, orphan){
    # Group query gene IDs by their genewise_class and get unique IDs
    q_id <- tapply(ortho$qgeneid, ortho$genewise_class, unique)
    n_q_id <- sapply(q_id, length)

    # Group subject gene IDs by their genewise_class and get unique IDs
    s_id <- tapply(ortho$sgeneid, ortho$genewise_class, unique)
    n_s_id <- sapply(s_id, length)

    # Get counts of orphan genes for query and subject
    n_orphan <- lapply(orphan, nrow)

    # Calculate total counts of query and subject genes, including orphan genes
    q_tot <- sum(n_q_id) + n_orphan$query
    s_tot <- sum(n_s_id) + n_orphan$subject

    # Create summary data frame
    out <- data.frame(
        query = c(q_tot, sum(n_q_id), n_q_id, n_orphan$query),
        subject = c(s_tot, sum(n_s_id), n_s_id, n_orphan$subject)
    )

    return(out)
}

#' Find the Best Pair of Orthologous Genes
#'
#' This function identifies the best pair of orthologous genes based on several criteria, including RBBH status, mutual e-value, mutual CI, and sequence length.
#'
#' @importFrom BiocGenerics width
.findBestPair <- function(object, ortho){
    # Generate a unique identifier for each pair
    pair_id <- paste(ortho$qgeneid, ortho$sgeneid, sep = "_")

    # Import and filter GFF files for query and subject
    q_gff <- .importAllGFF(object$query_gff)
    s_gff <- .importAllGFF(object$subject_gff)
    q_gff <- q_gff[q_gff$type %in% c("mRNA", "transcript")]
    s_gff <- s_gff[s_gff$type %in% c("mRNA", "transcript")]

    # Calculate lengths of sequences
    q_len <- width(q_gff[match(ortho$qseqid, q_gff$ID)])
    s_len <- width(s_gff[match(ortho$sseqid, s_gff$ID)])

    # Function to select the best pair based on criteria
    best_pair <- tapply(seq_along(pair_id), pair_id, function(index){
        # Check if RBBH status is TRUE
        check_rbbh <- ortho$rbbh[index] == "TRUE"
        if(sum(check_rbbh) == 1){
            return(index[check_rbbh])
        } else if(sum(check_rbbh) > 1){
            index <- index[check_rbbh]
        }

        # Check for minimum mutual e-value
        check_e <- ortho$mutual_e[index] == min(ortho$mutual_e[index])
        if(sum(check_e) == 1){
            return(index[check_e])
        } else if(sum(check_e) > 1){
            index <- index[check_e]
        }

        # Check for maximum mutual CI
        check_ci <- ortho$mutual_ci[index] == max(ortho$mutual_ci[index])
        if(sum(check_ci) == 1){
            return(index[check_ci])
        } else if(sum(check_ci) > 1){
            index <- index[check_ci]
        }

        # Check for maximum length if there are ties
        if(length(index) == 1){
            return(index)
        }

        len <- q_len[index] + s_len[index]
        check_len <- len == max(len)
        if(sum(check_len) == 1){
            return(index[check_len])
        } else {
            return(sample(index, 1)) # Randomly select one if there are still ties
        }
    })

    return(ortho[best_pair, ])
}


################################################################################
#' Split Genes Based on Orthologous Relationships
#'
#' This function processes genewise ortholog pairs to identify and handle split genes. It evaluates the genewise classes, replaces gene IDs, and generates orthologous gene groups.
#'
#' @param object A SynogDB object containing necessary HDF5 and GFF files.
#'
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
splitGenes <- function(object){
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))

    # Check if genewise ortholog info is available
    if(!H5Lexists(h5, "synog_gene/orthopairs")){
        stop("Run geneOrtho to obtain genewise ortholog info.")
    }

    # Load ortholog pairs
    ortho <- h5$synog_gene$orthopairs
    ortho$pair_id1 <- paste(ortho$qseqid, ortho$sgeneid, sep = "_")
    ortho$pair_id2 <- paste(ortho$qgeneid, ortho$sseqid, sep = "_")
    ortho$pair_id3 <- paste(ortho$qseqid, ortho$sseqid, sep = "_")

    # Prepare GFF data
    gff <- .prepGFF(object)

    # Evaluate genewise classes
    ortho$genewise_class <- .eval1toM(ortho = ortho, gff = gff)
    ortho$genewise_class <- .evalMto1(ortho = ortho, gff = gff)
    ortho$genewise_class <- .evalMtoM(ortho = ortho, gff = gff)

    # Replace gene IDs
    new_id <- .replaceIDs(ortho = ortho)
    ortho$qgeneid <- new_id$qgeneid
    ortho$sgeneid <- new_id$sgeneid
    ortho$pair_id <- paste(ortho$qgeneid, ortho$sgeneid, sep = "_")

    # Generate orthologous gene groups
    class_og <- .numberingOGgene(ortho = ortho)
    ortho$genewise_class <- class_og$class
    ortho$genewise_OG <- class_og$OG
    ortho <- ortho[order(ortho$genewise_OG, ortho$qgeneid, ortho$sgeneid), ]
    rownames(ortho) <- NULL

    # Load orphan data and generate summary
    orphan <- h5$synog_gene$orphan
    gene_summary <- .makeSummary(ortho = ortho, orphan = orphan)
    ortho$pair_id <- ortho$pair_id1 <- ortho$pair_id2 <- ortho$pair_id3 <- NULL

    # Save the processed data into new HDF5 groups
    .h5creategroup(object$h5,"synog_gene_split")
    .h5overwrite(obj = ortho, file = object$h5, "synog_gene_split/orthopairs")
    .h5overwrite(obj = gene_summary, file = object$h5, "synog_gene_split/summary")
    .h5overwrite(obj = orphan, file = object$h5, "synog_gene_split/orphan")
}

#' Prepare GFF Data for Gene Processing
#'
#' This function imports and filters GFF data to extract CDS features and ensures that the Parent attributes are properly unlisted.
#'
.prepGFF <- function(object){
    # Import GFF data for query and subject genomes
    query_gff <- .importAllGFF(object$query_gff)
    subject_gff <- .importAllGFF(object$subject_gff)

    # Filter GFF data to include only CDS features
    query_gff <- query_gff[query_gff$type == "CDS"]
    query_gff$Parent <- unlist(query_gff$Parent)
    subject_gff <- subject_gff[subject_gff$type == "CDS"]
    subject_gff$Parent <- unlist(subject_gff$Parent)

    # Return processed GFF data as a list
    return(list(query_gff = query_gff, subject_gff = subject_gff))
}
#' Evaluate 1-to-M Orthologs for Splits
#'
#' This function evaluates 1-to-M orthologs to identify cases where a 1-to-M relationship might actually be a result of gene splits.
#'
.eval1toM <- function(ortho, gff = gff){
    # Initialize output with current classifications
    out <- ortho$genewise_class

    # Subset the data to include only 1-to-M orthologs
    tmp <- subset(ortho, select = c(qgeneid, qseqid, sgeneid, pair_id1),
                  subset = genewise_class == "1toM")

    # Find 1-to-1 splits in the query genome
    og_11 <- .find1to1split(seqid = tmp$qseqid, geneid = tmp$sgeneid,
                            pair_id = tmp$pair_id1, gff = gff$query_gff)
    # Update the classification to "split_1toM" for identified splits
    out[ortho$pair_id1 %in% og_11] <- "split_1toM"

    # Process the remaining 1-to-M orthologs
    rest <- tmp[!tmp$pair_id1 %in% og_11, ]
    og_11 <- .pick1to1rest(rest = rest, pair_id = rest$pair_id1)
    out[ortho$pair_id1 %in% og_11] <- "split_1toM"

    return(out)
}
#' Evaluate M-to-1 Orthologs for Splits
#'
#' This function evaluates M-to-1 orthologs to identify cases where a M-to-1 relationship might actually be a result of gene splits.
#'
.evalMto1 <- function(ortho, gff = gff){
    # Initialize output with current classifications
    out <- ortho$genewise_class

    # Subset the data to include only M-to-1 orthologs
    tmp <- subset(ortho, select = c(qgeneid, sseqid, sgeneid, pair_id2),
                  subset = genewise_class == "Mto1")

    # Find 1-to-1 splits in the subject genome
    og_11 <- .find1to1split(seqid = tmp$sseqid, geneid = tmp$qgeneid,
                            pair_id = tmp$pair_id2, gff = gff$subject_gff)
    # Update the classification to "split_Mto1" for identified splits
    out[ortho$pair_id2 %in% og_11] <- "split_Mto1"

    # Process the remaining M-to-1 orthologs
    rest <- tmp[!tmp$pair_id2 %in% og_11, ]
    og_11 <- .pick1to1rest(rest = rest, pair_id = rest$pair_id2)
    out[ortho$pair_id2 %in% og_11] <- "split_Mto1"

    return(out)
}

#' Evaluate M-to-M Orthologs for Splits
#'
#' This function evaluates M-to-M orthologs to identify cases where an M-to-M relationship might actually be a result of gene splits.
#'
.evalMtoM <- function(ortho, gff = gff){
    # Initialize output with current classifications
    out <- ortho$genewise_class

    # Subset the data to include only M-to-M orthologs
    tmp <- subset(ortho, select = c(qseqid, qgeneid, sseqid, sgeneid,
                                    pair_id1, pair_id2, pair_id3),
                  subset = genewise_class == "MtoM")

    # Find 1-to-1 splits in the query and subject genomes
    og_11_1 <- .find1to1split(seqid = tmp$qseqid, geneid = tmp$sgeneid,
                              pair_id = tmp$pair_id3, gff = gff$query_gff)
    og_11_2 <- .find1to1split(seqid = tmp$sseqid, geneid = tmp$qgeneid,
                              pair_id = tmp$pair_id3, gff = gff$subject_gff)

    # Identify common 1-to-1 splits between query and subject genomes
    og_11 <- og_11_1[og_11_1 %in% og_11_2]
    out[ortho$pair_id3 %in% og_11] <- "split_MtoM"

    # Process the remaining M-to-M orthologs
    rest <- tmp[!tmp$pair_id3 %in% og_11, ]
    og_11 <- .pick1to1rest(rest = rest, pair_id = rest$pair_id3)
    out[ortho$pair_id3 %in% og_11] <- "split_MtoM"

    return(out)
}

#' Find 1-to-1 Gene Splits
#'
#' This function identifies cases where an M-to-M relationship might actually be a result of gene splits.
#'
#' @importFrom GenomicRanges findOverlaps
.find1to1split <- function(seqid, geneid, pair_id, gff){
    # Initialize output with input pair IDs
    out <- pair_id

    # Identify duplicated sequences and genes
    tx_dup <- duplicated(seqid)
    tx_dup <- seqid %in% seqid[tx_dup]
    gene_dup <- duplicated(geneid)
    gene_dup <- geneid %in% geneid[gene_dup]

    # Identify 1-to-1 sequence IDs and update output accordingly
    seqid_1to1 <- seqid[!(tx_dup | gene_dup)]
    out <- out[!(tx_dup | gene_dup)]

    # Filter GFF data for 1-to-1 sequence IDs
    gff_1to1 <- gff[gff$Parent %in% seqid_1to1]

    # Find overlapping genomic features within the same sequences
    ol <- findOverlaps(gff_1to1, gff_1to1)
    ol <- as.data.frame(ol)
    ol <- subset(ol, subset = queryHits != subjectHits)

    # Identify invalid sequence IDs based on overlaps
    hit_id <- unique(gff_1to1$Parent[ol$queryHits])
    seqid_1to1_candidate <- seqid_1to1[!seqid_1to1 %in% hit_id]
    seqid_1to1_rest <- seqid[!seqid %in% seqid_1to1_candidate]
    out <- out[seqid_1to1 %in% seqid_1to1_candidate]

    # Further filter GFF data for candidate and remaining sequences
    gff_1to1_candidate <- gff[gff$Parent %in% seqid_1to1_candidate]
    gff_1to1_rest <- gff[gff$Parent %in% seqid_1to1_rest]

    # Find overlaps between candidate and remaining sequences
    ol <- findOverlaps(gff_1to1_candidate, gff_1to1_rest)
    invalid_seqid <- gff_1to1_candidate$Parent[queryHits(ol)]

    # Update output by removing invalid sequence IDs
    out <- out[!seqid_1to1_candidate %in% invalid_seqid]
    return(out)
}

#' Pick 1-to-1 Gene Relationships from the Remaining Pairs
#'
#' This function identifies 1-to-1 gene relationships from the remaining pairs by removing pairs where either gene is duplicated.
#'
.pick1to1rest <- function(rest, pair_id){
    out <- pair_id

    # Identify duplicated query genes and subject genes
    q_dup <- duplicated(rest$qgeneid)
    q_dup <- rest$qgeneid %in% rest$qgeneid[q_dup]
    s_dup <- duplicated(rest$sgeneid)
    s_dup <- rest$sgeneid %in% rest$sgeneid[s_dup]

    # Filter out pairs where either query or subject gene is duplicated
    out <- out[!(q_dup | s_dup)]
    return(out)
}
#' Replace Gene IDs with Transcript IDs for Split Gene Relationships
#'
#' This function updates the gene IDs in the ortholog pairs with their respective transcript IDs for split gene relationships.
#'
.replaceIDs <- function(ortho){
    out <- list(qgeneid = ortho$qgeneid, sgeneid = ortho$sgeneid)

    # Replace query gene IDs with query transcript IDs for split_1toM relationships
    hit <- ortho$genewise_class == "split_1toM"
    out$qgeneid[hit] <- ortho$qseqid[hit]

    # Replace subject gene IDs with subject transcript IDs for split_Mto1 relationships
    hit <- ortho$genewise_class == "split_Mto1"
    out$sgeneid[hit] <- ortho$sseqid[hit]

    # Replace both query and subject gene IDs with their respective transcript IDs for split_MtoM relationships
    hit <- ortho$genewise_class == "split_MtoM"
    out$qgeneid[hit] <- ortho$qseqid[hit]
    out$sgeneid[hit] <- ortho$sseqid[hit]

    return(out)
}
