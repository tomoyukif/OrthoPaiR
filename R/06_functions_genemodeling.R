#' @export
#'
mapProt <- function(object, out_prefix, miniprot_bin, n_core, len_diff = 0.5){
    stopifnot(inherits(x = object, "SynogDB"))
    .mapEngine(object = object, subject_prot = object$subject_prot,
               query_prot = object$query_prot, out_prefix = out_prefix,
               miniprot_bin = miniprot_bin, n_core = n_core,
               overlap = TRUE, len_diff = len_diff)
}

#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
mapOrphan <- function(object, out_prefix, miniprot_bin, n_core, len_diff = 0.5){
    stopifnot(inherits(x = object, "SynogDB"))
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    if(!H5Lexists(h5, "synog_gene/orthopairs")){
        stop("Run geneOrtho to obtain genewise ortholog info.")
    }

    q_aa_fn <- paste0(out_prefix, "miniprot_query_prot.fa")
    s_aa_fn <- paste0(out_prefix, "miniprot_subject_prot.fa")
    .writeProtFASTA(id = h5$synog_tx$orphan$query$qseqid,
                    prot_fn = object$query_prot,
                    gff_fn = object$query_gff,
                    fn = q_aa_fn)

    .writeProtFASTA(id = h5$synog_tx$orphan$subject$sseqid,
                    prot_fn = object$subject_prot,
                    gff_fn = object$subject_gff,
                    fn = s_aa_fn)

    .mapEngine(object = object, subject_prot = s_aa_fn,
               query_prot = q_aa_fn, out_prefix = out_prefix,
               miniprot_bin = miniprot_bin, n_core = n_core,
               overlap = TRUE, len_diff = len_diff)
}

#' @importFrom rtracklayer export.gff3
.mapEngine <- function(object, subject_prot, query_prot,
                       out_prefix, miniprot_bin, n_core,
                       overlap, len_diff){
    s2q_out <- paste0(out_prefix, "query_miniprot_out")
    q2s_out <- paste0(out_prefix, "subject_miniprot_out")
    .miniprot(query_fn = subject_prot, genome_fn = object$query_genome,
              out_prefix = s2q_out, miniprot_bin = miniprot_bin,
              n_core = n_core)
    .miniprot(query_fn = query_prot, genome_fn = object$subject_genome,
              out_prefix = q2s_out, miniprot_bin = miniprot_bin,
              n_core = n_core)

    .h5creategroup(object$h5,"protmap")
    .h5overwrite(obj = paste0(q2s_out, ".gff"),
                 file = object$h5, "protmap/q2s_gff")
    .h5overwrite(obj = paste0(s2q_out, ".gff"),
                 file = object$h5, "protmap/s2q_gff")

    query_gff <- .importAllGFF(object$query_gff)
    subject_gff <- .importAllGFF(object$subject_gff)
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
    s2q_gff <- .mergeGFF(gff1 = s2q_gff, gff2 = query_gff)
    q2s_gff <- .mergeGFF(gff1 = q2s_gff, gff2 = subject_gff)

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
    export.gff3(s2q_gff, paste0(s2q_out, ".gff"))
    export.gff3(q2s_gff, paste0(q2s_out, ".gff"))
}

#' @importFrom Biostrings readAAStringSet writeXStringSet
.writeProtFASTA <- function(id, prot_fn, gff_fn, fn){
    aa <- readAAStringSet(prot_fn)
    gff <- .importAllGFF(gff_fn)
    tx_i <- gff$type %in% c("mRNA", "transcript")
    tx_hit <- gff$ID[tx_i] %in% id
    hit_tx <- gff$ID[tx_i][tx_hit]
    hit_aa <- aa[names(aa) %in% hit_tx]
    writeXStringSet(hit_aa, fn)
}

.miniprot <- function(query_fn, genome_fn, out_prefix,
                      miniprot_bin, n_core = 1){
    command <- miniprot_bin

    args <- paste(paste("-t", n_core,
                        "-d", paste0(out_prefix, ".mpi"),
                        genome_fn),
                  paste(miniprot_bin, "-t", n_core,
                        "--gff",
                        paste0(out_prefix, ".mpi"),
                        query_fn, ">", paste0(out_prefix, ".gff")),
                  sep = ";")
    system2(command = command, args = args)
}

#' @importFrom BioGenerics start
.orgGFF <- function(gff1, gff2, gff3, prefix, overlap = FALSE, len_diff){
    lv <- levels(gff1$type)
    lv[lv == "mRNA"] <- "transcript"
    levels(gff1$type) <- lv
    gff1 <- gff1[gff1$type %in% c("gene", "transcript",
                                  "five_prime_UTR", "CDS", "three_prime_UTR")]
    gff1 <- gff1[order(as.numeric(seqnames(gff1)), start(gff1))]
    gff2 <- gff2[order(as.numeric(seqnames(gff2)), start(gff2))]

    gff_valid <- .validTX(gff1 = gff1, gff2 = gff2, gff3 = gff3,
                          overlap = overlap, len_diff = len_diff)

    gff_tx <- .orgTX(gff = gff_valid, prefix = prefix)
    gff_cds <- .orgCDS(gff1 = gff1, gff_tx = gff_tx)
    gff_gene <- .orgGene(gff_tx = gff_tx)
    gff_tx$old_id <- gff_gene$old_id <- NULL
    out <- c(gff_gene, gff_tx, gff_cds)
    out <- out[order(as.numeric(seqnames(out)),
                     start(out),
                     as.numeric(out$type))]
    return(out)
}

#' @importFrom S4Vectors queryHits
#' @importFrom GenomicRanges findOverlaps
.validTX <- function(gff1, gff2, gff3, overlap, len_diff){
    if(overlap){
        gff1_block <- .getCDSblock(gff = gff1)
        gff1_block_uniq <- gff1_block[!duplicated(gff1_block)]

        gff2_block <- .getCDSblock(gff = gff2)
        gff1_block_uniq <- gff1_block_uniq[!gff1_block_uniq %in% gff2_block]
        out <- .getUniqTx(gff1 = gff1, gff1_block_uniq = gff1_block_uniq)

    } else {
        tx_i1 <- gff1$type == "transcript"
        tx_i2 <- gff2$type == "transcript"
        ol <- findOverlaps(gff1[tx_i1], gff2[tx_i2])
        hit <- unique(queryHits(ol))
        out <- gff1[tx_i1][!seq_len(sum(tx_i1)) %in% hit]
    }
    out <- .checkLength(gff = out, gff3 = gff3, len_diff = len_diff)
    return(out)
}

.getCDSblock <- function(gff){
    gff_cds_i <- gff$type == "CDS"
    gff_tx_i <-  gff$type == "transcript"
    gff_cds_start <- start(gff[gff_cds_i])
    gff_cds_end <- end(gff[gff_cds_i])
    gff_cds_exon <- paste(gff_cds_start, gff_cds_end, sep = "-")
    map_to_tx <- match(unlist(gff$Parent[gff_cds_i]), gff$ID[gff_tx_i])
    first_i <- !duplicated(map_to_tx)
    gff_cds_chr <- as.character(seqnames(gff[gff_cds_i]))
    gff_cds_exon[first_i] <- paste(gff_cds_chr[first_i],
                                   gff_cds_exon[first_i],
                                   sep = ":")
    out <- tapply(gff_cds_exon, map_to_tx, paste, collapse = ",")
    return(out)
}

.getUniqTx <- function(gff1, gff1_block_uniq){
    gff1_tx_i <-  gff1$type == "transcript"
    out_tx <- gff1[gff1_tx_i][as.numeric(names(gff1_block_uniq))]
    out_element <- unlist(gff1$Parent[!gff1_tx_i]) %in% out_tx$ID
    out_element <- gff1[!gff1_tx_i][out_element]
    out <- c(out_tx, out_element)
    return(out)
}

.checkLength <- function(gff, gff3, len_diff){
    gff_cds_len <- .getCDSlen(gff = gff)
    index <- as.numeric(names(gff_cds_len))
    gff_tx_i <- gff$type %in% "transcript"
    names(gff_cds_len) <- sub("\\s.+", "", gff$Target[gff_tx_i][index])

    gff3_cds_len <- .getCDSlen(gff = gff3)
    index <- as.numeric(names(gff3_cds_len))
    gff3_tx_i <- gff3$type %in% c("transcript", "mRNA")
    names(gff3_cds_len) <- gff3$ID[gff3_tx_i][index]

    id_hit <- match(names(gff_cds_len), names(gff3_cds_len))
    gff3_cds_len <- gff3_cds_len[id_hit]
    longer <- gff_cds_len
    is_gff3_longer <- longer < gff3_cds_len
    longer[is_gff3_longer] <- gff3_cds_len[is_gff3_longer]
    valid <- abs(gff_cds_len - gff3_cds_len) / longer <= len_diff
    valid_tx <- gff[gff_tx_i][as.vector(valid)]
    non_tx_gff <- gff[!gff_tx_i]
    out <- c(valid_tx,
             non_tx_gff[unlist(non_tx_gff$Parent) %in% valid_tx$ID])
    return(out)
}

#' @importFrom BiocGenerics start end
.getCDSlen <- function(gff){
    gff_cds_i <- gff$type == "CDS"
    gff_tx_i <-  gff$type == "transcript"
    gff_cds_start <- start(gff[gff_cds_i])
    gff_cds_end <- end(gff[gff_cds_i])
    gff_cds_len <- gff_cds_end - gff_cds_start
    map_to_tx <- match(unlist(gff$Parent[gff_cds_i]), gff$ID[gff_tx_i])
    out <- tapply(gff_cds_len, map_to_tx, sum)
    return(out)
}

#' @importFrom BiocGenerics start
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges reduce findOverlaps
#' @importFrom S4Vectors subjectHits
.orgTX <- function(gff, prefix){
    tx_i <- gff$type == "transcript"
    gff_tx <- gff[tx_i]
    gff_tx <- gff_tx[order(as.numeric(seqnames(gff_tx)), start(gff_tx))]
    gff_tx <- unique(gff_tx)
    uni_loci <- reduce(gff_tx)
    map_to_loci <- findOverlaps(gff_tx, uni_loci)
    gff_tx$gene_id <- paste0(prefix, "G",
                             sprintf("%05d", subjectHits(map_to_loci)))
    gff_tx$old_id <- gff_tx$ID
    gff_tx <- gff_tx[order(gff_tx$gene_id)]
    gff_tx$ID <- unlist(tapply(gff_tx$gene_id, gff_tx$gene_id, function(x){
        return(paste0(sub("G", "T", x[1]), ".", sprintf("%02d", seq_len(length(x)))))
    }))
    gff_tx$Parent <- lapply(gff_tx$gene_id, c)
    return(gff_tx)
}

.orgCDS <- function(gff1, gff_tx){
    cds_i <- gff1$type != "transcript"
    gff_cds <- gff1[cds_i][unlist(gff1$Parent[cds_i]) %in% gff_tx$old_id]
    gff_cds$gene_id <- gff_tx$gene_id[match(unlist(gff_cds$Parent),
                                            gff_tx$old_id)]
    gff_cds$Parent <- lapply(gff_tx$ID[match(unlist(gff_cds$Parent),
                                             gff_tx$old_id)],
                             c)
    gff_cds$ID <- paste0(unlist(gff_cds$Parent), ":CDS")
    return(gff_cds)
}

.orgGene <- function(gff_tx){
    gff_gene <- gff_tx[!duplicated(gff_tx$gene_id)]
    gff_gene$type <- "gene"
    gff_gene$ID <- gff_gene$gene_id
    gff_gene$Parent <- lapply(seq_along(gff_gene), function(i) character())
    return(gff_gene)
}

#' @importFrom BiocGenerics start
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges reduce findOverlaps
#' @importFrom S4Vectors subjectHits
.mergeGFF <- function(gff1, gff2 = NULL, ann_priority){
    gff <- c(gff1, gff2)
    gff <- gff[gff$type %in% c("gene", "transcript",
                               "five_prime_UTR", "CDS", "three_prime_UTR")]
    gene_i <- gff$type == "gene"
    gff_gene <- gff[gene_i]
    gff_gene <- gff_gene[order(as.numeric(seqnames(gff_gene)), start(gff_gene))]
    uni_loci <- reduce(gff_gene)
    map_to_loci <- findOverlaps(gff_gene, uni_loci)

    ol_list <- tapply(queryHits(map_to_loci), subjectHits(map_to_loci), c)
    ol_list <- ol_list[sapply(ol_list, length) > 1]
    id_map <- .getIDmap(gff_gene = gff_gene, ol_list = ol_list)
    gff <- .mapID(gff = gff, id_map = id_map)

    gff <- gff[order(as.numeric(seqnames(gff)), start(gff), as.numeric(gff$type))]
    return(gff)
}

.getIDmap <- function(gff_gene, ol_list){
    out <- lapply(ol_list, function(i){
        gene_id <- gff_gene$ID[i]
        set_gene_id <- gene_id[!grepl("^query_|^subject_", gene_id)][1]
        id_map <- cbind(gene_id, set_gene_id)
        return(id_map)
    })
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

.fixGFF <- function(gff){
    gff <- .setGeneID(gff = gff)
    gff <- .fixGFFrange(gff = gff)
    gff <- .fixGFFexon(gff = gff)
    gff <- .fixGFFphase(gff = gff)
    return(gff)
}

.setGeneID <- function(gff){
    tx_i <- gff$type %in% c("transcript", "mRNA")
    hit <- match(unlist(gff$Parent[tx_i]), gff$ID)
    gff$gene_id[tx_i] <- gff$ID[hit]

    element_i <- !gff$type %in% c("gene", "transcript", "mRNA")
    hit <- match(unlist(gff$Parent[element_i]), gff$ID[tx_i])
    gff$gene_id[element_i] <- gff$gene_id[tx_i][hit]
    return(gff)
}

#' @importFrom BiocGenerics start end
.fixGFFrange <- function(gff){
    # Fix the start and end position of each transcript to
    # cover whole ranges of member elements (CDS, exon, and UTRs)
    tx_i <- which(gff$type %in% c("transcript", "mRNA"))
    element_i <- !gff$type %in% c("gene", "transcript", "mRNA")
    tx_start <- start(gff[tx_i])
    tx_end <- end(gff[tx_i])
    min_start <- tapply(start(gff[element_i]),
                        unlist(gff$Parent[element_i]),
                        min)
    max_end <- tapply(end(gff[element_i]),
                      unlist(gff$Parent[element_i]),
                      max)
    hit <- match(gff$ID[tx_i], names(min_start))
    tx_start <- min_start[hit]
    tx_end <- max_end[hit]
    table(tx_start < tx_end)
    not_na <- !is.na(tx_start)
    start(gff[tx_i[not_na]]) <- tx_start[not_na]
    not_na <- !is.na(tx_end)
    end(gff[tx_i[not_na]]) <- tx_end[not_na]

    # Fix the start and end position of each gene to
    # cover whole ranges of member transcripts
    gene_i <- gff$type %in% "gene"
    tx_i <- gff$type %in% c("transcript", "mRNA")
    gene_start <- start(gff[gene_i])
    gene_end <- end(gff[gene_i])
    min_start <- tapply(start(gff[tx_i]),
                        gff$gene_id[tx_i],
                        min)
    max_end <- tapply(end(gff[tx_i]),
                      gff$gene_id[tx_i],
                      max)
    hit <- match(gff$gene_id[gene_i], names(min_start))
    gene_start <- min_start[hit]
    gene_end <- max_end[hit]
    table(gene_start < gene_end)
    start(gff[gene_i]) <- gene_start
    end(gff[gene_i]) <- gene_end

    return(gff)
}

#' @importFrom BiocGenerics start
.fixGFFexon <- function(gff){
    gff <- gff[gff$type != "exon"]
    gff_exon <- gff[gff$type %in% c("CDS", "five_prime_UTR", "three_prime_UTR")]
    gff_exon$type <- "exon"
    gff_exon$Name <- gff_exon$ID <- paste0(unlist(gff_exon$Parent), ":exon")
    gff_exon$score <- gff_exon$phase <- NA

    tx_i <- gff$type %in% c("transcript", "mRNA")
    non_cds <- gff[tx_i][!gff$ID[tx_i] %in% unlist(gff_exon$Parent)]
    non_cds$type <- "exon"
    non_cds$Parent <- lapply(non_cds$ID, c)
    non_cds$Name <- non_cds$ID <- paste0(non_cds$ID, ":exon")
    non_cds$score <- non_cds$phase <- NA
    gff_exon <- c(gff_exon, non_cds)

    out <- c(gff, gff_exon)
    out$type <- factor(out$type, levels = c("gene", "transcript", "mRNA",
                                            "five_prime_UTR", "exon",
                                            "CDS", "three_prime_UTR"))
    out$type <- droplevels(out$type)
    out <- out[order(as.numeric(seqnames(out)), start(out), as.numeric(out$type))]
    return(out)
}

#' @importFrom BiocGenerics start
.fixGFFphase <- function(gff){
    gff_cds <- gff[gff$type == "CDS"]
    gff_cds_plus <- .phasePlus(gff_cds = gff_cds)
    gff_cds_minus <- .phaseMinus(gff_cds = gff_cds)
    out <- c(gff[gff$type != "CDS"], gff_cds_plus, gff_cds_minus)
    out$type <- factor(out$type, levels = c("gene", "transcript", "mRNA",
                                            "five_prime_UTR", "exon",
                                            "CDS", "three_prime_UTR"))
    out$type <- droplevels(out$type)
    out <- out[order(as.numeric(seqnames(out)), start(out), as.numeric(out$type))]
    return(out)
}

#' @importFrom BiocGenerics start width
#' @importFrom GenomeInfoDb seqnames
.phasePlus <- function(gff_cds){
    gff_cds_plus <- gff_cds[as.character(strand(gff_cds)) == "+"]
    gff_cds_plus <- gff_cds_plus[order(as.numeric(seqnames(gff_cds_plus)),
                                       start(gff_cds_plus))]
    gff_cds_plus_parent <- unlist(gff_cds_plus$Parent)
    target_i <- which(!duplicated(gff_cds_plus_parent))
    gff_cds_plus$phase[target_i] <- 0
    next_phase <- (3 - (width(gff_cds_plus[target_i]) - gff_cds_plus$phase[target_i]) %% 3) %% 3
    names(next_phase) <- gff_cds_plus_parent[target_i]
    gff_cds_plus_parent[target_i] <- "NA"

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

#' @importFrom BiocGenerics end width
#' @importFrom GenomeInfoDb seqnames
.phaseMinus <- function(gff_cds){
    gff_cds_minus <- gff_cds[as.character(strand(gff_cds)) == "-"]
    gff_cds_minus <- gff_cds_minus[order(as.numeric(seqnames(gff_cds_minus)),
                                         end(gff_cds_minus), decreasing = TRUE)]
    gff_cds_minus_parent <- unlist(gff_cds_minus$Parent)
    target_i <- which(!duplicated(gff_cds_minus_parent))
    gff_cds_minus$phase[target_i] <- 0
    next_phase <- (3 - (width(gff_cds_minus[target_i]) - gff_cds_minus$phase[target_i]) %% 3) %% 3
    names(next_phase) <- gff_cds_minus_parent[target_i]
    gff_cds_minus_parent[target_i] <- "NA"

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
