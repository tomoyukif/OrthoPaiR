#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics start end
#' @importFrom S4Vectors mcols
#' @importFrom Biostrings writeXStringSet
#' @importFrom rhdf5 H5Fopen H5Fclose
#' 
#' @export
.getPromoterSeq <- function(hdf5_fn, genome, out_fn){
    h5 <- H5Fopen(hdf5_fn[i_comb_id])
    on.exit(H5Fclose(h5))
    
    
    gff <- gff[gff$type != "gene"]
    subset_gff <- gff[unlist(gff$Parent) %in% id]
    subset_gff <- subset_gff[subset_gff$type == "CDS"]
    subset_df <- mcols(subset_gff)
    subset_df$strand <- as.character(strand(subset_gff))
    subset_df$start <- start(subset_gff)
    subset_df$end <- end(subset_gff)
    start_pos <- tapply(subset_df$start,
                        unlist(subset_gff$Parent),
                        min)
    end_pos <- tapply(subset_df$end,
                      unlist(subset_gff$Parent),
                      max)
    df <- data.frame(id = names(start_pos),
                     start = start_pos,
                     end = end_pos)
    hit <- match(df$id, unlist(subset_gff$Parent))
    df$strand <- subset_df$strand[hit]
    df$gene_ID <- subset_df$gene_id[hit]
    df$chr <- as.character(seqnames(subset_gff)[hit])
    promoters_plus <- GRanges(seqnames = df$chr,
                              ranges = IRanges(end = df$start - 1,
                                               width = 3000),
                              strand = "+",
                              ID = df$id, gene_ID = df$gene_ID)
    promoters_minus <- GRanges(seqnames = df$chr,
                               ranges = IRanges(start = df$end + 1,
                                                width = 3000),
                               strand = "-",
                               ID = df$id, gene_ID = df$gene_ID)
    promoters <- promoters_plus
    promoters[df$strand == "-"] <- promoters_minus[df$strand == "-"]
    
    ol <- findOverlaps(promoters, subset_gff)
    ol <- as.data.frame(ol)
    ol$queryHits <- promoters$gene_ID[ol$queryHits]
    ol$subjectHits <- subset_gff$gene_id[ol$subjectHits]
    ol <- unique(ol[ol$queryHits != ol$subjectHits, ])
    
    hit <- match(ol$queryHits, df$gene_ID)
    ol$query_start <- df$start[hit]
    ol$query_end <- df$end[hit]
    ol$query_strand <- df$strand[hit]
    hit <- match(ol$subjectHits, df$gene_ID)
    ol$ol_start <- df$start[hit]
    ol$ol_end <- df$end[hit]
    ol$ol_strand <- df$strand[hit]
    ol <- ol[ol$query_strand == ol$ol_strand, ]
    
    ol_plus <- ol[ol$query_strand == "+", ]
    plus_target <- df$strand == "+" & df$gene_ID %in% ol_plus$queryHits
    hit <- match(promoters$gene_ID[plus_target], ol_plus$queryHits)
    plus_target_end <- end(promoters[plus_target])
    ol_plus_hit_end <- ol_plus$ol_end[hit]
    valid <- plus_target_end - ol_plus_hit_end < 3000 & plus_target_end - ol_plus_hit_end > 0
    start(promoters[plus_target][valid]) <- ol_plus$ol_end[hit][valid] + 1
    
    ol_minus <- ol[ol$query_strand == "-", ]
    minus_target <- df$strand == "-" & df$gene_ID %in% ol_minus$queryHits
    hit <- match(promoters$gene_ID[minus_target], ol_minus$queryHits)
    minus_target_start <- start(promoters[minus_target])
    ol_minus_hit_start <- ol_minus$ol_start[hit]
    valid <- ol_minus_hit_start - minus_target_start < 3000 & ol_minus_hit_start - minus_target_start > 0
    end(promoters[minus_target][valid]) <- ol_minus$ol_start[hit][valid] - 1
    
    invalid_start <- start(promoters) < 1
    start(promoters)[invalid_start] <- 1
    invalid_end <- end(promoters) > width(genome[seqnames(promoters)])
    end(promoters)[invalid_end] <- width(genome[seqnames(promoters[invalid_end])])
    promoter_seq <- genome[promoters]
    names(promoter_seq) <- paste(promoters$ID,
                                 as.character(seqnames(promoters)),
                                 start(promoters),
                                 end(promoters),
                                 as.character(strand(promoters)), sep = " ")
    
    writeXStringSet(promoter, filepath = out_fn)
    return(out_fn)
}


################################################################################
#' @importFrom Biostrings AAStringSet writeXStringSet
#' @importFrom mas msaClustalW
#' 
#' @export
doAAmsa <- function(index, og, aa1, aa2, out_dir){
    seq1 <- aa1[names(aa1) %in% og$query[index]]
    seq2 <- aa2[names(aa2) %in% og$subject[index]]
    queries <- c(seq1, seq2)
    queries <- AAStringSet(gsub("*", "", queries, fixed = TRUE))
    out <- msaClustalW(queries)
    out_fn <- file.path(out_dir,
                        paste(names(seq1), names(seq2), "aa_msa.fa", sep = "_"))
    writeXStringSet(AAStringSet(out), out_fn)
}

#' @importFrom Biostrings DNAStringSet writeXStringSet
#' @importFrom mas msaClustalW
#' 
#' @export
doCDSmsa <- function(index, og, cds1, cds2, out_dir){
    seq1 <- cds1[names(cds1) %in% og$query[index]]
    seq2 <- cds2[names(cds2) %in% og$subject[index]]
    queries <- c(seq1, seq2)
    queries <- DNAStringSet(gsub("*", "", queries, fixed = TRUE))
    out <- msaClustalW(queries)
    out_fn <- file.path(out_dir,
                        paste(names(seq1), names(seq2), "cds_msa.fa", sep = "_"))
    writeXStringSet(DNAStringSet(out), out_fn)
}

#' @importFrom Biostrings DNAStringSet writeXStringSet
#' @importFrom mas msaClustalW
#' 
#' @export
doPROMOTERmsa <- function(index, og, promoter1, promoter2, out_dir){
    seq1 <- promoter1[names(promoter1) %in% og$query[index]]
    seq2 <- promoter2[names(promoter2) %in% og$subject[index]]
    queries <- c(seq1, seq2)
    queries <- DNAStringSet(gsub("*", "", queries, fixed = TRUE))
    out <- msaClustalW(queries)
    out_fn <- file.path(out_dir,
                        paste(names(seq1), names(seq2), "promoter_msa.fa", sep = "_"))
    writeXStringSet(DNAStringSet(out), out_fn)
}
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics start end
#' @importFrom S4Vectors mcols
#' @importFrom Biostrings writeXStringSet
#' 
#' @export
evalAAmsa <- function(msa_fn){
    msa <- readAAMultipleAlignment(msa_fn)
    obj <- .makeMSAobj(msa = msa)
    
    del <- .evalDel(obj = obj, prefix = "p.")
    ins <- .evalIns(obj = obj, prefix = "p.")
    mis <- .evalMismatch(obj = obj, prefix = "p.")
    out <- cbind(seq1 = obj$id[1], seq2 = obj$id[2], 
                 aln_len = nrow(obj$df),
                 del, ins, mis)
    return(out)
}

evalCDSmsa <- function(msa_fn){
    msa <- readDNAMultipleAlignment(msa_fn)
    obj <- .makeMSAobj(msa = msa)
    
    del <- .evalDel(obj = obj, prefix = "c.")
    ins <- .evalIns(obj = obj, prefix = "c.")
    mis <- .evalMismatch(obj = obj, prefix = "c.")
    out <- cbind(seq1 = obj$id[1], seq2 = obj$id[2], 
                 aln_len = nrow(obj$df),
                 del, ins, mis)
    names(out) <- sub("AA_", "DNA_", names(out))
    return(out)
}

evalPROMOTERmsa <- function(msa_fn){
    msa <- readDNAMultipleAlignment(msa_fn)
    obj <- .makeMSAobj(msa = msa)
    
    del <- .evalDel(obj = obj, prefix = "g.")
    ins <- .evalIns(obj = obj, prefix = "g.")
    mis <- .evalMismatch(obj = obj, prefix = "g.")
    out <- cbind(seq1 = obj$id[1], seq2 = obj$id[2],
                 aln_len = nrow(obj$df),
                 del, ins, mis)
    names(out) <- sub("AA_", "DNA_", names(out))
    return(out)
}

.makeMSAobj <- function(msa){
    id <- rownames(msa)
    msa_mat <- as.matrix(msa)
    msa1 <- msa_mat[1, ]
    msa2 <- msa_mat[2, ]
    df <- data.frame(msa_pos = seq_along(msa1), 
                     seq1_pos = msa1 != "-",
                     seq2_pos = msa2 != "-" )
    df$seq1_pos[df$seq1_pos] <- seq_len(sum(df$seq1_pos))
    df$seq1_pos[df$seq1_pos == 0] <- NA
    df$seq2_pos[df$seq2_pos] <- seq_len(sum(df$seq2_pos))
    df$seq2_pos[df$seq2_pos == 0] <- NA
    out <- list(msa1 = msa1, msa2 = msa2, df = df, id = id)
    return(out)
}

#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics start end width
#' 
#' @export
.evalDel <- function(obj, prefix){
    out <- NULL
    del <- obj$msa1 != "-" & obj$msa2 == "-"
    del <- reduce(IRanges(start = which(del), width = 1))
    
    if(length(del) == 0){
        out <- NA
        l_del <- 0
        n_del <- 0
        
    } else {
        l_del <- sum(width(del), na.rm = TRUE)
        n_del <- length(del)
        
        head_del <- start(del) == 1
        if(any(head_del)){
            out <- c(out,
                     paste0(prefix, "start.", 
                            "del.", width(del[head_del])))
            del <- del[!head_del] 
        }
        tail_del <- end(del) == length(obj$msa2)
        if(any(tail_del)){
            out <- c(out, 
                     paste0(prefix, "end.", 
                            "del.", width(del[tail_del])))
            del <- del[!tail_del] 
        }
        
        if(length(del) > 0){
            start_hit <- match(start(del), obj$df$msa_pos)
            end_hit <- match(end(del), obj$df$msa_pos)
            out <- c(out, paste0(prefix, 
                                 obj$df$seq1_pos[start_hit], 
                                 "_",
                                 obj$df$seq1_pos[end_hit], 
                                 "del.",
                                 width(del)))
        }
    }
    out <- paste(out, collapse = ",")
    out <- data.frame(AA_Del = out, 
                      AA_Del_Length = l_del, 
                      AA_Del_Number = n_del)
    return(out)
}

#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics start end width
#' 
#' @export
.evalIns <- function(obj, prefix){
    out <- NULL
    ins <- obj$msa1 == "-" & obj$msa2 != "-"
    ins <- reduce(IRanges(start = which(ins), width = 1))
    if(length(ins) == 0){
        out <- NA
        l_ins <- 0
        n_ins <- 0
        
    } else {
        l_ins <- sum(width(ins), na.rm = TRUE)
        n_ins <- length(ins)
        
        head_ins <- start(ins) == 1
        if(any(head_ins)){
            out <- c(out, 
                     paste0(prefix, "start.", 
                            "ins.", width(ins[head_ins])))
            ins <- ins[!head_ins] 
        }
        tail_ins <- end(ins) == length(obj$msa1)
        if(any(tail_ins)){
            out <- c(out, 
                     paste0(prefix, "end.", 
                            "ins.", width(ins[tail_ins])))
            ins <- ins[!tail_ins] 
        }
        if(length(ins) > 0){
            start_hit <- match(start(ins), obj$df$msa_pos)
            end_hit <- match(end(ins), obj$df$msa_pos)
            out <- c(out, paste0(prefix, 
                                 obj$df$seq1_pos[start_hit - 1], 
                                 "_",
                                 obj$df$seq1_pos[end_hit + 1], 
                                 "ins.",
                                 width(ins)))
        }
    }
    out <- paste(out, collapse = ",")
    out <- data.frame(AA_Ins = out, 
                      AA_Ins_Length = l_ins, 
                      AA_Ins_Number = n_ins)
    return(out)
}

#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics start end
#' 
#' @export
.evalMismatch <- function(obj, prefix){
    missense <- obj$msa1 != obj$msa2 & obj$msa1 != "-" & obj$msa2 != "-"
    missense <- IRanges(start = which(missense), width = 1)
    if(length(missense) == 0){
        out <- NA
        n_mis <- 0
        
    } else {
        n_mis <- length(missense)
        
        start_hit <- match(start(missense), obj$df$msa_pos)
        end_hit <- match(end(missense), obj$df$msa_pos)
        out <- paste0(prefix, 
                      obj$df$seq1_pos[start_hit], 
                      obj$msa1[start(missense)], 
                      ">",
                      obj$msa2[start(missense)])
    }
    out <- paste(out, collapse = ",")
    out <- data.frame(AA_Mis = out, 
                      AA_Mis_Number = n_mis)
    return(out)
}
