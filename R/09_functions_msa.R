#' @export
#' 
#' 
compareOrthoSeq <- function(hdf5_fn, graph_df = NULL, n_threads = 1, verbose = TRUE){
    h5 <- H5Fopen(hdf5_fn)
    on.exit(H5Fclose(h5))
    
    check <- h5$data_type == "reorg_orthopair"
    if(verbose){
        message("The input data type: reorg_orthopair")
    }
    if(check){
        in_files <- .orgInputFiles(hdf5_fn = hdf5_fn, h5 = h5)
        if(verbose){
            message("The number of genome pairs: ", 
                    length(in_files$comb_id$pair_id))
        }
        out <- NULL
        if(!is.null(graph_df)){
            if(is.character(graph_df)){
                graph_df <- read.csv(file = graph_df)
            }
        }
        
        for(i in seq_along(in_files$comb_id$pair_id)){
            pair_id <- in_files$comb_id$pair_id[i]
            if(verbose){
                message("Running for the following genome pair: ", pair_id)
            }
            query_id <- in_files$comb_id$query[i]
            subject_id <- in_files$comb_id$subject[i]
            genome1_index <- which(in_files$genome_names %in% query_id)
            genome2_index <- which(in_files$genome_names %in% subject_id)
            
            if(!is.null(graph_df)){
                orthopair_gene <- data.frame(query_tx = graph_df[[paste0(query_id, "_TX")]], 
                                             subject_tx = graph_df[[paste0(subject_id, "_TX")]], 
                                             SOG = graph_df$SOG)
                orthopair_gene <- subset(orthopair_gene,
                                         subset = !is.na(query_tx) & !is.na(subject_tx))
                orthopair_gene <- unique(orthopair_gene)
                
            } else {
                orthopair_gene <- h5$orthopair_gene[[pair_id]]
            }
            out_i <- .compareSeqEngine(cds1 = in_files$cds[genome1_index],
                                       prot1 = in_files$prot[genome1_index],
                                       promoter1 = in_files$promoter[genome1_index],
                                       genome1 = in_files$genome[[genome1_index]],
                                       gff1 = in_files$gff[genome1_index],
                                       cds2 = in_files$cds[genome2_index],
                                       prot2 = in_files$prot[genome2_index],
                                       promoter2 = in_files$promoter[genome2_index],
                                       genome2 = in_files$genome[[genome2_index]],
                                       gff2 = in_files$gff[genome2_index],
                                       orthopair_gene = orthopair_gene,
                                       query_id = query_id,
                                       subject_id = subject_id,
                                       pair_id = pair_id,
                                       n_threads = n_threads,
                                       verbose = verbose)
            out <- c(out, list(out_i))
        }
        names(out) <- in_files$comb_id$pair_id
        
    } else {
        promoter1 <- h5$files$query_promoter
        promoter2 <- h5$files$subject_promoter
        promoter1[is.null(promoter1)] <- "null/path"
        promoter2[is.null(promoter2)] <- "null/path"
        out <- .compareSeqEngine(cds1 = h5$files$query_cds,
                                 prot1 = h5$files$query_prot,
                                 promoter1 = promoter1,
                                 genome1 = h5$files$query_genome,
                                 gff1 = h5$files$query_gff,
                                 cds2 = h5$files$subject_cds,
                                 prot2 = h5$files$subject_prot,
                                 promoter2 = promoter2,
                                 genome2 = h5$files$subject_genome,
                                 gff2 = h5$files$subject_gff,
                                 orthopair_gene = h5$orthopair_gene,
                                 query_id = "",
                                 subject_id = "",
                                 pair_id = "",
                                 n_threads = n_threads,
                                 verbose = verbose)
        out <- list(out)
        names(out) <- ""
    }
    return(out)
}

.orgInputFiles <- function(hdf5_fn, h5){
    genome_names <- names(h5$genome_fn)
    gff_list <- NULL
    cds_list <- NULL
    prot_list <- NULL
    promoter_list <- NULL
    genome_list <- h5$genome_fn
    for(i in seq_along(genome_names)){
        gff_list <- c(gff_list,
                      file.path(dirname(hdf5_fn), 
                                paste0(genome_names[i], "_orthopair.gff")))
        cds_list <- c(cds_list,
                      file.path(dirname(hdf5_fn), 
                                paste0(genome_names[i], "_orthopair.cds")))
        prot_list <- c(prot_list,
                       file.path(dirname(hdf5_fn), 
                                 paste0(genome_names[i], "_orthopair.prot")))
        promoter_list <- c(promoter_list,
                           file.path(dirname(hdf5_fn), 
                                     paste0(genome_names[i], "_orthopair.promoter")))}
    
    comb_id <- data.frame(pair_id = names(h5$orthopair_gene))
    comb_id$query <- sub("_.+", "", comb_id$pair_id)
    comb_id$subject <- sub(".+_", "", comb_id$pair_id)
    
    out <- list(gff = gff_list, 
                cds = cds_list,
                prot = prot_list,
                promoter = promoter_list,
                genome = genome_list,
                genome_names = genome_names,
                comb_id = comb_id)
    return(out)
}

#' @importFrom parallel mclapply
#' @importFrom Biostrings readAAStringSet
.compareSeqEngine <- function(cds1, prot1, promoter1, genome1, gff1, 
                              cds2, prot2, promoter2, genome2, gff2,
                              orthopair_gene,
                              query_id = "", subject_id = "", pair_id = "",
                              n_threads, verbose){
    
    if(!file.exists(promoter1)){
        if(verbose){
            message("Preparing a promoter FASTA file for the query genome.")
        }
        out_fn <- file.path(dirname(gff1), 
                            paste0(query_id, "_orthopair.promoter"))
        .getPromoterSeq(gff = gff1, genome = genome1, out_fn = out_fn)
    }
    
    if(!file.exists(promoter2)){
        if(verbose){
            message("Preparing a promoter FASTA file for the subject genome.")
        }
        out_fn <- file.path(dirname(gff2), 
                            paste0(subject_id, "_orthopair.promoter"))
        .getPromoterSeq(gff = gff2, genome = genome2, out_fn = out_fn)
    }
    
    file_list <- c(cds1, cds2, prot2, prot2, promoter1, promoter2)
    check <- file.exists(file_list)
    if(any(!check)){
        stop("Following file(s) not exist: \n", 
             paste(file_list[!check], collapse = "\n"), 
             call. = FALSE)
    }
    
    
    if(verbose){
        message("Performing pairwise protein sequence alignments.")
    }
    out_dir <- file.path(dirname(gff1), paste("aa_msa", pair_id, sep = "_"))
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    prot1 <- readAAStringSet(filepath = prot1)
    prot2 <- readAAStringSet(filepath = prot2)
    names(prot1) <- sub("\\s.+", "", names(prot1))
    names(prot2) <- sub("\\s.+", "", names(prot2))
    aa_out <- mclapply(X = seq_along(orthopair_gene$SOG),
                       mc.cores = n_threads,
                       FUN = .doMSA, 
                       seq1 = prot1,
                       seq2 = prot2,
                       orthopair_gene = orthopair_gene,
                       type = "aa",
                       out_dir = out_dir)
    aa_out <- do.call("rbind", aa_out)
    aa_out_fn <- file.path(out_dir, "msa_AA_summary.csv")
    write.csv(aa_out, aa_out_fn, row.names = FALSE)
    
    if(verbose){
        message("Performing pairwise CDS alignments.")
    }
    out_dir <- file.path(dirname(gff1), paste("cds_msa", pair_id, sep = "_"))
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    cds1 <- readDNAStringSet(filepath = cds1)
    cds2 <- readDNAStringSet(filepath = cds2)
    names(cds1) <- sub("\\s.+", "", names(cds1))
    names(cds2) <- sub("\\s.+", "", names(cds2))
    cds_out <- mclapply(X = seq_along(orthopair_gene$SOG),
                        mc.cores = n_threads,
                        FUN = .doMSA, 
                        seq1 = cds1,
                        seq2 = cds2,
                        orthopair_gene = orthopair_gene,
                        type = "cds",
                        out_dir = out_dir)
    cds_out <- do.call("rbind", cds_out)
    cds_out_fn <- file.path(out_dir, "msa_CDS_summary.csv")
    write.csv(cds_out, cds_out_fn, row.names = FALSE)
    
    if(verbose){
        message("Performing pairwise promoter sequence alignments.")
    }
    out_dir <- file.path(dirname(gff1), paste("promoter_msa", pair_id, sep = "_"))
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    promoter1 <- readDNAStringSet(filepath = promoter1)
    promoter2 <- readDNAStringSet(filepath = promoter2)
    names(promoter1) <- sub("\\s.+", "", names(promoter1))
    names(promoter2) <- sub("\\s.+", "", names(promoter2))
    promoter_out <- mclapply(X = seq_along(orthopair_gene$SOG),
                             mc.cores = n_threads,
                             FUN = .doMSA, 
                             seq1 = promoter1,
                             seq2 = promoter2,
                             orthopair_gene = orthopair_gene,
                             type = "promoter",
                             out_dir = out_dir)
    promoter_out <- do.call("rbind", promoter_out)
    promoter_out_fn <- file.path(out_dir, "msa_PROMOTER_summary.csv")
    write.csv(promoter_out, promoter_out_fn, row.names = FALSE)
    
    out <- c(aa = aa_out_fn, cds = cds_out_fn, promoter = promoter_out_fn)
    return(out)
}

#' @importFrom rtracklayer import.gff3
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics start end
#' @importFrom S4Vectors mcols
#' @import BSgenome
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @importFrom rhdf5 H5Fopen H5Fclose
#' 
#' @export
.getPromoterSeq <- function(gff, genome, out_fn){
    genome <- readDNAStringSet(genome)
    gff <- import.gff3(gff)
    gff <- gff[gff$type == "CDS"]
    subset_df <- mcols(gff)
    subset_df$strand <- as.character(strand(gff))
    subset_df$start <- start(gff)
    subset_df$end <- end(gff)
    start_pos <- tapply(subset_df$start,
                        unlist(gff$Parent),
                        min)
    end_pos <- tapply(subset_df$end,
                      unlist(gff$Parent),
                      max)
    df <- data.frame(id = names(start_pos),
                     start = start_pos,
                     end = end_pos)
    hit <- match(df$id, unlist(gff$Parent))
    df$strand <- subset_df$strand[hit]
    df$gene_ID <- subset_df$gene_id[hit]
    df$chr <- as.character(seqnames(gff)[hit])
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
    
    ol <- findOverlaps(promoters, gff)
    ol <- as.data.frame(ol)
    ol$queryHits <- promoters$gene_ID[ol$queryHits]
    ol$subjectHits <- gff$gene_id[ol$subjectHits]
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
    
    writeXStringSet(promoter_seq, filepath = out_fn)
    return(out_fn)
}

#' @importFrom pwalign pairwiseAlignment writePairwiseAlignments coverage
#' 
.doMSA <- function(index, seq1, seq2, orthopair_gene, type, out_dir){
    if(type == "aa"){
        prefix <- "p."
        sufix <- "aa_msa.fa"
        label <- "AA_"
        
    } else if(type == "cds"){
        prefix <- "c."
        sufix <- "cds_msa.fa"
        label <- "CDS_"
        
    } else if(type == "promoter"){
        prefix <- "g."
        sufix <- "promoter_msa.fa"
        label <- "PROMOTER_"
    }
    
    seq1_hit <- names(seq1) == orthopair_gene$query_tx[index]
    seq2_hit <- names(seq2) == orthopair_gene$subject_tx[index]
    
    if(sum(seq1_hit) == 0 | sum(seq2_hit) == 0){
        if(sum(seq1_hit) != 0){
            seq1 <- seq1[seq1_hit]
            seq1_len <- width(seq1)
        } else {
            seq1_len <- NA
        }
        if(sum(seq2_hit) != 0){
            seq2 <- seq2[seq2_hit]
            seq2_len <- width(seq2)
        } else {
            seq2_len <- NA
        }
        baseinfo <- cbind(seq1 = orthopair_gene$query_tx[index],
                     seq2 = orthopair_gene$subject_tx[index], 
                     seq1_len = seq1_len, seq2_len = seq2_len)
        
        if(type == "aa"){
            baseinfo <- cbind(baseinfo,
                         seq1_init = NA,
                         seq1_term = NA, 
                         seq1_trancated = NA, 
                         seq2_init = NA, 
                         seq2_term = NA, 
                         seq2_trancated = NA)
        }
        out <- cbind(aln_len = 0,
                     .evalDel(pwa = NA), 
                     .evalIns(pwa = NA), 
                     .evalMismatch(pwa = NA))
        
    } else {
        seq1 <- seq1[seq1_hit]
        seq2 <- seq2[seq2_hit]
        
        if(type == "aa"){
            seq1_eval <- .evalSeq(seq = seq1)
            seq2_eval <- .evalSeq(seq = seq2)
            seq1 <- seq1_eval$seq
            seq2 <- seq2_eval$seq
            baseinfo <- cbind(seq1 = names(seq1), seq2 = names(seq2), 
                              seq1_len = width(seq1), seq2_len = width(seq2), 
                              seq1_init = seq1_eval$init, 
                              seq1_term = seq1_eval$term, 
                              seq1_trancated = seq1_eval$trancated, 
                              seq2_init = seq2_eval$init, 
                              seq2_term = seq2_eval$term, 
                              seq2_trancated = seq2_eval$trancated)
        } else {
            baseinfo <- cbind(seq1 = names(seq1), seq2 = names(seq2), 
                              seq1_len = width(seq1), seq2_len = width(seq2))
        }
        
        if(width(seq1) == 0 | width(seq2) == 0){
            out <- cbind(aln_len = 0,
                         .evalDel(pwa = NA), 
                         .evalIns(pwa = NA), 
                         .evalMismatch(pwa = NA))
            
        } else {
            pwa <- pairwiseAlignment(pattern = seq1, subject = seq2)
            out_fn <- file.path(out_dir,
                                paste(names(seq1), names(seq2), sufix, sep = "_"))
            writePairwiseAlignments(pwa, out_fn)
            del <- .evalDel(pwa = pwa, prefix = prefix)
            ins <- .evalIns(pwa = pwa, prefix = prefix)
            mis <- .evalMismatch(pwa = pwa, prefix = prefix)
            out <- cbind(aln_len = length(coverage(pwa)),
                         del, ins, mis)
        }
    }
    
    names(out) <- paste0(label, names(out))
    out <- cbind(baseinfo, out)
    return(out)
}

#' @importFrom Biostrings gregexpr2 substr
.evalSeq <- function(seq){
    init <- substr(seq, 1, 1) == "M"
    term_hit <- gregexpr2(pattern = "*", text = seq)[[1]][1]
    term <- term_hit != -1
    len <- width(seq)
    if(term){
        trancated <- term_hit < len
        if(trancated){
            end_seq <- substr(seq, term_hit, len)
            n_end_seq <- nchar(end_seq)
            end_seq <- gregexpr2("*", end_seq)[[1]]
            only_asteriscs <- length(end_seq) == n_end_seq
            if(only_asteriscs){
                trancated <- FALSE
                len <- term_hit
            }
            seq <- substr(seq, 1, term_hit)
        }
    } else {
        trancated <- FALSE
    }
    out <- list(seq = seq, init = init, term = term, trancated = trancated)
    return(out)
}

#' @importFrom pwalign insertion aligned alignedSubject
#' @importFrom BiocGenerics start end width
#' 
.evalDel <- function(pwa, prefix){
    if(!inherits(x = pwa, "PairwiseAlignments")){
        out <- data.frame(Del = NA, 
                          Del_Length = NA, 
                          Del_Number = NA)
        return(out)
    }
    
    del <- insertion(x = pwa)[[1]]
    non_del_len <- c(start(del)[1] - 1, diff(start(del)))
    end_in_query <- cumsum(width(del) + non_del_len)
    start_in_query <- end_in_query - width(del) + 1
    end(del) <- end_in_query
    start_in_subject <- start(del)
    start(del) <- start_in_query
    
    out <- NULL
    l_del <- sum(width(del), na.rm = TRUE)
    n_del <- length(del)
    
    seq2_char <- as.character(alignedSubject(x = pwa))
    head_del <- regexpr(pattern = "^-+", text = seq2_char)
    if(head_del != -1){
        out <- c(out,
                 paste0(prefix, "start.", 
                        "del:", attributes(head_del)$match.length))
        l_del <- l_del + attributes(head_del)$match.length
        n_del < n_del + 1
    }
    
    tail_del <- regexpr(pattern = "-+$", text = seq2_char)
    if(tail_del != -1){
        out <- c(out, 
                 paste0(prefix, "end.", 
                        "del:", attributes(tail_del)$match.length))
        l_del <- l_del + attributes(tail_del)$match.length
        n_del < n_del + 1
    }
    
    if(length(del) > 0){
        out <- c(out, paste0(prefix, 
                             start(del), 
                             "_",
                             end(del), 
                             ".del.",
                             start_in_subject,
                             ":",
                             width(del)))
    }
    
    if(is.null(out)){
        out <- NA
    } else {
        out <- paste(out, collapse = ",")
    }
    out <- data.frame(Del = out, 
                      Del_Length = l_del, 
                      Del_Number = n_del)
    return(out)
}

#' @importFrom pwalign deletion aligned alignedPattern
#' @importFrom BiocGenerics start end width
#' 
.evalIns <- function(pwa, prefix){
    if(!inherits(x = pwa, "PairwiseAlignments")){
        out <- data.frame(Ins = NA,
                          Ins_Length = NA, 
                          Ins_Number = NA)
        return(out)
    }
    
    ins <- deletion(x = pwa)[[1]]
    non_ins_len <- c(start(ins)[1] - 1, diff(start(ins)))
    end_in_query <- cumsum(width(ins) + non_ins_len)
    start_in_query <- end_in_query - width(ins) + 1
    end(ins) <- end_in_query
    start_in_subject <- start(ins)
    start(ins) <- start_in_query
    
    out <- NULL
    l_ins <- sum(width(ins), na.rm = TRUE)
    n_ins <- length(ins)
    
    seq1_char <- as.character(alignedPattern(x = pwa))
    head_ins <- regexpr(pattern = "^-+", text = seq1_char)
    if(head_ins != -1){
        out <- c(out,
                 paste0(prefix, "start.", 
                        "ins:", attributes(head_ins)$match.length))
        l_ins <- l_ins + attributes(head_ins)$match.length
        n_ins <- n_ins + 1
    }
    
    tail_ins <- regexpr(pattern = "-+$", text = seq1_char)
    if(tail_ins != -1){
        out <- c(out, 
                 paste0(prefix, "end.", 
                        "ins:", attributes(tail_ins)$match.length))
        l_ins <- l_ins + attributes(tail_ins)$match.length
        n_ins <- n_ins + 1
    }
    
    if(length(ins) > 0){
        out <- c(out, paste0(prefix, 
                             start(ins), 
                             "_",
                             end(ins), 
                             ".ins.",
                             start_in_subject,
                             ":",
                             width(ins)))
    }
    
    if(is.null(out)){
        out <- NA
    } else {
        out <- paste(out, collapse = ",")
    }
    out <- data.frame(Ins = out, 
                      Ins_Length = l_ins, 
                      Ins_Number = n_ins)
    return(out)
}

#' @importFrom pwalign mismatchTable nmismatch
#' @importFrom BiocGenerics start end
.evalMismatch <- function(pwa, prefix){
    if(!inherits(x = pwa, "PairwiseAlignments")){
        out <- data.frame(Mis = NA, 
                          Mis_Number = NA)
        return(out)
    }
    
    missense <- mismatchTable(pwa)
    
    if(nrow(missense) == 0){
        out <- NA
        n_mis <- 0
        
    } else {
        n_mis <- nmismatch(pwa)
        
        out <- paste0(prefix, 
                      missense$SubjectStart, 
                      missense$SubjectSubstring, 
                      ">",
                      missense$PatternStart, 
                      missense$PatternSubstring)
    }
    out <- paste(out, collapse = ",")
    out <- data.frame(Mis = out, 
                      Mis_Number = n_mis)
    return(out)
}
