#' @export
#' 
#' 
compareOrthoSeq <- function(working_dir,
                            ortholog_source = c("orthopair", "reorg_out"),
                            n_threads = 1,
                            verbose = TRUE){
    ortholog_source <- match.arg(ortholog_source)
    out_dir <- file.path(working_dir, "msa_out")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    
    in_files <- .orgInputFiles(working_dir = working_dir)
    orthopair_gene_list <- .loadOrthoPairGeneList(working_dir = working_dir,
                                                  ortholog_source = ortholog_source)
    pair_id <- names(orthopair_gene_list)
    if(verbose){
        message("The ortholog source: ", ortholog_source)
        message("The number of genome pairs: ", length(pair_id))
    }
    
    out <- NULL
    for(i in seq_along(pair_id)){
        pair_id_i <- pair_id[i]
        if(verbose){
            message("Running for the following genome pair: ", pair_id_i)
        }
        genomes <- strsplit(pair_id_i, "_", fixed = TRUE)[[1]]
        genome1_id <- genomes[1]
        genome2_id <- genomes[2]
        genome1_index <- match(genome1_id, in_files$genome_ids)
        genome2_index <- match(genome2_id, in_files$genome_ids)
        if(is.na(genome1_index) || is.na(genome2_index)){
            warning("Skip pair ", pair_id_i, " because genome IDs were not found in working_dir/input.")
            next
        }
        required <- c(in_files$cds[genome1_index], in_files$cds[genome2_index],
                      in_files$prot[genome1_index], in_files$prot[genome2_index],
                      in_files$genome[genome1_index], in_files$genome[genome2_index],
                      in_files$gff[genome1_index], in_files$gff[genome2_index])
        if(any(is.na(required))){
            warning("Skip pair ", pair_id_i, " due to missing CDS/PROT/GENOME/GFF files in input folders.")
            next
        }
        
        orthopair_gene <- orthopair_gene_list[[pair_id_i]]
        out_i <- .compareSeqEngine(cds1 = as.character(in_files$cds[genome1_index]),
                                   prot1 = as.character(in_files$prot[genome1_index]),
                                   promoter1 = as.character(in_files$promoter[genome1_index]),
                                   genome1 = as.character(in_files$genome[genome1_index]),
                                   gff1 = as.character(in_files$gff[genome1_index]),
                                   cds2 = as.character(in_files$cds[genome2_index]),
                                   prot2 = as.character(in_files$prot[genome2_index]),
                                   promoter2 = as.character(in_files$promoter[genome2_index]),
                                   genome2 = as.character(in_files$genome[genome2_index]),
                                   gff2 = as.character(in_files$gff[genome2_index]),
                                   orthopair_gene = orthopair_gene,
                                   genome1_id = genome1_id,
                                   genome2_id = genome2_id,
                                   pair_id = pair_id_i,
                                   out_dir = out_dir,
                                   n_threads = n_threads,
                                   verbose = verbose)
        out <- c(out, list(out_i))
    }
    names(out) <- pair_id
    return(out)
}

.orgInputFiles <- function(working_dir){
    input_dir <- file.path(working_dir, "input")
    input_folders <- list.dirs(path = input_dir, recursive = FALSE, full.names = TRUE)
    if(length(input_folders) == 0){
        stop("No input folders found under: ", input_dir, call. = FALSE)
    }
    base <- basename(input_folders)
    genome_idx <- suppressWarnings(as.integer(sub("_.*", "", base)))
    valid <- !is.na(genome_idx)
    input_folders <- input_folders[valid]
    base <- base[valid]
    genome_ids <- as.character(genome_idx[valid] + 1000L)
    
    pick_existing <- function(candidates){
        hit <- candidates[file.exists(candidates)]
        if(length(hit) == 0) return(NA_character_)
        hit[1]
    }
    
    cds <- vapply(input_folders, function(d){
        pick_existing(file.path(d, c("cds.fa", "cds.fasta", "cds.fna")))
    }, character(1))
    prot <- vapply(input_folders, function(d){
        pick_existing(file.path(d, c("prot.fa", "prot.faa", "protein.fa", "protein.faa")))
    }, character(1))
    gff <- vapply(input_folders, function(d){
        pick_existing(file.path(d, c("gff3", "gff", "annotation.gff3", "annotation.gff")))
    }, character(1))
    promoter <- vapply(seq_along(input_folders), function(i){
        file.path(input_folders[i], paste0(base[i], "_orthopair.promoter"))
    }, character(1))
    genome <- vapply(input_folders, function(d){
        pick_existing(file.path(d, c("genome.fa", "genome.fna", "genome.fasta", "genome.fa.gz")))
    }, character(1))
    
    names(cds) <- names(prot) <- names(gff) <- names(promoter) <- names(genome) <- genome_ids
    list(gff = gff,
         cds = cds,
         prot = prot,
         promoter = promoter,
         genome = genome,
         genome_ids = genome_ids)
}

.loadOrthoPairGeneList <- function(working_dir, ortholog_source){
    normalize_pair_id <- function(x){
        x <- sub("\\.tsv$", "", basename(x))
        parts <- strsplit(x, "_", fixed = TRUE)[[1]]
        paste(sort(parts), collapse = "_")
    }
    
    load_from_orthopair <- function(){
        orthopair_dir <- file.path(working_dir, "orthopair")
        fns <- list.files(orthopair_dir, pattern = "\\.tsv$", full.names = TRUE)
        skip <- c("orthopair_pairwise_mutual_ci_stats.tsv", "orthopair_genome_mean_mutual_ci_matrix.tsv")
        fns <- fns[!basename(fns) %in% skip]
        out <- list()
        for(fn in fns){
            dt <- data.table::fread(fn, sep = "\t", header = TRUE)
            need <- c("genome1_tx", "genome2_tx")
            if(!all(need %in% names(dt))) next
            dt <- unique(dt[, c("genome1_tx", "genome2_tx"), with = FALSE])
            dt <- dt[!is.na(genome1_tx) & !is.na(genome2_tx)]
            if(nrow(dt) == 0) next
            dt$SOG <- seq_len(nrow(dt))
            out[[normalize_pair_id(fn)]] <- as.data.frame(dt[, .(genome1_tx, genome2_tx, SOG)])
        }
        out
    }
    
    if(ortholog_source == "orthopair"){
        out <- load_from_orthopair()
        if(length(out) == 0){
            stop("No usable pairwise ortholog TSV with genome1_tx/genome2_tx in working_dir/orthopair.")
        }
        return(out)
    }
    
    pairwise_dir <- file.path(working_dir, "reorg_out", "pairwise")
    pair_fns <- list.files(pairwise_dir, pattern = "\\.tsv$", full.names = TRUE)
    if(length(pair_fns) == 0){
        stop("No pairwise ortholog TSV files in: ", pairwise_dir, call. = FALSE)
    }
    orthopair_list <- load_from_orthopair()
    if(length(orthopair_list) == 0){
        stop("reorg_out source requires working_dir/orthopair/*.tsv to recover genome1_tx/genome2_tx.")
    }
    out <- list()
    for(fn in pair_fns){
        pair_id <- normalize_pair_id(fn)
        pair_gene <- data.table::fread(fn, sep = "\t", header = TRUE)
        if(!all(c("genome1_gene", "genome2_gene") %in% names(pair_gene))) next
        op_fn <- file.path(working_dir, "orthopair", paste0(pair_id, ".tsv"))
        if(!file.exists(op_fn)) next
        op_dt <- data.table::fread(op_fn, sep = "\t", header = TRUE)
        if(!all(c("genome1_gene", "genome2_gene", "genome1_tx", "genome2_tx") %in% names(op_dt))) next
        m <- merge(pair_gene, op_dt[, .(genome1_gene, genome2_gene, genome1_tx, genome2_tx)],
                   by = c("genome1_gene", "genome2_gene"), all.x = TRUE)
        m <- unique(m[, .(genome1_tx, genome2_tx)])
        m <- m[!is.na(genome1_tx) & !is.na(genome2_tx)]
        if(nrow(m) == 0) next
        m$SOG <- seq_len(nrow(m))
        out[[pair_id]] <- as.data.frame(m)
    }
    if(length(out) == 0){
        stop("Failed to build transcript-level pair list from reorg_out pairwise TSV files.")
    }
    out
}

#' @importFrom parallel mclapply
#' @importFrom Biostrings readAAStringSet AAStringSet
.compareSeqEngine <- function(cds1, prot1, promoter1, genome1, gff1, 
                              cds2, prot2, promoter2, genome2, gff2,
                              orthopair_gene,
                              genome1_id = "", genome2_id = "", pair_id = "",
                              out_dir,
                              n_threads, verbose){
    
    if(!file.exists(promoter1)){
        if(verbose){
            message("Preparing a promoter FASTA file for genome1.")
        }
        promoter1 <- file.path(dirname(gff1),
                            paste0(genome1_id, "_orthopair.promoter"))
        .getPromoterSeq(gff = gff1, genome = genome1, out_fn = promoter1)
    }
    
    if(!file.exists(promoter2)){
        if(verbose){
            message("Preparing a promoter FASTA file for genome2.")
        }
        promoter2 <- file.path(dirname(gff2), 
                            paste0(genome2_id, "_orthopair.promoter"))
        .getPromoterSeq(gff = gff2, genome = genome2, out_fn = promoter2)
    }
    
    file_list <- c(cds1, cds2, prot1, prot2, promoter1, promoter2)
    check <- file.exists(file_list)
    if(any(!check)){
        stop("Following file(s) not exist: \n", 
             paste(file_list[!check], collapse = "\n"), 
             call. = FALSE)
    }
    
    
    if(verbose){
        message("Performing pairwise protein sequence alignments.")
    }
    out_aa <- file.path(out_dir, paste("aa_msa", pair_id, sep = "_"))
    dir.create(out_aa, showWarnings = FALSE, recursive = TRUE)
    prot1 <- readAAStringSet(filepath = prot1)
    prot2 <- readAAStringSet(filepath = prot2)
    names(prot1) <- sub("\\s.+", "", names(prot1))
    names(prot2) <- sub("\\s.+", "", names(prot2))
    prot1 <- AAStringSet(sub("\\*.+", "", prot1))
    prot2 <- AAStringSet(sub("\\*.+", "", prot2))
    aa_out <- mclapply(X = seq_along(orthopair_gene$SOG),
                       mc.cores = n_threads,
                       FUN = .doMSA, 
                       seq1 = prot1,
                       seq2 = prot2,
                       orthopair_gene = orthopair_gene,
                       type = "aa",
                       out_dir = out_aa)
    aa_out <- do.call("rbind", aa_out)
    aa_out_fn <- file.path(out_dir, paste0("aa_msa_summary_", pair_id, ".csv"))
    write.csv(aa_out, aa_out_fn, row.names = FALSE)
    
    if(verbose){
        message("Performing pairwise CDS alignments.")
    }
    out_cds <- file.path(out_dir, paste("cds_msa", pair_id, sep = "_"))
    dir.create(out_cds, showWarnings = FALSE, recursive = TRUE)
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
                        out_dir = out_cds)
    cds_out <- do.call("rbind", cds_out)
    cds_out_fn <- file.path(out_dir, paste0("cds_msa_summary_", pair_id, ".csv"))
    write.csv(cds_out, cds_out_fn, row.names = FALSE)
    
    if(verbose){
        message("Performing pairwise promoter sequence alignments.")
    }
    out_promoter <- file.path(out_dir, paste("promoter_msa", pair_id, sep = "_"))
    dir.create(out_promoter, showWarnings = FALSE, recursive = TRUE)
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
                             out_dir = out_promoter)
    promoter_out <- do.call("rbind", promoter_out)
    promoter_out_fn <- file.path(out_dir, paste0("promoter_msa_summary_", pair_id, ".csv"))
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
#' 
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
    ol$genome1_start <- df$start[hit]
    ol$genome1_end <- df$end[hit]
    ol$genome1_strand <- df$strand[hit]
    hit <- match(ol$subjectHits, df$gene_ID)
    ol$ol_start <- df$start[hit]
    ol$ol_end <- df$end[hit]
    ol$ol_strand <- df$strand[hit]
    ol <- ol[ol$genome1_strand == ol$ol_strand, ]
    
    ol_plus <- ol[ol$genome1_strand == "+", ]
    plus_target <- df$strand == "+" & df$gene_ID %in% ol_plus$queryHits
    hit <- match(promoters$gene_ID[plus_target], ol_plus$queryHits)
    plus_target_end <- end(promoters[plus_target])
    ol_plus_hit_end <- ol_plus$ol_end[hit]
    valid <- plus_target_end - ol_plus_hit_end < 3000 & plus_target_end - ol_plus_hit_end > 0
    start(promoters[plus_target][valid]) <- ol_plus$ol_end[hit][valid] + 1
    
    ol_minus <- ol[ol$genome1_strand == "-", ]
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
    
    seq1_hit <- names(seq1) == orthopair_gene$genome1_tx[index]
    seq2_hit <- names(seq2) == orthopair_gene$genome2_tx[index]
    
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
        baseinfo <- cbind(seq1 = orthopair_gene$genome1_tx[index],
                     seq2 = orthopair_gene$genome2_tx[index], 
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
