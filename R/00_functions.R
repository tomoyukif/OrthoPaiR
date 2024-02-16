library(Biostrings)
library(ggplot2)
library(rtracklayer)
library(rBLAST)
library(reticulate)
library(parallel)
library(dplyr)
library(rhdf5)
library(BSgenome)
library(GenomicFeatures)

.h5overwrite <- function(obj, file, name){
    file <- H5Fopen(file)
    on.exit(H5Fclose(file))
    if(H5Lexists(h5loc = file, name = name)){
        h5delete(file = file, name = name)
    }
    h5write(obj = obj, file = file, name = name)
}

.h5creategroup <- function(file, name){
    file <- H5Fopen(file)
    on.exit(H5Fclose(file))
    if(!H5Lexists(h5loc = file, name = name)){
        h5createGroup(file, name)
    }
}

.importAllGFF <- function(fn){
    for(i in seq_along(fn)){
        if(i == 1){
            out <- import.gff(fn[i])
        } else {
            out <- c(out, import.gff(fn[i]))
        }
    }
    return(out)
}

################################################################################
makeSynogDB <- function(query_genome, subject_genome,
                        query_gff, subject_gff,
                        query_cds, subject_cds,
                        query_prot, subject_prot,
                        positive_list, negative_list,
                        hdf5_path = "./synog.h5"){
    out <- list(query_genome = query_genome, subject_genome = subject_genome,
                query_gff = query_gff, subject_gff = subject_gff,
                query_cds = query_cds, subject_cds = subject_cds,
                query_prot = query_prot, subject_prot = subject_prot, 
                positive_list = positive_list, negative_list = negative_list)
    
    for(i in seq_along(out)){
        if(!file.exists(out[[i]])){
            stop(out[[i]], " do not exists!")
        }
    }
    
    out$genome$query <- .genomeSummary(genome = query_genome)
    out$genome$subject <- .genomeSummary(genome = subject_genome)
    
    class(out) <- c(class(out), "SynogDB")
    out$h5 <- .makeHDF5(hdf5_path = hdf5_path)
    return(out)
}

.genomeSummary <- function(genome){
    genome <- readDNAStringSet(genome)
    out <- list(names = names(genome),
                length = width(genome))
    return(out)
}

.makeHDF5 <- function(hdf5_path){
    if(file.exists(hdf5_path)){
        message(hdf5_path, " already exists.")
    } else {
        h5createFile(file = hdf5_path)
    }
    return(hdf5_path)
}

################################################################################
runSibeliaz <- function(object, out_dir,
                        sibeliaz_bin = "sibeliaz",
                        maf2synteny_bin = "maf2synteny",
                        conda = "conda",
                        condaenv = NULL,
                        run_sibeliaz = TRUE){
    stopifnot(inherits(x = object, "SynogDB"))
    
    if(run_sibeliaz){
        q_fn <- .prepGenome(genome = object$query_genome, label = "query")
        s_fn <- .prepGenome(genome = object$subject_genome, label = "subject")
        
        if(!is.null(condaenv)){
            .condaExe(conda = conda, env = condaenv, command = sibeliaz_bin, 
                      args = paste("-o", out_dir, "-n", q_fn, s_fn))
        } else {
            system2(command = sibeliaz_bin, 
                    args = paste("-o", out_dir, "-n", q_fn, s_fn))
        }
        
        if(!is.null(condaenv)){
            .condaExe(conda = conda, env = condaenv, command = maf2synteny_bin, 
                      args = paste("-o", out_dir,
                                   file.path(out_dir, "blocks_coords.gff")))
        } else {
            system2(command = maf2synteny_bin,
                    args = paste("-o", out_dir, 
                                 file.path(out_dir, "blocks_coords.gff")))
        }
    }
    
    .h5creategroup(object$h5,"sibeliaz")
    .h5overwrite(obj = file.path(out_dir, "blocks_coords.gff"),
                 file = object$h5, "sibeliaz/raw")
    .h5overwrite(obj = file.path(out_dir, "5000/blocks_coords.txt"),
                 file = object$h5, "sibeliaz/blocks")
}

.prepGenome <- function(genome, label){
    dna <- readDNAStringSet(genome)
    names(dna) <- paste(label, names(dna), sep = "_")
    out <- tempfile(fileext = ".fa")
    writeXStringSet(dna, out)
    return(out)
}

.condaExe <- function(conda, env, command, args){
    system2(command = conda,
            args = paste("run -n", env, command, args))
}

################################################################################
# Define a function to load sibeliaz's raw output to data.frame
sibeliaRAW2DF <- function(object){
    stopifnot(inherits(x = object, "SynogDB"))
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    if(!H5Lexists(h5, "sibeliaz/raw")){
        stop("Run runSibeliaz to obtain LCB info.")
    }
    gff <- import.gff3(h5$sibeliaz$raw)
    
    df <- data.frame(chr = as.character(seqnames(gff)), 
                     start = start(gff),
                     end =  end(gff), 
                     ID = gff$ID)
    
    out <- tapply(df$ID, df$ID, function(id){
        id <- id[1]
        tmp <- subset(df, subset = ID == id)
        q <- grep("query", tmp$chr)
        if(length(q) == 0){
            return(NULL)
        }
        s <- grep("subject", tmp$chr)
        if(length(s) == 0){
            return(NULL)
        }
        
        out <- expand.grid(q, s)
        out <- cbind(tmp[out$Var1, ], tmp[out$Var2, ])
        return(out)
    })
    
    out <- do.call("rbind", out)
    names(out) <- c(paste(names(out)[1:4], "query", sep = "_"), 
                    paste(names(out)[5:8], "subject", sep = "_"))
    
    out$query_chr <- sub("query_", "", out$query_chr)
    out$subject_chr <- sub("subject_", "", out$subject_chr)
    rownames(out) <- NULL
    class(out) <- c(class(out), "RawLCB")
    .h5overwrite(obj = out, file = object$h5, "sibeliaz/lcb")
}

################################################################################
# Define a function to convert sibelia's output to data.frame
sibeliaLCB2DF <- function(object){
    stopifnot(inherits(x = object, "SynogDB"))
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    if(!H5Lexists(h5, "sibeliaz/blocks")){
        stop("Run runSibeliaz to obtain LCB info.")
    }
    if(!file.exists(h5$sibeliaz$blocks)){
        stop(h5$sibeliaz$blocks, " does not exists.")
    }
    
    sibelia <- scan(file = h5$sibeliaz$blocks, 
                    what = "character", sep = "\n")
    sep_lines <- grep("^-+", sibelia)
    
    header_lines <- 1:(sep_lines[1] - 1)
    header <- sibelia[header_lines]
    header <- do.call("rbind", strsplit(header, "\t"))
    colnames(header) <- header[1, ]
    header <- header[-1, ]
    header <- data.frame(header)
    header$Seq_id <- as.numeric(header$Seq_id)
    header$Size <- as.numeric(header$Size)
    
    block_i <- grep("^Block", sibelia)
    i_diff <- diff(c(block_i, length(sibelia) + 1))
    entry_n <- sapply(seq_along(i_diff), function(i){
        return(rep(i, i_diff[i] - 3))
    })
    entry_n <- unlist(entry_n)
    block_header_i <- block_i[entry_n]
    block_id <- as.numeric(sub(".*#", "", sibelia[block_header_i]))
    
    entries <- sibelia[-header_lines][!grepl("^Block|^Seq_id|^-+",
                                             sibelia[-header_lines])]
    entries <- do.call("rbind", strsplit(entries, "\t"))
    out <- data.frame(entries, block_id)
    names(out) <- c("chr", "strand", "start", "end", "length", "block_id")
    out$chr <- header$Description[match(out$chr, header$Seq_id)]
    out$start <- as.numeric(out$start)
    out$end <- as.numeric(out$end)
    out$length <- as.numeric(out$length)
    class(out) <- c(class(out), "LCB")
    .h5overwrite(obj = out, file = object$h5, "sibeliaz/lcb")
}

################################################################################
# Define a function to show LCB information
showLCB <- function(object){
    stopifnot(inherits(x = object, "SynogDB"))
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    if(!H5Lexists(h5, "sibeliaz/lcb")){
        stop("Run sibeliaLCB2DF to obtain LCB info.")
    }
    
    message("Chromosome IDs:")
    q_id <- grep("^query_", unique(h5$sibeliaz$lcb$chr), value = TRUE)
    q_id <- sort(sub("^query_", "", q_id))
    s_id <- grep("^subject_", unique(h5$sibeliaz$lcb$chr), value = TRUE)
    s_id <- sort(sub("^subject_", "", s_id))
    message("    Query genome:")
    print(q_id)
    message("    Subject genome:")
    print(s_id)
    message("No of entries: ")
    print(length(h5$sibeliaz$lcb$chr))
    message("No of LCBs: ")
    print(length(unique(h5$sibeliaz$lcb$block_id)))
}

################################################################################
# Define a function to classify LCBs
lcbClassify <- function(object){
    stopifnot(inherits(x = object, "SynogDB"))
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    if(!H5Lexists(h5, "sibeliaz/lcb")){
        stop("Run sibeliaLCB2DF to obtain LCB info.")
    }
    
    chr <- rep(0, length(h5$sibeliaz$lcb$chr))
    s_chr <- grepl("^subject_", h5$sibeliaz$lcb$chr)
    chr[s_chr] <- 1
    lcb_class <- tapply(chr, 
                        h5$sibeliaz$lcb$block_id, 
                        function(i){
                            if(sum(i) == 0){
                                out <- "q2q"
                            } else if(prod(i) == 0){
                                out <- "q2s"
                            } else {
                                out <- "s2s"
                            }
                        })
    hit <- match(h5$sibeliaz$lcb$block_id, as.numeric(names(lcb_class)))
    .h5overwrite(obj = lcb_class[hit], file = object$h5, "sibeliaz/lcb_class")
}

################################################################################
# Define a function to summarize LCB statistics
statsLCB <- function(object){
    stopifnot(inherits(x = object, "SynogDB"))
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    if(!H5Lexists(h5, "sibeliaz/lcb")){
        stop("Run sibeliaLCB2DF to obtain LCB info.")
    }
    if(!H5Lexists(h5, "sibeliaz/lcb_class")){
        stop("Run lcbClassify to obtain LCB classification info.")
    }
    
    if(is.null(h5$sibeliaz$lcb_class)){
        q_summary <- .calc_lcb_stats(object = h5$sibeliaz$lcb,
                                     target = "^query_")
        s_summary <- .calc_lcb_stats(object = h5$sibeliaz$lcb,
                                     target = "^subject_")
        out <- rbind(cbind(genome = "query", q_summary),
                     cbind(genome = "subject", s_summary))
        
    } else {
        q2q <- h5$sibeliaz$lcb_class == "q2q"
        q_q2q_summary <- .calc_lcb_stats(object = h5$sibeliaz$lcb[q2q, ],
                                         target = "^query_")
        s2s <- h5$sibeliaz$lcb_class == "s2s"
        s_s2s_summary <- .calc_lcb_stats(object = h5$sibeliaz$lcb[s2s, ],
                                         target = "^subject_")
        q2s <- h5$sibeliaz$lcb_class == "q2s"
        q_q2s_summary <- .calc_lcb_stats(object = h5$sibeliaz$lcb[q2s, ],
                                         target = "^query_")
        s_q2s_summary <- .calc_lcb_stats(object = h5$sibeliaz$lcb[q2s, ],
                                         target = "^subject_")
        out <- rbind(cbind(genome = "query", class = "Q2Q", q_q2q_summary),
                     cbind(genome = "subject", class = "S2S", s_s2s_summary),
                     cbind(genome = "query", class = "Q2S", q_q2s_summary),
                     cbind(genome = "subject", class = "Q2S", s_q2s_summary))
    }
    .h5overwrite(obj = out, file = object$h5, "sibeliaz/lcb_stats")
    print(out)
}

.calc_lcb_stats <- function(object, target){
    if(is.null(object$genome)){
        tmp_df <- subset(object, subset = grepl(target, chr))
    } else {
        tmp_df <- subset(object, subset = genome == target)    
    }
    
    out <- data.frame(total_length = sum(tmp_df$length),
                      mean_length = mean(tmp_df$length),
                      max_length = max(tmp_df$length),
                      min_length = min(tmp_df$length))
    
    minus <- tmp_df$strand == "-"
    tmp <- tmp_df$start[minus]
    tmp_df$start[minus] <- tmp_df$end[minus]
    tmp_df$end[minus] <- tmp
    tmp_df_nb_gr <- GRanges(seqnames = tmp_df$chr, 
                            ranges = IRanges(start = tmp_df$start, 
                                             end = tmp_df$end))
    out <- cbind(out,
                 data.frame(total_gap_length = sum(width(gaps(tmp_df_nb_gr))),
                            mean_gap_length = mean(width(gaps(tmp_df_nb_gr))),
                            max_gap_length = max(width(gaps(tmp_df_nb_gr))),
                            min_gap_length = min(width(gaps(tmp_df_nb_gr))))
    )
    return(out)
}

################################################################################
# Define a function to make a table that show all pairs of blocks between 
# the query genome and the subject genome
getLCBpairs <- function(object){
    stopifnot(inherits(x = object, "SynogDB"))
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    if(!H5Lexists(h5, "sibeliaz/lcb")){
        stop("Run sibeliaLCB2DF to obtain LCB info.")
    }
    if(!H5Lexists(h5, "sibeliaz/lcb_class")){
        stop("Run lcbClassify to obtain LCB classification info.")
    }
    
    q2s <- h5$sibeliaz$lcb_class == "q2s"
    lcb_q2s_id_tbl <- table(h5$sibeliaz$lcb$block_id[q2s])
    lcb_id_1to1 <- names(lcb_q2s_id_tbl)[lcb_q2s_id_tbl == 2]
    valid <- h5$sibeliaz$lcb$block_id %in% lcb_id_1to1
    
    lcb_pairs <- list(lcb_1to1 = .get1to1(h5 = h5, valid = valid),
                      lcb_non_1to1 = .getNon1to1(h5 = h5, valid = valid))
    
    .h5overwrite(obj = lcb_pairs, file = object$h5, "sibeliaz/lcb_pairs")
}


.get1to1 <- function(h5, valid){
    hit <- valid & h5$sibeliaz$lcb_class == "q2s"
    df <- h5$sibeliaz$lcb[hit, ]
    q_df <- subset(df, subset = grepl("^query_", chr), 
                   select = c(chr, start, end, block_id))
    s_df <- subset(df, subset = grepl("^subject_", chr), 
                   select = c(chr, start, end, block_id))
    q_df <- q_df[order(q_df$block_id), ]
    s_df <- s_df[order(s_df$block_id), ]
    q_df$chr <- sub("^query_", "", q_df$chr)
    s_df$chr <- sub("^subject_", "", s_df$chr)
    out <- cbind(subset(q_df, select = chr:end), 
                 subset(s_df, select = c(chr:end, block_id)))
    names(out) <- c("query_chr", "query_start", "query_end",  
                    "subject_chr", "subject_start", "subject_end", "block_id")
    return(out)
}

.getNon1to1 <- function(h5, valid){
    hit <- !valid & h5$sibeliaz$lcb_class == "q2s"
    df <- h5$sibeliaz$lcb[hit, ]
    q_df <- subset(df, subset = grepl("^query_", chr), 
                   select = c(chr, start, end, block_id))
    s_df <- subset(df, subset = grepl("^subject_", chr), 
                   select = c(chr, start, end, block_id))
    q_df$chr <- sub("^query_", "", q_df$chr)
    s_df$chr <- sub("^subject_", "", s_df$chr)
    out <- full_join(q_df, s_df, by = "block_id", relationship = "many-to-many")
    names(out) <- c("query_chr", "query_start", "query_end", "block_id",
                    "subject_chr", "subject_start", "subject_end")
    out <- out[, c(1:3, 5:7, 4)]
    return(out)
}

################################################################################
# Plot LCB
# plotRawLCBpairs <- function(object, 
#                             size = 0.5, 
#                             color = "darkgreen", 
#                             chr_border_color = "gray60",
#                             chr_len = FALSE){
#     object$query_chr <- .sortLabels(x = object$query_chr, rev = FALSE)
#     object$subject_chr <- .sortLabels(x = object$subject_chr, rev = TRUE)
#     
#     p <- ggplot(object) +
#         geom_point(aes(x = (query_start + query_end) / 2,
#                        y = (subject_start + subject_end) / 2),
#                    color = color, size = size) +
#         facet_grid(rows = vars(subject_chr),
#                    cols = vars(query_chr),
#                    scales = "free") +
#         xlab("Query") +
#         ylab("Subject") +
#         scale_y_continuous(expand = c(0, 0)) +
#         scale_x_continuous(expand = c(0, 0)) +
#         theme(axis.text = element_blank(),
#               axis.ticks = element_blank(),
#               strip.text.y.right = element_text(angle = -90),
#               panel.spacing = unit(0, "lines"),
#               panel.border = element_rect(color = "gray60",
#                                           fill = NA,
#                                           linewidth = 0.5),
#               legend.position = "none")
#     
#     p <- .addChrLen(p = p, chrLen = object$genome)
#     
#     return(p)
# }

plotLCBpairs <- function(object, 
                         class = "both", 
                         size = 0.5, 
                         color_1to1 = "darkgreen", 
                         color_non1to1 = "blue",
                         chr_border_color = "gray60"){
    stopifnot(inherits(x = object, "SynogDB"))
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    if(!H5Lexists(h5, "sibeliaz/lcb_pairs")){
        stop("Run getLCBpairs to obtain LCB pair info.")
    }
    
    if(class == "1to1"){
        df <- h5$sibeliaz$lcb_pairs$lcb_1to1
        df$class <- "1to1"
        color_value <- color_1to1
        color_break <- "1to1"
        
    } else if(class == "non1to1"){
        df <- h5$sibeliaz$lcb_pairs$lcb_non_1to1
        df$class <- "non1to1"
        color_value <- color_non1to1
        color_break <- "non1to1"
        
    } else {
        df1 <- h5$sibeliaz$lcb_pairs$lcb_1to1
        df1$class <- "1to1"
        df2 <- h5$sibeliaz$lcb_pairs$lcb_non_1to1
        df2$class <- "non1to1"
        df <- rbind(df1, df2)
        color_value <- c(color_1to1, color_non1to1)
        color_break <- c("1to1", "non1to1")
    }
    
    df$query_pos <- (df$query_start + df$query_end) / 2
    df$subject_pos <- (df$subject_start + df$subject_end) / 2
    df$query_chr <- .sortLabels(x = df$query_chr, rev = FALSE)
    df$subject_chr <- .sortLabels(x = df$subject_chr, rev = TRUE)
    
    p <- ggplot(df) +
        geom_point(aes(x = query_pos, y = subject_pos, color = class),
                   size = size) +
        scale_color_manual(values = color_value,
                           breaks = color_break) +
        facet_grid(rows = vars(subject_chr),
                   cols = vars(query_chr),
                   scales = "free") +
        xlab("Query") +
        ylab("Subject") +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              strip.text.y.right = element_text(angle = -90),
              panel.spacing = unit(0, "lines"),
              panel.border = element_rect(color = "gray60",
                                          fill = NA,
                                          linewidth = 0.5),
              legend.position = "none")
    
    p <- .addChrLen(p = p, chrLen = object$genome)
    
    return(p)
}

.sortLabels <- function(x, rev = FALSE){
    if(rev){
        lev <- sort(unique(x))
        rev_lev <- rev(lev)
        x <- factor(x = x, levels = rev_lev)
        
    } else {
        lev <- sort(unique(x))
        x <- factor(x = x, levels = lev)
    }
    return(x)
}

.addChrLen <- function(p, chrLen){
    dummy <- cbind(expand.grid(x = c(rep(0, length(chrLen$query$length)), 
                                     chrLen$query$length),
                               y = c(rep(0, length(chrLen$subject$length)), 
                                     chrLen$subject$length)),
                   expand.grid(query_chr = c(chrLen$query$names, 
                                             chrLen$query$names),
                               subject_chr = c(chrLen$subject$names, 
                                               chrLen$subject$names)))
    
    dummy$query_chr <- .sortLabels(x = dummy$query_chr, rev = TRUE)
    dummy$subject_chr <- .sortLabels(x = dummy$subject_chr, rev = FALSE)
    
    p <- p + 
        geom_point(data = dummy, mapping = aes(x = x, y = y), size = 0)
    return(p)
}


# ################################################################################
# # Deprecated
# ################################################################################
# # Define a function to remove ancient duplications
# rmAncientDup <- function(object){
#     stopifnot(inherits(x = object, "SynogDB"))
#     if(!inherits(x = object$sibeliaz$lcb, "LCB")){
#         stop("Invalid LCB info was given.")
#     }
#     
#     object$pairs$lcb_1to1 <- .order(object = object$pairs$lcb_1to1)
#     object$pairs$lcb_non1to1 <- .order(object = object$pairs$lcb_non1to1)
#     ad_candidate <- table(object$pairs$lcb_non1to1$block_id)
#     ad_candidate_id <- as.numeric(names(ad_candidate[ad_candidate > 2]))
#     ad_candidate_index <- object$pairs$lcb_non1to1$block_id %in% ad_candidate_id
#     df <- rbind(object$pairs$lcb_1to1, 
#                 object$pairs$lcb_non1to1[ad_candidate_index, ])
#     df <- .order(object = df)
#     df$index <- seq_along(df$block_id)
#     
#     ad_id <- sapply(ad_candidate_id, .evalAD, df = df)
#     out <- subset(df, subset = class == "non1to1")
#     out <- out[!out$index %in% unlist(ad_id), ]
#     object$pairs$lcb_non1to1 <- rbind(subset(out, select = -index), 
#                                       object$pairs$lcb_non1to1[!ad_candidate_index, ])
#     return(object)
# }
# 
# .order <- function(object, by = "query"){
#     if(by == "query"){
#         q_order <- order(object$query_chr, 
#                          object$query_start)
#         object <- object[q_order, ]
#     } else {
#         s_order <- order(object$subject_chr, 
#                          object$subject_start)
#         object <- object[s_order, ]
#     }
#     return(object)
# }
# 
# .evalAD <- function(df, id){
#     tmp <- subset(df, subset = block_id == id)
#     q_id <- paste(tmp$query_chr, 
#                   tmp$query_start,
#                   tmp$query_end, 
#                   sep = "_")
#     s_id <- paste(tmp$subject_chr, 
#                   tmp$subject_start,
#                   tmp$subject_end, 
#                   sep = "_")
#     
#     ad_list <- NULL
#     for(i in seq_along(tmp$block_id)){
#         s_id_i <- s_id[i]
#         s_id_hits <- s_id %in% s_id_i
#         s_id_hits[i] <- FALSE
#         s_id_hits <- which(s_id_hits)
#         q_id_hits <- q_id[s_id_hits]
#         
#         valid_route <- lapply(seq_along(q_id_hits),
#                               function(j){
#                                   q_tmp <- q_id[-i]
#                                   s_tmp <- s_id[-i]
#                                   s_tmp_hits <- s_tmp[q_tmp %in% q_id[i]]
#                                   s_end <- s_id[q_id %in% q_id_hits[j]]
#                                   hits <- s_end %in% s_tmp_hits
#                                   s_end_i <- which(q_id %in% q_id_hits[j])[hits]
#                                   out <- expand.grid(i,
#                                                      s_id_hits[j], 
#                                                      s_end_i,
#                                                      i)
#                                   return(out)
#                               })
#         ad_list <- rbind(ad_list, do.call("rbind", valid_route))
#     }
#     ad_id <- tmp$index[unique(unlist(ad_list))]
#     return(ad_id)
# }
# # Deprecated
# ################################################################################

################################################################################
# Define a function to execute BLAST search 
rbh <- function(object,
                db1 = NULL,
                db2 = NULL,
                n_threads = 1, 
                n_batch = NULL, 
                makedb = TRUE,
                max_target_seqs = 100000,
                pident = 90,
                qcovs = 0, 
                evalue = 1e-10
){
    stopifnot(inherits(x = object, "SynogDB"))
    
    fa1 <- .makeFASTA(cds_fn = object$query_cds)
    fa2 <- .makeFASTA(cds_fn = object$subject_cds)
    
    if(makedb | {is.null(db1) | is.null(db2)}){
        makeblastdb(file = fa1, dbtype = "nucl")
        makeblastdb(file = fa2, dbtype = "nucl")
        db1 <- fa1
        db2 <- fa2
    }
    
    fa1 <- readDNAStringSet(fa1)
    fa2 <- readDNAStringSet(fa2)
    db1 <- blast(db = db1, type = "blastn")
    db2 <- blast(db = db2, type = "blastn")
    
    options(scipen = 10^6)
    blast_out1 <- .blast_search(fa = fa1,
                                db = db2,
                                n_threads = n_threads, 
                                n_batch = n_batch,
                                max_target_seqs = max_target_seqs,
                                evalue = 1e-3,
                                task = "-task blastn")
    blast_out2 <- .blast_search(fa = fa2,
                                db = db1,
                                n_threads = n_threads, 
                                n_batch = n_batch,
                                max_target_seqs = max_target_seqs,
                                evalue = 1e-3,
                                task = "-task blastn")
    
    .h5creategroup(object$h5,"blast")
    .h5overwrite(obj = blast_out1, file = object$h5, "blast/blast_q2s")
    .h5overwrite(obj = blast_out2, file = object$h5, "blast/blast_s2q")
    
    # RBBH
    blast_best1 <- .blast_filter(blast_out = blast_out1, 
                                 pident = 0,
                                 qcovs = 0,
                                 evalue = 1e-3,
                                 best = TRUE)
    
    blast_best2 <- .blast_filter(blast_out = blast_out2, 
                                 pident = 0,
                                 qcovs = 0,
                                 evalue = 1e-3,
                                 best = TRUE)
    
    # RBH
    blast_out1 <- .blast_filter(blast_out = blast_out1, 
                                pident = pident,
                                qcovs = qcovs,
                                evalue = evalue,
                                best = FALSE)
    
    blast_out2 <- .blast_filter(blast_out = blast_out2, 
                                pident = pident,
                                qcovs = qcovs,
                                evalue = evalue,
                                best = FALSE)
    
    rbbh_out <- .find_reciprocal(df1 = blast_best1, df2 = blast_best2)
    rbh_out <- .find_reciprocal(df1 = blast_out1, df2 = blast_out2)
    options(scipen = 0)
    .h5overwrite(obj = rbbh_out, file = object$h5, "blast/rbbh")
    .h5overwrite(obj = rbh_out, file = object$h5, "blast/rbh")
}

.makeFASTA <- function(cds_fn){
    if(length(cds_fn) != 1){
        for(i in seq_along(cds_fn)){
            if(i == 1){
                cds <- readDNAStringSet(cds_fn[i])
            } else {
                cds <- c(cds, readDNAStringSet(cds_fn[i]))
            }
        }
        cds_fn <- sub("_cds.fa", "_merge_cds.fa", cds_fn[i])
        writeXStringSet(cds, cds_fn)
    }
    return(cds_fn)
}

.blast_search <- function(fa,
                          db, 
                          n_threads,
                          n_batch, 
                          max_target_seqs,
                          evalue,
                          task){
    blast_args <- paste(task,
                        "-best_hit_overhang 0.1",
                        "-best_hit_score_edge 0.1",
                        paste("-max_target_seqs",
                              max_target_seqs), 
                        paste("-evalue", evalue),
                        paste("-num_threads", n_threads))
    
    if(is.null(n_batch)){
        out <- predict(db, fa, 
                       silent = TRUE, 
                       BLAST_args = blast_args,
                       custom_format = "qseqid sseqid pident evalue qcovs")
        
    } else {
        batch <- split(seq_along(fa), cut(seq_along(fa), n_batch))
        out <- NULL
        for(i in batch){
            tmp <- predict(db, fa[i], 
                           silent = TRUE, 
                           BLAST_args = blast_args,
                           custom_format = "qseqid sseqid pident evalue qcovs")
            out <- rbind(out, tmp)
        }
    }
    return(out)
}

.blast_filter <- function(blast_out, 
                          pident = 0, 
                          qcovs = 0,
                          evalue = Inf, 
                          best){
    if(best){
        n_hit <- unlist(tapply(blast_out$qseqid, blast_out$qseqid, length))
        single_hit <- blast_out$qseqid %in% names(n_hit)[n_hit == 1]
        single_hit <- blast_out[single_hit, ]
        mult_hit <- blast_out$qseqid %in% names(n_hit)[n_hit > 1]
        mult_hit <- blast_out[mult_hit, ]
        
        filter <- tapply(mult_hit$evalue, mult_hit$qseqid, min)
        hit <- match(mult_hit$qseqid, names(filter))
        filter <- filter[hit]
        mult_hit <- subset(mult_hit, subset = evalue == filter)
        
        qcov_pident <- mult_hit$pident * mult_hit$qcovs * 1e-2
        
        filter <- tapply(qcov_pident, mult_hit$qseqid, max)
        hit <- match(mult_hit$qseqid, names(filter))
        filter <- filter[hit]
        mult_hit <- subset(mult_hit, subset = qcov_pident == filter)
        
        blast_out <- rbind(single_hit, mult_hit)
        
    } else {
        filter <- blast_out$pident >= pident & 
            blast_out$qcovs >= qcovs & 
            blast_out$evalue <= evalue
        blast_out <- subset(blast_out, subset = filter)
    }
    return(blast_out)
}

.find_reciprocal <- function(df1, df2){
    id1 <- paste(df1$qseqid, df1$sseqid, sep = "_")
    id2 <- paste(df2$sseqid, df2$qseqid, sep = "_")
    rhit <- id1 %in% id2
    out <- subset(df1, subset = rhit, select = c(qseqid, sseqid))
    out <- unique(out)
    
    out_id <- paste(out$qseqid, out$sseqid, sep = "_")
    out$q2s_pident <- df1$pident[match(out_id, id1)]
    out$q2s_qcovs <- df1$qcovs[match(out_id, id1)]
    out$q2s_evalue <- df1$evalue[match(out_id, id1)]
    out$s2q_pident <- df2$pident[match(out_id, id2)]
    out$s2q_qcovs <- df2$qcovs[match(out_id, id2)]
    out$s2q_evalue <- df2$evalue[match(out_id, id2)]
    
    return(out)
}

################################################################################
# Define a function to filter ortholog pairs based on LCBs
anchorOrtho <- function(object, 
                        non1to1 = TRUE){
    stopifnot(inherits(x = object, "SynogDB"))
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    if(!H5Lexists(h5, "sibeliaz/lcb_pairs")){
        stop("Run getLCBpairs to obtain LCB pair info.")
    }
    if(!H5Lexists(h5, "blast/rbbh")){
        stop("Run rbh with `best = TRUE` to obtain RBBH info.")
    }
    
    obj <- .makeObject(object = object, h5 = h5)
    
    obj <- .orthoIn1to1lcb(obj = obj)
    obj <- .orthoInNon1to1lcb(obj = obj)
    obj <- .orthoIn1to1cbi(obj = obj)
    
    .h5overwrite(obj = obj$out, file = object$h5, "anchor")
}

.makeObject <- function(object, h5){
    pairs <- list(lcb_1to1 = .order(df = h5$sibeliaz$lcb_pairs$lcb_1to1),
                  lcb_non_1to1 = .order(df = h5$sibeliaz$lcb_pairs$lcb_non_1to1))
    
    gr <- .df2gr(pairs = pairs)
    
    gap_gr <- list(query = .gapGR(gr = gr$query_1to1, 
                                  chrLen = object$genome$query),
                   subject = .gapGR(gr = gr$subject_1to1, 
                                    chrLen = object$genome$subject))
    
    query_gff <- .importAllGFF(object$query_gff)
    subject_gff <- .importAllGFF(object$subject_gff)
    gff <- .subsetGFF(query_gff = query_gff, subject_gff = subject_gff)
    
    rbbh <- h5$blast$rbbh
    rbbh$index <- seq_along(rbbh$qseqid)
    
    out <- list(gr = gr,
                rbbh = rbbh, 
                gap_gr = gap_gr, 
                gff = gff)
    return(out)
}

.order <- function(df, by = "query"){
    if(by == "query"){
        q_order <- order(df$query_chr, df$query_start)
        df <- df[q_order, ]
    } else {
        s_order <- order(df$subject_chr,
                         df$subject_start)
        df <- df[s_order, ]
    }
    return(df)
}

.df2gr <- function(pairs){
    gr <- list(query_1to1 = .makeGRanges(df = pairs$lcb_1to1, 
                                         genome = "query"),
               subject_1to1 = .makeGRanges(df = pairs$lcb_1to1, 
                                           genome = "subject"),
               query_non1to1 = .makeGRanges(df = pairs$lcb_non_1to1, 
                                            genome = "query"),
               subject_non1to1 = .makeGRanges(df = pairs$lcb_non_1to1, 
                                              genome = "subject"))
    return(gr)
}

.makeGRanges <- function(df, genome){
    if(genome == "query"){
        minus <- df$query_start > df$query_end
        tmp <- df$query_start[minus]
        df$query_start[minus] <- df$query_end[minus]
        df$query_end[minus] <- tmp
        gr <- GRanges(seqnames = df$query_chr,
                      ranges = IRanges(start = df$query_start, 
                                       end = df$query_end),
                      strand = "*")
        
    } else {
        minus <- df$subject_start > df$subject_end
        tmp <- df$subject_start[minus]
        df$subject_start[minus] <- df$subject_end[minus]
        df$subject_end[minus] <- tmp
        gr <- GRanges(seqnames = df$subject_chr,
                      ranges = IRanges(start = df$subject_start, 
                                       end = df$subject_end),
                      strand = "*")
    }
    return(gr)
}

.subsetGFF <- function(query_gff, subject_gff){
    query_gff <- query_gff[query_gff$type %in% c("transcript", "mRNA")]
    subject_gff <- subject_gff[subject_gff$type %in% c("transcript", "mRNA")]
    
    query_gff <- .orderGFF(gff = query_gff)
    subject_gff <- .orderGFF(gff = subject_gff)
    return(list(query = query_gff, subject = subject_gff))
}

.orderGFF <- function(gff){
    gff_order <- order(as.character(seqnames(gff)), start(gff))
    gff <- gff[gff_order]
    return(gff)
}

.gapGR <- function(gr, chrLen){
    n_gr <- length(gr)
    gap_gr <- data.frame(chr_start = as.character(seqnames(gr[-n_gr])),
                         chr_end = as.character(seqnames(gr[-1])),
                         start = end(gr[-n_gr]) + 1,
                         end = start(gr[-1]) - 1,
                         start_block = seq_along(gr)[-n_gr], 
                         end_block = seq_along(gr)[-1])
    
    gap_gr <- rbind(subset(gap_gr, subset = chr_start == chr_end), 
                    .gapGrStartEdge(gr = gr, gap_gr = gap_gr),
                    .gapGrEndEdge(gr = gr,
                                  gap_gr = gap_gr, 
                                  chrLen = chrLen,
                                  n_gr = n_gr))
    
    gap_gr <- .flipMinusStrand(df = gap_gr)
    
    gap_gr <- GRanges(seqnames = gap_gr$chr_start, 
                      ranges = IRanges(start = gap_gr$start,
                                       end = gap_gr$end), 
                      start_block = gap_gr$start_block,
                      end_block = gap_gr$end_block)
    
    gap_gr <- gap_gr[width(gap_gr) > 1]
    return(gap_gr)
}

.gapGrStartEdge <- function(gr, gap_gr){
    out <- subset(gap_gr, subset = chr_start != chr_end)
    out$chr_start <- out$chr_end
    out$start <- 1
    out$start_block <- NA
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

.gapGrEndEdge <- function(gr, gap_gr, chrLen, n_gr){
    out <- subset(gap_gr, subset = chr_start != chr_end)
    out$chr_end <- out$chr_start
    out$end <- chrLen$length[match(out$chr_start, chrLen$names)]
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

.flipMinusStrand <- function(df){
    minus <- df$start > df$end
    tmp <- df$start[minus]
    df$start[minus] <- df$end[minus]
    df$end[minus] <- tmp
    return(df)
}

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

.orthoInBlock <- function(gr_query, 
                          gr_subject, 
                          rbbh,
                          gff_query,
                          gff_subject,
                          interval = FALSE){
    q_g2b <-  .gene2block(gr = gr_query,
                          rbbh_id = rbbh$qseqid,
                          gff = gff_query)
    
    s_g2b <-  .gene2block(gr = gr_subject,
                          rbbh_id = rbbh$sseqid,
                          gff = gff_subject)
    if(interval){
        q_g2b <- .getInterval(g2b = q_g2b, gr = gr_query)
        s_g2b <- .getInterval(g2b = s_g2b, gr = gr_subject)
    }
    valid <- q_g2b == s_g2b
    out <- subset(rbbh, subset = valid)
    return(out)
}

.gene2block <- function(gr, rbbh_id, gff){
    ol <- findOverlaps(gff, gr)
    out <- rep(NA, length(gff))
    out[queryHits(ol)] <- subjectHits(ol)
    gene_index <- match(rbbh_id, gff$ID)
    return(out[gene_index])
}

.getInterval <- function(g2b, gr){
    not_na <- !is.na(g2b)
    start_int <- end_int <- g2b[not_na]
    start_int <- gr$start_block[start_int]
    end_int <- gr$end_block[end_int]
    out <- rep(NA, length(g2b))
    out[not_na] <- paste(start_int, end_int, sep = "_")
    return(out)
}

################################################################################
# Define a function to filter ortholog pairs based on gene synteny.
syntenyOrtho <- function(object,
                         omit_chr = "",
                         pident = 90,
                         evalue = 1e-50, 
                         qcovs = 50){
    stopifnot(inherits(x = object, "SynogDB"))
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    if(!H5Lexists(h5, "anchor")){
        stop("Run anchorOrtho to obtain ortholog anchors info.")
    }
    if(!H5Lexists(h5, "blast/rbh")){
        stop("Run rbh with `best = FALSE` to obtain RBH info.")
    }
    
    gff_ls <- .getGFFlist(object = object)
    anchor_ls <- .makeAnchorList(gff_ls = gff_ls, h5 = h5, omit_chr = omit_chr)
    rbh_dist <- .getDist(gff_ls = gff_ls, h5 = h5, anchor_ls = anchor_ls)
    tp_list <- read.csv(object$positive_list)
    fp_list <- read.csv(object$negative_list)
    metrics <- .getMetrics(rbh_dist = rbh_dist, 
                           tp_list = tp_list, 
                           fp_list = fp_list)
    rbh_dist$so_valid <- rbh_dist$dist <= metrics$so_threshold
    syn_og <- .makeSynogDF(rbh = rbh_dist, h5 = h5)
    syn_og <- .getMutualScore(ortho = syn_og)
    syn_og$class <- .classifySO(ortho = syn_og)
    
    rbbh_og <- .addRBBH(rbbh = h5$blast$rbbh, syn_og = syn_og)
    rbbh_og <- .filterOG(target = rbbh_og, ortho = rbh_dist, 
                         pident = 90, evalue = 1e-50, qcovs = 50)
    syn_og <- .mergeSynogDF(syn_og = syn_og, rbbh_og = rbbh_og, gff_ls = gff_ls)
    
    syn_og$class <- .classifySO(ortho = syn_og)
    orphan <- .getOrphan(gff_ls = gff_ls, syn_og = syn_og)
    out <- .makeOutput(syn_og = syn_og, orphan = orphan)
    .h5creategroup(object$h5,"synog_tx")
    .h5overwrite(obj = out$syn_og, file = object$h5, "synog_tx/orthopairs")
    .h5overwrite(obj = out$summary, file = object$h5, "synog_tx/summary")
    .h5overwrite(obj = out$orphan, file = object$h5, "synog_tx/orphan")
    .h5overwrite(obj = metrics$so_metrics, file = object$h5, "synog_tx/metrics")
    .h5overwrite(obj = metrics$so_threshold, file = object$h5, "synog_tx/so_threshold")
    invisible(metrics$metrics_plot)
}

.getGFFlist <- function(object){
    query_gff <- .importAllGFF(object$query_gff)
    subject_gff <- .importAllGFF(object$subject_gff)
    query_gff <- .orderGFF(gff = query_gff)
    subject_gff <- .orderGFF(gff = subject_gff)
    out <- list(query_gff =query_gff, subject_gff = subject_gff)
    return(out)
}

.makeAnchorList <- function(gff_ls, h5, omit_chr){
    anchor <- h5$anchor
    query_anchor <- .findNearestAnchor(gff = gff_ls$query_gff, 
                                       anchor = anchor$qseqid,
                                       omit_chr = omit_chr)
    subject_anchor <- .findNearestAnchor(gff = gff_ls$subject_gff, 
                                         anchor = anchor$sseqid,
                                         omit_chr = omit_chr)
    gene_id <- .setGeneID2Anchor(anchor = anchor, gff_ls = gff_ls)
    anchor <- cbind(anchor, gene_id)
    out <- .orgAnchor(anchor = anchor,
                      query_anchor = query_anchor,
                      subject_anchor = subject_anchor)
    return(out)
}

.findNearestAnchor <- function(gff, anchor, omit_chr){
    tx <- gff[gff$type %in% c("transcript", "mRNA")]
    anchor <- unlist(tx$Parent[match(anchor, tx$ID)])
    anchor <- na.omit(anchor)
    gff <- gff[gff$type == "gene"]
    gff <- gff[order(as.numeric(seqnames(gff)), start(gff))]
    gff <- data.frame(chr = as.character(seqnames(gff)), 
                      start = start(gff),
                      ID = gff$ID)
    gff <- gff[!grepl(omit_chr, gff$chr), ]
    gff$anchor <- FALSE
    gff$anchor[gff$ID %in% anchor] <- TRUE
    gff$nearest_anchor <- vapply(X = seq_along(gff$anchor), 
                                 FUN.VALUE = numeric(1), 
                                 FUN = .getNearestAnchorIndex, gff = gff)
    gff$location <- seq_along(gff$chr)
    return(gff)
}

.getNearestAnchorIndex <- function(index, gff){
    target <- gff$chr %in% gff$chr[index] & gff$anchor
    near <- which.min(abs(gff$start[target] - gff$start[index]))
    near <- which(target)[near]
    return(near)
}

.setGeneID2Anchor <- function(anchor, gff_ls){
    tx_i <- gff_ls$query_gff$type %in% c("transcript", "mRNA")
    tx <- gff_ls$query_gff[tx_i]
    hit <- match(anchor$qseqid, tx$ID)
    qgeneid <- unlist(tx$Parent[hit])
    qchr <- as.character(seqnames(tx[hit]))
    tx_i <- gff_ls$subject_gff$type %in% c("transcript", "mRNA")
    tx <- gff_ls$subject_gff[tx_i]
    hit <- match(anchor$sseqid, tx$ID)
    sgeneid <- unlist(tx$Parent[hit])
    schr <- as.character(seqnames(tx[hit]))
    out <- data.frame(qgeneid = qgeneid, sgeneid = sgeneid,
                      qchr = qchr, schr = schr)
    return(out)
}

.orgAnchor <- function(anchor, query_anchor, subject_anchor){
    anchor <- subset(anchor, select = qgeneid:schr)
    anchor <- unique(anchor)
    hit <- match(anchor$qgeneid, query_anchor$ID)
    anchor$query_anchor <- query_anchor$nearest_anchor[hit]
    tbl <- table(anchor$query_anchor)
    multiple <- anchor$query_anchor %in% names(tbl[tbl > 1])
    anchor$q2s_multiple <- FALSE
    anchor$q2s_multiple[multiple] <- TRUE
    hit <- match(anchor$sgeneid, subject_anchor$ID)
    anchor$subject_anchor <- subject_anchor$nearest_anchor[hit]
    tbl <- table(anchor$subject_anchor)
    multiple <- anchor$subject_anchor %in% names(tbl[tbl > 1])
    anchor$s2q_multiple <- FALSE
    anchor$s2q_multiple[multiple] <- TRUE
    out <- list(query = query_anchor, subject = subject_anchor, anchor = anchor)
    return(out)
}

.getDist <- function(gff_ls, h5, anchor_ls){
    rbh <- .extractRBH(h5 = h5, gff_ls = gff_ls)
    anchor2rbh <- .setAnchor2RBH(rbh = rbh, anchor_ls = anchor_ls)
    rbh <- cbind(rbh, anchor2rbh)
    q2s_dist <- .traceAnchorQ2S(rbh = rbh, anchor_ls = anchor_ls)
    s2q_dist <- .traceAnchorS2Q(rbh = rbh, anchor_ls = anchor_ls)
    rbh <- cbind(rbh, dist = rowMin(cbind(q2s_dist, s2q_dist)))
    return(rbh)
}

.extractRBH <- function(h5, gff_ls){
    out <- h5$blast$rbh
    q_hit <- match(out$qseqid, gff_ls$query_gff$ID)
    out$qgeneid <- gff_ls$query_gff$gene_id[q_hit]
    out$qchr <- as.character(seqnames(gff_ls$query_gff[q_hit]))
    s_hit <- match(out$sseqid, gff_ls$subject_gff$ID)
    out$sgeneid <- gff_ls$subject_gff$gene_id[s_hit]
    out$schr <- as.character(seqnames(gff_ls$subject_gff[s_hit]))
    out$pair_id <- paste(out$qgeneid, out$sgeneid, sep = "_")
    return(out)
}

.setAnchor2RBH <- function(rbh, anchor_ls){
    query_hit <- match(rbh$qgeneid, anchor_ls$query$ID)
    q2s_q_anchor <- anchor_ls$query$nearest_anchor[query_hit]
    s2q_q_location <- anchor_ls$query$location[query_hit]
    subject_hit <- match(rbh$sgeneid, anchor_ls$subject$ID)
    s2q_s_anchor <- anchor_ls$subject$nearest_anchor[subject_hit]
    q2s_s_location <- anchor_ls$subject$location[subject_hit]
    out <- data.frame(q2s_q_anchor = q2s_q_anchor, 
                      s2q_q_location = s2q_q_location,
                      s2q_s_anchor = s2q_s_anchor,
                      q2s_s_location = q2s_s_location)
    return(out)
}

.traceAnchorQ2S <- function(rbh, anchor_ls){
    q2s_q_anchor_s <- match(rbh$q2s_q_anchor, anchor_ls$anchor$query_anchor)
    q2s_q_anchor_s[is.na(rbh$q2s_q_anchor)] <- NA
    multiple_anchor <- anchor_ls$anchor$q2s_multiple[q2s_q_anchor_s]
    multiple_anchor[is.na(multiple_anchor)] <- FALSE
    single_projection <- anchor_ls$anchor$subject_anchor[q2s_q_anchor_s]
    single_projection[multiple_anchor] <- NA
    q2s_dist <- abs(rbh$q2s_s_location - single_projection)
    single_projection_chr <- anchor_ls$anchor$schr[q2s_q_anchor_s]
    q2s_dist_chr_unmatch <- rbh$schr != single_projection_chr
    q2s_dist[q2s_dist_chr_unmatch] <- Inf
    multiple_projection_rbh <- rbh[multiple_anchor, ]
    q2s_dist_multi <- vapply(X = seq_along(multiple_projection_rbh$q2s_q_anchor),
                             FUN.VALUE = numeric(1), FUN = .solveMultiProjectionQ2S, 
                             rbh = multiple_projection_rbh, anchor_ls = anchor_ls)
    q2s_dist[multiple_anchor] <- q2s_dist_multi
    q2s_dist[is.na(q2s_dist)] <- Inf
    return(q2s_dist)
}

.traceAnchorS2Q <- function(rbh, anchor_ls){
    s2q_s_anchor_q <- match(rbh$s2q_s_anchor, anchor_ls$anchor$subject_anchor)
    s2q_s_anchor_q[is.na(rbh$s2q_s_anchor)] <- NA
    multiple_anchor <- anchor_ls$anchor$s2q_multiple[s2q_s_anchor_q]
    multiple_anchor[is.na(multiple_anchor)] <- FALSE
    single_projection <- anchor_ls$anchor$query_anchor[s2q_s_anchor_q]
    single_projection[multiple_anchor] <- NA
    s2q_dist <- abs(rbh$s2q_q_location - single_projection)
    single_projection_chr <- anchor_ls$anchor$qchr[s2q_s_anchor_q]
    s2q_dist_chr_unmatch <- rbh$qchr != single_projection_chr
    s2q_dist[s2q_dist_chr_unmatch] <- Inf
    multiple_projection_rbh <- rbh[multiple_anchor, ]
    s2q_dist_multi <- vapply(X = seq_along(multiple_projection_rbh$s2q_s_anchor),
                             FUN.VALUE = numeric(1), FUN = .solveMultiProjectionS2Q, 
                             rbh = multiple_projection_rbh, anchor_ls = anchor_ls)
    s2q_dist[multiple_anchor] <- s2q_dist_multi
    s2q_dist[is.na(s2q_dist)] <- Inf
    return(s2q_dist)
}

.solveMultiProjectionQ2S <- function(index, rbh, anchor_ls){
    q_anchor <- rbh$q2s_q_anchor[index]
    hit <- which(anchor_ls$anchor$query_anchor == q_anchor)
    projection <- anchor_ls$anchor$subject_anchor[hit]
    q2s_dist <- abs(rbh$q2s_s_location[index] - projection)
    projection_chr <- anchor_ls$anchor$schr[hit]
    q2s_dist_chr_unmatch <- rbh$schr[index] != projection_chr
    q2s_dist[q2s_dist_chr_unmatch] <- Inf
    out <- min(q2s_dist)
    return(out)
}

.solveMultiProjectionS2Q <- function(index, rbh, anchor_ls){
    s_anchor <- rbh$s2q_s_anchor[index]
    hit <- which(anchor_ls$anchor$subject_anchor == s_anchor)
    projection <- anchor_ls$anchor$query_anchor[hit]
    s2q_dist <- abs(rbh$s2q_q_location[index] - projection)
    projection_chr <- anchor_ls$anchor$qchr[hit]
    s2q_dist_chr_unmatch <- rbh$qchr[index] != projection_chr
    s2q_dist[s2q_dist_chr_unmatch] <- Inf
    out <- min(s2q_dist)
    return(out)
}


.getMetrics <- function(rbh_dist, tp_list, fp_list){
    n_list <- unique(c(fp_list$query, fp_list$subject))
    n_hit <- rbh_dist$qgeneid %in% n_list | rbh_dist$sgeneid %in% n_list
    n_hit_index_dist <- rbh_dist$dist[n_hit]
    n_hit_index_dist <- tapply(n_hit_index_dist,
                               rbh_dist$qgeneid[n_hit],
                               min)
    n_tbl <- table(n_hit_index_dist)
    n_tbl <- n_tbl[names(n_tbl) != "Inf"]
    
    p_id <- paste(tp_list$query, tp_list$subject, sep = "_")
    p_id <- unique(p_id)
    p_hit <- which(rbh_dist$pair_id %in% p_id)
    p_hit_index_dist <- rbh_dist$dist[p_hit]
    p_hit_index_dist <- tapply(p_hit_index_dist,
                               rbh_dist$pair_id[p_hit],
                               min)
    p_tbl <- table(p_hit_index_dist)
    p_tbl <- p_tbl[names(p_tbl) != "Inf"]
    
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
    
    df <- rbind(data.frame(Metrics = "precision", 
                           Threshold = seq_along(precision),
                           Score = precision), 
                data.frame(Metrics = "recall", 
                           Threshold = seq_along(recall),
                           Score = recall),
                data.frame(Metrics = "specificity", 
                           Threshold = seq_along(specificity),
                           Score = specificity),
                data.frame(Metrics = "accuracy", 
                           Threshold = seq_along(accuracy),
                           Score = accuracy),
                data.frame(Metrics = "f_measure", 
                           Threshold = seq_along(f_measure),
                           Score = f_measure))
    p <- ggplot(data = df) +
        geom_line(aes(x = Threshold, y = Score, group = Metrics, color = Metrics)) +
        xlim(0, 100) +
        scale_color_manual(values = c("blue", 
                                      "darkgreen",
                                      "magenta", 
                                      "skyblue", 
                                      "darkorange"),
                           breaks = c("precision", 
                                      "recall",
                                      "specificity",
                                      "accuracy", 
                                      "f_measure"))
    out <- list()
    out$so_metrics <- data.frame(threshold = seq_along(precision),
                                 precision = precision,
                                 recall = recall,
                                 specificity = specificity,
                                 accuracy = accuracy,
                                 f_measure = f_measure)
    out$metrics_plot <- p
    out$so_threshold <- which.max(f_measure)
    return(out)
}

.makeSynogDF <- function(rbh, h5){
    out <- subset(rbh, subset = so_valid %in% TRUE, 
                  select = -c(so_valid, qchr, schr, q2s_q_anchor, 
                              s2q_q_location, s2q_s_anchor, q2s_s_location,
                              dist))
    out$syntenic <- TRUE
    out$rbbh <- FALSE
    out$pair_id <- paste(out$qseqid, out$sseqid, sep = "_")
    rbbh_id <- paste(h5$anchor$qseqid,
                     h5$anchor$sseqid, sep = "_")
    out$rbbh[out$pair_id %in% rbbh_id] <- TRUE
    return(out)
}

.mergeSynogDF <- function(syn_og, rbbh_og, gff_ls){
    if(nrow(rbbh_og) != 0){
        rbbh_og <- .setGeneIDsynog(rbbh_og = rbbh_og, gff_ls = gff_ls)
        rbbh_og$syntenic <- FALSE
        rbbh_og$rbbh <- TRUE
        rbbh_og$pair_id <- paste(rbbh_og$qseqid, rbbh_og$sseqid, sep = "_")
        rbbh_og <- .getMutualScore(ortho = rbbh_og)
        rbbh_og$class <- .classifySO(ortho = rbbh_og)
        syn_og <- rbind(syn_og, rbbh_og)
    }
    return(syn_og)
}

.makeOutput <- function(syn_og, orphan){
    out <- NULL
    out$syn_og <- syn_og
    out$syn_og$OG <- .numberingOG(ortho = syn_og)
    out$syn_og <- subset(out$syn_og, select = c(qseqid, sseqid, mutual_e,
                                                mutual_ci, syntenic, rbbh,
                                                class, OG, qgeneid, sgeneid))
    rownames(out$syn_og) <- NULL
    out$syn_og$syntenic <- as.character(out$syn_og$syntenic)
    out$syn_og$rbbh <- as.character(out$syn_og$rbbh)
    
    q_id <- tapply(out$syn_og$qseqid, out$syn_og$class, unique)
    n_q_id <- sapply(q_id, length)
    s_id <- tapply(out$syn_og$sseqid, out$syn_og$class, unique)
    n_s_id <- sapply(s_id, length)
    
    n_orphan <- lapply(orphan, nrow)
    
    q_tot <- sum(n_q_id) + n_orphan$query
    s_tot <- sum(n_s_id) + n_orphan$subject
    df <- data.frame(query = c(q_tot, sum(n_q_id), n_q_id, n_orphan$query), 
                     subject = c(s_tot, sum(n_s_id), n_s_id, n_orphan$subject))
    out$summary <- df
    out$orphan <- orphan
    rownames(out$orphan$query) <- NULL
    rownames(out$orphan$subject) <- NULL
    return(out)   
}

.classifySO <- function(ortho){
    tmp <- subset(ortho, select = c(qseqid, sseqid, pair_id))
    q_dup <- duplicated(tmp$qseqid)
    q_dup <- tmp$qseqid %in% tmp$qseqid[q_dup]
    s_dup <- duplicated(tmp$sseqid)
    s_dup <- tmp$sseqid %in% tmp$sseqid[s_dup]
    og_11 <- tmp[!(q_dup | s_dup), ]
    
    mult <- tmp[q_dup | s_dup, ]
    mult_q_dup <- mult$qseqid %in% mult$qseqid[duplicated(mult$qseqid)]
    mult_s_dup <- mult$sseqid %in% mult$sseqid[duplicated(mult$sseqid)]
    
    og_mm <- mult[mult_q_dup & mult_s_dup, ]
    og_mm_i <- mult$qseqid %in% og_mm$qseqid | mult$sseqid %in% og_mm$sseqid
    og_mm <- mult[og_mm_i, ]
    
    not_mm <- mult[!og_mm_i, ]
    not_mm_q_dup <- not_mm$qseqid %in% not_mm$qseqid[duplicated(not_mm$qseqid)]
    not_mm_s_dup <- not_mm$sseqid %in% not_mm$sseqid[duplicated(not_mm$sseqid)]
    og_1m <- not_mm[not_mm_q_dup, ]
    og_m1 <- not_mm[not_mm_s_dup, ]
    
    out <- og_11
    out$class <- "1to1"
    og_1m$class <-"1toM"
    out <- rbind(out, og_1m)
    og_m1$class <-"Mto1"
    out <- rbind(out, og_m1)
    og_mm$class <-"MtoM"
    out <- rbind(out, og_mm)
    
    hit <- match(ortho$pair_id, out$pair_id)
    return(out$class[hit])
}

.getOrphan <- function(gff_ls, syn_og){
    tx_i <- gff_ls$query_gff$type %in% c("transcript", "mRNA")
    tx <- gff_ls$query_gff[tx_i]
    hit <- tx$ID %in% syn_og$qseqid
    q_orphan <- data.frame(ID = tx$ID[!hit], gene_id = tx$gene_id[!hit])
    
    tx_i <- gff_ls$subject_gff$type %in% c("transcript", "mRNA")
    tx <- gff_ls$subject_gff[tx_i]
    hit <- tx$ID %in% syn_og$sseqid
    s_orphan <- data.frame(ID = tx$ID[!hit], gene_id = tx$gene_id[!hit])
    return(list(query = q_orphan, subject = s_orphan))
}

.getMutualScore <- function(ortho){
    q2s_covident <- ortho$q2s_pident * ortho$q2s_qcovs * 1e-2
    s2q_covident <- ortho$s2q_pident * ortho$s2q_qcovs * 1e-2
    mutual_e <- sqrt(ortho$q2s_evalue^2 + ortho$s2q_evalue^2)
    mutual_ci <- sqrt(q2s_covident^2 + s2q_covident^2)
    ortho$mutual_e <- mutual_e
    ortho$mutual_ci <- mutual_ci
    return(ortho)
}

# .find1to1 <- function(target){
#     og_target <- target[target$class == "MtoM", ]
#     score_index <- seq_along(og_target$mutual_e)
#     q2s <- tapply(score_index, og_target$qseqid, function(i){
#         min_e <- min(og_target$mutual_e[i])
#         min_e_i <- which(og_target$mutual_e[i] == min_e)
#         max_ci <- max(og_target$mutual_ci[i][min_e_i])
#         max_ci_i <- og_target$mutual_ci[i][min_e_i] == max_ci
#         best_i <- min_e_i[max_ci_i]
#         return(score_index[i][best_i])
#     })
#     q2s <- sapply(seq_along(q2s), function(i){
#         out <- q2s[[i]]
#         names(out) <- rep(names(q2s)[i], length(out))
#         return(out)
#     })
#     q2s <- unlist(q2s)
#     q2s <- data.frame(qseqid = names(q2s), sseqid = og_target$sseqid[q2s])
#     s2q <-  tapply(score_index, og_target$sseqid, function(i){
#         min_e <- min(og_target$mutual_e[i])
#         min_e_i <- which(og_target$mutual_e[i] == min_e)
#         max_ci <- max(og_target$mutual_ci[i][min_e_i])
#         max_ci_i <- og_target$mutual_ci[i][min_e_i] == max_ci
#         best_i <- min_e_i[max_ci_i]
#         return(score_index[i][best_i])
#     })
#     s2q <- sapply(seq_along(s2q), function(i){
#         out <- s2q[[i]]
#         names(out) <- rep(names(s2q)[i], length(out))
#         return(out)
#     })
#     s2q <- unlist(s2q)
#     s2q <- data.frame(qseqid = og_target$qseqid[s2q], sseqid = names(s2q))
#     add_1to1 <- unique(rbind(q2s, s2q))
#     add_1to1$pair_id <- paste(add_1to1$qseqid, add_1to1$sseqid, sep = "_")
#     
#     add_1to1_index <- which(og_target$pair_id %in% add_1to1$pair_id)
#     q_hit <- og_target$qseqid %in% og_target$qseqid[add_1to1_index]
#     s_hit <- og_target$sseqid %in% og_target$sseqid[add_1to1_index]
#     out <- og_target$pair_id[q_hit | s_hit]
#     out <- out[!out %in% add_1to1$pair_id]
#     return(out)
# }

.addRBBH <- function(rbbh, syn_og){
    candidate_rbbh <- !(rbbh$qseqid %in% syn_og$qseqid |
                            rbbh$sseqid %in% syn_og$sseqid)
    out <- rbbh[candidate_rbbh, ]
    return(out)
}

.filterOG <- function(target, ortho, pident, evalue, qcovs){
    ortho_id <- paste(ortho$qseqid, ortho$sseqid, sep = "_")
    target_id <- paste(target$qseqid, target$sseqid, sep = "_")
    target_metrics <- ortho[ortho_id %in% target_id, ]
    valid_og <- target_metrics$q2s_pident >= pident & 
        target_metrics$s2q_pident >= pident &
        target_metrics$q2s_qcovs >= qcovs &
        target_metrics$s2q_qcovs >= qcovs &
        target_metrics$q2s_evalue <= evalue &
        target_metrics$s2q_evalue <= evalue
    valid_og <- target_metrics[valid_og, ]
    valid_id <- paste(valid_og$qseqid, valid_og$sseqid, sep = "_")
    valid_og <- target[target_id %in% valid_id, ]
    return(valid_og)
}

.setGeneIDsynog <- function(rbbh_og, gff_ls){
    q_hit <- match(rbbh_og$qseqid, gff_ls$query_gff$ID)
    rbbh_og$qgeneid <- gff_ls$query_gff$gene_id[q_hit]
    s_hit <- match(rbbh_og$sseqid, gff_ls$subject_gff$ID)
    rbbh_og$sgeneid <- gff_ls$subject_gff$gene_id[s_hit]
    return(rbbh_og)
}

# .pickOmitOrtho <- function(syn_og, rbbh, ortho, omit, filt){
#     omit_ortho <- ortho[omit, ]
#     omit_ortho_id <- paste(omit_ortho$qseqid, omit_ortho$sseqid, sep = "_")
#     rbbh_id <- paste(rbbh$qseqid, rbbh$sseqid, sep = "_")
#     rbbh <- subset(rbbh,
#                    subset = rbbh_id %in% omit_ortho_id, 
#                    select = c(qseqid, sseqid))
#     rbbh_1to1 <- (rbbh = rbbh, syn_og = syn_og)
#     rbbh_1to1 <- .filterOG(ortho = ortho, target = rbbh_1to1, filt = filt)
#     out <- rbind(syn_og, rbbh_1to1)
#     return(out)
# }

.numberingOG <- function(ortho){
    tmp <- subset(ortho, select = c(qseqid, sseqid, pair_id))
    q_dup <- duplicated(tmp$qseqid)
    q_dup <- tmp$qseqid %in% tmp$qseqid[q_dup]
    s_dup <- duplicated(tmp$sseqid)
    s_dup <- tmp$sseqid %in% tmp$sseqid[s_dup]
    og_11 <- tmp[!(q_dup | s_dup), ]
    og_11 <- og_11[order(og_11$qseqid), ]
    og_11$OG <- seq_along(og_11$qseqid)
    
    mult <- tmp[q_dup | s_dup, ]
    mult_q_dup <- mult$qseqid %in% mult$qseqid[duplicated(mult$qseqid)]
    mult_s_dup <- mult$sseqid %in% mult$sseqid[duplicated(mult$sseqid)]
    
    og_mm <- mult[mult_q_dup & mult_s_dup, ]
    og_mm_i <- mult$qseqid %in% og_mm$qseqid | mult$sseqid %in% og_mm$sseqid
    og_mm <- mult[og_mm_i, ]
    
    not_mm <- mult[!og_mm_i, ]
    not_mm_q_dup <- not_mm$qseqid %in% not_mm$qseqid[duplicated(not_mm$qseqid)]
    not_mm_s_dup <- not_mm$sseqid %in% not_mm$sseqid[duplicated(not_mm$sseqid)]
    og_1m <- not_mm[not_mm_q_dup, ]
    og_m1 <- not_mm[not_mm_s_dup, ]
    
    og_1m <- og_1m[order(og_1m$qseqid), ]
    og_1m_index <- data.frame(id = unique(og_1m$qseqid), 
                              OG = seq_along(unique(og_1m$qseqid)))
    og_1m$OG <- og_1m_index$OG[match(og_1m$qseqid, og_1m_index$id)]
    
    og_m1 <- og_m1[order(og_m1$qseqid), ]
    og_m1_index <- data.frame(id = unique(og_m1$sseqid), 
                              OG = seq_along(unique(og_m1$sseqid)))
    og_m1$OG <- og_m1_index$OG[match(og_m1$sseqid, og_m1_index$id)]
    
    og_mm <- og_mm[order(og_mm$qseqid), ]
    og_mm$OG <- seq_along(og_mm$qseqid)
    og_mm_q_min <- tapply(og_mm$OG, og_mm$qseqid, min)
    og_mm$q_min <- og_mm_q_min[match(og_mm$qseqid, names(og_mm_q_min))]
    og_mm_s_min <- tapply(og_mm$q_min, og_mm$sseqid, min)
    og_mm$s_min <- og_mm_s_min[match(og_mm$sseqid, names(og_mm_s_min))]
    og_mm_q_min <- tapply(og_mm$s_min, og_mm$qseqid, min)
    og_mm$OG <- og_mm_q_min[match(og_mm$qseqid, names(og_mm_q_min))]
    og_mm <- subset(og_mm, select = -c(q_min, s_min))
    
    out <- og_11
    out$class <- "1to1"
    og_1m$OG <- og_1m$OG + max(out$OG)
    og_1m$class <-"1toM"
    out <- rbind(out, og_1m)
    og_m1$OG <- og_m1$OG + max(out$OG)
    og_m1$class <-"Mto1"
    out <- rbind(out, og_m1)
    og_mm$OG <- og_mm$OG + max(out$OG)
    og_mm$class <-"MtoM"
    out <- rbind(out, og_mm)
    
    hit <- match(ortho$pair_id, out$pair_id)
    return(out$OG[hit])
}

# ################################################################################
# # Define a function to pick up reliable ortholog pairs
# realibleOrtho <- function(rbbh,
#                           syn_ortho,
#                           add_rbh = TRUE,
#                           pident = 95,
#                           qcovs = 95){
#     
#     # ortho$q2s_ci <- ortho$q2s_pident * ortho$q2s_qcovs * 1e-2
#     # ortho$s2q_ci <- ortho$s2q_pident * ortho$s2q_qcovs * 1e-2
#     # ortho$score <- sqrt(ortho$q2s_ci^2 + ortho$s2q_ci^2)
#     # uni_score <- sort(unique(ortho$score), decreasing = TRUE)
#     # genome <- attributes(syn_ortho)$genome
#     # n <- syn_ortho$n[[genome$query]]
#     # n <- n[n %in% syn_ortho$orphan$query]
#     # n_hit <- which(ortho$qseqid %in% n)
#     # n_score <- ortho$score[n_hit]
#     # n_score <- tapply(n_score, ortho$qseqid[n_hit], min)
#     # 
#     # id <- paste(ortho$qseqid, ortho$sseqid, sep = "_")
#     # p <- syn_ortho$p
#     # p_id <- paste(p[[genome$query]], p[[genome$subject]], sep = "_")
#     # syn_id <- paste(syn_ortho$syn_og$qseqid, syn_ortho$syn_og$sseqid, sep = "_")
#     # p_id <- p_id[!p_id %in% syn_id]
#     # p_hit <- which(id %in% p_id)
#     # p_score <- ortho$score[p_hit]
#     # 
#     # fp <- sapply(uni_score, function(x){ sum(n_score >= x) })
#     # tp <- sapply(uni_score, function(x){ sum(p_score >= x) })
#     # fn <- length(p_id) - tp
#     # tn <- length(n) - fp
#     # precision <- tp / (tp + fp)
#     # recall <- tp / (tp + fn)
#     # f_measure <- (2 * recall * precision) / (recall + precision)
#     # f_measure[is.na(f_measure)] <- 0
#     # out <- ortho[ortho$score >= uni_score[which.max(f_measure)], ]
#     
#     orphan <- syn_ortho$orphan
#     q_id <- c(syn_ortho$syn_og$qseqid, rbbh_1to1$qseqid)
#     s_id <- c(syn_ortho$syn_og$sseqid, rbbh_1to1$sseqid)
#     orphan$query <- orphan$query[!orphan$query %in% q_id]
#     orphan$subject <- orphan$subject[!orphan$subject %in% s_id]
# 
#     if(add_rbh){
#         add_og <- .addRBH(ortho = ortho, orphan = orphan, te = syn_ortho$te)
#         add_og <- .filterOG(ortho = ortho, target = add_og, filt = filt)
#         syn_og <- rbind(syn_ortho$syn_og, rbbh_1to1, add_og)
# 
#     } else {
#         syn_og <- rbind(syn_ortho$syn_og, rbbh_1to1)
#     }
# 
#     syn_og <- .classifySO(ortho = syn_og)
#     syn_og_id <- paste(syn_og$qseqid, syn_og$sseqid, sep = "_")
#     rbbh_1to1_id <- paste(rbbh_1to1$qseqid, rbbh_1to1$sseqid, sep = "_")
#     rbbh_1to1_i <- syn_og_id %in% rbbh_1to1_id
#     syn_og$class[rbbh_1to1_i] <- paste0("rbbh_", syn_og$class[rbbh_1to1_i])
# 
#     if(add_rbh){
#         add_og_id <- paste(add_og$qseqid, add_og$sseqid, sep = "_")
#         add_og_id_i <- syn_og_id %in% add_og_id
#         syn_og$class[add_og_id_i] <- paste0("nonsyn_", syn_og$class[add_og_id_i])
#     }
# 
#     orphan <- syn_ortho$orphan
#     orphan$query <- orphan$query[!orphan$query %in% syn_og$qseqid]
#     orphan$subject <- orphan$subject[!orphan$subject %in% syn_og$sseqid]
# 
#     out <- .makeOutput(syn_og = syn_og, orphan = orphan, te = syn_ortho$te)
# 
#     return(out)
# }
# 
# .addRBH <- function(ortho, orphan, te){
#     q_hit_te <- ortho$qseqid %in% te$query
#     s_hit_te <- ortho$sseqid %in% te$subject
#     ortho <- subset(ortho, subset = !q_hit_te & !s_hit_te)
#     
#     q_orphan_hit <- ortho$qseqid %in% orphan$query
#     s_orphan_hit <- ortho$sseqid %in% orphan$subject
#     out <- ortho[q_orphan_hit & s_orphan_hit, ]
#     out$class <- "add"
#     
#     out <- .find1to1(ortho = ortho, target = out, class = "add")
#     out <- .classifySO(ortho = out)
#     out$class <- paste0("nonsyn_", out$class)
#     
#     return(out)
# }

geneOrtho <- function(object){
    stopifnot(inherits(x = object, "SynogDB"))
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    if(!H5Lexists(h5, "synog_tx/orthopairs")){
        stop("Run syntenyOrtho to obtain ortholog anchors info.")
    }
    ortho <- h5$synog_tx$orthopairs
    ortho <- .findBestPair(object = object, ortho = ortho)
    ortho$pair_id <- paste(ortho$qgeneid, ortho$sgeneid, sep = "_")
    ortho$genewise_class <- .classifySOgene(ortho = ortho)
    omit_id <- .find1to1gene(target = ortho)
    ortho <- subset(ortho, subset = !ortho$pair_id %in% omit_id)
    class_og <- .numberingOGgene(ortho = ortho)
    ortho$genewise_class <- class_og$class
    ortho$genewise_OG <- class_og$OG
    ortho <- ortho[order(ortho$genewise_OG, ortho$qgeneid, ortho$sgeneid), ]
    rownames(ortho) <- NULL
    orphan <- .getOrphanGene(ortho = ortho, tx_orphan = h5$synog_tx$orphan)
    gene_summary <- .makeSummary(ortho = ortho, orphan = orphan)
    ortho$pair_id <- NULL
    .h5creategroup(object$h5,"synog_gene")
    .h5overwrite(obj = ortho, file = object$h5, "synog_gene/orthopairs")
    .h5overwrite(obj = gene_summary, file = object$h5, "synog_gene/summary")
    .h5overwrite(obj = orphan, file = object$h5, "synog_gene/orphan")
}

.classifySOgene <- function(ortho){
    tmp <- subset(ortho, select = c(qgeneid, sgeneid, pair_id))
    tmp <- unique(tmp)
    q_dup <- duplicated(tmp$qgeneid)
    q_dup <- tmp$qgeneid %in% tmp$qgeneid[q_dup]
    s_dup <- duplicated(tmp$sgeneid)
    s_dup <- tmp$sgeneid %in% tmp$sgeneid[s_dup]
    og_11 <- tmp[!(q_dup | s_dup), ]
    
    mult <- tmp[q_dup | s_dup, ]
    mult_q_dup <- mult$qgeneid %in% mult$qgeneid[duplicated(mult$qgeneid)]
    mult_s_dup <- mult$sgeneid %in% mult$sgeneid[duplicated(mult$sgeneid)]
    
    og_mm <- mult[mult_q_dup & mult_s_dup, ]
    og_mm_i <- mult$qgeneid %in% og_mm$qgeneid | mult$sgeneid %in% og_mm$sgeneid
    og_mm <- mult[og_mm_i, ]
    
    not_mm <- mult[!og_mm_i, ]
    not_mm_q_dup <- not_mm$qgeneid %in% not_mm$qgeneid[duplicated(not_mm$qgeneid)]
    not_mm_s_dup <- not_mm$sgeneid %in% not_mm$sgeneid[duplicated(not_mm$sgeneid)]
    og_1m <- not_mm[not_mm_q_dup, ]
    og_m1 <- not_mm[not_mm_s_dup, ]
    
    out <- og_11
    out$class <- "1to1"
    og_1m$class <-"1toM"
    out <- rbind(out, og_1m)
    og_m1$class <-"Mto1"
    out <- rbind(out, og_m1)
    og_mm$class <-"MtoM"
    out <- rbind(out, og_mm)
    
    hit <- match(ortho$pair_id, out$pair_id)
    
    return(out$class[hit])
}

.find1to1gene <- function(target){
    og_target <- target[target$genewise_class == "MtoM", ]
    score_index <- seq_along(og_target$qgeneid)
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
    
    add_1to1 <- unique(rbind(q2s, s2q))
    add_1to1$pair_id <- paste(add_1to1$qgeneid, add_1to1$sgeneid, sep = "_")
    
    add_1to1_index <- which(og_target$pair_id %in% add_1to1$pair_id)
    q_hit <- og_target$qgeneid %in% og_target$qgeneid[add_1to1_index]
    s_hit <- og_target$sgeneid %in% og_target$sgeneid[add_1to1_index]
    out <- og_target$pair_id[q_hit | s_hit]
    out <- out[!out %in% add_1to1$pair_id]
    return(out)
}

.numberingOGgene <- function(ortho){
    tmp <- subset(ortho, select = c(qgeneid, sgeneid, pair_id))
    tmp <- unique(tmp)
    q_dup <- duplicated(tmp$qgeneid)
    q_dup <- tmp$qgeneid %in% tmp$qgeneid[q_dup]
    s_dup <- duplicated(tmp$sgeneid)
    s_dup <- tmp$sgeneid %in% tmp$sgeneid[s_dup]
    og_11 <- tmp[!(q_dup | s_dup), ]
    og_11$OG <- seq_along(og_11$qgeneid)
    
    mult <- tmp[q_dup | s_dup, ]
    mult_q_dup <- mult$qgeneid %in% mult$qgeneid[duplicated(mult$qgeneid)]
    mult_s_dup <- mult$sgeneid %in% mult$sgeneid[duplicated(mult$sgeneid)]
    
    og_mm <- mult[mult_q_dup & mult_s_dup, ]
    og_mm_i <- mult$qgeneid %in% og_mm$qgeneid | mult$sgeneid %in% og_mm$sgeneid
    og_mm <- mult[og_mm_i, ]
    
    not_mm <- mult[!og_mm_i, ]
    not_mm_q_dup <- not_mm$qgeneid %in% not_mm$qgeneid[duplicated(not_mm$qgeneid)]
    not_mm_s_dup <- not_mm$sgeneid %in% not_mm$sgeneid[duplicated(not_mm$sgeneid)]
    og_1m <- not_mm[not_mm_q_dup, ]
    og_m1 <- not_mm[not_mm_s_dup, ]
    
    og_1m <- og_1m[order(og_1m$qgeneid), ]
    og_1m_index <- data.frame(id = unique(og_1m$qgeneid), 
                              OG = seq_along(unique(og_1m$qgeneid)))
    og_1m$OG <- og_1m_index$OG[match(og_1m$qgeneid, og_1m_index$id)]
    
    og_m1 <- og_m1[order(og_m1$qgeneid), ]
    og_m1_index <- data.frame(id = unique(og_m1$sgeneid), 
                              OG = seq_along(unique(og_m1$sgeneid)))
    og_m1$OG <- og_m1_index$OG[match(og_m1$sgeneid, og_m1_index$id)]
    
    og_mm <- og_mm[order(og_mm$qgeneid), ]
    og_mm$OG <- seq_along(og_mm$qgeneid)
    og_mm_q_min <- tapply(og_mm$OG, og_mm$qgeneid, min)
    og_mm$q_min <- og_mm_q_min[match(og_mm$qgeneid, names(og_mm_q_min))]
    og_mm_s_min <- tapply(og_mm$q_min, og_mm$sgeneid, min)
    og_mm$s_min <- og_mm_s_min[match(og_mm$sgeneid, names(og_mm_s_min))]
    og_mm_q_min <- tapply(og_mm$s_min, og_mm$qgeneid, min)
    og_mm$OG <- og_mm_q_min[match(og_mm$qgeneid, names(og_mm_q_min))]
    og_mm <- subset(og_mm, select = -c(q_min, s_min))
    
    out <- og_11
    out$class <- "1to1"
    og_1m$class <-"1toM"
    og_1m$OG <- og_1m$OG + max(out$OG)
    out <- rbind(out, og_1m)
    og_m1$class <-"Mto1"
    og_m1$OG <- og_m1$OG + max(out$OG)
    out <- rbind(out, og_m1)
    og_mm$class <-"MtoM"
    og_mm$OG <- og_mm$OG + max(out$OG)
    out <- rbind(out, og_mm)
    
    hit <- match(ortho$pair_id, out$pair_id)
    
    return(list(class = out$class[hit], OG = out$OG[hit]))
}

.getOrphanGene <- function(ortho, tx_orphan){
    q_genes <- unique(tx_orphan$query$gene_id)
    q_orphan <- q_genes[!q_genes %in% ortho$qgeneid]
    q_orphan <- data.frame(qgeneid = q_orphan)
    s_genes <- unique(tx_orphan$subject$gene_id)
    s_orphan <- s_genes[!s_genes %in% ortho$sgeneid]
    s_orphan <- data.frame(sgeneid = s_orphan)
    return(list(query = q_orphan, subject = s_orphan))
}

.makeSummary <- function(ortho, orphan){
    q_id <- tapply(ortho$qgeneid, ortho$genewise_class, unique)
    n_q_id <- sapply(q_id, length)
    s_id <- tapply(ortho$sgeneid, ortho$genewise_class, unique)
    n_s_id <- sapply(s_id, length)
    
    n_orphan <- lapply(orphan, nrow)
    
    q_tot <- sum(n_q_id) + n_orphan$query
    s_tot <- sum(n_s_id) + n_orphan$subject
    out <- data.frame(query = c(q_tot, sum(n_q_id), n_q_id, n_orphan$query), 
                      subject = c(s_tot, sum(n_s_id), n_s_id, n_orphan$subject))
    return(out)   
}

.findBestPair <- function(object, ortho){
    pair_id <- paste(ortho$qgeneid, ortho$sgeneid, sep = "_")
    
    q_gff <- .importAllGFF(object$query_gff)
    s_gff <- .importAllGFF(object$subject_gff)
    q_gff <- q_gff[q_gff$type %in% c("mRNA", "transcript")]
    s_gff <- s_gff[s_gff$type %in% c("mRNA", "transcript")]
    q_len <- width(q_gff[match(ortho$qseqid, q_gff$ID)])
    s_len <- width(s_gff[match(ortho$sseqid, s_gff$ID)])
    
    best_pair <- tapply(seq_along(pair_id), pair_id, function(index){
        check_rbbh <- ortho$rbbh[index] == "TRUE"
        if(sum(check_rbbh) == 1){
            return(index[check_rbbh])
            
        } else if(sum(check_rbbh) > 1){
            index <- index[check_rbbh]
        }
        
        check_e <- ortho$mutual_e[index] == min(ortho$mutual_e[index])
        if(sum(check_e) == 1){
            return(index[check_e])
            
        } else if(sum(check_e) > 1){
            index <- index[check_e]
        }
        
        check_ci <- ortho$mutual_ci[index] == max(ortho$mutual_ci[index])
        if(sum(check_ci) == 1){
            return(index[check_ci])
            
        } else if(sum(check_ci) > 1){
            index <- index[check_ci]
        }
        
        if(length(index) == 1){
            return(index)
        }
        
        len <- q_len[index] + s_len[index]
        check_len <- len == max(len)
        if(sum(check_len) == 1){
            return(index[check_len])
            
        } else {
            return(sample(index, 1))
        }
    })
    return(ortho[best_pair, ])
}

################################################################################
summarySynog <- function(object, gene = FALSE, split = FALSE){
    stopifnot(inherits(x = object, "SynogDB"))
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    if(gene){
        if(split){
            if(!H5Lexists(h5, "synog_gene_split/summary")){
                stop("Run geneOrtho to obtain genewise ortholog info.")
            }
            df <- h5$synog_gene_split$summary
            
        } else {
            if(!H5Lexists(h5, "synog_gene/summary")){
                stop("Run geneOrtho to obtain genewise ortholog info.")
            }
            df <- h5$synog_gene$summary
        }
    } else {
        
        if(!H5Lexists(h5, "synog_tx/summary")){
            stop("Run syntenyOrtho to obtain ortholog info.")
        }
        df <- h5$synog_tx$summary
    }
    rownames(df) <- c("Total", 
                      "Classified", 
                      "1to1",
                      "1toM",
                      "Mto1",
                      "MtoM",
                      "Orphan")
    return(df)
}

################################################################################
mapProt <- function(object, out_prefix, miniprot_bin, n_core, len_diff = 0.5){
    stopifnot(inherits(x = object, "SynogDB"))
    .mapEngine(object = object, subject_prot = object$subject_prot,
               query_prot = object$query_prot, out_prefix = out_prefix,
               miniprot_bin = miniprot_bin, n_core = n_core,
               overlap = TRUE, len_diff = len_diff)
}

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

################################################################################
createFASTA <- function(object, out_dir){
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    q_cds <- .makeCDS(gff = as.vector(h5$protmap$s2q_gff), 
                      genome = object$query_genome)
    s_cds <- .makeCDS(gff = as.vector(h5$protmap$q2s_gff), 
                      genome = object$subject_genome)
    q_cds_fn <- sub("\\.gff", "_cds.fa", as.vector(h5$protmap$s2q_gff))
    s_cds_fn <- sub("\\.gff", "_cds.fa", as.vector(h5$protmap$q2s_gff))
    writeXStringSet(q_cds, q_cds_fn)
    writeXStringSet(s_cds, s_cds_fn)
    
    .h5overwrite(obj = q_cds_fn, file = object$h5, "protmap/s2q_cds")
    .h5overwrite(obj = s_cds_fn, file = object$h5, "protmap/q2s_cds")
}

.makeCDS <- function(gff, genome){
    txdb <- makeTxDbFromGFF(file = gff)
    genome <- readDNAStringSet(filepath = genome)
    cds_db <- cdsBy(x = txdb, by = "tx", use.names = TRUE)
    cds <- extractTranscriptSeqs(x = genome, transcripts = cds_db)
    cds <- cds[order(names(cds))]
    return(cds)
}

################################################################################
updateFiles <- function(object){
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    object$query_gff <- as.vector(h5$protmap$s2q_gff)
    object$subject_gff <- as.vector(h5$protmap$q2s_gff)
    object$query_cds <- as.vector(h5$protmap$s2q_cds)
    object$subject_cds <- as.vector(h5$protmap$q2s_cds)
    return(object)
}

################################################################################
getSynog <- function(object, gene = FALSE, split = FALSE){
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    if(gene){
        if(split){
            if(!H5Lexists(h5, "synog_gene_split/orthopairs")){
                stop("Run geneOrtho to obtain genewise ortholog info.")
            }
            out <- h5$synog_gene_split$orthopairs
            out <- out[, c(12:11, 6:5, 1:2, 9:10, 3:4)]
            names(out) <- c("OG", "class", "rbbh", "syntenic", 
                            "query_txID", "subject_txID", 
                            "query_geneID", "subject_geneID",
                            "mutual_e", "mutual_covidt")
            out <- out[order(out$OG), ]
            
        } else {
            if(!H5Lexists(h5, "synog_gene/orthopairs")){
                stop("Run geneOrtho to obtain genewise ortholog info.")
            }
            out <- h5$synog_gene$orthopairs
            out <- out[, c(12:11, 6:5, 1:2, 9:10, 3:4)]
            names(out) <- c("OG", "class", "rbbh", "syntenic", 
                            "query_txID", "subject_txID", 
                            "query_geneID", "subject_geneID",
                            "mutual_e", "mutual_covidt")
            out <- out[order(out$OG), ]
        }
        
    } else {
        if(!H5Lexists(h5, "synog_tx/orthopairs")){
            stop("Run syntenyOrtho to obtain ortholog info.")
        }
        out <- h5$synog_tx$orthopairs
        out <- out[, c(8:5, 1:2, 9:10, 3:4)]
        names(out) <- c("OG", "class", "rbbh", "syntenic", 
                        "query_txID", "subject_txID", 
                        "query_geneID", "subject_geneID",
                        "mutual_e", "mutual_covidt")
        out <- out[order(out$OG), ]
    }
    return(out)
}

getOrphan <- function(object, gene = FALSE, split = FALSE){
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    if(gene){
        if(split){
            if(!H5Lexists(h5, "synog_gene_split/orphan")){
                stop("Run geneOrtho to obtain genewise ortholog info.")
            }
            out <- h5$synog_gene_split$orphan
            names(out$query) <- "query_geneID"
            names(out$subject) <- "subject_geneID"
            
        } else {
            if(!H5Lexists(h5, "synog_gene/orphan")){
                stop("Run geneOrtho to obtain genewise ortholog info.")
            }
            out <- h5$synog_gene$orphan
            names(out$query) <- "query_geneID"
            names(out$subject) <- "subject_geneID"
        }
        
    } else {
        if(!H5Lexists(h5, "synog_tx/orphan")){
            stop("Run syntenyOrtho to obtain ortholog info.")
        }
        out <- h5$synog_tx$orphan
        names(out$query) <- "query_txID"
        names(out$subject) <- "subject_txID"
    }
    return(out)
}

################################################################################
splitGenes <- function(object){
    h5 <- H5Fopen(object$h5)
    on.exit(H5Fclose(h5))
    if(!H5Lexists(h5, "synog_gene/orthopairs")){
        stop("Run geneOrtho to obtain genewise ortholog info.")
    }
    
    ortho <- h5$synog_gene$orthopairs
    ortho$pair_id1 <- paste(ortho$qseqid, ortho$sgeneid, sep = "_")
    ortho$pair_id2 <- paste(ortho$qgeneid, ortho$sseqid, sep = "_")
    ortho$pair_id3 <- paste(ortho$qseqid, ortho$sseqid, sep = "_")
    
    gff <- .prepGFF(object)
    
    ortho$genewise_class <- .eval1toM(ortho = ortho, gff = gff)
    ortho$genewise_class <- .evalMto1(ortho = ortho, gff = gff)
    ortho$genewise_class <- .evalMtoM(ortho = ortho, gff = gff)
    
    new_id <- .replaceIDs(ortho = ortho)
    ortho$qgeneid <- new_id$qgeneid
    ortho$sgeneid <- new_id$sgeneid
    ortho$pair_id <- paste(ortho$qgeneid, ortho$sgeneid, sep = "_")
    
    class_og <- .numberingOGgene(ortho = ortho)
    ortho$genewise_class <- class_og$class
    ortho$genewise_OG <- class_og$OG
    ortho <- ortho[order(ortho$genewise_OG, ortho$qgeneid, ortho$sgeneid), ]
    rownames(ortho) <- NULL
    orphan <- h5$synog_gene$orphan
    gene_summary <- .makeSummary(ortho = ortho, orphan = orphan)
    ortho$pair_id <- ortho$pair_id1 <- ortho$pair_id2 <- ortho$pair_id3 <- NULL
    
    .h5creategroup(object$h5,"synog_gene_split")
    .h5overwrite(obj = ortho, file = object$h5, "synog_gene_split/orthopairs")
    .h5overwrite(obj = gene_summary, file = object$h5, "synog_gene_split/summary")
    .h5overwrite(obj = orphan, file = object$h5, "synog_gene_split/orphan")
}

.prepGFF <- function(object){
    query_gff <- .importAllGFF(object$query_gff)
    subject_gff <- .importAllGFF(object$subject_gff)
    query_gff <- query_gff[query_gff$type == "CDS"]
    query_gff$Parent <- unlist(query_gff$Parent)
    subject_gff <- subject_gff[subject_gff$type == "CDS"]
    subject_gff$Parent <- unlist(subject_gff$Parent)
    return(list(query_gff = query_gff, subject_gff = subject_gff))
}

.eval1toM <- function(ortho, gff = gff){
    out <- ortho$genewise_class
    tmp <- subset(ortho, select = c(qgeneid, qseqid, sgeneid, pair_id1), 
                  subset = genewise_class == "1toM")
    
    og_11 <- .find1to1split(seqid = tmp$qseqid, geneid = tmp$sgeneid, 
                            pair_id = tmp$pair_id1, gff = gff$query_gff)
    out[ortho$pair_id1 %in% og_11] <- "split_1toM"
    
    rest <- tmp[!tmp$pair_id1 %in% og_11, ]
    og_11 <- .pick1to1rest(rest = rest, pair_id = rest$pair_id1)
    out[ortho$pair_id1 %in% og_11] <- "split_1toM"
    
    return(out)
}

.evalMto1 <- function(ortho, gff = gff){
    out <- ortho$genewise_class
    tmp <- subset(ortho, select = c(qgeneid, sseqid, sgeneid, pair_id2), 
                  subset = genewise_class == "Mto1")
    
    og_11 <- .find1to1split(seqid = tmp$sseqid, geneid = tmp$qgeneid, 
                            pair_id = tmp$pair_id2, gff = gff$subject_gff)
    out[ortho$pair_id2 %in% og_11] <- "split_Mto1"
    
    rest <- tmp[!tmp$pair_id2 %in% og_11, ]
    og_11 <- .pick1to1rest(rest = rest, pair_id = rest$pair_id2)
    out[ortho$pair_id2 %in% og_11] <- "split_Mto1"
    
    return(out)
}

.evalMtoM <- function(ortho, gff = gff){
    out <- ortho$genewise_class
    tmp <- subset(ortho, select = c(qseqid, qgeneid, sseqid, sgeneid,
                                    pair_id1, pair_id2, pair_id3), 
                  subset = genewise_class == "MtoM")
    
    og_11_1 <- .find1to1split(seqid = tmp$qseqid, geneid = tmp$sgeneid, 
                              pair_id = tmp$pair_id3, gff = gff$query_gff)
    og_11_2 <- .find1to1split(seqid = tmp$sseqid, geneid = tmp$qgeneid, 
                              pair_id = tmp$pair_id3, gff = gff$subject_gff)
    
    og_11 <- og_11_1[og_11_1 %in% og_11_2]
    out[ortho$pair_id3 %in% og_11] <- "split_MtoM"
    
    rest <- tmp[!tmp$pair_id3 %in% og_11, ]
    og_11 <- .pick1to1rest(rest = rest, pair_id = rest$pair_id3)
    out[ortho$pair_id3 %in% og_11] <- "split_MtoM"
    
    return(out)
}

.find1to1split <- function(seqid, geneid, pair_id, gff){
    out <- pair_id
    tx_dup <- duplicated(seqid)
    tx_dup <- seqid %in% seqid[tx_dup]
    gene_dup <- duplicated(geneid)
    gene_dup <- geneid %in% geneid[gene_dup]
    seqid_1to1 <- seqid[!(tx_dup | gene_dup)]
    out <- out[!(tx_dup | gene_dup)]
    
    gff_1to1 <- gff[gff$Parent %in% seqid_1to1]
    ol <- findOverlaps(gff_1to1, gff_1to1)
    ol <- as.data.frame(ol)
    ol <- subset(ol, subset = queryHits != subjectHits)
    hit_id <- unique(gff_1to1$Parent[ol$queryHits])
    seqid_1to1_candidate <- seqid_1to1[!seqid_1to1 %in% hit_id]
    seqid_1to1_rest <- seqid[!seqid %in% seqid_1to1_candidate]
    out <- out[seqid_1to1 %in% seqid_1to1_candidate]
    
    gff_1to1_candidate <- gff[gff$Parent %in% seqid_1to1_candidate]
    gff_1to1_rest <- gff[gff$Parent %in% seqid_1to1_rest]
    ol <- findOverlaps(gff_1to1_candidate, gff_1to1_rest)
    invalid_seqid <- gff_1to1_candidate$Parent[queryHits(ol)]
    out <- out[!seqid_1to1_candidate %in% invalid_seqid]
    return(out)
}

.pick1to1rest <- function(rest, pair_id){
    out <- pair_id
    q_dup <- duplicated(rest$qgeneid)
    q_dup <- rest$qgeneid %in% rest$qgeneid[q_dup]
    s_dup <- duplicated(rest$sgeneid)
    s_dup <- rest$sgeneid %in% rest$sgeneid[s_dup]
    out <- out[!(q_dup | s_dup)]
    return(out)
}

.replaceIDs <- function(ortho){
    out <- list(qgeneid = ortho$qgeneid, sgeneid = ortho$sgeneid)
    hit <- ortho$genewise_class == "split_1toM"
    out$qgeneid[hit] <- ortho$qseqid[hit]
    
    hit <- ortho$genewise_class == "split_Mto1"
    out$sgeneid[hit] <- ortho$sseqid[hit]
    
    hit <- ortho$genewise_class == "split_MtoM"
    out$qgeneid[hit] <- ortho$qseqid[hit]
    out$sgeneid[hit] <- ortho$sseqid[hit]
    return(out)
}
