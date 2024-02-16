#'
#'
#' @export
#'
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

#' @importFrom Biostring readDNAStringSet writeXStringSet
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

#' Define a function to load sibeliaz's raw output to data.frame
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#' @importFrom rtracklayer import.gff3
#'
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
#' Define a function to convert sibelia's output to data.frame
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#'
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
#' Define a function to show LCB information
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#'
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
#' Define a function to classify LCBs
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#'
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
#' Define a function to summarize LCB statistics
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#'
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
#' Define a function to make a table that show all pairs of blocks between
#' the query genome and the subject genome
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#'
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
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#' @import ggplot2
#'
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
