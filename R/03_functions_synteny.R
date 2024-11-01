#' Run SibeliaZ
#'
#' This function runs SibeliaZ to detect synteny blocks between the query and subject genomes
#' and processes the output to store in the OrthoPairDB object.
#'
#' @param object A OrthoPairDB object.
#' @param out_dir Output directory for SibeliaZ results.
#' @param sibeliaz_bin Path to the SibeliaZ binary (default is "sibeliaz").
#' @param maf2synteny_bin Path to the maf2synteny binary (default is "maf2synteny").
#' @param conda Path to the conda executable (default is "conda").
#' @param condaenv Name of the conda environment to use (default is NULL).
#' @param run_sibeliaz Logical, whether to run SibeliaZ (default is TRUE).
#'
#' @export
#'
runSibeliaZ <- function(object, out_dir,
                        sibeliaz_bin = "sibeliaz",
                        maf2synteny_bin = "maf2synteny",
                        conda = "conda",
                        condaenv = NULL,
                        run_sibeliaz = TRUE){
    stopifnot(inherits(x = object, "OrthoPairDB"))

    dir.create(path = out_dir, showWarnings = FALSE, recursive = TRUE)

    if(run_sibeliaz){
        # Open the HDF5 file
        h5 <- H5Fopen(object$h5)
        # Ensure the HDF5 file is closed when the function exits
        on.exit(H5Fclose(h5))

        # Prepare query and subject genome files for SibeliaZ
        q_fn <- .prepGenome(genome = h5$files$query_genome, label = "query")
        s_fn <- .prepGenome(genome = h5$files$subject_genome, label = "subject")

        # Run SibeliaZ and maf2synteny with or without conda environment
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

    # Create HDF5 groups and overwrite with new results
    .h5creategroup(object$h5, "sibeliaz")
    .h5overwrite(obj = file.path(out_dir, "blocks_coords.gff"),
                 file = object$h5, "sibeliaz/raw")
    .h5overwrite(obj = file.path(out_dir, "5000/blocks_coords.txt"),
                 file = object$h5, "sibeliaz/blocks")
    
    .h5overwrite(obj = as.character(Sys.time()), file = object$h5, "timestamp/sibeliaz")
}

#' Prepare genome for SibeliaZ
#'
#' This function reads a genome file and modifies sequence names with a given label,
#' then writes the modified sequences to a temporary FASTA file.
#'
#' @param genome Path to the genome file.
#' @param label Label to prepend to sequence names.
#'
#' @return Path to the temporary FASTA file.
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#'
.prepGenome <- function(genome, label){
    dna <- readDNAStringSet(genome)
    names(dna) <- paste(label, names(dna), sep = "_")
    out <- tempfile(fileext = ".fa")
    writeXStringSet(dna, out)
    return(out)
}

#' Execute a command using conda environment
#'
#' This function executes a command within a specified conda environment.
#'
#' @param conda Path to the conda executable.
#' @param env Name of the conda environment to use.
#' @param command Command to execute.
#' @param args Arguments for the command.
.condaExe <- function(conda, env, command, args){
    system2(command = conda,
            args = paste("run -n", env, command, args))
}

#' Convert SibeliaZ raw output to a data.frame
#'
#' This function processes the raw output of SibeliaZ and converts it into a data.frame
#' containing information about Locally Collinear Blocks (LCBs).
#'
#' @param object A OrthoPairDB object.
#' @return A data.frame containing raw LCB information.
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#' @importFrom rtracklayer import.gff3
sibeliaRAW2DF <- function(object){
    # Check if the input object is of class "OrthoPairDB"
    stopifnot(inherits(x = object, "OrthoPairDB"))

    # Open the HDF5 file
    h5 <- H5Fopen(object$h5)
    # Ensure the HDF5 file is closed when the function exits
    on.exit(H5Fclose(h5))

    # Check if the "SibeliaZ/raw" group exists in the HDF5 file
    if(!H5Lexists(h5, "sibeliaz/raw")){
        stop("Run runSibeliaZ to obtain LCB info.")
    }

    # Import the GFF3 file from the HDF5 file
    gff <- import.gff3(h5$sibeliaz$raw)

    # Create a data frame from the imported GFF3 file
    df <- data.frame(chr = as.character(seqnames(gff)),
                     start = start(gff),
                     end = end(gff),
                     ID = gff$ID)

    # Process the data frame to create pairs of query and subject chromosomes
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

    # Combine the processed data into a single data frame
    out <- do.call("rbind", out)
    names(out) <- c(paste(names(out)[1:4], "query", sep = "_"),
                    paste(names(out)[5:8], "subject", sep = "_"))

    # Remove the "query_" and "subject_" prefixes from chromosome names
    out$query_chr <- sub("query_", "", out$query_chr)
    out$subject_chr <- sub("subject_", "", out$subject_chr)
    rownames(out) <- NULL

    # Set the class of the output data frame
    class(out) <- c(class(out), "RawLCB")

    # Overwrite the "sibeliaz/lcb" group in the HDF5 file with the new data frame
    .h5overwrite(obj = out, file = object$h5, "sibeliaz/lcb")
}


#' Convert SibeliaZ output to a data.frame
#'
#' This function processes the output of SibeliaZ and converts it into a data.frame
#' containing information about Locally Collinear Blocks (LCBs).
#'
#' @param object A OrthoPairDB object.
#' @return A data.frame containing LCB information.
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
sibeliaLCB2DF <- function(object){
    # Check if the input object is of class "OrthoPairDB"
    stopifnot(inherits(x = object, "OrthoPairDB"))

    # Open the HDF5 file
    h5 <- H5Fopen(object$h5)
    # Ensure the HDF5 file is closed when the function exits
    on.exit(H5Fclose(h5))

    # Check if the "sibeliaz/blocks" group exists in the HDF5 file
    if(!H5Lexists(h5, "sibeliaz/blocks")){
        stop("Run runSibeliaZ to obtain LCB info.")
    }

    # Check if the blocks file exists within the HDF5 file
    if(!file.exists(h5$sibeliaz$blocks)){
        stop(h5$sibeliaz$blocks, " does not exists.")
    }

    # Read the SibeliaZ blocks file
    sibelia <- scan(file = h5$sibeliaz$blocks, what = "character", sep = "\n")
    # Identify lines separating different sections in the blocks file
    sep_lines <- grep("^-+", sibelia)

    # Extract header lines and process them into a data frame
    header_lines <- 1:(sep_lines[1] - 1)
    header <- sibelia[header_lines]
    header <- do.call("rbind", strsplit(header, "\t"))
    colnames(header) <- header[1, ]
    header <- header[-1, ]
    header <- data.frame(header)
    header$Seq_id <- as.numeric(header$Seq_id)
    header$Size <- as.numeric(header$Size)

    # Identify lines containing block headers and calculate their lengths
    block_i <- grep("^Block", sibelia)
    i_diff <- diff(c(block_i, length(sibelia) + 1))
    entry_n <- sapply(seq_along(i_diff), function(i){
        return(rep(i, i_diff[i] - 3))
    })
    entry_n <- unlist(entry_n)
    block_header_i <- block_i[entry_n]
    block_id <- as.numeric(sub(".*#", "", sibelia[block_header_i]))

    # Extract entries and combine them into a data frame
    entries <- sibelia[-header_lines][!grepl("^Block|^Seq_id|^-+", sibelia[-header_lines])]
    entries <- do.call("rbind", strsplit(entries, "\t"))
    out <- data.frame(entries, block_id)
    names(out) <- c("chr", "strand", "start", "end", "length", "block_id")
    out$chr <- header$Description[match(out$chr, header$Seq_id)]
    out$start <- as.numeric(out$start)
    out$end <- as.numeric(out$end)
    out$length <- as.numeric(out$length)

    # Set the class of the output data frame
    class(out) <- c(class(out), "LCB")

    # Overwrite the "sibeliaz/lcb" group in the HDF5 file with the new data frame
    .h5overwrite(obj = out, file = object$h5, "sibeliaz/lcb")
}

#' Show LCB information
#'
#' This function displays information about the Locally Collinear Blocks (LCBs) in the OrthoPairDB object.
#'
#' @param object A OrthoPairDB object.
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
showLCB <- function(object){
    # Check if the input object is of class "OrthoPairDB"
    stopifnot(inherits(x = object, "OrthoPairDB"))

    # Open the HDF5 file
    h5 <- H5Fopen(object$h5)
    # Ensure the HDF5 file is closed when the function exits
    on.exit(H5Fclose(h5))

    # Check if the "sibeliaz/lcb" group exists in the HDF5 file
    if(!H5Lexists(h5, "sibeliaz/lcb")){
        stop("Run sibeliaLCB2DF to obtain LCB info.")
    }

    # Display chromosome IDs for query and subject genomes
    message("Chromosome IDs:")
    q_id <- grep("^query_", unique(h5$sibeliaz$lcb$chr), value = TRUE)
    q_id <- sort(sub("^query_", "", q_id))
    s_id <- grep("^subject_", unique(h5$sibeliaz$lcb$chr), value = TRUE)
    s_id <- sort(sub("^subject_", "", s_id))

    message("    Query genome:")
    print(q_id)

    message("    Subject genome:")
    print(s_id)

    # Display the number of entries and LCBs
    message("No of entries: ")
    print(length(h5$sibeliaz$lcb$chr))

    message("No of LCBs: ")
    print(length(unique(h5$sibeliaz$lcb$block_id)))
}

#' Classify LCBs
#'
#' This function classifies Locally Collinear Blocks (LCBs) based on their genome origin (query or subject).
#'
#' @param object A OrthoPairDB object.
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
lcbClassify <- function(object){
    # Check if the input object is of class "OrthoPairDB"
    stopifnot(inherits(x = object, "OrthoPairDB"))

    # Open the HDF5 file
    h5 <- H5Fopen(object$h5)
    # Ensure the HDF5 file is closed when the function exits
    on.exit(H5Fclose(h5))

    # Check if the "sibeliaz/lcb" group exists in the HDF5 file
    if(!H5Lexists(h5, "sibeliaz/lcb")){
        stop("Run sibeliaLCB2DF to obtain LCB info.")
    }

    # Initialize a vector to classify the chromosomes
    chr <- rep(0, length(h5$sibeliaz$lcb$chr))
    s_chr <- grepl("^subject_", h5$sibeliaz$lcb$chr)
    chr[s_chr] <- 1

    # Classify LCBs based on the presence of query and subject chromosomes
    lcb_class <- tapply(chr, h5$sibeliaz$lcb$block_id, function(i){
        if(sum(i) == 0){
            out <- "q2q"  # Only query chromosomes present
        } else if(prod(i) == 0){
            out <- "q2s"  # Both query and subject chromosomes present
        } else {
            out <- "s2s"  # Only subject chromosomes present
        }
    })

    # Match the block IDs and overwrite the HDF5 file with the LCB classifications
    hit <- match(h5$sibeliaz$lcb$block_id, as.numeric(names(lcb_class)))
    .h5overwrite(obj = lcb_class[hit], file = object$h5, "sibeliaz/lcb_class")
}

#' Summarize LCB statistics
#'
#' This function calculates and summarizes statistics for LCBs (Locally Collinear Blocks) in the OrthoPairDB object.
#'
#' @param object A OrthoPairDB object.
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#'
statsLCB <- function(object){
    # Check if the input object is of class "OrthoPairDB"
    stopifnot(inherits(x = object, "OrthoPairDB"))

    # Open the HDF5 file
    h5 <- H5Fopen(object$h5)
    # Ensure the HDF5 file is closed when the function exits
    on.exit(H5Fclose(h5))

    # Check if the "sibeliaz/lcb" group exists in the HDF5 file
    if(!H5Lexists(h5, "sibeliaz/lcb")){
        stop("Run sibeliaLCB2DF to obtain LCB info.")
    }

    # Check if the "sibeliaz/lcb_class" group exists in the HDF5 file
    if(!H5Lexists(h5, "sibeliaz/lcb_class")){
        stop("Run lcbClassify to obtain LCB classification info.")
    }

    # Calculate statistics based on the presence of LCB classification
    if(is.null(h5$sibeliaz$lcb_class)){
        # Calculate statistics for query and subject genomes separately
        q_summary <- .calc_lcb_stats(object = h5$sibeliaz$lcb, target = "^query_")
        s_summary <- .calc_lcb_stats(object = h5$sibeliaz$lcb, target = "^subject_")
        out <- rbind(cbind(genome = "query", q_summary),
                     cbind(genome = "subject", s_summary))

    } else {
        # Calculate statistics for different classes of LCBs
        q2q <- h5$sibeliaz$lcb_class == "q2q"
        q_q2q_summary <- .calc_lcb_stats(object = h5$sibeliaz$lcb[q2q, ], target = "^query_")
        s2s <- h5$sibeliaz$lcb_class == "s2s"
        s_s2s_summary <- .calc_lcb_stats(object = h5$sibeliaz$lcb[s2s, ], target = "^subject_")
        q2s <- h5$sibeliaz$lcb_class == "q2s"
        q_q2s_summary <- .calc_lcb_stats(object = h5$sibeliaz$lcb[q2s, ], target = "^query_")
        s_q2s_summary <- .calc_lcb_stats(object = h5$sibeliaz$lcb[q2s, ], target = "^subject_")
        out <- rbind(cbind(genome = "query", class = "Q2Q", q_q2q_summary),
                     cbind(genome = "subject", class = "S2S", s_s2s_summary),
                     cbind(genome = "query", class = "Q2S", q_q2s_summary),
                     cbind(genome = "subject", class = "Q2S", s_q2s_summary))
    }

    # Overwrite the "sibeliaz/lcb_stats" group in the HDF5 file with the new statistics
    .h5overwrite(obj = out, file = object$h5, "sibeliaz/lcb_stats")

    # Print the summarized statistics
    print(out)
}

#' Calculate LCB statistics
#'
#' This function calculates statistics for a given set of LCBs (Locally Collinear Blocks).
#'
#' @param object A data.frame containing LCB information.
#' @param target A regex pattern to match the target genome (e.g., "^query_").
#' @return A data.frame with LCB statistics.
.calc_lcb_stats <- function(object, target){
    # Subset the LCB data frame based on the target genome
    if(is.null(object$genome)){
        tmp_df <- subset(object, subset = grepl(target, chr))
    } else {
        tmp_df <- subset(object, subset = genome == target)
    }

    # Calculate basic statistics
    out <- data.frame(total_length = sum(tmp_df$length),
                      mean_length = mean(tmp_df$length),
                      max_length = max(tmp_df$length),
                      min_length = min(tmp_df$length))

    # Adjust for reverse strand entries
    minus <- tmp_df$strand == "-"
    tmp <- tmp_df$start[minus]
    tmp_df$start[minus] <- tmp_df$end[minus]
    tmp_df$end[minus] <- tmp

    # Create a GRanges object and calculate gap statistics
    tmp_df_nb_gr <- GRanges(seqnames = tmp_df$chr,
                            ranges = IRanges(start = tmp_df$start,
                                             end = tmp_df$end))
    out <- cbind(out,
                 data.frame(total_gap_length = sum(width(gaps(tmp_df_nb_gr))),
                            mean_gap_length = mean(width(gaps(tmp_df_nb_gr))),
                            max_gap_length = max(width(gaps(tmp_df_nb_gr))),
                            min_gap_length = min(width(gaps(tmp_df_nb_gr)))))
    return(out)
}

#' Get pairs of blocks between query and subject genomes
#'
#' This function generates a table showing all pairs of blocks between the query and subject genomes.
#'
#' @param object A OrthoPairDB object.
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
getLCBpairs <- function(object){
    # Check if the input object is of class "OrthoPairDB"
    stopifnot(inherits(x = object, "OrthoPairDB"))

    # Open the HDF5 file
    h5 <- H5Fopen(object$h5)
    # Ensure the HDF5 file is closed when the function exits
    on.exit(H5Fclose(h5))

    # Check if the "sibeliaz/lcb" group exists in the HDF5 file
    if(!H5Lexists(h5, "sibeliaz/lcb")){
        stop("Run sibeliaLCB2DF to obtain LCB info.")
    }

    # Check if the "sibeliaz/lcb_class" group exists in the HDF5 file
    if(!H5Lexists(h5, "sibeliaz/lcb_class")){
        stop("Run lcbClassify to obtain LCB classification info.")
    }

    # Identify LCBs classified as "q2s" (query to subject)
    q2s <- h5$sibeliaz$lcb_class == "q2s"
    lcb_q2s_id_tbl <- table(h5$sibeliaz$lcb$block_id[q2s])
    lcb_id_1to1 <- names(lcb_q2s_id_tbl)[lcb_q2s_id_tbl == 2]
    valid <- h5$sibeliaz$lcb$block_id %in% lcb_id_1to1

    # Create lists of 1-to-1 and non-1-to-1 LCB pairs
    lcb_pairs <- list(lcb_1to1 = .get1to1(h5 = h5, valid = valid),
                      lcb_non_1to1 = .getNon1to1(h5 = h5, valid = valid))

    # Overwrite the "sibeliaz/lcb_pairs" group in the HDF5 file with the new LCB pairs
    .h5overwrite(obj = lcb_pairs, file = object$h5, "sibeliaz/lcb_pairs")
    
    .h5overwrite(obj = as.character(Sys.time()), file = object$h5, "timestamp/lcbpair")
}

#' Get 1-to-1 LCB pairs
#'
#' This function retrieves 1-to-1 Locally Collinear Block (LCB) pairs from the HDF5 file.
#'
#' @param h5 An HDF5 file handle.
#' @param valid A logical vector indicating valid LCBs.
#'
#' @return A data.frame containing 1-to-1 LCB pairs.
.get1to1 <- function(h5, valid){
    # Identify the 1-to-1 LCB pairs where the LCB class is "q2s" (query to subject)
    hit <- valid & h5$sibeliaz$lcb_class == "q2s"

    # Subset the LCB data for the identified 1-to-1 pairs
    df <- h5$sibeliaz$lcb[hit, ]

    # Separate the data into query and subject subsets
    q_df <- subset(df, subset = grepl("^query_", chr), select = c(chr, start, end, block_id))
    s_df <- subset(df, subset = grepl("^subject_", chr), select = c(chr, start, end, block_id))

    # Order the subsets by block ID
    q_df <- q_df[order(q_df$block_id), ]
    s_df <- s_df[order(s_df$block_id), ]

    # Remove the "query_" and "subject_" prefixes from chromosome names
    q_df$chr <- sub("^query_", "", q_df$chr)
    s_df$chr <- sub("^subject_", "", s_df$chr)

    # Combine the query and subject subsets into a single data frame
    out <- cbind(subset(q_df, select = chr:end),
                 subset(s_df, select = c(chr:end, block_id)))

    # Rename the columns of the combined data frame
    names(out) <- c("query_chr", "query_start", "query_end",
                    "subject_chr", "subject_start", "subject_end", "block_id")

    return(out)
}

#' Get non-1-to-1 LCB pairs
#'
#' This function retrieves non-1-to-1 Locally Collinear Block (LCB) pairs from the HDF5 file.
#'
#' @param h5 An HDF5 file handle.
#' @param valid A logical vector indicating valid LCBs.
#'
#' @return A data.frame containing non-1-to-1 LCB pairs.
#'
#' @importFrom dplyr full_join
#'
.getNon1to1 <- function(h5, valid){
    # Identify the non-1-to-1 LCB pairs where the LCB class is "q2s" (query to subject)
    hit <- !valid & h5$sibeliaz$lcb_class == "q2s"

    # Subset the LCB data for the identified non-1-to-1 pairs
    df <- h5$sibeliaz$lcb[hit, ]

    # Separate the data into query and subject subsets
    q_df <- subset(df, subset = grepl("^query_", chr), select = c(chr, start, end, block_id))
    s_df <- subset(df, subset = grepl("^subject_", chr), select = c(chr, start, end, block_id))

    # Remove the "query_" and "subject_" prefixes from chromosome names
    q_df$chr <- sub("^query_", "", q_df$chr)
    s_df$chr <- sub("^subject_", "", s_df$chr)

    # Perform a full join on the query and subject subsets based on block ID
    out <- full_join(q_df, s_df, by = "block_id", relationship = "many-to-many")

    # Rename the columns of the combined data frame
    names(out) <- c("query_chr", "query_start", "query_end", "block_id",
                    "subject_chr", "subject_start", "subject_end")

    # Reorder the columns
    out <- out[, c(1:3, 5:7, 4)]

    return(out)
}

#'
#' Plot LCB pairs
#'
#' This function plots pairs of blocks between the query genome and the subject genome.
#'
#' @param object A OrthoPairDB object.
#' @param class The class of LCB pairs to plot ("1to1", "non1to1", or "both").
#' @param size Point size for the plot (default is 0.5).
#' @param color_1to1 Color for 1-to-1 pairs (default is "darkgreen").
#' @param color_non1to1 Color for non-1-to-1 pairs (default is "blue").
#' @param chr_border_color Color for chromosome borders (default is "gray60").
#' @return A ggplot object.
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#' @import ggplot2
plotLCBpairs <- function(object,
                         class = "both",
                         size = 0.5,
                         color_1to1 = "darkgreen",
                         color_non1to1 = "blue",
                         chr_border_color = "gray60"){
    # Check if the input object is of class "OrthoPairDB"
    stopifnot(inherits(x = object, "OrthoPairDB"))

    # Open the HDF5 file
    h5 <- H5Fopen(object$h5)
    # Ensure the HDF5 file is closed when the function exits
    on.exit(H5Fclose(h5))

    # Check if the "sibeliaz/lcb_pairs" group exists in the HDF5 file
    if(!H5Lexists(h5, "sibeliaz/lcb_pairs")){
        stop("Run getLCBpairs to obtain LCB pair info.")
    }

    # Select the appropriate LCB pairs based on the class parameter
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

    # Calculate the midpoint positions for query and subject
    df$query_pos <- (df$query_start + df$query_end) / 2
    df$subject_pos <- (df$subject_start + df$subject_end) / 2

    # Sort the labels for query and subject chromosomes
    df$query_chr <- .sortLabels(x = df$query_chr, rev = FALSE)
    df$subject_chr <- .sortLabels(x = df$subject_chr, rev = TRUE)

    # Create the ggplot object
    p <- ggplot(df) +
        geom_point(aes(x = query_pos, y = subject_pos, color = class), size = size) +
        scale_color_manual(values = color_value, breaks = color_break) +
        facet_grid(rows = vars(subject_chr), cols = vars(query_chr), scales = "free") +
        xlab("Query") +
        ylab("Subject") +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              strip.text.y.right = element_text(angle = -90),
              panel.spacing = unit(0, "lines"),
              panel.border = element_rect(color = chr_border_color, fill = NA, linewidth = 0.5),
              legend.position = "none")

    # Add chromosome lengths to the plot
    p <- .addChrLen(p = p, chrLen = object$genome)

    return(p)
}

#' Sort labels
#'
#' This function sorts labels, optionally in reverse order.
#'
#' @param x A vector of labels to sort.
#' @param rev Logical, whether to reverse the sort order (default is FALSE).
#'
#' @return A factor with sorted labels.
.sortLabels <- function(x, rev = FALSE){
    if(rev){
        # Sort the unique labels and then reverse the order
        lev <- sort(unique(x))
        rev_lev <- rev(lev)
        x <- factor(x = x, levels = rev_lev)

    } else {
        # Sort the unique labels
        lev <- sort(unique(x))
        x <- factor(x = x, levels = lev)
    }
    return(x)
}

#' Add chromosome lengths to plot
#'
#' This function adds chromosome lengths to a ggplot object.
#'
#' @param p A ggplot object.
#' @param chrLen A list containing chromosome lengths.
#'
#' @return A ggplot object with chromosome lengths added.
.addChrLen <- function(p, chrLen){
    # Create a dummy data frame with expanded grid of chromosome lengths
    dummy <- cbind(expand.grid(x = c(rep(0, length(chrLen$query$length)),
                                     chrLen$query$length),
                               y = c(rep(0, length(chrLen$subject$length)),
                                     chrLen$subject$length)),
                   expand.grid(query_chr = c(chrLen$query$names,
                                             chrLen$query$names),
                               subject_chr = c(chrLen$subject$names,
                                               chrLen$subject$names)))

    # Sort the query and subject chromosome labels
    dummy$query_chr <- .sortLabels(x = dummy$query_chr, rev = TRUE)
    dummy$subject_chr <- .sortLabels(x = dummy$subject_chr, rev = FALSE)

    # Add points to the plot with the dummy data frame
    p <- p +
        geom_point(data = dummy, mapping = aes(x = x, y = y), size = 0)

    return(p)
}
