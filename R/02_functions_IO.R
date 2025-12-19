#' Create a OrthoPairDB object
#'
#' This function creates a OrthoPairDB object containing genome, GFF, CDS, and protein information
#' for query and subject genomes. It also initializes an HDF5 file to store related data.
#'
#' @param query_genome Path to the query genome file.
#' @param subject_genome Path to the subject genome file.
#' @param query_gff Path to the query GFF file.
#' @param subject_gff Path to the subject GFF file.
#' @param query_cds Path to the query CDS file.
#' @param subject_cds Path to the subject CDS file.
#' @param query_prot Path to the query protein file.
#' @param subject_prot Path to the subject protein file.
#' @param hdf5_path Path to the HDF5 file (default is "./orthopair.h5").
#'
#' @return A OrthoPairDB object.
#' @export
#'
makeOrthoPairDB <- function(input_files, overwrite = FALSE, verbose = TRUE){
    out <- NULL
    out$h5 <- .makeHDF5(hdf5_path = input_files$hdf5_path, 
                        overwrite = overwrite, verbose = verbose)
    
    out$resume <- .checkResumePoint(hdf5_path = out$h5)
    
    .h5overwrite(obj = input_files$pair_id,
                 file = out$h5, "pair_id")
    .h5creategroup(out$h5,"files")
    .h5overwrite(obj = input_files$query_genome,
                 file = out$h5, "files/query_genome")
    .h5overwrite(obj = input_files$subject_genome,
                 file = out$h5, "files/subject_genome")
    .h5overwrite(obj = input_files$query_prot,
                 file = out$h5, "files/query_prot")
    .h5overwrite(obj = input_files$subject_prot,
                 file = out$h5, "files/subject_prot")
    .h5overwrite(obj = input_files$query_gff,
                 file = out$h5, "files/query_gff")
    .h5overwrite(obj = input_files$subject_gff,
                 file = out$h5, "files/subject_gff")
    .h5overwrite(obj = input_files$query_cds,
                 file = out$h5, "files/query_cds")
    .h5overwrite(obj = input_files$subject_cds,
                 file = out$h5, "files/subject_cds")
    .h5overwrite(obj = input_files$query_h5,
                 file = out$h5, "files/query_h5")
    .h5overwrite(obj = input_files$subject_h5,
                 file = out$h5, "files/subject_h5")
    
    # Assign class and initialize HDF5 file
    class(out) <- c(class(out), "OrthoPairDB")
    
    .h5creategroup(out$h5, "timestamp")
    .h5overwrite(obj = as.character(Sys.time()),
                 file = out$h5, "timestamp/makedb")
    
    return(out)
}

#' @importFrom Biostrings readDNAStringSet
#' @importFrom BiocGenerics width
.genomeSummary <- function(genome){
    # Read the genome file as a DNAStringSet object
    genome <- readDNAStringSet(genome)
    
    # Create a summary list with sequence names and lengths
    out <- list(names = names(genome),
                length = width(genome))
    return(out)
}

#' Summarize OrthoPairDB object
#'
#' This function provides a summary of the OrthoPairDB object, optionally for gene or transcript information.
#'
#' @param object A OrthoPairDB object.
#' @param hdf5_fn A path to a hdf5 file storing OrthoPaiR results.
#'
#' @return A summary dataframe.
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#'
#' @export
#'
summaryOrthoPair <- function(object = NULL, hdf5_fn = NULL){
    if(is.null(object)){
        h5 <- H5Fopen(hdf5_fn)
        
    } else {
        h5 <- H5Fopen(object$h5)
    }
    on.exit(H5Fclose(h5))
    
    if(!H5Lexists(h5, "orthopair")){
        stop("No genewise ortholog info.")
    }
    
    check <- h5$data_type == "reorg_orthopair"
    if(check){
        stop("This function only accepts the output object or hdf5_fn for a pair of genomes.\n", 
             "You specified a hdf5_fn for the output from the orgMPgenes() function.", 
             call. = FALSE)
    }
    
    orthopair <- h5$orthopair[[1]]
    query <- subset(orthopair, subset = !duplicated(orthopair$query_gene))
    query_summary <- table(factor(query$class, levels = c("1to1", "1toM", "Mto1", "MtoM")))
    subject <- subset(orthopair, subset = !duplicated(orthopair$subject_gene))
    subject_summary <- table(factor(subject$class, levels = c("1to1", "1toM", "Mto1", "MtoM")))
    query_orphan <- length(h5$orphan[[1]])
    subject_orphan <- length(h5$orphan[[2]])
    query_total <- sum(query_summary) + query_orphan
    subject_total <- sum(subject_summary) + subject_orphan
    files <- h5$files
    query_name <- basename(sub("/input.h5", "", files$query_h5))
    subject_name <- basename(sub("/input.h5", "", files$subject_h5))
    df <- data.frame(Query = c(query_name,
                               query_total,
                               query_total - query_orphan,
                               query_summary,
                               query_orphan),
                     Subject = c(subject_name,
                                 subject_total,
                                 subject_total - subject_orphan,
                                 subject_summary,
                                 subject_orphan))
    rownames(df) <- c("Name",
                      "Total",
                      "Classified",
                      "1to1",
                      "1toM",
                      "Mto1",
                      "MtoM",
                      "Orphan")
    return(df)
}

#' Get OrthoPair information
#'
#' This function retrieves ortholog pairs information from the OrthoPairDB object.
#'
#' @param object A OrthoPairDB object.
#' @param hdf5_fn A path to a hdf5 file storing OrthoPaiR results.
#'
#' @return A dataframe containing ortholog pairs.
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#' @export
getOrthoPair <- function(object = NULL,
                         hdf5_fn = NULL, 
                         score = FALSE,
                         loc = FALSE){
    if(is.null(object)){
        h5 <- H5Fopen(hdf5_fn)
        
    } else {
        h5 <- H5Fopen(object$h5)
    }
    on.exit(H5Fclose(h5))
    
    if(!H5Lexists(h5, "orthopair")){
        stop("Run syntenicOrtho() to obtain genewise ortholog info.")
    }
    
    target_col <- c("original_query_gene", "original_subject_gene", 
                    "query_gene", "query_tx", 
                    "subject_gene", "subject_tx", 
                    "query_synteny_block", "subject_synteny_block",
                    "class", "SOG", 
                    "query_collapse", "subject_collapse")
    if(score){
        target_col <- c(target_col, 
                        "pident", "qcovs_q2s", "qcovs_s2q", 
                        "ci_q2s", "ci_s2q", "mutual_ci")
    }
    
    if(loc){
        target_col <- c(target_col, 
                        "query_chr", "query_start",
                        "query_end", "query_strand",
                        "subject_chr", "subject_start",
                        "subject_end", "subject_strand")
    }
    
    out <- NULL
    name <- names(h5$orthopair)
    for(i in seq_along(name)){
        col_names <- names(h5$orthopair[[i]])
        out_i <- h5$orthopair[[i]][, target_col]
        
        if(all(out_i$query_collapse == 0)){
            out_i$query_collapse <- NA
        }
        if(all(out_i$subject_collapse == 0)){
            out_i$subject_collapse <- NA
        }
        out <- c(out, list(out_i))
    }
    names(out) <- name
    return(out)
}

#' Get orphan information
#'
#' This function retrieves orphan genes or transcripts information from the OrthoPairDB object.
#'
#' @param object A OrthoPairDB object.
#' @return A dataframe containing orphan genes or transcripts.
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#' @export
#'
getOrphan <- function(object = NULL, hdf5_fn = NULL, loc = FALSE){
    if(is.null(object)){
        h5 <- H5Fopen(hdf5_fn)
        
    } else {
        h5 <- H5Fopen(object$h5)
    }
    on.exit(H5Fclose(h5))
    
    check <- h5$data_type == "reorg_orthopair"
    if(check){
        out <- h5$orphan_gene
        
    } else {
        if(!H5Lexists(h5, "orphan")){
            stop("Run geneOrtho to obtain genewise ortholog info.")
        }
        out <- list(query = unique(h5$orphan[[1]]),
                    subject = unique(h5$orphan[[2]]))
    }
    return(out)
}


#' Get Meta data of pairing
#'
#' This function retrieves meta data from the OrthoPairDB object.
#'
#' @param object A OrthoPairDB object.
#' @param hdf5_fn A path to a hdf5 file storing OrthoPaiR results.
#'
#' @return A dataframe containing ortholog pairs.
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#' @export
getMeta <- function(object = NULL, hdf5_fn = NULL){
    if(is.null(object)){
        h5 <- H5Fopen(hdf5_fn)
        
    } else {
        h5 <- H5Fopen(object$h5)
    }
    on.exit(H5Fclose(h5))
    
    if(!H5Lexists(h5, "orthopair")){
        stop("Run syntenicOrtho() to obtain genewise ortholog info.")
    }
    if(H5Lexists(h5loc = h5, name = "orphan_gene")){
        genomes <- names(h5$orphan_gene)
    }
    if(H5Lexists(h5loc = h5, name = "orphan")){
        genomes <- names(h5$orphan)
    }
    if(H5Lexists(h5loc = h5, name = "orthopair")){
        pair_id <- names(h5$orthopair)
    }
    if(H5Lexists(h5loc = h5, name = "orthopair_gene")){
        pair_id <- names(h5$orthopair_gene)
    }
    if(H5Lexists(h5loc = h5, name = "files")){
        files <- names(h5$files)
    } else {
        files <- NULL
    }
    if(H5Lexists(h5loc = h5, name = "timestamp")){
        timestamp <- h5$timestamp
    } else {
        timestamp <- NULL
    }
    out <- c(list(genomes = genomes),
             list(pair_id = pair_id),
             list(files = files),
             list(timestamp = timestamp))
    return(out)
}
