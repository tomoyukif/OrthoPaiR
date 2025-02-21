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
makeOrthoPairDB <- function(query_genome, subject_genome,
                            query_gff, subject_gff,
                            query_cds, subject_cds,
                            query_prot, subject_prot,
                            hdf5_path = "./orthopair.h5",
                            miniprot_out_dir ="./miniprot",
                            overwrite = FALSE,
                            resume = FALSE,
                            module = NULL,
                            param_list = NULL){
    
    # Create a list containing all input file paths
    files <- list(query_genome = query_genome,
                  subject_genome = subject_genome,
                  query_gff = query_gff,
                  subject_gff = subject_gff,
                  query_cds = query_cds,
                  subject_cds = subject_cds,
                  query_prot = query_prot,
                  subject_prot = subject_prot)
    
    out <- NULL
    out$h5 <- .makeHDF5(hdf5_path = hdf5_path, overwrite = overwrite)
    
    if(is.null(files$query_genome) | is.null(files$subject_genome)){
        query_genome <- "no_query_genome"
        subject_genome <- "no_subject_genome"
        no_genome <- TRUE
        
    } else {
        # Summarize the query and subject genomes
        out$genome$query <- .genomeSummary(genome = query_genome)
        out$genome$subject <- .genomeSummary(genome = subject_genome)
        no_genome <- FALSE
    }
    
    if(is.null(files$query_prot) | is.null(files$subject_prot)){
        query_prot <- "no_query_prot"
        subject_prot <- "no_subject_prot"
        no_prot <- TRUE
    } else {
        no_prot <- FALSE
    }
    
    out$resume <- .checkResumePoint(hdf5_path = out$h5,
                                    resume = resume,
                                    module = module,
                                    no_genome = no_genome,
                                    no_prot = no_prot)
    
    .h5creategroup(out$h5,"files")
    .h5overwrite(obj = query_genome,
                 file = out$h5, "files/query_genome")
    .h5overwrite(obj = subject_genome,
                 file = out$h5, "files/subject_genome")
    .h5overwrite(obj = query_prot,
                 file = out$h5, "files/query_prot")
    .h5overwrite(obj = subject_prot,
                 file = out$h5, "files/subject_prot")
    .h5overwrite(obj = query_gff,
                 file = out$h5, "files/query_gff")
    .h5overwrite(obj = subject_gff,
                 file = out$h5, "files/subject_gff")
    .h5overwrite(obj = query_cds,
                 file = out$h5, "files/query_cds")
    .h5overwrite(obj = subject_cds,
                 file = out$h5, "files/subject_cds")
    
    if(out$resume$set_mp){
        message("Use gene models including miniprot predicted genes.")
        .h5overwrite(obj = file.path(miniprot_out_dir, "miniprot_merge_query.gff"),
                     file = out$h5, "files/query_gff")
        .h5overwrite(obj = file.path(miniprot_out_dir, "miniprot_merge_subject.gff"),
                     file = out$h5, "files/subject_gff")
        .h5overwrite(obj = file.path(miniprot_out_dir, "miniprot_merge_query.cds"),
                     file = out$h5, "files/query_cds")
        .h5overwrite(obj = file.path(miniprot_out_dir, "miniprot_merge_subject.cds"),
                     file = out$h5, "files/subject_cds")
    }
    
    out$param_list <- .validateParamList(param_list = param_list)
    
    # Assign class and initialize HDF5 file
    class(out) <- c(class(out), "OrthoPairDB")
    
    .h5creategroup(out$h5, "timestamp")
    .h5overwrite(obj = as.character(Sys.time()), file = out$h5, "timestamp/makedb")
    
    return(out)
}

.validateParamList <- function(param_list){
    out <- list(len_diff = 0.2,
                pident = 0,
                qcovs = 0,
                evalue = 1e-4,
                rbbh_mci_threshold = 0.1,
                rbh_mci_threshold = 0.4)
    if(!is.null(param_list)){
        check <- names(param_list) %in% names(out)
        if(!all(check)){
            stop("The param_list object must be a named list with the following names:",
                 "\n'len_diff', 'pident', 'qcovs', 'evalue', 'rbbh_mci_threshold', 'rbh_mci_threshold'.")
        }
        if(!is.null(param_list$len_diff)){
            out$len_diff <- param_list$len_diff
        }
        if(!is.null(param_list$pident)){
            out$pident <- param_list$pident
        }
        if(!is.null(param_list$qcovs)){
            out$qcovs <- param_list$qcovs
        }
        if(!is.null(param_list$evalue)){
            out$evalue <- param_list$evalue
        }
        if(!is.null(param_list$rbbh_mci_threshold)){
            out$rbbh_mci_threshold <- param_list$rbbh_mci_threshold
        }
        if(!is.null(param_list$rbh_mci_threshold)){
            out$rbh_mci_threshold <- param_list$rbh_mci_threshold
        }
    }
    message("Use following parameters:")
    message("len_diff = ", out$len_diff)
    message("pident = ", out$pident)
    message("qcovs = ", out$qcovs)
    message("evalue = ", out$evalue)
    message("rbbh_mci_threshold = ", out$rbbh_mci_threshold)
    message("rbh_mci_threshold = ", out$rbh_mci_threshold)
    return(out)
}

#' Summarize genome information
#'
#' This function reads a genome file and returns a summary containing the names and lengths of the sequences.
#'
#' @param genome Path to the genome file.
#'
#' @return A list containing genome names and lengths.
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
summaryOrthoPair <- function(object = NULL, hdf5_fn = NULL, gene = FALSE){
    if(is.null(object)){
        h5 <- H5Fopen(hdf5_fn)
        
    } else {
        h5 <- H5Fopen(object$h5)
    }
    
    if(!H5Lexists(h5, "orthopair_gene")){
        stop("No genewise ortholog info.")
    }
    
    check <- h5$data_type == "reorg_orthopair"
    if(check){
        stop("This function only accepts the output object or hdf5_fn for a pair of genomes.\n", 
             "You specified a hdf5_fn for the output from the orgMPgenes() function.", 
             call. = FALSE)
    }
    
    orthopair <- h5$orthopair_gene
    query <- subset(orthopair, subset = !duplicated(orthopair$query_gene))
    query_summary <- table(query$class)
    subject <- subset(orthopair, subset = !duplicated(orthopair$subject_gene))
    subject_summary <- table(subject$class)
    query_orphan <- length(h5$orphan_query)
    subject_orphan <- length(h5$orphan_subject)
    query_total <- sum(query_summary) + query_orphan
    subject_total <- sum(subject_summary) + subject_orphan
    df <- data.frame(Query = c(query_total,
                               query_total - query_orphan,
                               query_summary,
                               query_orphan),
                     Subject = c(subject_total,
                                 subject_total - subject_orphan,
                                 subject_summary,
                                 subject_orphan))
    rownames(df) <- c("Total",
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
#' @param gene Logical, whether to get genewise orthology information (default is FALSE).
#' @param split Logical, whether to get split genewise orthology information (default is FALSE).
#'
#' @return A dataframe containing ortholog pairs.
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#' @export
getOrthoPair <- function(object = NULL, hdf5_fn = NULL, gene = FALSE){
    if(is.null(object)){
        h5 <- H5Fopen(hdf5_fn)
        
    } else {
        h5 <- H5Fopen(object$h5)
    }
    on.exit(H5Fclose(h5))
    
    check <- h5$data_type == "reorg_orthopair"
    
    if(gene){
        if(!H5Lexists(h5, "orthopair_gene")){
            stop("Run syntenicOrtho() to obtain genewise ortholog info.")
        }
        if(check){
            out <- NULL
            name <- names(h5$orthopair_gene)
            for(i in seq_along(name)){
                out <- c(out, list(h5$orthopair_gene[[i]]))
            }
            names(out) <- name
            
        } else {
            out <- h5$orthopair_gene
        }
        
    } else {
        if(!H5Lexists(h5, "orthopair_tx")){
            stop("Run syntenicOrtho() to obtain genewise ortholog info.")
        }
        
        if(check){
            out <- NULL
            name <- names(h5$orthopair_tx)
            for(i in seq_along(name)){
                out <- c(out, list(h5$orthopair_tx[[i]]))
            }
            names(out) <- name
            
        } else {
            out <- h5$orthopair_tx
        }
    }
    return(out)
}

#' Get orphan information
#'
#' This function retrieves orphan genes or transcripts information from the OrthoPairDB object.
#'
#' @param object A OrthoPairDB object.
#' @param gene Logical, whether to get orphan gene information.
#' Orphan transcript information will be returned if FALSE. (default is FALSE).
#' @param split Logical, whether to obtain split orphan gene information (default is FALSE).
#'
#' @return A dataframe containing orphan genes or transcripts.
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#' @export
#'
getOrphan <- function(object = NULL, hdf5_fn = NULL){
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
        if(!H5Lexists(h5, "orphan_query")){
            stop("Run geneOrtho to obtain genewise ortholog info.")
        }
        out <- list(query = unique(h5$orphan_query),
                    subject = unique(h5$orphan_subject))
    }
    return(out)
}
