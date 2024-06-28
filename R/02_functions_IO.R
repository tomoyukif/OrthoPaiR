#' Create a SynogDB object
#'
#' This function creates a SynogDB object containing genome, GFF, CDS, and protein information
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
#' @param hdf5_path Path to the HDF5 file (default is "./synog.h5").
#'
#' @return A SynogDB object.
#' @export
#'
makeSynogDB <- function(query_genome, subject_genome,
                        query_gff, subject_gff,
                        query_cds, subject_cds,
                        query_prot, subject_prot,
                        hdf5_path = "./synog.h5"){
    # Create a list containing all input file paths
    out <- list(query_genome = query_genome,
                subject_genome = subject_genome,
                query_gff = query_gff,
                subject_gff = subject_gff,
                query_cds = query_cds,
                subject_cds = subject_cds,
                query_prot = query_prot,
                subject_prot = subject_prot)

    # Check if all input files exist
    for(i in seq_along(out)){
        if(!file.exists(out[[i]])){
            stop(out[[i]], " do not exists!")
        }
    }

    # Summarize the query and subject genomes
    out$genome$query <- .genomeSummary(genome = query_genome)
    out$genome$subject <- .genomeSummary(genome = subject_genome)

    # Assign class and initialize HDF5 file
    class(out) <- c(class(out), "SynogDB")
    out$h5 <- .makeHDF5(hdf5_path = hdf5_path)
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

#' Summarize SynogDB object
#'
#' This function provides a summary of the SynogDB object, optionally for gene or transcript information.
#'
#' @param object A SynogDB object.
#' @param h5_fn A path to a hdf5 file storing Synog results.
#' @param gene Logical, whether to summarize genewise orthology information (default is FALSE).
#' @param split Logical, whether to summarize split genewise orthology information (default is FALSE).
#'
#' @return A summary dataframe.
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#'
#' @export
#'
summarySynog <- function(object = NULL, h5_fn = NULL, gene = FALSE, split = FALSE){
    if(is.null(object)){
        h5 <- H5Fopen(h5_fn)

    } else {
        h5 <- H5Fopen(object$h5)
    }

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

#' Get Synog information
#'
#' This function retrieves ortholog pairs information from the SynogDB object.
#'
#' @param object A SynogDB object.
#' @param h5_fn A path to a hdf5 file storing Synog results.
#' @param gene Logical, whether to get genewise orthology information (default is FALSE).
#' @param split Logical, whether to get split genewise orthology information (default is FALSE).
#'
#' @return A dataframe containing ortholog pairs.
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#' @export
getSynog <- function(object = NULL, h5_fn = NULL, gene = FALSE, split = FALSE){
    if(is.null(object)){
        h5 <- H5Fopen(h5_fn)

    } else {
        h5 <- H5Fopen(object$h5)
    }
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

#' Get orphan information
#'
#' This function retrieves orphan genes or transcripts information from the SynogDB object.
#'
#' @param object A SynogDB object.
#' @param gene Logical, whether to get orphan gene information.
#' Orphan transcript information will be returned if FALSE. (default is FALSE).
#' @param split Logical, whether to obtain split orphan gene information (default is FALSE).
#'
#' @return A dataframe containing orphan genes or transcripts.
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#' @export
#'
getOrphan <- function(object = NULL, h5_fn = NULL, gene = FALSE, split = FALSE){
    if(is.null(object)){
        h5 <- H5Fopen(h5_fn)

    } else {
        h5 <- H5Fopen(object$h5)
    }
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
