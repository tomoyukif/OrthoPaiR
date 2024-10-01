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
                        overwrite = FALSE){

    # Create a list containing all input file paths
    files <- list(query_genome = query_genome,
                  subject_genome = subject_genome,
                  query_gff = query_gff,
                  subject_gff = subject_gff,
                  query_cds = query_cds,
                  subject_cds = subject_cds,
                  query_prot = query_prot,
                  subject_prot = subject_prot)

    # Check if all input files exist
    for(i in seq_along(files)){
        if(!file.exists(files[[i]])){
            stop(files[[i]], " do not exists!")
        }
    }

    out <- NULL
    out$h5 <- .makeHDF5(hdf5_path = hdf5_path, overwrite = overwrite)

    .h5creategroup(out$h5,"files")
    .h5overwrite(obj = query_genome,
                 file = out$h5, "files/query_genome")
    .h5overwrite(obj = subject_genome,
                 file = out$h5, "files/subject_genome")
    .h5overwrite(obj = query_gff,
                 file = out$h5, "files/query_gff")
    .h5overwrite(obj = subject_gff,
                 file = out$h5, "files/subject_gff")
    .h5overwrite(obj = query_cds,
                 file = out$h5, "files/query_cds")
    .h5overwrite(obj = subject_cds,
                 file = out$h5, "files/subject_cds")
    .h5overwrite(obj = query_prot,
                 file = out$h5, "files/query_prot")
    .h5overwrite(obj = subject_prot,
                 file = out$h5, "files/subject_prot")

    # Summarize the query and subject genomes
    out$genome$query <- .genomeSummary(genome = query_genome)
    out$genome$subject <- .genomeSummary(genome = subject_genome)

    # Assign class and initialize HDF5 file
    class(out) <- c(class(out), "OrthoPairDB")
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

        if(!any(check)){
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
        gene_list <- h5$gene_list
        genomes <- unique(gene_list$genome)
        genome_comb <- names(h5$orthopair_gene)
        out <- NULL
        for(i in seq_along(genomes)){
            genome_out <- subset(gene_list,
                                 subset = genome == genomes[i],
                                 select = gene)
            genome_out <- unlist(genome_out)
            hit <- grep(genomes[i], genome_comb)
            for(j in hit){
                gene_in_j <- h5$orthopair_gene[[genome_comb[j]]]
                genome_out <- genome_out[!genome_out %in% gene_in_j$query_gene]
                genome_out <- genome_out[!genome_out %in% gene_in_j$subject_gene]
            }
            out <- c(out, list(unique(genome_out)))
        }
        names(out) <- genomes
        
    } else {
        if(!H5Lexists(h5, "orphan_query")){
            stop("Run geneOrtho to obtain genewise ortholog info.")
        }
        out <- list(query = unique(h5$orphan_query),
                    subject = unique(h5$orphan_subject))
    }
    return(out)
}
