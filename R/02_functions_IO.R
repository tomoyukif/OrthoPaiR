#' @export
makeSynogDB <- function(query_genome, subject_genome,
                        query_gff, subject_gff,
                        query_cds, subject_cds,
                        query_prot, subject_prot,
                        positive_list, negative_list,
                        hdf5_path = "./synog.h5"){
    out <- list(query_genome = query_genome,
                subject_genome = subject_genome,
                query_gff = query_gff,
                subject_gff = subject_gff,
                query_cds = query_cds,
                subject_cds = subject_cds,
                query_prot = query_prot,
                subject_prot = subject_prot,
                positive_list = positive_list,
                negative_list = negative_list)

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

#' @importFrom Biostrings readDNAStringSet
#' @importFrom BiocGenerics width
.genomeSummary <- function(genome){
    genome <- readDNAStringSet(genome)
    out <- list(names = names(genome),
                length = width(genome))
    return(out)
}


#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
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
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#' @importFrom Biostrings writeXStringSet
#'
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

#' @importFrom Biostrings readDNAStringSet
#' @importFrom GenomicFeatures makeTxDbFromGFF cdsBy extractTranscriptSeqs
.makeCDS <- function(gff, genome){
    txdb <- makeTxDbFromGFF(file = gff)
    genome <- readDNAStringSet(filepath = genome)
    cds_db <- cdsBy(x = txdb, by = "tx", use.names = TRUE)
    cds <- extractTranscriptSeqs(x = genome, transcripts = cds_db)
    cds <- cds[order(names(cds))]
    return(cds)
}

################################################################################
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose
#'
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
#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#'
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

#' @export
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
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

