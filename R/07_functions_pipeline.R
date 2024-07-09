#' Function to go through the Synog pipeline
#'
#' This is a wrapper function to execute the series of functions to go through
#' the Synog pipeline.
#'
#' @export
#'
#' @importFrom parallel detectCores
runSynog <- function(query_genome,
                     subject_genome,
                     query_gff,
                     subject_gff,
                     query_cds,
                     subject_cds,
                     query_prot,
                     subject_prot,
                     hdf5_path = "./synog.h5",
                     sibeliaz_out_dir = "./sibeliaz_out",
                     sibeliaz_bin = "sibeliaz",
                     maf2synteny_bin = "maf2synteny",
                     conda = "conda",
                     sibeliaz_condaenv = "sibeliaz",
                     miniprot_bin = "miniprot",
                     miniprot_condaenv = "miniprot",
                     miniprot_out_dir = "./miniprot_out/",
                     n_threads = NULL,
                     omit_chr = NULL,
                     verbose = TRUE){

    if(is.null(n_threads)){
        core <- detectCores()
        n_threads <- core - 1
    }

    if(verbose){
        message("Preparing a SynogDB class object.")
    }

    object <- makeSynogDB(query_genome = query_genome,
                          subject_genome = subject_genome,
                          query_gff = query_gff,
                          subject_gff = subject_gff,
                          query_cds = query_cds,
                          subject_cds = subject_cds,
                          query_prot = query_prot,
                          subject_prot = subject_prot,
                          hdf5_path = hdf5_path,
                          overwrite = FALSE)

    if(verbose){
        message("Running SibeliaZ for local collinear block (LCB) detection.")
    }

    runSibeliaZ(object = object,
                out_dir = sibeliaz_out_dir,
                sibeliaz_bin = sibeliaz_bin,
                maf2synteny_bin = maf2synteny_bin,
                conda = conda,
                condaenv = sibeliaz_condaenv,
                run_sibeliaz = TRUE)

    if(verbose){
        message("Summarize the results for the LCB detection.")
    }

    sibeliaLCB2DF(object = object)
    lcbClassify(object = object)
    getLCBpairs(object = object)
#
#     if(verbose){
#         message("Performinig reciprocal BLAST.")
#     }
#
#     rbh(object = object, n_threads = n_threads)
#
#     if(verbose){
#         message("Find anchor orthologs.")
#     }
#
#     anchorOrtho(object = object)
#
#     if(verbose){
#         message("Pairing orthologs.")
#     }
#
#     syntenyOrtho(object = object, omit_chr = omit_chr)
#
#     if(verbose){
#         message("Sammarize genewise orthology information.")
#     }
#
#     geneOrtho(object = object)

    if(verbose){
        message("Gene modeling by Miniprot to find missing genes.")
    }

    mapProt(object = object,
            miniprot_bin = miniprot_bin,
            conda = conda,
            condaenv = miniprot_condaenv,
            n_threads = n_threads,
            out_dir = miniprot_out_dir)

    if(verbose){
        message("Performinig reciprocal BLAST.")
    }

    rbh(object = object, n_threads = n_threads)

    # if(verbose){
    #     message("Find anchor orthologs.")
    # }
    #
    # anchorOrtho(object = object)

    if(verbose){
        message("Pairing orthologs.")
    }
    syntenicOrtho(object = object)
    #
    # syntenyOrtho(object = object, omit_chr = omit_chr)
    #
    # if(verbose){
    #     message("Sammarize genewise orthology information.")
    # }
    #
    # geneOrtho(object = object)
    #
    # if(verbose){
    #     message("Split chemric annotations based on the ortholog pairing result.")
    # }
    #
    # splitGenes(object = object)
    return(object)
}
