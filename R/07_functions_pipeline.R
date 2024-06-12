#' Function to go through the Synog pipeline
#'
#' This is a wrapper function to execute the series of functions to go through
#' the Synog pipeline.
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
                     hdf5_fn = "./synog.h5",
                     sibeliaz_out_dir = "./sibeliaz_out",
                     sibeliaz_bin = "sibeliaz",
                     maf2synteny_bin = "maf2synteny",
                     conda = "conda",
                     sibeliaz_condaenv = "sibeliaz",
                     miniprot_bin = "miniprot",
                     miniprot_condaenv = "miniprot",
                     miniprot_out_dir = "./miniprot_out",
                     n_threads = NULL,
                     omit_chr = NULL){

    if(is.null(n_threads)){
        core <- detectCores()
        n_threads <- core - 1
    }

    object <- makeSynogDB(query_genome = query_genome,
                          subject_genome = subject_genome,
                          query_gff = query_gff,
                          subject_gff = subject_gff,
                          query_cds = query_cds,
                          subject_cds = subject_cds,
                          query_prot = query_prot,
                          subject_prot = subject_prot,
                          hdf5_path = hdf5_path)

    runSibeliaZ(object = object,
                out_dir = sibeliaz_out_dir,
                sibeliaz_bin = sibeliaz_bin,
                maf2synteny_bin = maf2synteny_bin,
                conda = conda,
                condaenv = sibeliaz_condaenv,
                run_sibeliaz = TRUE)

    sibeliaLCB2DF(object = object)
    lcbClassify(object = object)
    getLCBpairs(object = object)

    rbh(object = object, n_threads = n_threads)

    anchorOrtho(object = object)

    syntenyOrtho(object = object, omit_chr = omit_chr)

    geneOrtho(object = object)

    mapProt(object = object,
            miniprot_bin = miniprot_bin,
            conda = conda,
            condaenv = miniprot_condaenv,
            n_core = n_threads,
            out_prefix = out_prefix)

    createFASTA(object = object)
    object <- updateFiles(object = object)
    rbh(object = object, n_threads = n_threads)
    anchorOrtho(object = object)
    syntenyOrtho(object = object, omit_chr = omit_chr)
    geneOrtho(object = object)
    splitGenes(object = object)
    return(object)
}
