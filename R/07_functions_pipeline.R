#' Function to go through the OrthoPair pipeline for all combinations of genomes
#'
#' This is a wrapper function to execute the series of functions to go through
#' the OrthoPair pipeline for all combinations of genomes.
#'
#' @export
#'
#' @importFrom parallel detectCores
orthopair <- function(in_list,
                      hdf5_out_dir = "./hdf5",
                      sibeliaz_out_dir = "./sibeliaz_out",
                      sibeliaz_bin = "sibeliaz",
                      maf2synteny_bin = "maf2synteny",
                      conda = "conda",
                      sibeliaz_condaenv = "sibeliaz",
                      miniprot_bin = "miniprot",
                      miniprot_condaenv = "miniprot",
                      miniprot_out_dir = "./miniprot_out/",
                      n_threads = NULL,
                      one_to_many = FALSE,
                      verbose = TRUE,
                      overwrite = FALSE){
    
    if(is.null(n_threads)){
        core <- detectCores()
        n_threads <- core - 1
    }
    
    dir.create(path = hdf5_out_dir, showWarnings = FALSE, recursive = TRUE)
    
    .validateInput(in_list = in_list)
    
    pairwise_input <- .prepPairs(in_list = in_list,
                                 hdf5_out_dir = hdf5_out_dir,
                                 sibeliaz_out_dir = sibeliaz_out_dir,
                                 miniprot_out_dir = miniprot_out_dir,
                                 one_to_many = one_to_many)
    hdf5_fn <- NULL
    for(i in seq_along(pairwise_input)){
        object <- .runOrthoPair(query_genome = pairwise_input[[i]]$query_genome,
                                subject_genome = pairwise_input[[i]]$subject_genome,
                                query_gff = pairwise_input[[i]]$query_gff,
                                subject_gff = pairwise_input[[i]]$subject_gff,
                                query_cds = pairwise_input[[i]]$query_cds,
                                subject_cds = pairwise_input[[i]]$subject_cds,
                                query_prot = pairwise_input[[i]]$query_prot,
                                subject_prot = pairwise_input[[i]]$subject_prot,
                                hdf5_path = pairwise_input[[i]]$hdf5_path,
                                sibeliaz_out_dir = pairwise_input[[i]]$sibeliaz_out_dir,
                                sibeliaz_bin = sibeliaz_bin,
                                maf2synteny_bin = maf2synteny_bin,
                                conda = conda,
                                sibeliaz_condaenv = sibeliaz_condaenv,
                                miniprot_bin = miniprot_bin,
                                miniprot_condaenv = miniprot_condaenv,
                                miniprot_out_dir = pairwise_input[[i]]$miniprot_out_dir,
                                n_threads = n_threads,
                                overwrite = overwrite)
        hdf5_fn <- c(hdf5_fn, pairwise_input[[i]]$hdf5_path)
    }
    names(hdf5_fn) <- names(pairwise_input)
    
    hdf5_fn <- reorgOrthopiars(hdf5_fn = hdf5_fn,
                               out_dir = hdf5_out_dir,
                               out_fn = "reorg_orthopair.h5",
                               makeFASTA = TRUE,
                               overwrite = TRUE)
    
    graph <- makeOrthoGraph(hdf5_fn = hdf5_fn)
    df <- graph2df(graph = graph)
    orthopair_fn <- file.path(hdf5_out_dir, "orthopair_list.csv")
    write.csv(x = df, file = orthopair_fn, row.names = FALSE)
    
    orphan <- getOrphan(hdf5_fn = hdf5_fn)
    orphan <- lapply(seq_along(orphan), function(i){
        i_out <- data.frame(orphan[[i]])
        fill_na <- matrix(data = NA, nrow = nrow(i_out), ncol = length(orphan) - 1)
        i_out <- cbind(i_out, fill_na)
        names(i_out) <- names(orphan)
        return(i_out)
    })
    orphan <- do.call("rbind", orphan)
    orphan_fn <- file.path(hdf5_out_dir, "orphan_list.csv")
    write.csv(x = orphan, file = orphan_fn, row.names = FALSE)
    
    out <- c(hdf5_fn = hdf5_fn, orthopair_fn = orthopair_fn, orphan_fn = orphan_fn)
    return(out)
}

.prepPairs <- function(in_list,
                       hdf5_out_dir,
                       sibeliaz_out_dir,
                       miniprot_out_dir,
                       one_to_many = FALSE){
    in_list <- as.data.frame(x = in_list)
    combs <- combn(x = seq_len(nrow(in_list)), m = 2)
    if(one_to_many){
        combs <- combs[, combs[1, ] == 1]
    }
    comb_names <- matrix(data = in_list$name[combs], nrow = 2)
    comb_id <- apply(X = comb_names, MARGIN = 2, FUN = paste, collapse = "_")
    
    out <- NULL
    for(i in seq_along(comb_id)){
        prefix <- comb_id[i]
        i_out <- list(query_genome = in_list$genome[combs[1, i]],
                      subject_genome = in_list$genome[combs[2, i]],
                      query_gff = in_list$gff[combs[1, i]],
                      subject_gff = in_list$gff[combs[2, i]],
                      query_cds = in_list$cds[combs[1, i]],
                      subject_cds = in_list$cds[combs[2, i]],
                      query_prot = in_list$prot[combs[1, i]],
                      subject_prot = in_list$prot[combs[2, i]],
                      hdf5_path = file.path(hdf5_out_dir, paste0(prefix, ".h5")),
                      sibeliaz_out_dir = file.path(sibeliaz_out_dir, prefix),
                      miniprot_out_dir = file.path(miniprot_out_dir, prefix))
        out <- c(out, list(i_out))
    }
    names(out) <- comb_id
    return(out)
}

#' @importFrom parallel detectCores
.runOrthoPair <- function(query_genome,
                          subject_genome,
                          query_gff,
                          subject_gff,
                          query_cds,
                          subject_cds,
                          query_prot,
                          subject_prot,
                          hdf5_path = "./orthopair.h5",
                          sibeliaz_out_dir = "./sibeliaz_out",
                          sibeliaz_bin = "sibeliaz",
                          maf2synteny_bin = "maf2synteny",
                          conda = "conda",
                          sibeliaz_condaenv = "sibeliaz",
                          miniprot_bin = "miniprot",
                          miniprot_condaenv = "miniprot",
                          miniprot_out_dir = "./miniprot_out/",
                          n_threads = NULL,
                          verbose = TRUE,
                          overwrite = FALSE){
    
    if(is.null(n_threads)){
        core <- detectCores()
        n_threads <- core - 1
    }
    
    if(verbose){
        message("Preparing a OrthoPairDB class object.")
    }
    
    object <- makeOrthoPairDB(query_genome = query_genome,
                              subject_genome = subject_genome,
                              query_gff = query_gff,
                              subject_gff = subject_gff,
                              query_cds = query_cds,
                              subject_cds = subject_cds,
                              query_prot = query_prot,
                              subject_prot = subject_prot,
                              hdf5_path = hdf5_path,
                              overwrite = overwrite)
    
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
    
    if(verbose){
        message("Pairing orthologs.")
    }
    syntenicOrtho(object = object)
    return(object)
}

#' @export
orgInputFiles <- function(object = NULL, name, genome, gff, cds, prot){
    .validateInput(name = name,
                   genome = genome,
                   gff = gff,
                   cds = cds,
                   prot = prot)
    
    if(is.null(object)){
        object <- list(name = name,
                       genome = genome,
                       gff = gff,
                       cds = cds,
                       prot = prot)
        class(object) <- c(class(object), "PairwiseOrthoPairInput")
        
    } else {
        if(!inherits(x = object, what = "PairwiseOrthoPairInput")){
            stop("The input object must be a PairwiseOrthoPairInput class object.",
                 call. = FALSE)
        }
        object$name <- c(object$name, name)
        object$genome <- c(object$genome, genome)
        object$gff <- c(object$gff, gff)
        object$cds <- c(object$cds, cds)
        object$prot <- c(object$prot, prot)
    }
    return(object)
}

.validateInput <- function(in_list = NULL, name, genome, gff, cds, prot){
    if(is.null(in_list)){
        check <- is.na(name) | name == ""
        if(check){
            stop("Empty name is not allowed.", call. = FALSE)
        }
        
        check <- !file.exists(genome)
        if(check){
            stop('The file "genome" does not exist.', call. = FALSE)
        }
        
        check <- !file.exists(gff)
        if(check){
            stop('The file "gff" does not exist.', call. = FALSE)
        }
        
        check <- !file.exists(cds)
        if(check){
            stop('The file "cds" does not exist.', call. = FALSE)
        }
        
        check <- !file.exists(prot)
        if(check){
            stop('The file "prot" does not exist.', call. = FALSE)
        }
        
    } else {
        check <- !c("name", "genome", "gff", "cds", "prot") %in% names(in_list)
        if(any(check)){
            stop('The "in_list" object should be a named list ',
                 'containing the following elements:',
                 '"name", "genome", "gff", "cds", "prot"', call. = FALSE)
        }
        for(i in seq_along(in_list$name)){
            .validateInput(name = in_list$name[i],
                           genome = in_list$genome[i],
                           gff = in_list$gff[i],
                           cds = in_list$cds[i],
                           prot = in_list$prot[i])
        }
    }
}
