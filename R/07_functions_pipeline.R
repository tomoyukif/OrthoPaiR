#' Function to go through the OrthoPair pipeline for all combinations of genomes
#'
#' This is a wrapper function to execute the series of functions to go through
#' the OrthoPair pipeline for all combinations of genomes.
#'
#' @export
#'
#' @importFrom parallel detectCores
#' @importFrom igraph write_graph read_graph
orthopair <- function(in_list,
                      working_dir = ".",
                      conda = "conda",
                      miniprot_bin = "miniprot",
                      miniprot_condaenv = "miniprot",
                      diamond_exec_path = "diamond",
                      n_threads = NULL,
                      target_pair = NULL,
                      use_prot = FALSE,
                      verbose = TRUE,
                      overwrite = FALSE,
                      resume = FALSE,
                      module = NULL,
                      orthopair = TRUE,
                      reorg = TRUE,
                      makegraph = TRUE,
                      output_table = FALSE){
    if(!inherits(x = in_list, what = "PairwiseOrthoPairInput")){
        stop("The input object must be a PairwiseOrthoPairInput class object.",
             call. = FALSE)
    }
    
    if(is.null(n_threads)){
        core <- detectCores()
        n_threads <- core - 1
    }
    hdf5_out_dir <- file.path(working_dir, "hdf5_out")
    miniprot_out_dir <- file.path(working_dir, "miniprot_out")
    diamond_out_dir <- file.path(working_dir, "diamond_out")
    dir.create(path = hdf5_out_dir, showWarnings = FALSE, recursive = TRUE)
    
    if(use_prot & !file.exists(diamond_exec_path)){
       stop("Diamond executable binary file was not found at ", 
            diamond_exec_path, 
            call. = FALSE) 
    }
    
    pairwise_input <- .prepPairs(in_list = in_list,
                                 hdf5_out_dir = hdf5_out_dir,
                                 working_dir = working_dir,
                                 miniprot_out_dir = miniprot_out_dir,
                                 diamond_out_dir = diamond_out_dir,
                                 target_pair = target_pair)
    hdf5_fn <- NULL
    on.exit({return(hdf5_fn)})
    error <- FALSE
    for(i in seq_along(pairwise_input)){
        if(!orthopair){
            hdf5_fn <- c(hdf5_fn, pairwise_input[[i]]$hdf5_path)
            next
        }
        message("Start proccessing the ", names(pairwise_input)[i], " pair")
        
        # For debug #######################################################
        # query_genome = pairwise_input[[i]]$query_genome
        # subject_genome = pairwise_input[[i]]$subject_genome
        # query_gff = pairwise_input[[i]]$query_gff
        # subject_gff = pairwise_input[[i]]$subject_gff
        # query_cds = pairwise_input[[i]]$query_cds
        # subject_cds = pairwise_input[[i]]$subject_cds
        # query_prot = pairwise_input[[i]]$query_prot
        # subject_prot = pairwise_input[[i]]$subject_prot
        # hdf5_path = pairwise_input[[i]]$hdf5_path
        # miniprot_out_dir = pairwise_input[[i]]$miniprot_out_dir
        #############################################################
        object <- try(expr = {.runOrthoPair(query_genome = pairwise_input[[i]]$query_genome,
                                            subject_genome = pairwise_input[[i]]$subject_genome,
                                            query_gff = pairwise_input[[i]]$query_gff,
                                            subject_gff = pairwise_input[[i]]$subject_gff,
                                            query_cds = pairwise_input[[i]]$query_cds,
                                            subject_cds = pairwise_input[[i]]$subject_cds,
                                            query_prot = pairwise_input[[i]]$query_prot,
                                            subject_prot = pairwise_input[[i]]$subject_prot,
                                            hdf5_path = pairwise_input[[i]]$hdf5_path,
                                            conda = conda,
                                            miniprot_bin = miniprot_bin,
                                            miniprot_condaenv = miniprot_condaenv,
                                            miniprot_out_dir = pairwise_input[[i]]$miniprot_out_dir,
                                            diamond_out_dir = pairwise_input[[i]]$diamond_out_dir,
                                            diamond_exec_path = diamond_exec_path,
                                            n_threads = n_threads,
                                            overwrite = overwrite,
                                            use_prot = use_prot,
                                            resume = resume,
                                            module = module)},
                      silent = TRUE)
        if(inherits(x = object, what = "try-error")){
            hdf5_fn <- c(hdf5_fn, object)
            print(object)
            message("\n")
            error <- TRUE
            
        } else {
            hdf5_fn <- c(hdf5_fn, pairwise_input[[i]]$hdf5_path)
        }
    }
    names(hdf5_fn) <- names(pairwise_input)
    
    if(error){
        message("Error occured in process(es).", 
                "\nStop the pipeline and output the file paths for no-error",
                " processes and error messages for error processes.")
        on.exit()
        return(hdf5_fn)
    }
    
    if(reorg){
        message("Start reorganizing orthopairs ...")
        hdf5_fn <- reorgOrthopiars(hdf5_fn = hdf5_fn,
                                   out_dir = file.path(working_dir, "reorg_out"),
                                   out_fn = "reorg_orthopair.h5",
                                   makeFASTA = TRUE,
                                   overwrite = overwrite)
        
    }
    
    if(makegraph){
        message("Summarizing orthopairs into graphs ...")
        hdf5_fn <- file.path(working_dir, "reorg_out", "reorg_orthopair.h5")
        graph <- makeOrthoGraph(hdf5_fn = hdf5_fn)
        graph_fn <- file.path(working_dir, "orthopair.graphml")
        write_graph(graph = graph, file = graph_fn, format = "graphml")
        
    } else {
        graph_fn <- file.path(working_dir, "orthopair.graphml")
        
        if(output_table){
            if(file.exists(graph_fn)){
                hdf5_fn <- file.path(working_dir, "reorg_out", "reorg_orthopair.h5")
                graph <- read_graph(file = graph_fn, format = "graphml")
                
            } else {
                message("You set 'makegraph = FALSE' and 'output_table = TRUE'.", 
                        "\nHowever, the graph file ", graph_fn, 
                        ", which is required to create a table output, is not found",
                        "The table creation step will be skipped.")
                output_table <- FALSE
            }
        } else {
            graph_fn <- NULL
        }
    }
    
    if(output_table){
        message("Summarizing orthopairs into a spreadsheat ...")
        orthopair_fn <- file.path(working_dir, "orthopair_list.csv")
        unlink(x = orthopair_fn, force = TRUE)
        graph2df(hdf5_fn = hdf5_fn, graph = graph, orthopair_fn = orthopair_fn)
        
        orphan <- getOrphan(hdf5_fn = hdf5_fn)
        orphan <- lapply(seq_along(orphan), function(i){
            i_orphan <- orphan[[i]]
            i_out <- matrix(data = NA, nrow = length(i_orphan), ncol = length(orphan))
            i_out[, i] <- i_orphan
            i_out <- data.frame(i_out)
            names(i_out) <- names(orphan)
            return(i_out)
        })
        orphan <- do.call("rbind", orphan)
        orphan_fn <- file.path(working_dir, "orphan_list.csv")
        write.csv(x = orphan, file = orphan_fn, row.names = FALSE)
    } else {
        orthopair_fn <- orphan_fn <- NULL
    }
    
    on.exit()
    out <- c(hdf5_fn = hdf5_fn, graph_fn = graph_fn, 
             orthopair_fn = orthopair_fn, orphan_fn = orphan_fn)
    
    message("Finished all processes in the pipeline!")
    return(out)
}

.prepPairs <- function(in_list,
                       hdf5_out_dir,
                       working_dir,
                       miniprot_out_dir,
                       diamond_out_dir,
                       target_pair = NULL){
    in_list <- as.data.frame(x = in_list)
    combs <- combn(x = seq_len(nrow(in_list)), m = 2)
    comb_names <- matrix(data = in_list$name[combs], nrow = 2)
    comb_id <- apply(X = comb_names, MARGIN = 2, FUN = paste, collapse = "_")
    if(!is.null(target_pair)){
        target_pair <- apply(X = comb_names, 
                             MARGIN = 2,
                             FUN = function(x){
                                 any(apply(X = target_pair,
                                           MARGIN = 1, 
                                           FUN = function(y){
                                               return(all(x %in% y))
                                           }))
                             })
        comb_id <- comb_id[target_pair]
        combs <- combs[, target_pair]
    }
    
    out <- NULL
    for(i in seq_along(comb_id)){
        prefix <- comb_id[i]
        input_dir <- file.path(working_dir, "input", prefix)
        dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
        
        fn_list <- list(query_genome = in_list$genome[combs[1, i]],
                        subject_genome = in_list$genome[combs[2, i]],
                        query_gff = in_list$gff[combs[1, i]],
                        subject_gff = in_list$gff[combs[2, i]],
                        query_cds = in_list$cds[combs[1, i]],
                        subject_cds = in_list$cds[combs[2, i]],
                        query_prot = in_list$prot[combs[1, i]],
                        subject_prot = in_list$prot[combs[2, i]])
        
        fn_list <- .prepInputFiles(fn_list = fn_list, input_dir = input_dir)
        
        fn_list$hdf5_path <- file.path(hdf5_out_dir, paste0(prefix, ".h5"))
        fn_list$miniprot_out_dir <- file.path(miniprot_out_dir, prefix)
        fn_list$diamond_out_dir <- file.path(diamond_out_dir, prefix)
        out <- c(out, list(fn_list))
    }
    names(out) <- comb_id
    return(out)
}

.prepInputFiles <- function(fn_list, input_dir){
    for(i in seq_along(fn_list)){
        if(is.null(fn_list[[i]])){
            next
        }
        target <- names(fn_list)[i]
        if(!is.null(fn_list[[i]])){
            if(grepl("gff", target)){
                new_fn <- paste0(target, ".gff")
                
            } else {
                new_fn <- paste0(target, ".fa")
            }
            new_path <- file.path(input_dir, new_fn)
            unlink(new_path, force = TRUE)
            full_path <- sub("\\./", "", fn_list[[i]])
            if(!grepl("^\\/|^~\\/|^[^\\/]:", full_path)){
                full_path <- file.path(getwd(), full_path)
            }
            
            check <- file.symlink(full_path, new_path)
            if(check){
                fn_list[[i]] <- new_path
                
            } else {
                
                
                stop(paste("Faild to create a symlink to ",
                           full_path), 
                     call. = FALSE)
            }
        }
    }
    return(fn_list)
}

#' @importFrom parallel detectCores
.runOrthoPair <- function(query_genome = "",
                          subject_genome = "",
                          query_gff = "",
                          subject_gff = "",
                          query_cds = "",
                          subject_cds = "",
                          query_prot = "",
                          subject_prot = "",
                          hdf5_path = "./orthopair.h5",
                          conda = "conda",
                          miniprot_bin = "miniprot",
                          miniprot_condaenv = "miniprot",
                          miniprot_out_dir = "",
                          diamond_out_dir = "",
                          diamond_exec_path = NULL,
                          n_threads = NULL,
                          verbose = TRUE,
                          overwrite = FALSE,
                          use_prot = FALSE,
                          resume = FALSE,
                          module = NULL){
    
    if(is.null(n_threads)){
        core <- detectCores()
        n_threads <- core - 1
    }
    
    if(resume){
        if(overwrite){
            stop("Set 'overwrite = FALSE', if you wanna resume the process.", 
                 call. = FALSE)
        }
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
                              miniprot_out_dir = miniprot_out_dir,
                              overwrite = overwrite,
                              resume = resume,
                              module = module)
    
    if(object$resume$miniprot){
        if(verbose){
            message("Gene modeling by Miniprot to find missing genes.")
        }
        mapProt(object = object,
                miniprot_bin = miniprot_bin,
                conda = conda,
                condaenv = miniprot_condaenv,
                n_threads = n_threads,
                out_dir = miniprot_out_dir,
                len_diff = 0.2)
        
    } else {
        if(verbose){
            message("skip gene modeling.")
        }
    }
    
    if(object$resume$blast){
        if(verbose){
            message("Performinig reciprocal BLAST.")
        }
        rbh(object = object,
            use_prot = use_prot,
            n_threads = n_threads, 
            diamond_exec_path = diamond_exec_path,
            diamond_out_dir = diamond_out_dir)
        
    } else {
        if(verbose){
            message("skip performinig reciprocal BLAST.")
        }
    }
    
    if(object$resume$pairing){
        if(verbose){
            message("Pairing orthologs.")
        }
        syntenicOrtho(object = object)
        
    } else {
        if(verbose){
            message("skip pairing orthologs.")
        }
    }
    
    message("Finish processing.\n")
    return(object)
}

#' @importFrom rhdf5 H5Lexists
.checkResumePoint <- function(hdf5_path, files, resume, module, no_genome, no_prot){
    out <- list(miniprot = TRUE,
                blast = TRUE,
                pairing = TRUE,
                set_mp = FALSE)
    h5 <- H5Fopen(hdf5_path)
    # Ensure the HDF5 file is closed when the function exits
    on.exit(H5Fclose(h5))
    if(!is.null(module)){
        check <- names(module) %in% names(out)
        if(!all(check)){
            stop("The module object must be a named list with the following names:",
                 "\n'miniprot', 'blast', 'pairing'.")
        }
    }
    
    if(!is.null(module$miniprot)){
        out$miniprot <- module$miniprot
        if(out$miniprot){
            module$pairing <- module$blast <- TRUE
        }
    }
    if(!is.null(module$blast)){
        out$blast <- module$blast
        if(out$blast){
            module$pairing <- TRUE
        }
    }
    if(!is.null(module$pairing)){
        out$pairing <- module$pairing
    }
    
    if(resume){
        if(H5Lexists(h5, "timestamp/miniprot")){
            out$miniprot <- FALSE
            out$set_mp <- TRUE
        }
        if(H5Lexists(h5, "timestamp/blast")){
            out$blast <- FALSE
        }
        if(H5Lexists(h5, "timestamp/pairing")){
            out$pairing <- FALSE
        }
    }
    
    if(no_genome | no_prot){
        out$miniprot <- FALSE
    }
    
    return(out)
}

#' @export
orgInputFiles <- function(object = NULL, name, genome = NULL, gff, cds, prot = NULL){
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
        name_hit <- object$name %in% name
        if(any(name_hit)){
            message(name, " exsits in the input file list.")
            while(TRUE){
                check <- readline(prompt = "Overwrite?(y/n): ")
                if(check == "y"){
                    break
                    
                } else if(check == "n"){
                    stop("Use another name.", call. = FALSE)
                }
            }
            object$genome[name_hit] <- genome
            object$gff[name_hit] <- gff
            object$cds[name_hit] <- cds
            object$prot[name_hit] <- prot
            
        } else {
            object$name <- c(object$name, name)
            object$genome <- c(object$genome, genome)
            object$gff <- c(object$gff, gff)
            object$cds <- c(object$cds, cds)
            object$prot <- c(object$prot, prot)
        }
    }
    return(object)
}

#' @importFrom GenomeInfoDb seqlevels
#' @importFrom Biostrings readDNAStringSet
#' @importFrom rtracklayer import.gff
.validateInput <- function(name, genome, gff, cds, prot){
    check <- is.na(name) | name == ""
    if(check){
        stop("Empty name is not allowed.", call. = FALSE)
    }
    
    if(!is.null(genome)){
        check <- !file.exists(genome)
        if(check){
            stop('The file "genome" does not exist.', call. = FALSE)
        }
    }
    
    check <- !file.exists(gff)
    if(check){
        stop('The file "gff" does not exist.', call. = FALSE)
    }
    
    check <- !file.exists(cds)
    if(check){
        stop('The file "cds" does not exist.', call. = FALSE)
    }
    
    if(!is.null(prot)){
        check <- !file.exists(prot)
        if(check){
            stop('The file "prot" does not exist.', call. = FALSE)
        }
    }
    
    if(!is.null(genome)){
        genome <- readDNAStringSet(filepath = genome)
        genome_seq_name <- names(genome)
    }
    gff <- import.gff(con = gff)
    cds <- readDNAStringSet(filepath = cds)
    if(!is.null(prot)){
        prot <- readAAStringSet(filepath = prot)
    }
    gff_seq_lev <- seqlevels(gff)
    
    check <- .checkGFFstr(gff = gff)
    if(check$check){
        stop(paste0("In input data validation for ", name), "\n",
             check$msg,
             call. = FALSE)
    }
    
    check <- is.null(gff$gene_id)
    if(check){
        stop(paste0("In input data validation for ", name),
             "\nThe specified GFF does not have the 'gene_id' column.",
             call. = FALSE)
    }
    
    check <- any(is.na(gff$gene_id))
    if(check){
        stop(paste0("In input data validation for ", name),
             "\nThe specified GFF contains NA in the 'gene_id' column.",
             call. = FALSE)
    }
    
    if(!is.null(genome)){
        check <-  gff_seq_lev %in% genome_seq_name
        if(!all(check)){
            stop(paste0("In input data validation for ", name),
                 "\nFollowing sequence names appeared in the GFF but not in the genome:\n",
                 paste(head(gff_seq_lev[!gff_seq_lev %in% genome_seq_name]),
                       collapse = ", "), call. = FALSE)
        }
    }
    cds_names <- names(cds)
    check <- cds_names %in% gff$ID
    if(!all(check)){
        stop(paste0("In input data validation for ", name),
             "\nFollowing sequence names appeared in the CDS but not in the GFF:\n",
             paste(head(cds_names[!cds_names %in% gff$ID]),
                   collapse = ", "), call. = FALSE)
    }
    
    if(!is.null(prot)){
        prot_names <- names(prot)
        check <- prot_names %in% gff$ID
        if(!all(check)){
            stop(paste0("In input data validation for ", name),
                 "\nFollowing sequence names appeared in the protein but not in the GFF:\n",
                 paste(head(prot_names[!prot_names %in% gff$ID]),
                       collapse = ", "), call. = FALSE)
        }
    }
}

.checkGFFstr <- function(gff){
    check <- any(is.na(gff$type)) | any(gff$type == "")
    if(check){
        out <- list(check = TRUE, msg = "NA or missing values in the type column.")
        return(out)
    }
    
    gene <- gff[gff$type %in% "gene"]
    tx <- gff[gff$type %in% c("mRNA", "transcript")]
    tx_p <- unlist(tx$Parent)
    check <- !all(tx_p %in% gene$ID)
    if(check){
        out <- list(check = TRUE, 
                    msg = "The Parent column of the transcripts contains values that do not match the values in the ID column of the genes.")
        return(out)
    }
    check <- !all(gene$ID %in% tx_p)
    if(check){
        out <- list(check = TRUE, 
                    msg = "The ID column of the genes contains values that do not match the values in the Parent column of the transcripts.")
        return(out)
    }
    
    out <- list(check = FALSE)
    return(out)
}
