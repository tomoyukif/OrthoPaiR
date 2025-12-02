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
                      miniprot_path = "",
                      blast_path = "",
                      diamond_path = "",
                      target_pair = NULL,
                      use_prot = FALSE,
                      run_miniprot = FALSE,
                      orthopair = TRUE,
                      reorg = TRUE,
                      makegraph = TRUE,
                      output_table = FALSE,
                      n_threads = NULL,
                      overwrite = FALSE,
                      verbose = TRUE){
    if(!inherits(x = in_list, what = "OrthoPairInput")){
        stop("The input object must be a OrthoPairInput class object.",
             call. = FALSE)
    }
    
    if(is.null(n_threads)){
        core <- detectCores()
        n_threads <- core - 1
    }
    hdf5_out_dir <- file.path(working_dir, "hdf5_out")
    dir.create(path = hdf5_out_dir, showWarnings = FALSE, recursive = TRUE)
    input_dir <- file.path(working_dir, "input")
    dir.create(path = input_dir, showWarnings = FALSE, recursive = TRUE)
    
    if(run_miniprot){
        if(verbose){
            message("Start missing ORF complementation...")
        }
        in_list <- mapProt(in_list = in_list,
                           out_dir = input_dir, 
                           overwrite = overwrite,
                           miniprot_path = miniprot_path,
                           n_threads = n_threads)
    }
    
    if(verbose){
        message("Start organizing input data...")
    }
    in_list <- .orgInput(in_list = in_list, 
                         out_dir = input_dir,
                         use_prot = use_prot,
                         blast_path = blast_path,
                         diamond_path = diamond_path,
                         overwrite = overwrite)
    
    pairwise_input <- .prepPairs(in_list = in_list,
                                 hdf5_out_dir = hdf5_out_dir,
                                 working_dir = working_dir,
                                 target_pair = target_pair)
    hdf5_fn <- NULL
    on.exit({return(hdf5_fn)})
    error <- FALSE
    for(i in seq_along(pairwise_input)){
        if(!orthopair){
            hdf5_fn <- c(hdf5_fn, pairwise_input[[i]]$hdf5_path)
            next
        }
        if(verbose){
            message("Start proccessing the ", names(pairwise_input)[i], " pair")
        }
        
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
        object <- try(expr = {.runOrthoPair(input_files = pairwise_input[[i]],
                                            use_prot = use_prot,
                                            blast_path = blast_path,
                                            diamond_path = diamond_path,
                                            overwrite = overwrite,
                                            n_threads = n_threads,
                                            verbose = verbose)},
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
                "\nStop the pipeline and output the file paths for non-error",
                " processes and error messages for error processes.")
        on.exit()
        return(hdf5_fn)
    }
    
    if(reorg){
        if(verbose){
            message("Start reorganizing orthopairs ...")
        }
        hdf5_fn <- reorgOrthopiars(hdf5_fn = hdf5_fn,
                                   out_dir = file.path(working_dir, "reorg_out"),
                                   out_fn = "reorg_orthopair.h5",
                                   makeFASTA = TRUE,
                                   overwrite = overwrite)
        
    }
    
    if(makegraph){
        if(verbose){
            message("Summarizing orthopairs into graphs ...")
        }
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
                if(verbose){
                    message("You set 'makegraph = FALSE' and 'output_table = TRUE'.", 
                            "\nHowever, the graph file ", graph_fn, 
                            ", which is required to create a table output, is not found",
                            "The table creation step will be skipped.")
                }
                output_table <- FALSE
            }
        } else {
            graph_fn <- NULL
        }
    }
    
    if(output_table){
        if(verbose){
            message("Summarizing orthopairs into a spreadsheat ...")
        }
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
    
    if(verbose){
        message("Finished all processes in the pipeline!")
    }
    return(out)
}

.orgInput <- function(in_list, out_dir, use_prot, blast_path, diamond_path, overwrite){
    in_list$h5 <- rep(NA, length(in_list$name))
    
    for(i in seq_along(in_list$name)){
        out_dir_i <- file.path(out_dir, in_list$name[i])
        dir.create(out_dir_i, showWarnings = FALSE, recursive = TRUE)
        h5_fn <- .makeHDF5(hdf5_path = file.path(out_dir_i, "input.h5"), 
                           overwrite = overwrite, verbose = FALSE)
        in_list$h5[i] <- h5_fn
        new_cds_fn <- file.path(out_dir_i, "cds.fa")
        
        if(overwrite){
            gff_fn <-  in_list$gff[i]
            gff <- import.gff3(gff_fn)
            gff_df <- as.data.frame(gff)
            gff_df <- subset(gff_df, select = -c(source, score, phase, Name, width))
            gff_df$type <- as.character(gff_df$type)
            gff_df$seqnames <- as.character(gff_df$seqnames)
            gff_df <- gff_df[order(gff_df$seqnames, gff_df$start), ]
            tx_i <- gff_df$type %in% c("mRNA", "transcript")
            gene_i <- gff_df$type %in% "gene"
            gff_df$gene_index <- gff_df$tx_index <- NA
            gff_df$tx_index[tx_i] <- seq_len(sum(tx_i))
            gff_df$gene_index[gene_i] <- seq_len(sum(gene_i))
            hit <- match(gff_df$gene_id[tx_i], gff_df$gene_id[gene_i])
            gff_df$gene_index[tx_i] <- gff_df$gene_index[gene_i][hit]
            gff_df$Parent[gene_i] <- ""
            gff_df$Parent <- unlist(gff_df$Parent)
            .h5overwrite(obj = gff_df, file = h5_fn, "gff")
        }
        
        new_cds_fn <- file.path(out_dir_i, "cds.fa")
        if(overwrite){
            cds <- readDNAStringSet(in_list$cds[i])
            rep_tx_id <- .getRepTx(gff_cds = gff[gff$type == "CDS"])
            rep_cds <- cds[names(cds) %in% rep_tx_id]
            hit <- match(names(rep_cds), gff_df$ID)
            names(rep_cds) <- gff_df$tx_index[hit]
            writeXStringSet(rep_cds, new_cds_fn)
            rep_cds <- data.frame(id = names(rep_cds), seq = as.character(rep_cds))
            rownames(rep_cds) <- NULL
            .h5overwrite(obj = rep_cds, file = h5_fn, "cds")
            .makeblastdb(blast_path = blast_path, fn = new_cds_fn)
        }
        .h5overwrite(obj = new_cds_fn, file = h5_fn, "cds_fn")
        
        if(!is.na(in_list$prot[i]) & use_prot){
            new_prot_fn <- file.path(out_dir_i, "prot.fa")
            if(overwrite){
                prot <- readDNAStringSet(in_list$prot[i])
                rep_prot <- prot[names(prot) %in% rep_tx_id]
                hit <- match(names(rep_prot), gff_df$ID)
                names(rep_prot) <- gff_df$tx_index[hit]
                writeXStringSet(rep_prot, new_prot_fn)
                .makediamonddb(diamond_path = diamond_path, fn = new_prot_fn)
            }
            .h5overwrite(obj = new_prot_fn, file = h5_fn, "prot_fn")
        }
    }
    return(in_list)
}

.getRepTx <- function(gff_cds){
    gff_cds$Parent <- unlist(gff_cds$Parent)
    exon_edges <- paste(start(gff_cds), end(gff_cds), sep = "-")
    exon_str <- tapply(exon_edges, gff_cds$Parent, paste, collapse = "_")
    exon_chr <- tapply(as.character(seqnames(gff_cds)), gff_cds$Parent, unique)
    exon_str <- paste(exon_chr, exon_str, sep = "_")
    rep_tx <- names(exon_chr)[!duplicated(exon_str)]
    return(rep_tx)
}

.makeblastdb <- function(blast_path, fn){
    blast_args <- paste(paste("-in", fn),
                        paste("-dbtype nucl"))
    check <- try({
        system2(command = file.path(blast_path, "makeblastdb"),
                args = blast_args, 
                stdout = FALSE)
    }, silent = TRUE)
    if(inherits(check, "try-error")){
        stop("Error in makeblastdb of BLAST.\n", check)
    }
}

.makediamonddb <- function(diamond_path, fn){
    diamond_args <- paste("makedb",
                          paste("--in", fn),
                          paste("--db", paste0(fn, ".dmnd")))
    check <- try({
        system2(command = file.path(diamond_path, "diamond"),
                args = diamond_args, 
                stdout = FALSE)
    }, silent = TRUE)
    if(inherits(check, "try-error")){
        stop("Error in makedb of DIAMOND.\n", check)
    }
}

.prepPairs <- function(in_list,
                       hdf5_out_dir,
                       working_dir,
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
        # input_dir <- file.path(working_dir, "input", prefix)
        # dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
        
        fn_list <- list(pair_id = prefix,
                        query_genome = in_list$genome[combs[1, i]],
                        subject_genome = in_list$genome[combs[2, i]],
                        query_gff = in_list$gff[combs[1, i]],
                        subject_gff = in_list$gff[combs[2, i]],
                        query_cds = in_list$cds[combs[1, i]],
                        subject_cds = in_list$cds[combs[2, i]],
                        query_prot = in_list$prot[combs[1, i]],
                        subject_prot = in_list$prot[combs[2, i]],
                        query_h5 = in_list$h5[combs[1, i]],
                        subject_h5 = in_list$h5[combs[2, i]])
        
        # fn_list <- .prepInputFiles(fn_list = fn_list, input_dir = input_dir)
        
        fn_list$hdf5_path <- file.path(hdf5_out_dir, paste0(prefix, ".h5"))
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
        if(!is.na(fn_list[[i]])){
            if(grepl("gff", target)){
                new_fn <- paste0(target, ".gff")
                
            } else if(grepl("h5", target)){
                new_fn <- paste0(target, ".h5")
                
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
.runOrthoPair <- function(input_files,
                          use_prot = FALSE,
                          blast_path = "",
                          diamond_path = "",
                          overwrite = FALSE,
                          n_threads = NULL,
                          verbose = TRUE){
    
    if(is.null(n_threads)){
        core <- detectCores()
        n_threads <- core - 1
    }
    
    if(verbose){
        message("Preparing a OrthoPairDB class object.")
    }
    
    object <- makeOrthoPairDB(input_files, 
                              overwrite = overwrite, 
                              verbose = verbose)
    
    if(object$resume$blast){
        if(verbose){
            message("Performinig reciprocal BLAST.")
        }
        rbh(object = object,
            use_prot = use_prot,
            n_threads = n_threads, 
            blast_path = blast_path,
            diamond_path = diamond_path)
        
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
.checkResumePoint <- function(hdf5_path, files){
    out <- list(blast = TRUE, pairing = TRUE)
    h5 <- H5Fopen(hdf5_path)
    on.exit(H5Fclose(h5))
    
    if(H5Lexists(h5, "timestamp/blast")){
        out$blast <- FALSE
    }
    if(H5Lexists(h5, "timestamp/pairing")){
        out$pairing <- FALSE
    }
    return(out)
}

#' @export
orgInputFiles <- function(object = NULL, name, genome = NA, gff, cds, prot = NA){
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
        class(object) <- c(class(object), "OrthoPairInput")
        
    } else {
        if(!inherits(x = object, what = "OrthoPairInput")){
            stop("The input object must be a OrthoPairInput class object.",
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
    
    if(all(!is.na(genome))){
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
    
    if(all(!is.na(prot))){
        check <- !file.exists(prot)
        if(check){
            stop('The file "prot" does not exist.', call. = FALSE)
        }
    }
    
    if(any(!is.na(genome))){
        genome <- readDNAStringSet(filepath = genome)
        genome_seq_name <- names(genome)
    }
    gff <- import.gff(con = gff)
    cds <- readDNAStringSet(filepath = cds)
    if(!is.na(prot)){
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
    
    if(all(!is.na(genome))){
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
    
    if(all(!is.na(prot))){
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
