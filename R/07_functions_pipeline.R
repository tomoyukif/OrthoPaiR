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
                      rename = TRUE,
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
    
    if(verbose){
        message("Start organizing input data...")
    }
    in_list <- .orgInput(in_list = in_list,
                         input_dir = input_dir,
                         use_prot = use_prot,
                         run_miniprot = run_miniprot,
                         blast_path = blast_path,
                         diamond_path = diamond_path,
                         overwrite = overwrite,
                         n_threads = n_threads,
                         verbose = verbose)
    saveRDS(object = in_list, file = file.path(input_dir, "in_list.rds"))
    
    if(run_miniprot){
        if(verbose){
            message("Start missing ORF complementation...")
        }
        in_list <- mapProt(in_list = in_list,
                           input_dir = input_dir, 
                           overwrite = overwrite,
                           miniprot_path = miniprot_path,
                           blast_path = blast_path, 
                           diamond_path = diamond_path, 
                           use_prot = use_prot,
                           n_threads = n_threads,
                           verbose = verbose)
    }
    saveRDS(object = in_list, file = file.path(input_dir, "in_list.rds"))
    
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
        reorg_hdf5_fn <- reorgOrthopiars(hdf5_fn = hdf5_fn,
                                   out_dir = file.path(working_dir, "reorg_out"),
                                   out_fn = "reorg_orthopair.h5",
                                   rename = rename,
                                   overwrite = overwrite,
                                   verbose = verbose)
    }
    
    if(makegraph){
        if(verbose){
            message("Summarizing orthopairs into graphs ...")
        }
        reorg_hdf5_fn <- file.path(working_dir, "reorg_out", "reorg_orthopair.h5")
        graph <- makeOrthoGraph(hdf5_fn = reorg_hdf5_fn)
        graph_fn <- file.path(working_dir, "orthopair.graphml")
        write_graph(graph = graph, file = graph_fn, format = "graphml")
        
    } else {
        graph_fn <- file.path(working_dir, "orthopair.graphml")
        
        if(output_table){
            if(file.exists(graph_fn)){
                reorg_hdf5_fn <- file.path(working_dir, "reorg_out", "reorg_orthopair.h5")
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
        graph2df(hdf5_fn = reorg_hdf5_fn, graph = graph, orthopair_fn = orthopair_fn)
        
        orphan <- getOrphan(hdf5_fn = reorg_hdf5_fn)
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
    out <- c(hdf5_fn = reorg_hdf5_fn, graph_fn = graph_fn, 
             orthopair_fn = orthopair_fn, orphan_fn = orphan_fn)
    
    if(verbose){
        message("Finished all processes in the pipeline!")
    }
    return(out)
}

.orgInput <- function(in_list, 
                      input_dir,
                      use_prot, 
                      run_miniprot,
                      blast_path,
                      diamond_path, 
                      overwrite,
                      n_threads,
                      verbose){
    ma <- .mem_available_bytes()
    ma_par_cpu <- ma / n_threads
    if(ma_par_cpu < 3e09){
        n_threads <- floor(ma / 3e09)
        if(n_threads < 1){
            n_threads <- 1
        }
    }
    
    out_list <- mclapply(X = seq_along(in_list$name),
                         mc.cores = n_threads, 
                         FUN = .orgEngine,
                         in_list = in_list,
                         out_dir = input_dir,
                         use_prot = use_prot, 
                         run_miniprot = run_miniprot,
                         blast_path = blast_path,
                         diamond_path = diamond_path, 
                         overwrite = overwrite)
    out_list <- do.call(rbind, out_list)
    in_list$gff <- out_list$gff_fn
    in_list$gff_df <- out_list$gff_df_fn
    in_list$cds <- out_list$cds_fn
    in_list$prot <- out_list$prot_fn
    return(in_list)
}

.orgEngine <- function(i, 
                       in_list, 
                       out_dir,
                       use_prot, 
                       run_miniprot,
                       blast_path,
                       diamond_path, 
                       overwrite){
    out_dir_i <- file.path(out_dir, in_list$name[i])
    dir.create(out_dir_i, showWarnings = FALSE, recursive = TRUE)
    
    if(overwrite){
        org_gff <- .orgGFF(gff_fn = in_list$gff[i], out_dir_i = out_dir_i)
    }
    
    if(overwrite){
        org_cds <- .orgCDS(org_gff = org_gff,
                           genome = in_list$genome[i],
                           cds = in_list$cds[i],
                           blast_path = blast_path,
                           out_dir_i = out_dir_i)
    }
    
    if(overwrite & {use_prot | run_miniprot}){
        .orgProt(org_cds = org_cds,
                 prot = in_list$prot[i],
                 gff_df = org_gff$gff_df,
                 use_prot = use_prot,
                 diamond_path = diamond_path,
                 out_dir_i = out_dir_i)
        
    }
    return(data.frame(gff_fn = file.path(out_dir_i, "gff.rds"),
                      gff_df_fn = file.path(out_dir_i, "gff_df.rds"),
                      cds_fn = file.path(out_dir_i, "cds.fa"),
                      prot_fn = file.path(out_dir_i, "prot.fa")))
}

.orgGFF <- function(gff_fn, out_dir_i){
    gff <- import.gff3(gff_fn)
    entry_type <-  c("gene", "transcript", "mRNA", "CDS")
    gff <- gff[gff$type %in% entry_type]
    gff$Name <- gff$ID
    gff <- .setGeneID(gff = gff)
    gff <- .orderGFF(gff = gff)
    gff <- .setIndex(gff = gff)
    
    gff_df <- as.data.frame(gff)
    gff_df <- subset(gff_df, select = c(seqnames, start, end,
                                        strand, type, ID, 
                                        Parent, gene_id,
                                        tx_index, gene_index))
    gff_fn <- file.path(out_dir_i, "gff.rds")
    saveRDS(object = gff, file = gff_fn)
    gff_df_fn <- file.path(out_dir_i, "gff_df.rds")
    saveRDS(object = gff_df, file = gff_df_fn)
    out <- list(gff = gff, gff_df = gff_df)
    return(out)
}

.setIndex <- function(gff){
    tx_i <- gff$type %in% c("mRNA", "transcript")
    gene_i <- gff$type %in% "gene"
    gff$gene_index <- gff$tx_index <- NA
    gff$tx_index[tx_i] <- seq_len(sum(tx_i))
    gff$gene_index[gene_i] <- seq_len(sum(gene_i))
    hit <- match(gff$gene_id[tx_i], gff$gene_id[gene_i])
    gff$gene_index[tx_i] <- gff$gene_index[gene_i][hit]
    gff$Parent[gene_i] <- ""
    gff$Parent <- unlist(gff$Parent)
    hit <- match(gff$Parent[!{gene_i | tx_i}], gff$ID[tx_i])
    gff$tx_index[!{gene_i | tx_i}] <- gff$tx_index[tx_i][hit]
    gff$gene_index[!{gene_i | tx_i}] <- gff$gene_index[tx_i][hit]
    return(gff)
}

.orgCDS <- function(org_gff, genome, cds, blast_path, out_dir_i){
    if(is.na(cds)){
        cds <- .makeCDS(gff = org_gff$gff, genome = genome)
        
    } else {
        cds <- readDNAStringSet(cds)
    }
    gff_cds <- org_gff$gff[org_gff$gff$type %in% "CDS"]
    check <- names(cds) %in% unlist(gff_cds$Parent)
    if(!all(check)){
        stop("The following transcript names does not appear in the GFF file:",
             "\n",
             paste(names(cds)[!check], collapse = ",\n"))
    }
    
    rep_tx_id <- .getRepTx(gff_cds = gff_cds)
    rep_cds <- cds[names(cds) %in% rep_tx_id]
    hit <- match(names(rep_cds), org_gff$gff_df$ID)
    names(rep_cds) <- org_gff$gff_df$tx_index[hit]
    cds_fn <- file.path(out_dir_i, "cds.fa")
    writeXStringSet(rep_cds, cds_fn)
    .makeblastdb(blast_path = blast_path, fn = cds_fn)
    out <- list(cds = cds, rep_tx_id = rep_tx_id)
    return(out)
}

.orgProt <- function(org_cds, prot, gff_df, use_prot, diamond_path, out_dir_i){
    if(is.na(prot)){
        prot <- translate(org_cds$cds, if.fuzzy.codon = "solve")
        
    } else {
        prot <- readDNAStringSet(prot)
    }
    
    check <- names(prot) %in% names(org_cds$cds)
    if(!all(check)){
        stop("The following protein names does not appear in the GFF file:",
             "\n",
             paste(names(prot)[!check], collapse = ",\n"))
    }
    
    rep_prot <- prot[names(prot) %in% org_cds$rep_tx_id]
    hit <- match(names(rep_prot), gff_df$ID)
    names(rep_prot) <- gff_df$tx_index[hit]
    prot_fn <- file.path(out_dir_i, "prot.fa")
    writeXStringSet(rep_prot, prot_fn)
    if(use_prot){
        .makediamonddb(diamond_path = diamond_path, fn = prot_fn)
    }
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
        fn_list <- list(pair_id = prefix,
                        query_genome = in_list$genome[combs[1, i]],
                        subject_genome = in_list$genome[combs[2, i]],
                        query_gff = in_list$gff[combs[1, i]],
                        subject_gff = in_list$gff[combs[2, i]],
                        query_gff_df = in_list$gff_df[combs[1, i]],
                        subject_gff_df = in_list$gff_df[combs[2, i]],
                        query_cds = in_list$cds[combs[1, i]],
                        subject_cds = in_list$cds[combs[2, i]],
                        query_prot = in_list$prot[combs[1, i]],
                        subject_prot = in_list$prot[combs[2, i]])
        
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
orgInputFiles <- function(object = NULL,
                          name,
                          genome = NULL,
                          gff,
                          cds = NULL, 
                          prot = NULL,
                          validation = FALSE){
    if(length(gff) > 1){
        message("Multiple input files were specified.")
        if(is.null(genome)){
            genome <- rep(NA, length(gff))
        }
        if(is.null(cds)){
            cds <- rep(NA, length(gff))
        }
        if(is.null(prot)){
            prot <- rep(NA, length(gff))
        }
        for(i in seq_along(gff)){
            message("Processing input ", i, ": ", name[i], "..." )
            object <- orgInputFiles(object = object,
                                    name = name[i],
                                    genome = genome[i], 
                                    gff = gff[i], 
                                    cds = cds[i],
                                    prot = prot[i],
                                    validation = validation)
        }
        return(object)
        
    } else {
        if(is.null(genome)){
            genome <- NA
        }
        if(is.null(cds)){
            cds <- NA
        }
        if(is.null(prot)){
            prot <- NA
        }
    }
    
    if(validation){
        .validateInput(name = name,
                       genome = genome,
                       gff = gff,
                       cds = cds,
                       prot = prot)
    }
    
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
    
    # Confirm the existence of the genome FASTA file
    if(!is.na(genome)){
        check <- !file.exists(genome)
        if(check){
            stop('The file "genome" does not exist.', call. = FALSE)
        }
    }
    
    # Confirm the existence of the annotation GFF file
    check <- !file.exists(gff)
    if(check){
        stop('The file "gff" does not exist.', call. = FALSE)
    }
    
    # Confirm the existence of the CDS FASTA file
    if(!is.na(cds)){
        check <- !file.exists(cds)
        if(check){
            stop('The file "cds" does not exist.', call. = FALSE)
        }
    }
    
    # Confirm the existence of the protein FASTA file
    if(!is.na(prot)){
        check <- !file.exists(prot)
        if(check){
            stop('The file "prot" does not exist.', call. = FALSE)
        }
    }
    
    gff <- import.gff(con = gff)
    check <- .checkGFFstr(gff = gff)
    if(check$check){
        stop(paste0("In input data validation for ", name), "\n",
             check$msg,
             call. = FALSE)
    }
    
    if(!is.na(genome)){
        gff_seq_lev <- seqlevels(gff)
        genome <- readDNAStringSet(filepath = genome)
        genome_seq_name <- names(genome)
        check <-  gff_seq_lev %in% genome_seq_name
        if(!all(check)){
            stop(paste0("In input data validation for ", name),
                 "\nFollowing sequence names appeared in the GFF but not in the genome:\n",
                 paste(head(gff_seq_lev[!gff_seq_lev %in% genome_seq_name]),
                       collapse = ", "), call. = FALSE)
        }
        
    } else {
        if(is.na(cds)){
            stop("the cds argument must be specified when the genome argument is missing.",
                 call. = FALSE)
        }
    }
    
    if(!is.na(cds)){
        cds <- readDNAStringSet(filepath = cds)
        cds_names <- names(cds)
        check <- cds_names %in% gff$ID
        if(!all(check)){
            stop(paste0("In input data validation for ", name),
                 "\nFollowing sequence names appeared in the CDS but not in the GFF:\n",
                 paste(head(cds_names[!cds_names %in% gff$ID]),
                       collapse = ", "), call. = FALSE)
        }
    }
    
    if(!is.na(prot)){
        prot <- readAAStringSet(filepath = prot)
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
