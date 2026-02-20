#' Make input object and validate working directory
#'
#' This function organize input files and structure working directory with
#' confirming previous outputs
#'
#' @export
#'
init_opr <- function(object, working_dir, overwrite, n_threads){
    input <- .checkInput(object, working_dir)
    input <- .assignIndex(input = input, 
                          working_dir = working_dir,
                          overwrite = overwrite)
    input <- .orgInput(input = input,
                       working_dir = working_dir,
                       overwrite = overwrite,
                       n_threads = n_threads,
                       verbose = verbose)
    out <- c(input = list(input),
             working_dir = list(working_dir))
    class(out) <- c(class(out), "initOPR")
    return(out)
}

.checkInput <- function(input, working_dir){
    if(is.character(input)){
        input <- read.csv(file = input, stringsAsFactors = FALSE)
    }
    if(!is.data.frame(input) & !is.list(input)){
        stop("Input must be a path to a CSV file, a data.frame, or a list.",
             call. = FALSE)
    }
    input <- as.data.frame(input)
    .checkFilesExist(input = input)
    
    n_genome <- nrow(input)
    n_genome_limit <- 8999
    if(n_genome > n_genome_limit){
        stop("Max number of genomes to be processed must be less than ", 
             n_genome_limit,
             call. = FALSE)
    }
    
    dir.create(working_dir, showWarnings = FALSE, recursive = TRUE)
    if(!dir.exists(working_dir)){
        stop("Failed to create the working directory at ", working_dir, 
             call. = FALSE)
    }
    return(input)
}

.checkFilesExist <- function(input){
    sapply(seq_len(nrow(input)), function(i){
        if(!is.null(input$genome[i]) & !is.na(input$genome[i])){
            if(!file.exists(input$genome[i])){
                stop(input$genome[i], " does not exist!", 
                     call. = FALSE)
            }
        }
        if(!is.null(input$gff[i]) & !is.na(input$gff[i])){
            if(!file.exists(input$gff[i])){
                stop(input$gff[i], " does not exist!", 
                     call. = FALSE)
            }
        }
        if(!is.null(input$cds[i]) & !is.na(input$cds[i])){
            if(!file.exists(input$cds[i])){
                stop(input$cds[i], " does not exist!", 
                     call. = FALSE)
            }
        }
        if(!is.null(input$prot[i]) & !is.na(input$prot[i])){
            if(!file.exists(input$prot[i])){
                stop(input$prot[i], " does not exist!", 
                     call. = FALSE)
            }
        }
    })
}

.assignIndex <- function(input, working_dir, overwrite){
    input_dir <- file.path(working_dir, "input")
    
    if(dir.exists(input_dir) & overwrite){
        base_name <- list.files(input_dir)
        base_name[!grepl("_", base_name)] <- NA
        index_name <- do.call(rbind, strsplit(base_name, "_"))
        hit <- match(input$name, index_name[, 2])
        index <- as.integer(index_name[hit, 1])
        index_offset <- max(index, na.rm = TRUE)
        index[is.na(index)] <- seq_len(sum(is.na(index))) + index_offset
        
    } else {
        index <- seq_along(input$name)
    }
    
    input$index <- as.integer(index)
    return(input)
}

.orgInput <- function(input, 
                      working_dir,
                      overwrite,
                      n_threads,
                      verbose){
    n_threads <- .getThreadsLimit(n_threads = n_threads, min_ma = 2e09)
    
    input_list <- mclapply(X = seq_along(input$name),
                           mc.cores = n_threads, 
                           FUN = .orgEngine,
                           input = input,
                           working_dir = working_dir,
                           overwrite = overwrite)
    input_list <- do.call(rbind, input_list)
    input$gff <- input_list$gff_fn
    input$gff_df <- input_list$gff_df_fn
    input$cds <- input_list$cds_fn
    input$prot <- input_list$prot_fn
    return(input)
}


.orgEngine <- function(i, 
                       input, 
                       working_dir,
                       overwrite){
    out_dir_i <- file.path(working_dir, "input", 
                           paste(input$index[i], input$name[i], sep = "_"))
    dir.create(out_dir_i, showWarnings = FALSE, recursive = TRUE)
    
    org_gff <- .orgGFF(gff = input$gff[i],
                       out_dir_i = out_dir_i,
                       index = i,
                       overwrite = overwrite)
    
    org_cds <- .orgCDS(org_gff = org_gff,
                       genome = input$genome[i],
                       cds = input$cds[i],
                       out_dir_i = out_dir_i,
                       index = i,
                       overwrite = overwrite)
    
    out <- data.frame(gff_fn = file.path(out_dir_i, "gff.rds"),
                      gff_df_fn = file.path(out_dir_i, "gff_df.rds"),
                      cds_fn = file.path(out_dir_i, "cds.fa"))
    return(out)
}

.getThreadsLimit <- function(n_threads, min_ma){
    ma <- .mem_available_bytes()
    ma_par_cpu <- ma / n_threads
    if(ma_par_cpu < min_ma){
        n_threads <- floor(ma / min_ma)
        if(n_threads < 1){
            n_threads <- 1
        }
    }
    return(n_threads)
}

.orgGFF <- function(gff, out_dir_i, index, overwrite){
    gff_fn <- file.path(out_dir_i, "gff.rds")
    gff_df_fn <- file.path(out_dir_i, "gff_df.rds")
    
    if(!overwrite){
        if(file.exists(gff_df_fn)){
            gff <- readRDS(file = gff_fn)
            gff_df <- readRDS(file = gff_df_fn)
            out <- list(gff = gff, gff_df = gff_df)
            return(out)
        }
    }
    
    entry_type <-  c("gene", "transcript", "mRNA", "CDS")
    gff <- import.gff3(gff, feature.type = entry_type)
    gff <- subset(gff, select = c(source:ID, Parent))
    gff$Name <- gff$ID
    gff <- .setGeneID(gff = gff)
    gff <- .fixGFF(gff = gff)
    gff <- .orderGFF(gff = gff)
    gff <- .setIndex(gff = gff, index = index)
    
    gff_df <- as.data.table(gff)
    cds_df <- gff_df[type %in% "CDS"][, .(seqnames, start, end, gene_id)]
    gene_df <- gff_df[type %in% "gene", .(seqnames, start, end, strand, gene_id, gene_index)]
    tx_df <- gff_df[type %in% c("transcript", "mRNA")][, .(transcript_id = ID, gene_id, tx_index)]
    gff_df <- merge(tx_df, gene_df, by = "gene_id", all.x = TRUE, allow.cartesian = TRUE)
    saveRDS(object = list(gff_df, cds_df), file = gff_df_fn)
    out <- list(gff = gff, gff_df = gff_df)
    return(out)
}

.setIndex <- function(gff, index){
    tx_i <- gff$type %in% c("mRNA", "transcript")
    gene_i <- gff$type %in% "gene"
    gff$gene_index <- gff$tx_index <- NA
    n_tx <- sum(tx_i)
    n_gene <- sum(gene_i)
    n_index_limit <- 999999
    if(n_tx > n_index_limit){
        stop("Max number of transcript models per chromosome must be less than ", 
             n_index_limit,
             call. = FALSE)
    }
    gff$tx_index[tx_i] <- (index + 1000L) * 1e6L + seq_len(n_tx)
    gff$gene_index[gene_i] <- (index + 1000L) * 1e6L + seq_len(n_gene)
    hit <- match(gff$gene_id[tx_i], gff$gene_id[gene_i])
    gff$gene_index[tx_i] <- gff$gene_index[gene_i][hit]
    gff$Parent[gene_i] <- ""
    gff$Parent <- unlist(gff$Parent)
    hit <- match(gff$Parent[!{gene_i | tx_i}], gff$ID[tx_i])
    gff$tx_index[!{gene_i | tx_i}] <- gff$tx_index[tx_i][hit]
    gff$gene_index[!{gene_i | tx_i}] <- gff$gene_index[tx_i][hit]
    return(gff)
}

.orgCDS <- function(org_gff, genome, cds, out_dir_i, index, overwrite){
    cds_fn <- file.path(out_dir_i, "cds.fa")
    
    if(!overwrite){
        if(file.exists(taxid_map_fn)){
            cds <- readDNAStringSet(cds_fn)
            rep_tx_id <- read.table(file = taxid_map_fn, 
                                    header = FALSE,
                                    sep = "\t")
            out <- list(cds = cds, rep_tx_id = rep_tx_id)
            return(out)
        }
    }
    
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
    hit <- match(names(rep_cds), org_gff$gff_df$transcript_id)
    names(rep_cds) <- org_gff$gff_df$tx_index[hit]
    writeXStringSet(rep_cds, cds_fn)
    invisible(rep_cds)
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
