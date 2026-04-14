#' Make input object and validate working directory
#'
#' This function organize input files and structure working directory with
#' confirming previous outputs
#'
#' @export
#'
reformatFiles <- function(object, working_dir, overwrite, n_threads){
    input <- .checkInput(object, working_dir)
    input <- .assignIndex(input = input, 
                          working_dir = working_dir,
                          overwrite = overwrite)
    .orgInput(input = input,
              working_dir = working_dir,
              overwrite = overwrite,
              n_threads = n_threads,
              verbose = verbose)
    invisible(TRUE)
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
    cds_df <- gff_df[type %in% "CDS"][, .(seqnames, start, end, gene_id, transcript_id = unlist(Parent))]
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

.mem_available_bytes <- function() {
    os <- Sys.info()[["sysname"]]
    
    ## -------- Linux --------
    if (os == "Linux") {
        mi <- tryCatch(readLines("/proc/meminfo", warn = FALSE), error = function(e) NULL)
        if (!is.null(mi)) {
            ma <- grep("^MemAvailable:", mi, value = TRUE)
            if (length(ma)) {
                kb <- as.numeric(gsub("[^0-9]", "", ma))
                return(kb * 1024)
            }
        }
    }
    
    ## -------- macOS --------
    if (os == "Darwin") {
        # ページサイズ
        ps <- tryCatch(
            as.numeric(system("sysctl -n hw.pagesize", intern = TRUE)),
            error = function(e) NA_real_
        )
        
        vm <- tryCatch(system("vm_stat", intern = TRUE), error = function(e) NULL)
        if (!is.na(ps) && !is.null(vm)) {
            get <- function(name) {
                x <- grep(paste0("^", name), vm, value = TRUE)
                if (length(x)) as.numeric(gsub("[^0-9]", "", x)) else 0
            }
            free     <- get("Pages free")
            inactive <- get("Pages inactive")
            speculative <- get("Pages speculative")
            
            # macOS では inactive + speculative も実質使える
            pages <- free + inactive + speculative
            return(pages * ps)
        }
    }
    
    ## -------- Windows --------
    if (os == "Windows") {
        out <- tryCatch(
            system("wmic OS get FreePhysicalMemory /Value", intern = TRUE),
            error = function(e) NULL
        )
        if (!is.null(out)) {
            x <- grep("FreePhysicalMemory=", out, value = TRUE)
            if (length(x)) {
                kb <- as.numeric(gsub("[^0-9]", "", x))
                return(kb * 1024)
            }
        }
    }
    
    NA_real_
}

.setGeneID <- function(gff){
    gff$gene_id <- NA
    
    gene_i <- gff$type %in% "gene"
    gff$gene_id[gene_i] <- gff$ID[gene_i]
    
    # For transcripts or mRNA, set gene_id to their ID
    tx_i <- gff$type %in% c("transcript", "mRNA")
    tx_p <- unlist(gff$Parent[tx_i])
    if(length(tx_p) > sum(tx_i)){
        tx_p <- sapply(gff$Parent[tx_i], "[", 1)
        gff$Parent[tx_i] <- lapply(tx_p, c)
    }
    hit <- match(tx_p, gff$ID)
    gff$gene_id[tx_i] <- gff$ID[hit]
    
    # For other elements, set gene_id based on their parent transcript/mRNA
    element_i <- !gff$type %in% c("gene", "transcript", "mRNA")
    element_p <- unlist(gff$Parent[element_i])
    if(length(element_p) > sum(element_i)){
        element_p <- sapply(gff$Parent[element_i], "[", 1)
        gff$Parent[element_i] <- lapply(element_p, c)
    }
    hit <- match(element_p, gff$ID[tx_i])
    gff$gene_id[element_i] <- gff$gene_id[tx_i][hit]
    return(gff)
}

#' @importFrom rtracklayer start end start<- end<-

.fixGFF <- function(gff){
    entry_type <-  c("gene", "transcript", "mRNA", "CDS")
    gff <- gff[gff$type %in% entry_type]
    gff <- .setGeneID(gff = gff)
    gff <- .fixGFFrange(gff = gff)
    gff <- .fixGFFphase(gff = gff)
    return(gff)
}

.fixGFFrange <- function(gff){
    # Fix the start and end positions of each transcript to
    # cover the whole range of member elements (CDS, exon, and UTRs)
    tx_i <- which(gff$type %in% c("transcript", "mRNA"))
    element_i <- !gff$type %in% c("gene", "transcript", "mRNA")
    
    # Determine minimum start and maximum end positions for member elements
    min_start <- tapply(start(gff[element_i]), unlist(gff$Parent[element_i]), min)
    max_end <- tapply(end(gff[element_i]), unlist(gff$Parent[element_i]), max)
    
    # Match and update transcript start and end positions
    hit <- match(gff$ID[tx_i], names(min_start))
    tx_start <- min_start[hit]
    tx_end <- max_end[hit]
    not_na_start <- !is.na(tx_start)
    not_na_end <- !is.na(tx_end)
    start(gff[tx_i[not_na_start]]) <- tx_start[not_na_start]
    end(gff[tx_i[not_na_end]]) <- tx_end[not_na_end]
    
    # Fix the start and end positions of each gene to
    # cover the whole range of member transcripts
    gene_i <- which(gff$type == "gene")
    tx_i <- which(gff$type %in% c("transcript", "mRNA"))
    
    # Determine minimum start and maximum end positions for transcripts
    min_start <- tapply(start(gff[tx_i]), gff$gene_id[tx_i], min)
    max_end <- tapply(end(gff[tx_i]), gff$gene_id[tx_i], max)
    
    # Match and update gene start and end positions
    hit <- match(gff$gene_id[gene_i], names(min_start))
    gene_start <- min_start[hit]
    gene_end <- max_end[hit]
    start(gff[gene_i][!is.na(hit)]) <- gene_start[!is.na(hit)]
    end(gff[gene_i][!is.na(hit)]) <- gene_end[!is.na(hit)]
    
    return(gff)
}

.fixGFFphase <- function(gff){
    # Extract CDS features
    gff_cds <- gff[gff$type == "CDS"]
    
    # Adjust phase for plus strand
    gff_cds_plus <- .phasePlus(gff_cds = gff_cds)
    
    # Adjust phase for minus strand
    gff_cds_minus <- .phaseMinus(gff_cds = gff_cds)
    
    # Combine non-CDS features with adjusted CDS features
    out <- c(gff[gff$type != "CDS"], gff_cds_plus, gff_cds_minus)
    
    # Set factor levels for feature types
    out$type <- factor(out$type, levels = c("gene", "transcript", "mRNA", "CDS"))
    out$type <- droplevels(out$type)
    out$phase <- as.integer(out$phase)
    
    # Order the features by sequence name, start position, and type
    out <- out[order(as.numeric(seqnames(out)), start(out), as.numeric(out$type))]
    return(out)
}

#' @importFrom BiocGenerics start width
#' @importFrom GenomeInfoDb seqnames

.phasePlus <- function(gff_cds){
    # Filter CDS features on the plus strand
    gff_cds_plus <- gff_cds[as.character(strand(gff_cds)) == "+"]
    gff_cds_plus <- gff_cds_plus[order(as.numeric(seqnames(gff_cds_plus)), start(gff_cds_plus))]
    
    # Extract parent IDs for CDS features
    gff_cds_plus_parent <- unlist(gff_cds_plus$Parent)
    
    # Initialize phase for the first CDS in each transcript
    target_i <- which(!duplicated(gff_cds_plus_parent))
    gff_cds_plus$phase[target_i] <- 0
    
    # Calculate the phase for the next CDS feature
    next_phase <- (3 - (width(gff_cds_plus[target_i]) - gff_cds_plus$phase[target_i]) %% 3) %% 3
    names(next_phase) <- gff_cds_plus_parent[target_i]
    gff_cds_plus_parent[target_i] <- "NA"
    
    # Iterate over the remaining CDS features to set their phases
    while(TRUE){
        target_i <- which(!duplicated(gff_cds_plus_parent))[-1]
        if(length(target_i) == 0){
            break
        }
        next_phase <- next_phase[names(next_phase) %in% gff_cds_plus_parent[target_i]]
        gff_cds_plus$phase[target_i] <- next_phase
        next_phase <- (3 - (width(gff_cds_plus[target_i]) - gff_cds_plus$phase[target_i]) %% 3) %% 3
        names(next_phase) <- gff_cds_plus_parent[target_i]
        gff_cds_plus_parent[target_i] <- "NA"
    }
    
    return(gff_cds_plus)
}

#' @importFrom BiocGenerics end width
#' @importFrom GenomeInfoDb seqnames

.phaseMinus <- function(gff_cds){
    # Filter CDS features on the minus strand
    gff_cds_minus <- gff_cds[as.character(strand(gff_cds)) == "-"]
    gff_cds_minus <- gff_cds_minus[order(as.numeric(seqnames(gff_cds_minus)),
                                         end(gff_cds_minus), decreasing = TRUE)]
    
    # Extract parent IDs for CDS features
    gff_cds_minus_parent <- unlist(gff_cds_minus$Parent)
    
    # Initialize phase for the first CDS in each transcript
    target_i <- which(!duplicated(gff_cds_minus_parent))
    gff_cds_minus$phase[target_i] <- 0
    
    # Calculate the phase for the next CDS feature
    next_phase <- (3 - (width(gff_cds_minus[target_i]) - gff_cds_minus$phase[target_i]) %% 3) %% 3
    names(next_phase) <- gff_cds_minus_parent[target_i]
    gff_cds_minus_parent[target_i] <- "NA"
    
    # Iterate over the remaining CDS features to set their phases
    while(TRUE){
        target_i <- which(!duplicated(gff_cds_minus_parent))[-1]
        if(length(target_i) == 0){
            break
        }
        next_phase <- next_phase[names(next_phase) %in% gff_cds_minus_parent[target_i]]
        gff_cds_minus$phase[target_i] <- next_phase
        next_phase <- (3 - (width(gff_cds_minus[target_i]) - gff_cds_minus$phase[target_i]) %% 3) %% 3
        names(next_phase) <- gff_cds_minus_parent[target_i]
        gff_cds_minus_parent[target_i] <- "NA"
    }
    
    return(gff_cds_minus)
}

#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#' @importFrom Biostrings writeXStringSet translate

.orderGFF <- function(gff){
    # Order GFF data by chromosome and start position
    gff_order <- order(as.character(seqnames(gff)), start(gff))
    gff <- gff[gff_order]
    return(gff)
}

.makeCDS <- function(gff, genome){
    # Create a TxDb object from the GFF file
    if(inherits(gff, "GRanges")){
        txdb <- suppressMessages({makeTxDbFromGRanges(gr = gff)})
        
    } else {
        txdb <- suppressMessages({makeTxDbFromGFF(file = gff)})
    }
    
    # Read the genome file as a DNAStringSet object
    genome <- readDNAStringSet(filepath = genome)
    names(genome) <- sub("\\s.*", "", names(genome))
    
    # Extract CDS sequences from the TxDb object
    cds_db <- suppressMessages({cdsBy(x = txdb, by = "tx", use.names = TRUE)})
    
    cds <- suppressMessages({extractTranscriptSeqs(x = genome, transcripts = cds_db)})
    
    # Order CDS sequences by their names
    cds <- cds[order(names(cds))]
    return(cds)
}

#' Build/append input file definitions
#'
#' Create an `OrthoPairInput` object that stores per-genome input file paths.
#' This helper is kept for backward compatibility with previous OrthoPaiR workflows.
#'
#' @param object Existing `OrthoPairInput` object or `NULL`.
#' @param name Genome/sample name.
#' @param genome Genome FASTA path (optional; use `NA` when unavailable).
#' @param gff GFF3 path.
#' @param cds CDS FASTA path (optional; use `NA` when unavailable).
#' @param prot Protein FASTA path (optional; use `NA` when unavailable).
#' @param validation Logical; run lightweight file/format checks.
#'
#' @return `OrthoPairInput` object (list with `name`, `genome`, `gff`, `cds`, `prot`).
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
        if(is.null(genome)) genome <- rep(NA, length(gff))
        if(is.null(cds)) cds <- rep(NA, length(gff))
        if(is.null(prot)) prot <- rep(NA, length(gff))
        for(i in seq_along(gff)){
            message("Processing input ", i, ": ", name[i], "...")
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
        if(is.null(genome)) genome <- NA
        if(is.null(cds)) cds <- NA
        if(is.null(prot)) prot <- NA
    }
    
    if(validation){
        .validateInput(name = name, genome = genome, gff = gff, cds = cds, prot = prot)
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
            stop("The input object must be an OrthoPairInput class object.", call. = FALSE)
        }
        name_hit <- object$name %in% name
        if(any(name_hit)){
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

.validateInput <- function(name, genome, gff, cds, prot){
    if(is.na(name) || name == ""){
        stop("Empty name is not allowed.", call. = FALSE)
    }
    if(!is.na(genome) && !file.exists(genome)){
        stop('The file "genome" does not exist.', call. = FALSE)
    }
    if(!file.exists(gff)){
        stop('The file "gff" does not exist.', call. = FALSE)
    }
    if(!is.na(cds) && !file.exists(cds)){
        stop('The file "cds" does not exist.', call. = FALSE)
    }
    if(!is.na(prot) && !file.exists(prot)){
        stop('The file "prot" does not exist.', call. = FALSE)
    }
    
    gff_obj <- import.gff3(gff)
    chk <- .checkGFFstr(gff_obj)
    if(isTRUE(chk$check)){
        stop(paste0("In input data validation for ", name, "\n", chk$msg), call. = FALSE)
    }
    
    if(!is.na(genome)){
        gff_seq_lev <- seqlevels(gff_obj)
        genome_obj <- readDNAStringSet(filepath = genome)
        genome_seq_name <- names(genome_obj)
        check <- gff_seq_lev %in% genome_seq_name
        if(!all(check)){
            stop(paste0("In input data validation for ", name),
                 "\nFollowing sequence names appeared in the GFF but not in the genome:\n",
                 paste(head(gff_seq_lev[!check]), collapse = ", "),
                 call. = FALSE)
        }
    } else if(is.na(cds)){
        stop("the cds argument must be specified when the genome argument is missing.", call. = FALSE)
    }
}

.checkGFFstr <- function(gff){
    check <- any(is.na(gff$type)) || any(gff$type == "")
    if(check){
        return(list(check = TRUE, msg = "NA or missing values in the type column."))
    }
    gene <- gff[gff$type %in% "gene"]
    tx <- gff[gff$type %in% c("mRNA", "transcript")]
    tx_p <- unlist(tx$Parent)
    if(!all(tx_p %in% gene$ID)){
        return(list(check = TRUE,
                    msg = "The Parent column of transcripts contains values not present in gene IDs."))
    }
    if(!all(gene$ID %in% tx_p)){
        return(list(check = TRUE,
                    msg = "The ID column of genes contains values not present in transcript Parent values."))
    }
    list(check = FALSE)
}
