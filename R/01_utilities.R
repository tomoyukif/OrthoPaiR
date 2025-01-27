#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists h5delete h5write
# Function to overwrite an HDF5 dataset
.h5overwrite <- function(obj, file, name){
    # Open the HDF5 file
    file <- H5Fopen(file)
    
    # Ensure the file is closed when the function exits
    on.exit(H5Fclose(file))
    
    # Check if the specified name exists in the file
    if(H5Lexists(h5loc = file, name = name)){
        # Delete the existing dataset if it exists
        h5delete(file = file, name = name)
    }
    
    # Write the new object to the HDF5 file
    h5write(obj = obj, file = file, name = name)
}

#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists h5createGroup
# Function to create a group in an HDF5 file
.h5creategroup <- function(file, name){
    # Open the HDF5 file
    file <- H5Fopen(file)
    
    # Ensure the file is closed when the function exits
    on.exit(H5Fclose(file))
    
    # Check if the specified group name does not exist in the file
    if(!H5Lexists(h5loc = file, name = name)){
        # Create a new group in the HDF5 file
        h5createGroup(file, name)
    }
}

#' @importFrom rtracklayer import.gff
# Function to import multiple GFF files
.importAllGFF <- function(fn){
    for(i in seq_along(fn)){
        if(i == 1){
            # Import the first GFF file
            out <- import.gff(fn[i])
            
        } else {
            # Concatenate subsequent GFF files
            out <- c(out, import.gff(fn[i]))
        }
    }
    # Return the concatenated GFF data
    return(out)
}

#' @importFrom rhdf5 h5createFile
# Function to create an HDF5 file
.makeHDF5 <- function(hdf5_path, overwrite){
    if(file.exists(hdf5_path)){
        # Notify if the file already exists
        message(hdf5_path, " already exists.")
        
        if(overwrite){
            message("Since user specified 'overwrite = TRUE', the exsiting HDF5 will be overwritten.")
            
            unlink(x = hdf5_path)
            
            # Create a new HDF5 file
            h5createFile(file = hdf5_path)
        }
        
    } else {
        # Create a new HDF5 file
        h5createFile(file = hdf5_path)
    }
    
    # Return the HDF5 file path
    return(hdf5_path)
}


#' @importFrom Biostrings readDNAStringSet writeXStringSet translate
#' @importFrom rtracklayer import.gff export.gff3
#' @importFrom GenomeInfoDb seqlevels seqnames seqlevels<- dropSeqlevels
#' @importFrom GenomicFeatures cdsBy extractTranscriptSeqs
#' @importFrom txdbmaker makeTxDbFromGFF
#' @importFrom S4Vectors mcols mcols<-
#' @import BSgenome
#' 
#' @export
fixInfiles <- function(genome, gff){
    new_genome_fn <- paste0(sub(pattern = "\\.[^\\.]*$", replacement = "", genome), ".renamed.fa")
    new_gff_fn <- paste0(sub(pattern = "\\.[^\\.]*$", replacement = "", gff), ".renamed.gff")
    new_cds_fn <- paste0(sub(pattern = "\\.[^\\.]*$", replacement = "", new_gff_fn), ".cds.fa")
    new_prot_fn <- paste0(sub(pattern = "\\.[^\\.]*$", replacement = "", new_gff_fn), ".prot.fa")
    
    genome <- readDNAStringSet(genome)
    old_genome_names <- genome_names <- names(genome)
    message("Sequence names in the genome FASTA:\n",
            paste(genome_names, collapse = "\n"))
    chr_prefix <- NULL
    check <- readline(prompt = "Rename genome sequence(s)? (y/n): ")
    if(check != "n"){
        while(TRUE){
            while(TRUE){
                pattern <- readline(prompt = "Set pattern argument: ")
                replacement <- readline(prompt = "Set replacement argument: ")
                if(pattern != ""){
                    genome_names <- sub(pattern, replacement, genome_names)
                }
                message("Sequence names in the genome FASTA:\n",
                        paste(genome_names, collapse = "\n"))
                check <- readline(prompt = "Replace more? (y/n): ")
                if(check == "n"){
                    break
                }
            }
            
            check <- readline(prompt = "Pad numbers? (y/n): ")
            if(check != "n"){
                digits <- readline(prompt = "How many digits?: ")
                if(digits != ""){
                    is_numeric <- !is.na(suppressWarnings(as.numeric(genome_names)))
                    chr_prefix <- readline(prompt = "Set prefix for chromosomes, e.g. chr: ")
                    genome_names[is_numeric] <- paste0(chr_prefix, sprintf(paste0("%0", digits, "d"),
                                                                           as.numeric(genome_names[is_numeric])))
                }
                message("Sequence names in the genome FASTA:\n",
                        paste(genome_names, collapse = "\n"))
            }
            
            check <- readline(prompt = "Replace other sequence names? (y/n): ")
            if(check == "n"){
                names(genome) <- genome_names
                break
            }
        }
    }
    genome_names_full <- genome_names
    check <- readline(prompt = "Remove sequence(s)? (y/n): ")
    if(check != "n"){
        while(TRUE){
            pattern <- readline(prompt = "Set pattern argument: ")
            if(pattern != ""){
                genome <- genome[!grepl(pattern, genome_names)]
            }
            genome_names <- names(genome)
            message("Sequence names in the genome FASTA:\n",
                    paste(genome_names, collapse = "\n"))
            check <- readline(prompt = "Remove more? (y/n): ")
            if(check == "n"){
                names(genome) <- genome_names
                break
            }
        }
    }
    
    message("Saving new genome FASTA file.")
    writeXStringSet(x = genome, new_genome_fn)
    
    gff <- import.gff(gff)
    gff$Name <- gff$ID
    gff_levels <- seqlevels(gff)
    while(TRUE){
        if(all(gff_levels %in% old_genome_names)){
            break
            
        } else {
            message("Old genome sequece names does not match seqlevels in GFF.")
            message("Replacement is required.")
            message("Sequence names in the old genome FASTA:\n",
                    paste(old_genome_names, collapse = "\n"))
            message("Sequence levels in the GFF:\n",
                    paste(gff_levels, collapse = "\n"))
            pattern <- readline(prompt = "Set pattern argument: ")
            replacement <- readline(prompt = "Set replacement argument: ")
            if(pattern != ""){
                old_genome_names <- sub(pattern, replacement, old_genome_names)
            }
        }
    }
    
    hit <- match(gff_levels, old_genome_names)
    seqlevels(gff) <- genome_names_full[hit]
    genome_names <- names(genome)
    gff_levels <- seqlevels(gff)
    pattern <- gff_levels[!gff_levels %in% genome_names]
    gff <- dropSeqlevels(gff, pattern, pruning.mode = "coarse")
    
    gff_levels <- seqlevels(gff)
    message("Sequence levels in the GFF:\n",
            paste(gff_levels, collapse = "\n"))
    check <- readline(prompt = "Drop sequence level(s)? (y/n): ")
    if(check != "n"){
        while(TRUE){
            pattern <- readline(prompt = "Set pattern argument: ")
            if(pattern != ""){
                pattern <- gff_levels[grepl(pattern, gff_levels)]
                gff <- dropSeqlevels(gff, pattern, pruning.mode = "coarse")
            }
            message("Sequence levels in the GFF:\n",
                    paste(seqlevels(gff), collapse = "\n"))
            check <- readline(prompt = "Drop more? (y/n): ")
            if(check == "n"){
                break
            }
        }
    }
    
    message("\nDrop GFF entries of unnecessary type(s).\n")
    gene <- gff[grep("^gene$", gff$type)]
    tx <- gff[gff$type %in% c("mRNA", "transcript")]
    tx_p <- sapply(tx$Parent, function(x){if(length(x) == 0) {return(NA)} else {return(x[1])}})
    hit_tx <- tx[tx_p %in% gene$ID]
    gene <- gene[gene$ID %in% tx_p]
    element <- gff[gff$type %in% c("exon", "CDS")]
    element_p <- sapply(element$Parent, function(x){if(length(x) == 0) {return(NA)} else {return(x[1])}})
    hit_element <- element[element_p %in% hit_tx$ID]
    check <- any(sapply(hit_element$Parent, length) > 1)
    if(check){
        hit_element$Parent <- lapply(hit_element$Parent, "[", 1)
    }
    gff <- c(gene, hit_tx, hit_element)
    gff <- gff[order(as.numeric(seqnames(gff)), start(gff), as.numeric(gff$type))]
    Sys.sleep(3)
    
    check <- any(is.na(gff$ID))
    if(check){
        gff$ID[gff$type %in% c("CDS", "exon")] <- paste(unlist(gff$Parent[gff$type %in% c("CDS", "exon")]), 
                                                        gff$type[gff$type %in% c("CDS", "exon")], 
                                                        sep = ":")
    }
    message("Gene IDs in the GFF:\n",
            paste(head(gff$ID[gff$type == "gene"]), collapse = "\n"))
    check <- readline(prompt = "Replace IDs? (y/n): ")
    if(check != "n"){
        id_prefix <- readline(prompt = "Set ID prefix, e.g. Os or Osat for Oryza sativa: ")
        
        while(TRUE){
            if(id_prefix == ""){
                message("Need to set ID prefix.\n")
                check <- readline(prompt = "Do not replave IDs? (y/n): ")
                if(check == "n"){
                    break
                    
                } else {
                    id_prefix <- readline(prompt = "Set ID prefix, e.g. Os_ or Osat_ for Oryza sativa: ")
                }
            } else {
                break
            }
        }
        
        if(id_prefix != ""){
            message("Replaced IDs will '[ID prefix][chr number]G[serial number]'\n",
                    "The chr numbers will be otbained from the sequence levels of the GFF",
                    "Please specify the prefix of the sequence levels of chromosomes",
                    " that will be removed to obtain numbers only.",
                    "Sequence levels in the GFF:\n",
                    paste(seqlevels(gff), collapse = "\n"))
            chr_prefix <- readline(prompt = "Set chr prefix, e.g. chr: ")
            message("Two zeros will be automatically added at the end of serial numbers.")
            digits_fmt <- readline(prompt = "How many digits for serial numbers? (default = 5): ")
            if(digits_fmt == ""){
                digits_fmt <- "%05d00"
            } else {
                digits_fmt <- paste0("%0", digits_fmt, "d00")
            }
            
            gff <- .renameGFF(gff = gff,
                              id_prefix = id_prefix,
                              chr_prefix = chr_prefix,
                              digits_fmt = digits_fmt)
        }
    }
    
    message(prompt = "\nChecking the 'gene_id' column in the GFF.\n")
    need_fix_gene_id <- FALSE
    if(is.null(gff$gene_id)){
        need_fix_gene_id <- TRUE
        
    } else {
        if(any(is.na(gff$gene_id))){
            need_fix_gene_id <- TRUE
        }
    }
    if(need_fix_gene_id){
        gff <- .setGeneID(gff)
    }
    
    message(prompt = "\nChecking the 'Parent' column in the GFF.\n")
    non_gene_i <- !gff$type %in% c("gene")
    n_entry <- sapply(gff$Parent[non_gene_i], length)
    if(any(n_entry > 1)){
        message(prompt = "\nSome entires in the Parent column have multuple IDs.\n")
        message(prompt = "\nEach entry in the Parent column must be a single value to process in OrthoPaiR.\n")
        message(prompt = "\nThe first ID in each entry usually is the ID that is required to process in OrthoPaiR.\n")
        while(TRUE){
            check <- readline(prompt = "Are you sure to take only the first ID for each entry? (y/n): ")
            if(check == "y"){
                non_gene_i <- !gff$type %in% c("gene")
                non_gene_p <- sapply(gff$Parent[non_gene_i], "[", 1)
                gff$Parent[non_gene_i] <- lapply(non_gene_p, c)
                hit <- match(non_gene_p, gff$ID)
                gff$gene_id[non_gene_i] <- gff$ID[hit]
                break
                
            } else if(check == "n"){
                message(prompt = "\nOrthoPaiR may produce an error during the process.\n")
                break
            }
            
        }
    }
    
    message(prompt = "\nDrop unnecessary GFF annotation(s).\n")
    gff_mcols <- mcols(gff)
    if(is.null(gff_mcols$oldID)){
        gff_mcols <- subset(gff_mcols, select = c(source:ID, Name, Parent, gene_id))
        
    } else {
        gff_mcols <- subset(gff_mcols, select = c(source:ID, Name, Parent, oldID, gene_id))
    }
    mcols(gff) <- gff_mcols
    Sys.sleep(3)
    
    message("Saving new GFF file.")
    export.gff3(gff, new_gff_fn)
    
    message("Generate CDS and protein FASTA files.")
    txdb <- makeTxDbFromGFF(file = new_gff_fn, format = "gff3")
    genome <- readDNAStringSet(filepath = new_genome_fn)
    cds_db <- cdsBy(x = txdb, by = "tx", use.names = TRUE)
    cds_seqs <- extractTranscriptSeqs(x = genome, transcripts = cds_db)
    cds_seqs <- cds_seqs[order(names(cds_seqs))]
    prot_seqs <- translate(x = cds_seqs, no.init.codon = TRUE, if.fuzzy.codon = "X")
    prot_seqs <- prot_seqs[order(names(prot_seqs))]
    writeXStringSet(x = cds_seqs, filepath = new_cds_fn)
    writeXStringSet(x = prot_seqs, filepath = new_prot_fn)
    
    out <- list(genome = new_genome_fn, gff = new_gff_fn, cds = new_cds_fn, prot = new_prot_fn)
    return(out)
}

#' @importFrom GenomeInfoDb seqlevels seqnames seqlevels<-
#' @importFrom rtracklayer start 
.renameGFF <- function(gff, id_prefix, chr_prefix, digits_fmt){
    seqlevels(gff) <- sort(seqlevels(gff))
    gff <- gff[order(as.numeric(seqnames(gff)), start(gff), as.numeric(gff$type))]
    gene_i <- gff$type %in% "gene"
    transcript_i <- gff$type %in% c("transcript", "mRNA")
    element_i <- !gene_i & !transcript_i
    
    chr <- seqnames(gff[gene_i])
    id_table <- tapply(seq_along(gff[gene_i]), chr, function(x){
        chr_gene <- gff[gene_i][x]
        this_chr <- unique(seqnames(chr_gene))
        prefix <- paste0(id_prefix, sub(chr_prefix, "", this_chr), "G")
        new_id <- paste0(prefix, sprintf(fmt = digits_fmt,
                                         seq_along(chr_gene)))
        return(data.frame(old = chr_gene$ID, new = new_id))
    })
    id_table <- do.call("rbind", id_table)
    gff$oldID <- gff$ID
    gff[gene_i]$oldID <- id_table$old[match(gff[gene_i]$ID, id_table$old)]
    gff[gene_i]$ID <- id_table$new[match(gff[gene_i]$ID, id_table$old)]
    hit <- match(unlist(gff$Parent[transcript_i]), gff[gene_i]$oldID)
    gff$Parent[transcript_i] <- lapply(gff[gene_i]$ID[hit],
                                       function(x){return(x)})
    
    p_id <- unlist(gff$Parent[transcript_i])
    n_tx_max <- max(table(p_id))
    valid <- rep(TRUE, length(p_id))
    for(i in seq_len(n_tx_max)){
        uniq_p_id <- unique(p_id[valid])
        hit <- match(uniq_p_id, p_id[valid])
        gff$ID[transcript_i][valid][hit] <- paste(sub("G", "T",
                                                      p_id[valid][hit]),
                                                  sprintf("%02d", i),
                                                  sep = ".")
        valid[valid][hit] <- FALSE
    }
    
    hit <- match(unlist(gff$Parent[element_i]), gff[transcript_i]$oldID)
    gff$Parent[element_i] <- lapply(gff[transcript_i]$ID[hit],
                                    function(x){return(x)})
    gff$ID[element_i] <- paste(unlist(gff$Parent[element_i]),
                               gff$type[element_i],
                               sep = ":")
    
    gff$Name <- gff$ID
    gff$gene_id <- sub("T", "G", sub("\\..+", "", gff$ID))
    gff <- gff[order(as.numeric(seqnames(gff)), start(gff), as.numeric(gff$type))]
    return(gff)
}

#' @importFrom rtracklayer import.gff export.gff3
#' @importFrom S4Vectors mcols mcols<-
#'
#' @export
formatGFF <- function(gff, suffix = ".reformat.gff"){
  new_gff_fn <- paste0(sub(pattern = "\\.gff3?\\.?.?$", replacement = "", gff), suffix)
  
  gff <- import.gff(gff)
  gff$Name <- gff$ID
  
  message("\nDrop GFF entries of unnecessary type(s).\n")
  gene <- gff[grep("^gene$", gff$type)]
  tx <- gff[gff$type %in% c("mRNA", "transcript")]
  tx_p <- sapply(tx$Parent, function(x){if(length(x) == 0) {return(NA)} else {return(x[1])}})
  hit_tx <- tx[tx_p %in% gene$ID]
  gene <- gene[gene$ID %in% tx_p]
  element <- gff[gff$type %in% c("exon", "CDS")]
  element_p <- sapply(element$Parent, function(x){if(length(x) == 0) {return(NA)} else {return(x[1])}})
  hit_element <- element[element_p %in% hit_tx$ID]
  check <- any(sapply(hit_element$Parent, length) > 1)
  if(check){
    hit_element$Parent <- lapply(hit_element$Parent, "[", 1)
  }
  gff <- c(gene, hit_tx, hit_element)
  gff <- gff[order(as.numeric(seqnames(gff)), start(gff), as.numeric(gff$type))]
  Sys.sleep(3)
  
  check <- any(is.na(gff$ID))
  if(check){
    gff$ID[gff$type %in% c("CDS", "exon")] <- paste(unlist(gff$Parent[gff$type %in% c("CDS", "exon")]),
                                                    gff$type[gff$type %in% c("CDS", "exon")],
                                                    sep = ":")
  }
  
  message(prompt = "\nChecking the 'gene_id' column in the GFF.\n")
  need_fix_gene_id <- FALSE
  if(is.null(gff$gene_id)){
    need_fix_gene_id <- TRUE
    
  } else {
    if(any(is.na(gff$gene_id))){
      need_fix_gene_id <- TRUE
    }
  }
  if(need_fix_gene_id){
    gff <- .setGeneID(gff)
  }
  
  message(prompt = "\nChecking the 'Parent' column in the GFF.\n")
  non_gene_i <- !gff$type %in% c("gene")
  n_entry <- sapply(gff$Parent[non_gene_i], length)
  if(any(n_entry > 1)){
    while(TRUE){
      non_gene_i <- !gff$type %in% c("gene")
      non_gene_p <- sapply(gff$Parent[non_gene_i], "[", 1)
      gff$Parent[non_gene_i] <- lapply(non_gene_p, c)
      hit <- match(non_gene_p, gff$ID)
      gff$gene_id[non_gene_i] <- gff$ID[hit]
    }
  }
  
  message(prompt = "\nDrop unnecessary GFF annotation(s).\n")
  gff_mcols <- mcols(gff)
  if(is.null(gff_mcols$oldID)){
    gff_mcols <- subset(gff_mcols, select = c(source:ID, Name, Parent, gene_id))
    
  } else {
    gff_mcols <- subset(gff_mcols, select = c(source:ID, Name, Parent, oldID, gene_id))
  }
  mcols(gff) <- gff_mcols
  Sys.sleep(3)
  
  message("Saving new GFF file.")
  export.gff3(gff, new_gff_fn)
}
