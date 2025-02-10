#' @importFrom rtracklayer export.gff3
#' @importFrom Biostrings writeXStringSet
#' @export
reorgOrthopiars <- function(hdf5_fn,
                            out_dir = "./",
                            out_fn = "reorg_orthopair.h5",
                            makeFASTA = TRUE,
                            overwrite = TRUE,
                            reorg = TRUE){
    dir.create(path = out_dir, showWarnings = FALSE, recursive = TRUE)
    
    genomewise_list <- .getGenomewiseList(hdf5_fn = hdf5_fn)
    
    id2id_list <- NULL
    gene_list <- NULL
    genome_fn <- NULL
    for(i in seq_along(genomewise_list)){
        target_data <- genomewise_list[[i]]
        target_genome <- names(genomewise_list)[i]
        if(reorg){
            message("reorganizing gene models in ", target_genome, " ...")
            data <- .importData(hdf5_fn = hdf5_fn, target_data = target_data, index = 1)
            
            ref_gff <- data$gff
            ref_prefix <- paste(target_genome, 
                                target_data$pair_genome[1], 
                                sep = "_")
            ref_gff <- .renameMP(gff = ref_gff, prefix = ref_prefix)
            ref_gff$oldGeneID <- ref_gff$gene_id
            ref_gff <- .setSplitGeneGFF(gff = ref_gff, 
                                        df = data$orthopair, 
                                        target_data = target_data[1, ],
                                        prefix = ref_prefix)
            
            if(nrow(target_data) == 1){
                .outputGFFdata(gff = ref_gff,
                               prefix = target_genome,
                               genome_fn = data$genome_fn,
                               out_dir = out_dir,
                               makeFASTA = makeFASTA)
            } else {
                id2id_df <- NULL
                for(j in seq_along(target_data$pair_genome)[-1]){
                    data_j <- .importData(hdf5_fn = hdf5_fn, target_data = target_data, index = j)
                    ref_prefix <- paste(target_genome, 
                                        target_data$pair_genome[j], 
                                        sep = "_")
                    ref_gff <- .setSplitGeneGFF(gff = ref_gff, 
                                                df = data_j$orthopair, 
                                                target_data = target_data[j, ],
                                                prefix = ref_prefix)
                    .importData(hdf5_fn = hdf5_fn, target_data = target_data, index = 1)
                    mp_gff <- data_j$gff
                    mp_gff <- .renameMP(gff = mp_gff, prefix = ref_prefix)
                    mp_gff <- mp_gff[mp_gff$source %in% "miniprot"]
                    mp_gff$oldGeneID <- mp_gff$gene_id
                    
                    mp_gff <- .setSplitGeneGFF(gff = mp_gff, 
                                               df = data_j$orthopair, 
                                               target_data = target_data[j, ],
                                               prefix = ref_prefix)
                    
                    id2id <- .findMiniprotTxOverlaps(ref_gff = ref_gff,
                                                     mp_gff = mp_gff)
                    
                    ref_gff <- .rerogRefGFF(ref_gff = ref_gff,
                                            mp_gff = mp_gff,
                                            id2id = id2id)
                    
                    id2id_df <- rbind(id2id_df, id2id$identical, id2id$id2id)
                    if(!is.null(id2id_df)){
                        hit <- match(id2id_df$ref_tx_id, ref_gff$ID)
                        id2id_df$ref_gene_id <- ref_gff$gene_id[hit]
                    }
                }
                id2id_list <- c(id2id_list, list(id2id_df))
                .outputGFFdata(gff = ref_gff,
                               id2id_df = id2id_df,
                               prefix = target_genome,
                               genome_fn = data$genome_fn,
                               out_dir = out_dir,
                               makeFASTA = makeFASTA)
            }
            
        } else {
            message("skip gene model reorganization for ",  target_genome, ".")
            data <- .importData(hdf5_fn = hdf5_fn, target_data = target_data, index = 1, reorg = reorg)
            ref_gff <- data$gff
        }
        
        if(is.null(ref_gff$oldGeneID)){
            ref_gff$oldGeneID <- NA
        }
        genome_fn <- c(genome_fn, data$genome_fn)
        gene_list <- rbind(gene_list, 
                           data.frame(genome = target_genome, 
                                      gene = ref_gff$gene_id[ref_gff$type %in% c("transcript", "mRNA")], 
                                      old_gene_id = ref_gff$oldGeneID[ref_gff$type %in% c("transcript", "mRNA")],
                                      tx_id = ref_gff$ID[ref_gff$type %in% c("transcript", "mRNA")]))
    }
    
    names(genome_fn) <- names(genomewise_list)
    if(!is.null(id2id_list)){
        names(id2id_list) <- names(genomewise_list)
    }
    orthopair_list <- .renameOrthoPair(id2id_list = id2id_list,
                                       gene_list = gene_list,
                                       hdf5_fn = hdf5_fn,
                                       genomewise_list = genomewise_list)
    
    out <- .outputOrthoPairData(orthopair_list = orthopair_list,
                                gene_list = gene_list,
                                genome_fn = genome_fn,
                                out_dir = out_dir,
                                out_fn = out_fn,
                                overwrite = overwrite)
    
    return(out)
}

.getGenomewiseList <- function(hdf5_fn){
    comb_id <- names(hdf5_fn)
    comb_id <- strsplit(x = comb_id, split = "_")
    comb_id <- do.call("rbind", comb_id)
    if(any(comb_id %in% c("", NA))){
        stop("The input hdf5_fn vector was not properly named.\n",
             "Each name assgined to each element of the vector should indicate",
             "a pair of genomes separated by an underscore.",
             call. = FALSE)
    }
    is_file_exist <- file.exists(hdf5_fn)
    if(any(!is_file_exist)){
        stop("Follwoing file(s) not exist: \n", 
             paste(hdf5_fn[!is_file_exist], collapse = "\n"), 
             call. = FALSE)
    }
    genome_id <- unique(as.vector(comb_id))
    out <- NULL
    for(i in seq_along(genome_id)){
        hit <- which(x = comb_id == genome_id[i], arr.ind = TRUE)
        i_target <- c("query", "subject")[hit[, 2]]
        hit[, 2] <- abs(hit[, 2] - 3)
        i_pair_genome <- comb_id[hit]
        out <- c(out, list(data.frame(comb_id = hit[, 1],
                                      target = i_target,
                                      pair_genome = i_pair_genome)))
    }
    names(out) <- genome_id
    return(out)
}

#' @importFrom rtracklayer import.gff3
#' @importFrom rhdf5 H5Fopen H5Fclose
#' @importFrom Biostrings readDNAStringSet readAAStringSet
.importData <- function(hdf5_fn, target_data, index, reorg){
    i_comb_id <- target_data$comb_id[index]
    i_target_data <- target_data$target[index]
    h5 <- H5Fopen(hdf5_fn[i_comb_id])
    on.exit(H5Fclose(h5))
    orthopair <- NULL
    
    if(index == 1){
        if(i_target_data == "query"){
            genome_fn <- h5$files$query_genome[1]
            
        } else {
            genome_fn <- h5$files$subject_genome[1]
        }
        
    } else {
        genome_fn <- NULL
    }
    
    if(i_target_data == "query"){
        gff <- import.gff3(h5$files$query_gff[1])
        
    } else {
        gff <- import.gff3(h5$files$subject_gff[1])
    }
    
    if(reorg){
        orthopair <- h5$orthopair_gene
        if(i_target_data == "query"){
            gff <- .getMiniprotGFF(gff = gff,
                                   gene = i_orthopair$query_gene,
                                   tx = i_orthopair$query_tx)
            
        } else {
            gff <- .getMiniprotGFF(gff = gff,
                                   gene = i_orthopair$subject_gene,
                                   tx = i_orthopair$subject_tx)
        }
    }
    
    out <- list(gff = gff, orthopair = orthopair, genome_fn = genome_fn)
    
    return(out)
}

.getMiniprotGFF <- function(gff, gene, tx){
    non_mp_gff <- gff[!gff$source %in% "miniprot"]
    gff <- gff[gff$source %in% "miniprot"]
    gff_id <- sub("\\:.+", "", gff$ID)
    hit <- gff_id %in% c(sub("\\:.+", "", gene), tx)
    gff <- c(non_mp_gff, gff[hit])
    return(gff)
}

#' @importFrom Biostrings translate
.outputGFFdata <- function(gff, id2id_df = NULL, prefix, genome_fn, out_dir, makeFASTA){
    message("preparing the GFF file for ", prefix, " ...")
    
    out_fn <- paste0(prefix, "_orthopair.gff")
    out_gff_fn <- file.path(out_dir, out_fn)
    if(!is.null(id2id_df)){
        gff <- .renameMPgff(gff = gff, id2id_df = id2id_df)
    }
    gff <- .setGeneID(gff = gff)
    gff <- .fixGFFrange(gff = gff)
    gff <- .fixGFFexon(gff = gff)
    gff <- .fixGFFphase(gff = gff)
    message("exporting the GFF file for ", prefix, " ...")
    export.gff3(object = gff, con = out_gff_fn)
    
    if(makeFASTA){
        if(genome_fn %in% c("no_query_genome", "no_subject_genome")){
            message("Skip CDS and protein FASTA generation for reorganized GFF becuase of no genome FASTA input.")
            
        } else {
            out_cds <- .makeCDS(gff = out_gff_fn, genome = genome_fn)
            out_cds_fn <- sub("\\.gff", ".cds", out_gff_fn)
            writeXStringSet(out_cds, out_cds_fn)
            
            out_prot <- translate(x = out_cds, no.init.codon = TRUE, if.fuzzy.codon = "X")
            out_prot_fn <- sub("\\.gff", ".prot", out_gff_fn)
            writeXStringSet(out_prot, out_prot_fn)
        }
    }
}

.renameMP <- function(gff, prefix){
    mp_entries <- grepl("^MP[0-9]+", gff$ID)
    gff$oldID <- NA
    gff$oldID[mp_entries] <- gff$ID[mp_entries]
    gff$ID[mp_entries] <- paste(prefix, gff$ID[mp_entries], sep = "_")
    not_gene <- gff$type != "gene"
    not_gene_parent <- unlist(gff$Parent[not_gene])
    mp_entries <- grepl("^MP[0-9]+", not_gene_parent)
    not_gene_parent[mp_entries] <- paste(prefix, not_gene_parent[mp_entries], sep = "_")
    gff$Parent[not_gene] <- lapply(not_gene_parent, c)
    mp_entries <- grepl("^MP[0-9]+", gff$gene_id)
    gff$gene_id[mp_entries] <- paste(prefix, gff$gene_id[mp_entries], sep = "_")
    gff$Name <- gff$ID
    return(gff)
}

.setSplitGeneGFF <- function(gff, df, target_data, prefix){
    if(target_data$target == "query"){
        target_gene <- "query_gene"
        target_tx <- "query_tx"
        
    } else {
        target_gene <- "subject_gene"
        target_tx <- "subject_tx"
    }
    
    splited_gene <- subset(df, subset = grepl("split_gene", df[[target_gene]]))
    if(nrow(splited_gene) == 0){
        return(gff)
    }
    mp <- grep("^MP[0-9]+", splited_gene[[target_tx]])
    splited_gene[[target_tx]][mp] <- paste(prefix, splited_gene[[target_tx]][mp], sep = "_")
    mp <- grep("^MP[0-9]+", splited_gene[[target_gene]])
    splited_gene[[target_gene]][mp] <- paste(prefix, splited_gene[[target_gene]][mp], sep = "_")
    
    splited_gene <- splited_gene[!splited_gene[[target_gene]] %in% gff$ID, ]
    splited_gene_gene <- splited_gene_tx <- gff[gff$ID %in% splited_gene[[target_tx]]]
    if(length(splited_gene_tx) == 0){
        return(gff)
    }
    splited_gene_tx$Parent <- lapply(paste0(splited_gene_tx$ID, ":split_gene"), c)
    splited_gene_tx$gene_id <- paste0(splited_gene_tx$ID, ":split_gene")
    splited_gene_gene$type <- "gene"
    splited_gene_gene$ID <- paste0(splited_gene_gene$ID, ":split_gene")
    splited_gene_gene$gene_id <- splited_gene_gene$Name <- splited_gene_gene$ID
    splited_gene_gene$Parent <- lapply(rep("", length(splited_gene_gene$Parent)), c)
    element_i <- !gff$type %in% c("gene", "transcript", "mRNA")
    splited_gene_element <- gff[element_i][unlist(gff$Parent[element_i]) %in% splited_gene[[target_tx]]]
    splited_gene_element$gene_id <- paste0(unlist(splited_gene_element$Parent), ":split_gene")
    splited_gff <- c(splited_gene_gene, splited_gene_tx, splited_gene_element)
    
    gff <- gff[!gff$ID %in% splited_gene[[target_tx]]]
    element_i <- !gff$type %in% c("gene", "transcript", "mRNA")
    gff_element <- gff[element_i]
    gff_element <- gff_element[!unlist(gff_element$Parent) %in% splited_gene[[target_tx]]]
    gff <- gff[!element_i]
    gff <- c(gff, gff_element)
    
    tx_i <- gff$type %in% c("transcript", "mRNA")
    no_tx_retained <- splited_gff$oldGeneID[!splited_gff$oldGeneID %in% gff$gene_id[tx_i]]
    if(length(no_tx_retained) > 0){
        gff <- gff[!gff$gene_id %in% no_tx_retained]
    }
    gff <- c(gff, splited_gff)
    
    return(gff)
}

.findMiniprotTxOverlaps <- function(ref_gff, mp_gff){
    idt_tx <- .findIdenticalMiniprotTx(ref_gff = ref_gff, mp_gff = mp_gff)
    ol_tx <- .findMiniprotTxOnGeneLoci(ref_gff = ref_gff,
                                       rest_gff = idt_tx$rest_mp)
    out <- c(ol_tx, identical = list(idt_tx$id2id))
    return(out)
}

.findIdenticalMiniprotTx <- function(ref_gff, mp_gff){
    mp_gff_block <- .getCDSblock(gff = mp_gff)
    
    # Get CDS blocks from gff2 and remove overlaps from gff1_block_uniq
    ref_gff_block <- .getCDSblock(gff = ref_gff)
    
    # Get unique transcripts from gff1
    out <- .mapIdenticalTxID(mp_gff = mp_gff, ref_gff = ref_gff,
                             mp_gff_block = mp_gff_block,
                             ref_gff_block = ref_gff_block)
    return(out)
}

#' @importFrom BiocGenerics start end
.getCDSblock <- function(gff){
    # Identify CDS and transcript indices
    gff_cds_i <- gff$type == "CDS"
    gff_tx_i <- gff$type %in% c("transcript", "mRNA")
    
    # Extract start and end positions of CDS
    gff_cds_start <- start(gff[gff_cds_i])
    gff_cds_end <- end(gff[gff_cds_i])
    
    # Create exon strings combining start and end positions
    gff_cds_exon <- paste(gff_cds_start, gff_cds_end, sep = "-")
    
    # Map CDS to their parent transcripts
    map_to_tx <- match(unlist(gff$Parent[gff_cds_i]), gff$ID[gff_tx_i])
    
    # Identify the first occurrence of each transcript
    first_i <- !duplicated(map_to_tx)
    
    # Get the chromosome names for CDS
    gff_cds_chr <- as.character(seqnames(gff[gff_cds_i]))
    
    # Add chromosome information to the first occurrence of each transcript's CDS block
    gff_cds_exon[first_i] <- paste(gff_cds_chr[first_i], gff_cds_exon[first_i], sep = ":")
    
    # Concatenate CDS blocks for each transcript
    out <- tapply(gff_cds_exon, map_to_tx, paste, collapse = ",")
    
    return(out)
}

#' @importFrom BiocGenerics unlist
.mapIdenticalTxID <- function(mp_gff, ref_gff, mp_gff_block, ref_gff_block){
    # Identify transcript indices
    ref_gff_tx_i <- ref_gff$type %in% c("transcript", "mRNA")
    mp_gff_tx_i <- mp_gff$type %in% c("transcript", "mRNA")
    
    hit <- match(mp_gff_block, ref_gff_block)
    if(all(is.na(hit))){
        id2id <- NULL
        out <- list(rest_mp = mp_gff, id2id = id2id)
        return(out)
    }
    mp_hit_tx_i <- as.numeric(names(mp_gff_block[!is.na(hit)]))
    ref_hit_tx_i <- as.numeric(names(ref_gff_block[na.omit(hit)]))
    mp_hit_tx <- mp_gff$ID[mp_gff_tx_i][mp_hit_tx_i]
    ref_hit_tx <- ref_gff$ID[ref_gff_tx_i][ref_hit_tx_i]
    hit <- match(ref_hit_tx, ref_gff$ID)
    id2id <- data.frame(mp_tx_id = mp_hit_tx,
                        ref_tx_id = ref_hit_tx,
                        mp_gene_id = NA,
                        ref_gene_id = ref_gff$gene_id[hit])
    
    retained_tx <- mp_gff[mp_gff_tx_i][!mp_gff$ID[mp_gff_tx_i] %in% mp_hit_tx]
    retained_gene <- mp_gff[mp_gff$ID %in% unlist(retained_tx$Parent)]
    element_i <- !mp_gff$type %in% c("gene", "mRNA", "transcript")
    retained_element <- mp_gff[element_i][!unlist(mp_gff$Parent[element_i]) %in% mp_hit_tx]
    out <- list(rest_mp = c(retained_gene, retained_tx, retained_element),
                id2id = id2id)
    return(out)
}

.findMiniprotTxOnGeneLoci <- function(ref_gff, rest_gff){
    # Find overlapping Miniprot Tx on original Tx
    ref_gff_cds_i <- ref_gff$type %in% c("CDS")
    rest_gff_cds_i <- rest_gff$type %in% c("CDS")
    rest_gff_tx_i <- rest_gff$type %in% c("transcript", "mRNA")
    ol <- findOverlaps(rest_gff[rest_gff_cds_i], ref_gff[ref_gff_cds_i])
    ol <- as.data.frame(ol)
    if(nrow(ol) == 0){
        out <- list(non_ol_mp_gff = rest_gff,
                    id2id = NULL,
                    id2id_split_mp = NULL,
                    id2id_split_ref = NULL)
        return(out)
    }
    ol$query_gene <- rest_gff$gene_id[rest_gff_cds_i][ol$queryHits]
    ol$queryHits <- unlist(rest_gff$Parent[rest_gff_cds_i][ol$queryHits])
    tx_id <- rest_gff$ID[rest_gff_tx_i]
    non_ol_mp_tx <- tx_id[!tx_id %in% ol$queryHits]
    ol$subject_gene <- ref_gff$gene_id[ref_gff_cds_i][ol$subjectHits]
    ol$subjectHits <- unlist(ref_gff$Parent[ref_gff_cds_i][ol$subjectHits])
    ol <- unique(ol)
    
    # Filter out chimeric Miniprot Tx overlapping more than one gene loci
    check_chimeric_in_query <- tapply(X = ol$subject_gene, INDEX = ol$query_gene, FUN = unique)
    check_chimeric_in_query <- sapply(check_chimeric_in_query, length)
    chimeric_gene_in_query <- names(check_chimeric_in_query[check_chimeric_in_query > 1])
    chimeric_in_query_ol <- subset(ol, subset = query_gene %in% chimeric_gene_in_query)
    names(chimeric_in_query_ol) <- c("mp_tx_id", "ref_tx_id", "mp_gene_id", "ref_gene_id")
    check_chimeric_in_subject <- tapply(X = ol$query_gene, INDEX = ol$subject_gene, FUN = unique)
    check_chimeric_in_subject <- sapply(check_chimeric_in_subject, length)
    chimeric_gene_in_subject <- names(check_chimeric_in_subject[check_chimeric_in_subject > 1])
    chimeric_in_subject_ol <- subset(ol, subset = subject_gene %in% chimeric_gene_in_subject)
    names(chimeric_in_subject_ol) <- c("mp_tx_id", "ref_tx_id", "mp_gene_id", "ref_gene_id")
    
    id2id <- subset(ol,
                    subset = !query_gene %in% chimeric_gene_in_query &
                        !subject_gene %in% chimeric_gene_in_subject)
    names(id2id) <- c("mp_tx_id", "ref_tx_id", "mp_gene_id", "ref_gene_id")
    
    retained_tx <- rest_gff[rest_gff_tx_i][rest_gff$ID[rest_gff_tx_i] %in% non_ol_mp_tx]
    retained_gene <- rest_gff[rest_gff$ID %in% unlist(retained_tx$Parent)]
    element_i <- !rest_gff$type %in% c("gene", "mRNA", "transcript")
    retained_element <- rest_gff[element_i][unlist(rest_gff$Parent[element_i]) %in% non_ol_mp_tx]
    non_ol_mp_gff <- c(retained_gene, retained_tx, retained_element)
    
    out <- list(non_ol_mp_gff = non_ol_mp_gff,
                id2id = id2id,
                id2id_split_mp = chimeric_in_query_ol,
                id2id_split_ref = chimeric_in_subject_ol)
    return(out)
}

.rerogRefGFF <- function(ref_gff, mp_gff, id2id){
    ref_gff <- .mergeOverlappedTx(ref_gff = ref_gff,
                                  mp_gff = mp_gff,
                                  id2id = id2id)
    ref_gff <- .mergeChimericTx(ref_gff = ref_gff,
                                mp_gff = mp_gff,
                                id2id = id2id)
    ref_gff <- c(ref_gff, id2id$non_ol_mp_gff)
    return(ref_gff)
}

.mergeOverlappedTx <- function(ref_gff, mp_gff, id2id){
    tx_id <- id2id$id2id$mp_tx_id
    merging_tx <- mp_gff[mp_gff$ID %in% tx_id]
    element_i <- !mp_gff$type %in% c("gene", "mRNA", "transcript")
    merging_element <- mp_gff[element_i][unlist(mp_gff$Parent[element_i]) %in% tx_id]
    hit <- match(merging_tx$ID, tx_id)
    merging_tx$Parent <- lapply(id2id$id2id$ref_gene_id[hit], c)
    merging_tx$gene_id <- id2id$id2id$ref_gene_id[hit]
    hit <- match(unlist(merging_element$Parent), tx_id)
    merging_element$gene_id <- id2id$id2id$ref_gene_id[hit]
    merging_gff <- c(merging_tx, merging_element)
    ref_gff <- c(ref_gff, merging_gff)
    return(ref_gff)
}

.mergeChimericTx <- function(ref_gff, mp_gff, id2id){
    tx_id <- c(id2id$id2id_split_mp$mp_tx_id, id2id$id2id_split_ref$mp_tx_id)
    merging_tx <- mp_gff[mp_gff$ID %in% tx_id]
    element_i <- !mp_gff$type %in% c("gene", "mRNA", "transcript")
    merging_element <- mp_gff[element_i][unlist(mp_gff$Parent[element_i]) %in% tx_id]
    gene_i <- mp_gff$type %in% "gene"
    merging_gene <- mp_gff[gene_i][mp_gff$ID[gene_i] %in% unlist(merging_tx$Parent)]
    merging_gff <- c(merging_gene, merging_tx, merging_element)
    ref_gff <- c(ref_gff, merging_gff)
    return(ref_gff)
}

.renameOrthoPair <- function(id2id_list, gene_list, hdf5_fn, genomewise_list){
    genomewise_list <- lapply(seq_along(genomewise_list), function(i){
        out <- data.frame(genomewise_list[[i]],
                          target_genome = names(genomewise_list[i]))
        return(out)
    })
    genomewise_list <- do.call("rbind", genomewise_list)
    genomewise_list <- subset(x = genomewise_list, subset = !duplicated(comb_id))
    genomewise_list <- genomewise_list[order(genomewise_list$comb_id), ]
    out_tx <- NULL
    out_gene <- NULL
    for(i in genomewise_list$comb_id){
        h5 <- H5Fopen(name = hdf5_fn[i])
        on.exit(H5Fclose(h5))
        orthopair_tx <- h5$orthopair_tx
        orthopair_tx <- .renameMPorthopair(df = orthopair_tx, 
                                           genome = genomewise_list[i, ])
        orthopair_gene <- h5$orthopair_gene
        orthopair_gene <- .renameMPorthopair(df = orthopair_gene, genome = genomewise_list[i, ])
        
        if(!is.null(id2id_list)){
            orthopair_tx <- .replaceOrthoPairID(df = orthopair_tx,
                                                id2id_list = id2id_list,
                                                genome = genomewise_list[i, ])
            orthopair_gene <- .replaceOrthoPairID(df = orthopair_gene,
                                                  id2id_list = id2id_list,
                                                  genome = genomewise_list[i, ])
        }
        
        orthopair_gene <- .replaceSplitGeneID(df = orthopair_gene, 
                                              gene_list = gene_list)
        out_tx <- c(out_tx, list(orthopair_tx))
        out_gene <- c(out_gene, list(orthopair_gene))
    }
    out <- list(tx = out_tx, gene = out_gene, genomewise_list = genomewise_list)
    return(out)
}

.renameMPorthopair <- function(df, genome){
    if(genome$target == "query"){
        df$query_tx <- .replaceMPid(id = df$query_tx,
                                    prefix = paste(genome$target_genome,
                                                   genome$pair_genome, 
                                                   sep = "_"))
        df$query_gene <- .replaceMPid(id = df$query_gene,
                                      prefix = paste(genome$target_genome,
                                                     genome$pair_genome, 
                                                     sep = "_"))
        df$subject_tx <- .replaceMPid(id = df$subject_tx,
                                      prefix = paste(genome$pair_genome,
                                                     genome$target_genome, 
                                                     sep = "_"))
        df$subject_gene <- .replaceMPid(id = df$subject_gene,
                                        prefix = paste(genome$pair_genome,
                                                       genome$target_genome, 
                                                       sep = "_"))
        
    } else {
        df$query_tx <- .replaceMPid(id = df$query_tx,
                                    prefix = paste(genome$pair_genome,
                                                   genome$target_genome, 
                                                   sep = "_"))
        df$query_gene <- .replaceMPid(id = df$query_gene,
                                      prefix = paste(genome$pair_genome,
                                                     genome$target_genome, 
                                                     sep = "_"))
        df$subject_tx <- .replaceMPid(id = df$subject_tx,
                                      prefix = paste(genome$target_genome,
                                                     genome$pair_genome, 
                                                     sep = "_"))
        df$subject_gene <- .replaceMPid(id = df$subject_gene,
                                        prefix = paste(genome$target_genome,
                                                       genome$pair_genome, 
                                                       sep = "_"))
    }
    return(df)
}

.replaceMPid <- function(id, prefix){
    mp_id <- grepl("^MP[0-9]+", id)
    id[mp_id] <- paste(prefix, id[mp_id], sep = "_")
    return(id)
}

.replaceOrthoPairID <- function(df, id2id_list, genome){
    if(genome$target == "query"){
        if(!is.null(id2id_list[[genome$target_genome]])){
            id2id <- subset(id2id_list[[genome$target_genome]], subset = is.na(mp_gene_id))
            hit <- match(df$query_tx, id2id$mp_tx_id)
            if(any(!is.na(hit))){
                df$query_tx[!is.na(hit)] <- id2id$ref_tx_id[na.omit(hit)]
                df$query_gene[!is.na(hit)] <- id2id$ref_gene_id[na.omit(hit)]
            }
        }
        
        if(!is.null(id2id_list[[genome$pair_genome]])){
            id2id <- subset(id2id_list[[genome$pair_genome]], subset = is.na(mp_gene_id))
            hit <- match(df$subject_tx, id2id$mp_tx_id)
            if(any(!is.na(hit))){
                df$subject_tx[!is.na(hit)] <- id2id$ref_tx_id[na.omit(hit)]
                df$subject_gene[!is.na(hit)] <- id2id$ref_gene_id[na.omit(hit)]
            }
        }
    } else {
        if(!is.null(id2id_list[[genome$target_genome]])){
            id2id <- subset(id2id_list[[genome$target_genome]], subset = is.na(mp_gene_id))
            hit <- match(df$subject_tx, id2id$mp_tx_id)
            if(any(!is.na(hit))){
                df$subject_tx[!is.na(hit)] <- id2id$ref_tx_id[na.omit(hit)]
                df$subject_gene[!is.na(hit)] <- id2id$ref_gene_id[na.omit(hit)]
            }
        }
        
        if(!is.null(id2id_list[[genome$pair_genome]])){
            id2id <- subset(id2id_list[[genome$pair_genome]], subset = is.na(mp_gene_id))
            hit <- match(df$query_tx, id2id$mp_tx_id)
            if(any(!is.na(hit))){
                df$query_tx[!is.na(hit)] <- id2id$ref_tx_id[na.omit(hit)]
                df$query_gene[!is.na(hit)] <- id2id$ref_gene_id[na.omit(hit)]
            }
        }
    }
    
    if(genome$target == "query"){
        if(!is.null(id2id_list[[genome$target_genome]])){
            id2id <- subset(id2id_list[[genome$target_genome]], subset = !is.na(mp_gene_id))
            hit <- match(df$query_tx, id2id$mp_tx_id)
            if(any(!is.na(hit))){
                df$query_gene[!is.na(hit)] <- id2id$ref_gene_id[na.omit(hit)]
            }
        }
        
        if(!is.null(id2id_list[[genome$pair_genome]])){
            id2id <- subset(id2id_list[[genome$pair_genome]], subset = !is.na(mp_gene_id))
            hit <- match(df$subject_tx, id2id$mp_tx_id)
            if(any(!is.na(hit))){
                df$subject_gene[!is.na(hit)] <- id2id$ref_gene_id[na.omit(hit)]
            }
        }
    } else {
        if(!is.null(id2id_list[[genome$pair_genome]])){
            id2id <- subset(id2id_list[[genome$pair_genome]], subset = !is.na(mp_gene_id))
            hit <- match(df$query_tx, id2id$mp_tx_id)
            if(any(!is.na(hit))){
                df$query_gene[!is.na(hit)] <- id2id$ref_gene_id[na.omit(hit)]
            }
        }
        
        if(!is.null(id2id_list[[genome$target_genome]])){
            id2id <- subset(id2id_list[[genome$target_genome]], subset = !is.na(mp_gene_id))
            hit <- match(df$subject_tx, id2id$mp_tx_id)
            if(any(!is.na(hit))){
                df$subject_gene[!is.na(hit)] <- id2id$ref_gene_id[na.omit(hit)]
            }
        }
    }
    
    return(df)
}

.replaceSplitGeneID <- function(df, gene_list){
    hit <- match(df$query_tx, gene_list$tx_id)
    if(any(!is.na(hit))){
        df$query_gene[!is.na(hit)] <- gene_list$gene[na.omit(hit)]
    }
    
    hit <- match(df$subject_tx, gene_list$tx_id)
    if(any(!is.na(hit))){
        df$subject_gene[!is.na(hit)] <- gene_list$gene[na.omit(hit)]
    }
    return(df)
}

.outputOrthoPairData <- function(orthopair_list, gene_list, genome_fn, 
                                 out_dir, out_fn, overwrite){
    hdf5_path <- file.path(out_dir, out_fn)
    hdf5_fn <- .makeHDF5(hdf5_path = hdf5_path, overwrite = overwrite)
    
    .h5creategroup(hdf5_fn, "orthopair_tx")
    .h5creategroup(hdf5_fn, "orthopair_gene")
    for(i in seq_along(orthopair_list$genomewise_list$comb_id)){
        if(orthopair_list$genomewise_list$target[i] == "query"){
            data_name <- paste(orthopair_list$genomewise_list$target_genome[i],
                               orthopair_list$genomewise_list$pair_genome[i],
                               sep = "_")
            
        } else {
            data_name <- paste(orthopair_list$genomewise_list$pair_genome[i],
                               orthopair_list$genomewise_list$target_genome[i],
                               sep = "_")
        }
        .h5overwrite(obj = orthopair_list$tx[[i]],
                     file = hdf5_fn,
                     name = paste0("orthopair_tx/", data_name))
        .h5overwrite(obj = orthopair_list$gene[[i]],
                     file = hdf5_fn,
                     name = paste0("orthopair_gene/", data_name))
    }
    
    orphan_list <- .orgOrphan(hdf5_fn = hdf5_fn, 
                              orthopair_list = orthopair_list,
                              gene_list = gene_list)
    
    .h5creategroup(hdf5_fn, "genome_fn")
    for(i in seq_along(genome_fn)){
        .h5overwrite(obj = genome_fn[i],
                     file = hdf5_fn,
                     name = paste0("genome_fn/", names(genome_fn)[i]))
    }
    
    .h5overwrite(obj = gene_list,
                 file = hdf5_fn,
                 name = "gene_list")
    
    .h5overwrite(obj = orphan_list,
                 file = hdf5_fn,
                 name = "orphan_list")
    
    .h5overwrite(obj = "reorg_orthopair",
                 file = hdf5_fn,
                 name = "data_type")
    return(hdf5_fn)
}

.renameMPgff <- function(gff, id2id_df){
    query_mp <- grepl("_MP[0-9]+", id2id_df$query_tx)
    subject_mp <- grepl("_MP[0-9]+", id2id_df$subject_tx)
    id2id_df <- subset(id2id_df, subset = query_mp & subject_mp)
    hit <- match(gff$ID, id2id_df$mp_tx_id)
    gff$ID[!is.na(hit)] <- id2id_df$ref_tx_id[na.omit(hit)]
    hit <- match(gff$ID, id2id_df$mp_gene_id)
    gff$ID[!is.na(hit)] <- id2id_df$ref_gene_id[na.omit(hit)]
    
    tx_i <- gff$type %in% c("transcript", "mRNA")
    hit <- match(unlist(gff$Parent[tx_i]), id2id_df$mp_gene_id)
    gff$Parent[tx_i][!is.na(hit)] <- lapply(id2id_df$ref_gene_id[na.omit(hit)], c)
    
    element_i <- !gff$type %in% c("gene", "transcript", "mRNA")
    hit <- match(unlist(gff$Parent[element_i]), id2id_df$mp_tx_id)
    gff$Parent[element_i][!is.na(hit)] <- lapply(id2id_df$ref_tx_id[na.omit(hit)], c)
    gff$ID[element_i] <- paste(unlist(gff$Parent[element_i]),
                               gff$type[element_i],
                               sep = ":")
    return(gff)
}

.orgOrphan <- function(hdf5_fn, orthopair_list, gene_list){
    .h5creategroup(hdf5_fn, "orphan_tx")
    .h5creategroup(hdf5_fn, "orphan_gene")
    genomes <- unique(gene_list$genome)
    out <- NULL
    for(i in seq_along(genomes)){
        hit <- grepl(genomes[i], orthopair_list$genomewise_list$pair_genome) | 
            grepl(genomes[i], orthopair_list$genomewise_list$target_genome)
        hit <- which(hit)
        hit_tx <- NULL
        for(j in hit){
            pairs_in_j <- orthopair_list$gene[[j]]
            check <- genomes[i] == orthopair_list$genomewise_list$target_genome[j]
            if(check){
                hit_tx <- c(hit_tx, pairs_in_j$query_tx)
            } else {
                hit_tx <- c(hit_tx, pairs_in_j$subject_tx)
            }
            hit_tx <- unique(hit_tx)
        }
        genome_out <- subset(gene_list,
                             subset = genome == genomes[i])
        orphan_tx <- unique(genome_out$tx_id[!genome_out$tx_id %in% hit_tx])
        hit_old_gene_id <- genome_out$old_gene_id[genome_out$tx_id %in% hit_tx]
        orphan_gene <- unique(genome_out$old_gene_id[!genome_out$old_gene_id %in% hit_old_gene_id])
        if(length(orphan_gene) > 0){
            out <- rbind(out, data.frame(genome = genomes[i], 
                                         gene = unique(orphan_gene)))
        }
        
        .h5overwrite(obj = orphan_tx,
                     file = hdf5_fn,
                     name = paste0("orphan_tx/", genomes[i]))
        .h5overwrite(obj = orphan_gene,
                     file = hdf5_fn,
                     name = paste0("orphan_gene/", genomes[i]))
    }
    return(out)
}

################################################################################
#' @importFrom igraph make_empty_graph add_edges add_vertices V<-
#' @export
#' 
makeOrthoGraph <- function(hdf5_fn){
    h5 <- H5Fopen(hdf5_fn)
    on.exit(H5Fclose(h5))
    
    check <- h5$data_type != "reorg_orthopair"
    if(check){
        stop("This function onyl accepts an output hdf5 file created",
             " by the reorgOrthopiars() function.")
    }
    n_files <- length(h5$orthopair_gene)
    orphan_gene_id <- unlist(h5$orphan_gene)
    graph <- make_empty_graph(n = 0, directed = FALSE)
    for(i in seq_len(n_files)){
        orthopair_gene <- h5$orthopair_gene[[i]]
        edges <- subset(orthopair_gene, select = c(query_gene, subject_gene))
        edges <- as.vector(t(edges))
        check <- !edges %in% V(graph)$name
        if(any(check)){
            vertices <- unique(edges[check])
            graph <- add_vertices(graph = graph, nv = length(vertices), name = vertices)
        }
        graph <- add_edges(graph = graph,
                           edges = edges, 
                           attr = list(mutual_ci = orthopair_gene$mutual_ci,
                                       class = orthopair_gene$class))
    }
    graph <- add_vertices(graph = graph, nv = length(orphan_gene_id), name = orphan_gene_id)
    return(graph)
}

#' @importFrom igraph components ego V subgraph vcount
#' @export
graph2df <- function(hdf5_fn, graph, orthopair_fn){
    h5 <- H5Fopen(hdf5_fn)
    on.exit(H5Fclose(h5))
    gene_list <- h5$gene_list
    gene_list$assign <- FALSE
    genomes <- sort(unique(gene_list$genome))
    graph_grp <- .groupGraph(graph = graph)
    
    count <- 0
    full <- TRUE
    while(TRUE){
        n_len <- length(genomes) - count
        count <- count + 1
        if(count == 2){
            full <- FALSE
        }
        if(n_len == 1){
            break
        }
        if(vcount(graph) == 0){
            break
        }
        
        full_con <- .orgFullConnected(graph = graph,
                                      genomes = genomes,
                                      gene_list = gene_list,
                                      n_len = n_len,
                                      full = full)
        if(!is.null(full_con$full_con)){
            .writeEntries(x = full_con$full_con, file = orthopair_fn)
        }
        
        if(is.null(full_con$rest_genes)){
            full_mem <- list(rest_genes = NULL)
            
        } else {
            full_mem <- .orgFullMembered(graph = graph,
                                         rest_genes = full_con$rest_genes, 
                                         genomes = genomes, 
                                         gene_list = gene_list,
                                         n_len = n_len,
                                         full = full)
            if(!is.null(full_mem$full_mem)){
                .writeEntries(x = full_mem$full_mem, file = orthopair_fn)
            }
        }
        
        if(is.null(full_mem$rest_genes)){
            full_chain <- list(full_chain = NULL)
            
        } else {
            full_chain <- .orgFullChained(graph = graph, 
                                          rest_genes = full_mem$rest_genes, 
                                          genomes = genomes, 
                                          gene_list = gene_list,
                                          n_len = n_len,
                                          full = full)
            
            if(!is.null(full_chain$full_chain)){
                .writeEntries(x = full_chain$full_chain, file = orthopair_fn)
            }
        }
        
        out <- rbind(full_con$full_con, 
                     full_mem$full_mem,
                     full_chain$full_chain)
        check <- all(full_chain$rest_genes %in% unlist(out))
        if(check){
            break
        }
        graph_ego <- ego(graph = graph, order = 1, nodes = full_chain$rest_genes)
        graph_ego <- lapply(graph_ego, names)
        graph <- subgraph(graph = graph, vids = unlist(graph_ego))
    }
    out <- read.csv(file = orthopair_fn)
    out$SOG <- NA
    i <- 1
    while(TRUE){
        if(all(!is.na(out$SOG))){
            break
        }
        hit <- match(out[, i], graph_grp$gene_id)
        out$SOG[!is.na(hit)] <- graph_grp$SOG[na.omit(hit)]
        i <- i + 1
    }
    out <- out[order(out$SOG), ]
    out <- .setTxID(out = out, gene_list = gene_list)
    write.csv(out, file = orthopair_fn, row.names = FALSE)
    invisible(orthopair_fn)
}

.writeEntries <- function(x, file){
    if(file.exists(file)){
        write.table(x = x, 
                    file = file,
                    append = TRUE, 
                    quote = FALSE,
                    sep = ",", 
                    col.names = FALSE,
                    row.names = FALSE)
        
    } else {
        write.table(x = x, 
                    file = file,
                    append = FALSE, 
                    quote = FALSE,
                    sep = ",", 
                    col.names = TRUE,
                    row.names = FALSE)
    }
}

.is_mult_mp_in_grp <- function(x){
    mp <- grep("_MP[0-9]+", x, value = TRUE)
    out <- FALSE
    if(length(mp) > 1){
        check <- length(unique(sub("_.+", "", mp))) > 1
        if(check){
            out <- TRUE
        }
    }
    return(out)
}

.groupGraph <- function(graph){
    grp <- split(V(graph)$name, components(graph)$membership)
    mult_mp_in_grp <- sapply(grp, .is_mult_mp_in_grp)
    
    n_gene_in_grp <- lapply(grp, length)
    out <- sapply(seq_along(n_gene_in_grp), function(i){
        return(rep(i, n_gene_in_grp[[i]]))
    })
    out <- data.frame(gene_id = unlist(grp), SOG = unlist(out))
    hit <- match(out$SOG, as.numeric(names(mult_mp_in_grp)))
    out$mult_mp_in_grp <- mult_mp_in_grp[hit]
    return(out)
}

#' @importFrom igraph max_cliques
.getGraphCliques <- function(graph){
    graph_cliques <- max_cliques(graph = graph)
    graph_cliques <- lapply(graph_cliques, attributes)
    graph_cliques <- lapply(graph_cliques, "[[", "names")
    n_cliques <- sapply(graph_cliques, length)
    out <- list(cliques = graph_cliques, n = n_cliques)
    return(out)
}

.orgFullConnected <- function(graph, genomes, gene_list, n_len, full){
    graph_cliques <- .getGraphCliques(graph = graph)
    is_full_con <- graph_cliques$n == n_len
    if(sum(is_full_con) == 0){
        out <- list(rest_genes = unlist(graph_cliques$cliques[!is_full_con]),
                    full_con = NULL)
        return(out)
    }
    full_con <- graph_cliques$cliques[is_full_con]
    full_con <- data.frame(gene_id = unlist(full_con),
                           grp_id = rep(seq_along(full_con),
                                        each = n_len))
    hit <- match(full_con$gene_id, gene_list$gene)
    full_con$genome <- gene_list$genome[hit]
    full_df <- data.frame(id = paste(rep(unique(full_con$grp_id),
                                         each = length(genomes)),
                                     genomes,
                                     sep = "_"))
    full_con$id <- paste(full_con$grp_id, full_con$genome, sep = "_")
    full_con <- left_join(x = full_df, y = full_con, by = "id")
    
    full_con <- matrix(data = full_con$gene_id, 
                       ncol = length(genomes),
                       byrow = TRUE)
    full_con <- as.data.frame(full_con)
    names(full_con) <- genomes
    
    if(full){
        full_con$class <- "full_connected"
        
    } else {
        full_con$class <- "partial_connected"
    }
    
    rest_genes <- unlist(graph_cliques$cliques[!is_full_con])
    out <- list(rest_genes = rest_genes,
                full_con = full_con)
    return(out)
}

.orgFullMembered <- function(graph, rest_genes, genomes, gene_list, n_len, full){
    rest_ego <- ego(graph = graph, order = 1, nodes = rest_genes)
    rest_ego <- lapply(rest_ego, names)
    n_rest_ego <- sapply(rest_ego, length)
    ego_id <- sapply(seq_along(n_rest_ego), function(i){
        return(rep(i, n_rest_ego[[i]]))
    })
    rest_ego <- data.frame(gene_id = unlist(rest_ego), ego_id = unlist(ego_id))
    hit <- match(rest_ego$gene_id, gene_list$gene)
    rest_ego$genome <- gene_list$genome[hit]
    order_rest_ego <- order(rest_ego$ego_id, rest_ego$genome)
    rest_ego <- rest_ego[order_rest_ego, ]
    n_rest_ego_nodes <- table(rest_ego$ego_id)
    full_mem_candidate <- as.numeric(n_rest_ego_nodes) == n_len
    if(sum(full_mem_candidate) == 0){ return(list(rest_genes = rest_genes)) }
    full_mem_candidate <- rest_ego[rest_ego$ego_id %in% as.numeric(names(n_rest_ego_nodes[full_mem_candidate])), ]
    
    full_ego_id <- rep(unique(full_mem_candidate$ego_id), each = length(genomes))
    full_df <- data.frame(id = paste(full_ego_id, genomes, sep = "_"),
                          full_ego_id = full_ego_id)
    full_mem_candidate$id <- paste(full_mem_candidate$ego_id, 
                                   full_mem_candidate$genome, sep = "_")
    full_mem_candidate <- left_join(x = full_df, y = full_mem_candidate, by = "id")
    
    n_entry_full_mem_candidate <- table(full_mem_candidate$full_ego_id)
    non_multientries <- n_entry_full_mem_candidate == length(genomes)
    if(sum(non_multientries) == 0){ return(list(rest_genes = rest_genes)) }
    valid_ego_id <- names(n_entry_full_mem_candidate)[non_multientries]
    full_mem_candidate <- subset(full_mem_candidate, 
                                 subset = full_ego_id %in% valid_ego_id)
    
    is_na_full_mem_candidate <- tapply(full_mem_candidate$ego_id, 
                                       full_mem_candidate$full_ego_id,
                                       is.na)
    n_entry_full_mem_candidate <- sapply(is_na_full_mem_candidate, sum)
    n_entry_full_mem_candidate <- length(genomes) - n_entry_full_mem_candidate
    full_mem <- n_entry_full_mem_candidate == n_len
    if(sum(full_mem) == 0){ return(list(rest_genes = rest_genes)) }
    valid_ego_id <- names(n_entry_full_mem_candidate[full_mem])
    full_mem <- full_mem_candidate[full_mem_candidate$full_ego_id %in% as.numeric(valid_ego_id), ]
    full_mem <- matrix(data = full_mem$gene_id, 
                       ncol = length(genomes),
                       byrow = TRUE)
    full_mem <- unique(full_mem)
    full_mem <- as.data.frame(full_mem)
    names(full_mem) <- genomes
    
    if(full){
        full_mem$class <- "full_membered"
        
    } else {
        full_mem$class <- "partial_membered"
    }
    
    rest_genes <- rest_genes[!rest_genes %in% unlist(full_mem)]
    out <- list(rest_genes = rest_genes, full_mem = full_mem)
    return(out)
}

.orgFullChained <- function(graph, rest_genes, genomes, gene_list, n_len, full){
    rest_ego <- ego(graph = graph, order = length(genomes), nodes = rest_genes)
    rest_ego <- lapply(rest_ego, names)
    n_rest_ego <- sapply(rest_ego, length)
    ego_id <- sapply(seq_along(n_rest_ego), function(i){
        return(rep(i, n_rest_ego[[i]]))
    })
    rest_ego <- data.frame(gene_id = unlist(rest_ego), ego_id = unlist(ego_id))
    hit <- match(rest_ego$gene_id, gene_list$gene)
    rest_ego$genome <- gene_list$genome[hit]
    order_rest_ego <- order(rest_ego$ego_id, rest_ego$genome)
    rest_ego <- rest_ego[order_rest_ego, ]
    n_rest_ego_nodes <- table(rest_ego$ego_id)
    full_chain_candidate <- as.numeric(n_rest_ego_nodes) == n_len
    if(sum(full_chain_candidate) == 0){ return(list(rest_genes = rest_genes)) }
    full_chain_candidate <- rest_ego[rest_ego$ego_id %in% as.numeric(names(n_rest_ego_nodes[full_chain_candidate])), ]
    
    full_ego_id <- rep(unique(full_chain_candidate$ego_id), each = length(genomes))
    full_df <- data.frame(id = paste(full_ego_id, genomes, sep = "_"),
                          full_ego_id = full_ego_id)
    full_chain_candidate$id <- paste(full_chain_candidate$ego_id, 
                                     full_chain_candidate$genome, sep = "_")
    full_chain_candidate <- left_join(x = full_df, y = full_chain_candidate, by = "id")
    
    n_entry_full_chain_candidate <- table(full_chain_candidate$full_ego_id)
    non_multientries <- n_entry_full_chain_candidate == length(genomes)
    if(sum(non_multientries) == 0){ return(list(rest_genes = rest_genes)) }
    valid_ego_id <- names(n_entry_full_chain_candidate)[non_multientries]
    full_chain_candidate <- subset(full_chain_candidate, 
                                   subset = full_ego_id %in% valid_ego_id)
    
    is_na_full_chain_candidate <- tapply(full_chain_candidate$ego_id, 
                                         full_chain_candidate$full_ego_id,
                                         is.na)
    n_entry_full_chain_candidate <- sapply(is_na_full_chain_candidate, sum)
    n_entry_full_chain_candidate <- length(genomes) - n_entry_full_chain_candidate
    full_chain <- n_entry_full_chain_candidate == n_len
    if(sum(full_chain) == 0){ return(list(rest_genes = rest_genes)) }
    valid_ego_id <- names(n_entry_full_chain_candidate[full_chain])
    
    full_chain <- full_chain_candidate[full_chain_candidate$full_ego_id %in% as.numeric(valid_ego_id), ]
    full_chain <- matrix(data = full_chain$gene_id, 
                         ncol = length(genomes),
                         byrow = TRUE)
    full_chain <- unique(full_chain)
    full_chain <- as.data.frame(full_chain)
    names(full_chain) <- genomes
    
    if(full){
        full_chain$class <- "full_chained"
        
    } else {
        full_chain$class <- "partial_chained"
    }
    
    rest_genes <- rest_genes[!rest_genes %in% unlist(full_chain)]
    out <- list(rest_genes = rest_genes, full_chain = full_chain)
    return(out)
}

.setTxID <- function(out = out, gene_list = gene_list){
    tx_out <- subset(out, select = -c(class, SOG))
    names(tx_out) <- paste(names(tx_out), "TX", sep = "_")
    for(i in seq_len(ncol(tx_out))){
        hit <- match(tx_out[, i], gene_list$gene)
        tx_out[, i] <- gene_list$tx_id[hit]
    }
    out <- cbind(out, tx_out)
}