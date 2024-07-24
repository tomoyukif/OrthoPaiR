#' @importFrom rtracklayer export.gff3
#' @importFrom Biostrings writeXStringSet
orgMPgenes <- function(hdf5_fn,
                       out_dir = "./",
                       out_fn = "reorg_synog.h5",
                       makeFASTA = TRUE,
                       overwrite = TRUE){
    dir.create(path = out_dir, showWarnings = FALSE, recursive = TRUE)

    genomewise_list <- .getGenomewiseList(hdf5_fn = hdf5_fn)

    id2id_list <- NULL
    for(i in seq_along(genomewise_list)){
        target_data <- genomewise_list[[i]]
        data <- .importData(hdf5_fn = hdf5_fn, target_data = target_data)

        ref_gff <- data$gff[[1]]
        ref_gff <- .renameMP(gff = ref_gff,
                             pair_genome = target_data$pair_genome[1])

        if(nrow(target_data) == 1){
            .outputGFFdata(gff = ref_gff,
                           prefix = names(genomewise_list)[1],
                           genome_fn = data$genome_fn,
                           out_dir = out_dir,
                           makeFASTA = makeFASTA)
            next
        }

        id2id_df <- NULL
        for(j in seq_along(data$gff)[-1]){
            mp_gff <- data$gff[[j]]
            mp_gff <- .renameMP(gff = mp_gff,
                                pair_genome = target_data$pair_genome[j])
            mp_gff <- mp_gff[mp_gff$source %in% "miniprot"]

            id2id <- .findMiniprotTxOverlaps(ref_gff = ref_gff,
                                             mp_gff = mp_gff)

            ref_gff <- .rerogRefGFF(ref_gff = ref_gff,
                                    mp_gff = mp_gff,
                                    id2id = id2id)

            id2id_df <- rbind(id2id_df, id2id$identical, id2id$id2id)
        }
        id2id_list <- c(id2id_list, list(id2id_df))
        .outputGFFdata(gff = ref_gff,
                       prefix = names(genomewise_list)[i],
                       genome_fn = data$genome_fn,
                       out_dir = out_dir,
                       makeFASTA = makeFASTA)
    }
    names(id2id_list) <- names(genomewise_list)

    synog_list <- .renameSynog(id2id_list = id2id_list,
                               hdf5_fn = hdf5_fn,
                               genomewise_list = genomewise_list)

    out <- .outputSynogData(synog_list = synog_list,
                            out_dir = out_dir,
                            out_fn = out_fn,
                            overwrite = overwrite)

    return(out)
}

.getGenomewiseList <- function(hdf5_fn){
    comb_id <- names(hdf5_fn)
    comb_id <- strsplit(x = comb_id, split = "_")
    comb_id <- do.call("rbind", comb_id)
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
.importData <- function(hdf5_fn, target_data){
    for(i in seq_along(target_data$comb_id)){
        i_comb_id <- target_data$comb_id[i]
        i_target_data <- target_data$target[i]
        h5 <- H5Fopen(hdf5_fn[i_comb_id])
        i_synog <- h5$synog_gene

        if(i_target_data == "query"){
            i_gff <- import.gff3(h5$files$query_gff[1])
            i_gff <- .getMiniprotGFF(gff = i_gff,
                                     gene = i_synog$query_gene,
                                     tx = i_synog$query_tx)
            i_cds <- readDNAStringSet(h5$files$query_cds[1])
            i_cds <- i_cds[names(i_cds) %in% i_gff$ID]

        } else {
            i_gff <- import.gff3(h5$files$subject_gff[1])
            i_gff <- .getMiniprotGFF(gff = i_gff,
                                     gene = i_synog$subject_gene,
                                     tx = i_synog$subject_tx)
            i_cds <- readDNAStringSet(h5$files$subject_cds[1])
            i_cds <- i_cds[names(i_cds) %in% i_gff$ID]
        }

        if(i == 1){
            gff <- list(i_gff)
            cds <- list(i_cds)
            synog <- list(i_synog)

            if(i_target_data == "query"){
                genome_fn <- h5$files$query_genome[1]

            } else {
                genome_fn <- h5$files$subject_genome[1]
            }

        } else {
            gff <- c(gff, list(i_gff))
            synog <- c(synog, list(i_synog))
            cds <- c(cds, list(i_cds))
        }

        H5Fclose(h5)
    }
    out <- list(gff = gff, synog = synog, cds = cds, genome_fn = genome_fn)
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
.outputGFFdata <- function(gff, prefix, genome_fn, out_dir, makeFASTA){
    out_fn <- paste0(prefix, "_synog.gff")
    out_gff_fn <- file.path(out_dir, out_fn)
    export.gff3(object = gff, con = out_gff_fn)

    if(makeFASTA){
        out_cds <- .makeCDS(gff = out_gff_fn, genome = genome_fn)
        out_cds_fn <- sub("\\.gff", ".cds", out_gff_fn)
        writeXStringSet(out_cds, out_cds_fn)

        out_prot <- translate(x = out_cds, no.init.codon = TRUE, if.fuzzy.codon = "X")
        out_prot_fn <- sub("\\.gff", ".prot", out_gff_fn)
        writeXStringSet(out_prot, out_prot_fn)
    }
}

.renameMP <- function(gff, pair_genome){
    mp_entries <- grepl("^MP", gff$ID)
    gff$oldID <- NA
    gff$oldID[mp_entries] <- gff$ID[mp_entries]
    gff$ID[mp_entries] <- paste(pair_genome, gff$ID[mp_entries], sep = "_")
    not_gene <- gff$type != "gene"
    not_gene_parent <- unlist(gff$Parent[not_gene])
    mp_entries <- grepl("^MP", not_gene_parent)
    not_gene_parent[mp_entries] <- paste(pair_genome, not_gene_parent[mp_entries], sep = "_")
    gff$Parent[not_gene] <- lapply(not_gene_parent, c)
    mp_entries <- grepl("^MP", gff$gene_id)
    gff$gene_id[mp_entries] <- paste(pair_genome, gff$gene_id[mp_entries], sep = "_")
    gff$Name <- gff$ID
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
    ref_gff2 <- .mergeChimericTx(ref_gff = ref_gff,
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

.renameSynog <- function(id2id_list, hdf5_fn, genomewise_list){
    genomewise_list <- lapply(seq_along(genomewise_list), function(i){
        out <- data.frame(genomewise_list[[i]],
                          target_genome = names(genomewise_list[i]))
        return(out)
    })
    genomewise_list <- do.call("rbind", genomewise_list)
    genomewise_list <- subset(x = genomewise_list, subset = !duplicated(comb_id))
    out_tx <- NULL
    out_gene <- NULL
    for(i in genomewise_list$comb_id){
        h5 <- H5Fopen(name = hdf5_fn[i])
        synog_tx <- h5$synog_tx
        synog_tx <- .renameMPsynog(df = synog_tx, genome = genomewise_list[i, ])

        synog_tx <- .replaceSynogID(df = synog_tx,
                                    id2id_list = id2id_list,
                                    genome = genomewise_list[i, ])
        out_tx <- c(out_tx, list(synog_tx))

        synog_gene <- h5$synog_gene
        synog_gene <- .renameMPsynog(df = synog_gene, genome = genomewise_list[i, ])

        synog_gene <- .replaceSynogID(df = synog_gene,
                                      id2id_list = id2id_list,
                                      genome = genomewise_list[i, ])
        out_gene <- c(out_gene, list(synog_gene))
    }
    out <- list(tx = out_tx, gene = out_gene, genomewise_list = genomewise_list)
    return(out)
}

.renameMPsynog <- function(df, genome){
    if(genome$target == "query"){
        df$query_tx <- .replaceMPid(id = df$query_tx,
                                    prefix = genome$pair_genome)
        df$query_gene <- .replaceMPid(id = df$query_gene,
                                      prefix = genome$pair_genome)
        df$subject_tx <- .replaceMPid(id = df$subject_tx,
                                      prefix = genome$target_genome)
        df$subject_gene <- .replaceMPid(id = df$subject_gene,
                                        prefix = genome$target_genome)

    } else {
        df$query_tx <- .replaceMPid(id = df$query_tx,
                                    prefix = genome$target_genome)
        df$query_gene <- .replaceMPid(id = df$query_gene,
                                      prefix = genome$target_genome)
        df$subject_tx <- .replaceMPid(id = df$subject_tx,
                                      prefix = genome$pair_genome)
        df$subject_gene <- .replaceMPid(id = df$subject_gene,
                                        prefix = genome$pair_genome)
    }
    return(df)
}

.replaceMPid <- function(id, prefix){
    mp_id <- grepl("^MP", id)
    id[mp_id] <- paste(prefix, id[mp_id], sep = "_")
    return(id)
}

.replaceSynogID <- function(df, id2id_list, genome){
    if(genome$target == "query"){
        id2id <- subset(id2id_list[[genome$target_genome]], subset = is.na(mp_gene_id))
        hit <- match(df$query_tx, id2id$mp_tx_id)
        if(any(!is.na(hit))){
            df$query_tx[!is.na(hit)] <- id2id$ref_tx_id[na.omit(hit)]
            df$query_gene[!is.na(hit)] <- id2id$ref_gene_id[na.omit(hit)]
        }

    } else {
        id2id <- subset(id2id_list[[genome$pair_genome]], subset = is.na(mp_gene_id))
        hit <- match(df$subject_tx, id2id$mp_tx_id)
        if(any(!is.na(hit))){
            df$subject_tx[!is.na(hit)] <- id2id$ref_tx_id[na.omit(hit)]
            df$subject_gene[!is.na(hit)] <- id2id$ref_gene_id[na.omit(hit)]
        }
    }

    id2id <- subset(id2id_list[[genome$target_genome]], subset = !is.na(mp_gene_id))
    if(genome$target == "query"){
        hit <- match(df$query_tx, id2id$mp_tx_id)
        if(any(!is.na(hit))){
            df$query_gene[!is.na(hit)] <- id2id$ref_gene_id[na.omit(hit)]
        }

    } else {
        id2id <- subset(id2id_list[[genome$pair_genome]], subset = !is.na(mp_gene_id))
        hit <- match(df$subject_tx, id2id$mp_tx_id)
        if(any(!is.na(hit))){
            df$subject_gene[!is.na(hit)] <- id2id$ref_gene_id[na.omit(hit)]
        }
    }

    return(df)
}

.outputSynogData <- function(synog_list, out_dir, out_fn, overwrite){
    hdf5_path <- file.path(out_dir, out_fn)
    h5_fn <- .makeHDF5(hdf5_path = hdf5_path, overwrite = overwrite)

    .h5creategroup(h5_fn,"synog_tx")
    .h5creategroup(h5_fn,"synog_gene")
    for(i in seq_along(synog_list$genomewise_list$comb_id)){
        data_name <- paste(synog_list$genomewise_list$target_genome[i],
                           synog_list$genomewise_list$pair_genome[i],
                           sep = "_")
        .h5overwrite(obj = synog_list$tx[[i]],
                     file = h5_fn,
                     name = paste0("synog_tx/", data_name))
        .h5overwrite(obj = synog_list$gene[[i]],
                     file = h5_fn,
                     name = paste0("synog_gene/", data_name))
    }
    return(h5_fn)
}
