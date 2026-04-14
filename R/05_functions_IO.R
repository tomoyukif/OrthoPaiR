#' Get genome ID mapping from input folders
#'
#' List sample names and assigned IDs from `working_dir/input/<index>_<sample>/`.
#'
#' @param working_dir Working directory.
#'
#' @return data.frame with `sample`, `id_1digit`, `id_4digit`, and `input_dir`.
#' @export
getGenomeID <- function(working_dir) {
    input_dir <- file.path(working_dir, "input")
    dirs <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)
    if(!length(dirs)) {
        return(data.frame(sample = character(0),
                          id_1digit = integer(0),
                          id_4digit = character(0),
                          input_dir = character(0),
                          stringsAsFactors = FALSE))
    }
    bn <- basename(dirs)
    id1 <- suppressWarnings(as.integer(sub("_.*", "", bn)))
    sample <- sub("^[0-9]+_", "", bn)
    keep <- !is.na(id1)
    out <- data.frame(sample = sample[keep],
                      id_1digit = id1[keep],
                      id_4digit = as.character(id1[keep] + 1000L),
                      input_dir = dirs[keep],
                      stringsAsFactors = FALSE)
    out[order(out$id_1digit), , drop = FALSE]
}

#' Get ortholog table for a genome pair
#'
#' Read pairwise ortholog table from `working_dir/orthopair/<pair>.tsv`.
#'
#' @param working_dir Working directory.
#' @param pair Genome pair (e.g. `c(1001, 1002)` or `"1001_1002"`).
#' @param score Logical; include score-related columns when available.
#' @param loc Logical; include location-related columns when available.
#'
#' @return data.frame of ortholog pairs.
#' @export
getOrthoPair <- function(working_dir, pair, score = FALSE, loc = FALSE) {
    pair_id <- .normalize_pair_id(pair)
    fn <- file.path(working_dir, "orthopair", paste0(pair_id, ".tsv"))
    if(!file.exists(fn)) {
        stop("ortholog pair file not found: ", fn, call. = FALSE)
    }
    dt <- data.table::fread(fn, sep = "\t", header = TRUE)
    base_col <- c("genome1_original_gene", "genome2_original_gene",
                  "original_genome1_gene", "original_genome2_gene",
                  "genome1_gene", "genome1_tx",
                  "genome2_gene", "genome2_tx",
                  "genome1_synteny_block", "genome2_synteny_block",
                  "class", "SOG", "genome1_collapse", "genome2_collapse")
    score_col <- c("pident", "q2s_qcovs", "s2q_qcovs",
                   "q2s_ci", "s2q_ci", "mutual_ci",
                   "genome1_is_anchor", "genome2_is_anchor", "is_anchor_pair")
    loc_col <- c("genome1_chr", "genome1_start", "genome1_end", "genome1_strand",
                 "genome2_chr", "genome2_start", "genome2_end", "genome2_strand")
    target <- base_col
    if(isTRUE(score)) target <- c(target, score_col)
    if(isTRUE(loc)) target <- c(target, loc_col)
    target <- target[target %in% names(dt)]
    out <- as.data.frame(dt[, ..target], stringsAsFactors = FALSE)
    if("genome1_collapse" %in% names(out) && all(out$genome1_collapse == 0, na.rm = TRUE)) {
        out$genome1_collapse <- NA
    }
    if("genome2_collapse" %in% names(out) && all(out$genome2_collapse == 0, na.rm = TRUE)) {
        out$genome2_collapse <- NA
    }
    out
}

#' Get orphan genes for a genome pair
#'
#' Read orphan list from `working_dir/orphan/<pair>.tsv`.
#'
#' @param working_dir Working directory.
#' @param pair Genome pair (e.g. `c(1001, 1002)` or `"1001_1002"`).
#'
#' @return data.frame with columns `genome` and `gene`.
#' @export
getOrphan <- function(working_dir, pair) {
    pair_id <- .normalize_pair_id(pair)
    fn <- file.path(working_dir, "orphan", paste0(pair_id, ".tsv"))
    if(!file.exists(fn)) {
        stop("orphan file not found: ", fn, call. = FALSE)
    }
    as.data.frame(data.table::fread(fn, sep = "\t", header = TRUE), stringsAsFactors = FALSE)
}

#' Summarize ortholog class counts for one pair
#'
#' @param working_dir Working directory.
#' @param pair Genome pair (e.g. `c(1001, 1002)` or `"1001_1002"`).
#'
#' @return data.frame summary by genome side and class.
#' @export
summaryOrthoPair <- function(working_dir, pair) {
    pair_id <- .normalize_pair_id(pair)
    parts <- strsplit(pair_id, "_", fixed = TRUE)[[1L]]
    op <- getOrthoPair(working_dir = working_dir, pair = pair, score = FALSE, loc = FALSE)
    orphan <- getOrphan(working_dir = working_dir, pair = pair)
    
    cls <- c("1to1", "1toM", "Mto1", "MtoM")
    q <- op[!duplicated(op$genome1_gene), , drop = FALSE]
    s <- op[!duplicated(op$genome2_gene), , drop = FALSE]
    q_tbl <- table(factor(q$class, levels = cls))
    s_tbl <- table(factor(s$class, levels = cls))
    q_orphan <- sum(orphan$genome == parts[1L], na.rm = TRUE)
    s_orphan <- sum(orphan$genome == parts[2L], na.rm = TRUE)
    q_total <- sum(q_tbl) + q_orphan
    s_total <- sum(s_tbl) + s_orphan
    
    q_counts <- c(q_total, q_total - q_orphan, as.integer(q_tbl), q_orphan)
    s_counts <- c(s_total, s_total - s_orphan, as.integer(s_tbl), s_orphan)
    q_ratio <- if(q_total > 0) q_counts / q_total else rep(NA_real_, length(q_counts))
    s_ratio <- if(s_total > 0) s_counts / s_total else rep(NA_real_, length(s_counts))
    
    df <- data.frame(
        genome1 = c(parts[1L], q_counts, round(q_ratio, 4)),
        genome2 = c(parts[2L], s_counts, round(s_ratio, 4)),
        stringsAsFactors = FALSE
    )
    rownames(df) <- c(
        "Name", "Total", "Classified", "1to1", "1toM", "Mto1", "MtoM", "Orphan",
        "Total_ratio", "Classified_ratio", "1to1_ratio", "1toM_ratio", "Mto1_ratio", "MtoM_ratio", "Orphan_ratio"
    )
    df
}

.normalize_pair_id <- function(pair) {
    if(length(pair) == 1L) {
        x <- as.character(pair)
        if(grepl("_", x)) {
            p <- strsplit(x, "_", fixed = TRUE)[[1L]]
        } else {
            stop("`pair` must be like c(1001,1002) or \"1001_1002\".", call. = FALSE)
        }
    } else if(length(pair) == 2L) {
        p <- as.character(pair)
    } else {
        stop("`pair` must have two genome IDs.", call. = FALSE)
    }
    p <- sort(p)
    paste(p, collapse = "_")
}
