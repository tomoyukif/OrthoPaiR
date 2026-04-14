#' Create BLAST databases and perform reciprocal BLAST
#'
#' This function prepare BLAST databases and perform BLAST searches reciprocally
#' on given input datasets to create the list of reciprocal BLAST hits.
#'
#' @export
#'
rbh <- function(working_dir, 
                blast_path, 
                target_pair = NULL, 
                n_threads = 1, 
                overwrite = FALSE){
    blast_dir <- file.path(working_dir, "blast")
    dir.create(path = blast_dir, showWarnings = FALSE, recursive = TRUE)
    
    cds_files <- list.files(path = file.path(working_dir, "input"),
                            "cds.fa", 
                            full.names = TRUE,
                            recursive = TRUE)
    dir_names <- list.dirs(path = file.path(working_dir, "input"),
                           recursive = FALSE)
    sample_names <- sub(".+input\\/[0-9]+_", "", dir_names)
    input_list <- list(cds = cds_files, names = sample_names)
    
    tp <- .normalize_target_pair(target_pair = target_pair, 
                                 sample_names = input_list$names)
    .validate_target_pair_samples(target_pair = tp,
                                  sample_names = input_list$names)
    
    cds_sizes <- .cds_file_sizes(input_list)
    scheme <- .plan_blast_scheme_from_target_pair(target_pair = tp,
                                                  cds_sizes = cds_sizes)
    
    db_map <- .makeBlastDB_from_scheme(scheme = scheme,
                                       input_list = input_list,
                                       blast_dir = blast_dir,
                                       blast_path = blast_path,
                                       overwrite = overwrite)
    
    check_blast_out <- .checkBLASTout_target(blast_jobs = scheme$blast_jobs,
                                             db_map = db_map,
                                             blast_dir = blast_dir,
                                             overwrite = overwrite)
    
    job_assign <- .jobAssign(check_blast_out = check_blast_out,
                             input_list = input_list,
                             n_threads = n_threads,
                             min_threads = 8L)
    
    if(nrow(job_assign$fasta_chunk) > 0L) {
        mclapply(X = unique(job_assign$fasta_chunk$chunk),
                 mc.cores = job_assign$n_parallel_jobs,
                 FUN = .blastn_search,
                 input_list = input_list,
                 blast_dir = blast_dir,
                 blast_path = blast_path,
                 fasta_chunk = job_assign$fasta_chunk,
                 threads_per_job = job_assign$threads_per_job,
                 index_offset = check_blast_out$blast_out_index_offset)
    }
    
    rbh_dir <- file.path(working_dir, "rbh")
    dir.create(rbh_dir, showWarnings = FALSE, recursive = TRUE)
    tmp_dir <- file.path(rbh_dir, "tmp")
    dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
    sorttmp_dir <- file.path(rbh_dir, "sort_tmp")
    dir.create(sorttmp_dir, showWarnings = FALSE, recursive = TRUE)
    out_rbh_fn <- file.path(rbh_dir, "all.rbh.tsv")
    blast_out_list <- list.files(file.path(working_dir, "blast"), "blast.out$",
                                 full.names = TRUE,
                                 recursive = TRUE)
    rbh_extract(blast_out_list = blast_out_list,
                out_rbh_fn = out_rbh_fn,
                strategy = "auto",
                n_threads = n_threads,
                tmpdir = tmp_dir,
                sort_tmp = sorttmp_dir,
                final_sort_by_mutual_ci = TRUE,
                keep_intermediate = TRUE)
    
    .splitRBHbyGenomePair(rbh_fn = out_rbh_fn,
                          rbh_dir = rbh_dir,
                          genome_width = 4L)
    invisible(TRUE)
}

.normalize_target_pair <- function(target_pair, sample_names) {
    if(is.null(target_pair)) {
        message("Process all possible pairs.")
        target_pair <- t(combn(x = sample_names, m = 2))
    }
    if(is.data.frame(target_pair)) {
        target_pair <- as.matrix(target_pair)
    }
    if(!is.matrix(target_pair)) {
        stop("`target_pair` must be a matrix (N x 2) or a data.frame convertible to a matrix.", call. = FALSE)
    }
    if(ncol(target_pair) != 2L) {
        stop("`target_pair` must have exactly 2 columns.", call. = FALSE)
    }
    storage.mode(target_pair) <- "character"
    target_pair <- trimws(target_pair)
    if(anyNA(target_pair) || any(!nzchar(target_pair))) {
        stop("`target_pair` must not contain NA or empty strings.", call. = FALSE)
    }
    keep <- target_pair[, 1] != target_pair[, 2]
    target_pair <- target_pair[keep, , drop = FALSE]
    if(!nrow(target_pair)) {
        stop("`target_pair` contains no valid pairs after removing self-pairs.", call. = FALSE)
    }
    return(target_pair)
}

.validate_target_pair_samples <- function(target_pair, sample_names) {
    samples <- unique(c(target_pair[, 1], target_pair[, 2]))
    missing <- setdiff(samples, sample_names)
    if(length(missing)) {
        stop(
            "The following samples in `target_pair` were not found in `working_dir/input/*_<sample>/`:\n",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }
    invisible(TRUE)
}

.cost_db <- function(sum_db_bytes) {
    as.numeric(sum_db_bytes)
}

.cost_blast <- function(sum_query_bytes, sum_db_bytes) {
    as.numeric(sum_query_bytes) * as.numeric(sum_db_bytes)
}

.connected_components_undirected <- function(vertices, edges) {
    vertices <- sort(unique(as.character(vertices)))
    if(!length(vertices)) {
        return(data.frame(vertex = character(), comp_id = integer(), stringsAsFactors = FALSE))
    }
    parent <- seq_along(vertices)
    names(parent) <- vertices
    
    find_root <- function(v) {
        i <- match(v, vertices)
        while(parent[[i]] != i) {
            parent[[i]] <- parent[[parent[[i]]]]
            i <- parent[[i]]
        }
        i
    }
    
    union <- function(v1, v2) {
        r1 <- find_root(v1)
        r2 <- find_root(v2)
        if(r1 != r2) {
            parent[[r2]] <<- r1
        }
    }
    
    if(nrow(edges)) {
        for(i in seq_len(nrow(edges))) {
            union(as.character(edges$u[i]), as.character(edges$v[i]))
        }
    }
    
    roots <- vapply(vertices, function(v) vertices[[find_root(v)]], character(1L))
    comp_levels <- unique(roots)
    comp_id <- match(roots, comp_levels)
    data.frame(vertex = vertices, comp_id = as.integer(comp_id), stringsAsFactors = FALSE)
}

.plan_blast_scheme_from_target_pair <- function(target_pair, cds_sizes) {
    a <- as.character(target_pair[, 1])
    b <- as.character(target_pair[, 2])
    lo <- pmin(a, b)
    hi <- pmax(a, b)
    edges <- unique(data.frame(u = lo, v = hi, stringsAsFactors = FALSE))
    
    comps <- .connected_components_undirected(vertices = unique(c(edges$u, edges$v)), edges = edges)
    comp_ids <- unique(comps$comp_id)
    
    db_sets <- list()
    blast_jobs <- list()
    comp_plan <- list()
    
    for(cid in comp_ids) {
        V <- comps$vertex[comps$comp_id == cid]
        V <- sort(unique(as.character(V)))
        if(length(V) < 2L) next
        
        edges_c <- edges[edges$u %in% V & edges$v %in% V, , drop = FALSE]
        if(!nrow(edges_c)) next
        
        sumV <- sum(as.numeric(cds_sizes[V]), na.rm = TRUE)
        costA <- .cost_db(sumV) + .cost_blast(sumV, sumV)
        best <- list(type = "all_in_one", cost = costA, hub = NA_character_)
        
        for(hub in V) {
            rest <- setdiff(V, hub)
            if(!length(rest)) next
            
            edges_rest <- edges_c[edges_c$u %in% rest & edges_c$v %in% rest, , drop = FALSE]
            if(nrow(edges_rest) > 0L) next
            
            sumHub <- as.numeric(cds_sizes[[hub]])
            sumRest <- sum(as.numeric(cds_sizes[rest]), na.rm = TRUE)
            costB <- .cost_db(sumHub) + .cost_db(sumRest) +
                .cost_blast(sumHub, sumRest) + .cost_blast(sumRest, sumHub)
            if(is.finite(costB) && costB < best$cost) {
                best <- list(type = "hub_rest", cost = costB, hub = hub)
            }
        }
        
        if(identical(best$type, "all_in_one")) {
            db_sets[[length(db_sets) + 1L]] <- V
            db_id <- length(db_sets)
            blast_jobs[[length(blast_jobs) + 1L]] <- list(query = V, db_set_id = db_id)
        } else {
            hub <- best$hub
            rest <- setdiff(V, hub)
            db_sets[[length(db_sets) + 1L]] <- hub
            db_hub_id <- length(db_sets)
            db_sets[[length(db_sets) + 1L]] <- rest
            db_rest_id <- length(db_sets)
            blast_jobs[[length(blast_jobs) + 1L]] <- list(query = hub, db_set_id = db_rest_id)
            blast_jobs[[length(blast_jobs) + 1L]] <- list(query = rest, db_set_id = db_hub_id)
        }
        
        comp_plan[[as.character(cid)]] <- best
    }
    
    if(!length(db_sets) || !length(blast_jobs)) {
        stop("No BLAST jobs planned from `target_pair` (check pairs).", call. = FALSE)
    }
    
    list(db_sets = db_sets, blast_jobs = blast_jobs, comp_plan = comp_plan)
}

.splitRBHbyGenomePair <- function(rbh_fn, rbh_dir, genome_width = 4L){
    if(!file.exists(rbh_fn)) {
        warning("RBH file not found: ", rbh_fn, ". Skipping split by genome pairs.")
        return(invisible(NULL))
    }
    
    # Create rbh output directory
    dir.create(path = rbh_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Delete any existing .rbh files in the output directory before writing new ones
    existing_rbh <- list.files(path = rbh_dir, pattern = "\\.rbh$", full.names = TRUE)
    if(length(existing_rbh) > 0L) {
        unlink(existing_rbh, force = TRUE)
    }
    
    # Use awk for fast single-pass splitting (much faster than R line-by-line processing)
    # Columns: query_tx, subject_tx, pident, q2s_qcovs, s2q_qcovs, q2s_ci, s2q_ci, mutual_ci
    sys <- Sys.info()[["sysname"]]
    if(!sys %in% c("Linux", "Darwin")) {
        stop("Linux/macOS only (requires awk)")
    }
    
    if(Sys.which("awk") == "") {
        stop("Required command not found in PATH: awk")
    }
    
    # Normalize paths
    rbh_fn <- normalizePath(rbh_fn, mustWork = TRUE)
    rbh_dir <- normalizePath(rbh_dir, mustWork = TRUE)
    
    # awk script: extract genome IDs from first genome_width chars, create pair ID, write to file
    awk_script <- sprintf(
        'BEGIN{FS="\\t"; OFS="\\t"; gw=%d; outdir="%s"}
NF >= 2 {
  gA = substr($1, 1, gw);
  gB = substr($2, 1, gw);
  if(gA <= gB) {
    pair_id = gA "_" gB;
  } else {
    pair_id = gB "_" gA;
  }
  fn = sprintf("%%s/%%s.rbh", outdir, pair_id);
  print $0 >> fn;
}
', as.integer(genome_width), gsub('"', '\\"', rbh_dir))
    
    awk_fn <- tempfile("split_rbh_", tmpdir = dirname(rbh_fn), fileext = ".awk")
    writeLines(awk_script, awk_fn)
    on.exit(unlink(awk_fn, force = TRUE), add = TRUE)
    
    # Run awk in single pass
    rc <- system2("awk", c("-f", awk_fn, rbh_fn), stdout = TRUE, stderr = TRUE)
    if(!is.null(attr(rc, "status")) && attr(rc, "status") != 0L) {
        stop("awk split failed:\n", paste(rc, collapse = "\n"))
    }
    
    invisible(TRUE)
}

.makeBlastDB_from_scheme <- function(scheme,
                                     input_list,
                                     blast_dir,
                                     blast_path,
                                     overwrite) {
    stopifnot(is.list(scheme), !is.null(scheme$db_sets))
    dir.create(blast_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Existing DB definitions
    existing_list_files <- list.files(blast_dir, pattern = "_blastdb\\.list$", full.names = TRUE, recursive = TRUE)
    existing_sets <- list()
    if(length(existing_list_files)) {
        existing_sets <- lapply(existing_list_files, function(fn) {
            x <- try(fread(file = fn, sep = "\t", header = FALSE, stringsAsFactors = FALSE), silent = TRUE)
            if(inherits(x, "try-error")) return(character(0))
            sort(unique(as.character(unlist(x))))
        })
        names(existing_sets) <- existing_list_files
    }
    
    # index offset
    idx_offset <- 0L
    if(length(existing_list_files)) {
        idx <- suppressWarnings(as.integer(sub("_.+", "", basename(existing_list_files))))
        idx_offset <- max(idx, na.rm = TRUE)
        if(!is.finite(idx_offset)) idx_offset <- 0L
    }
    
    db_map <- vector("list", length(scheme$db_sets))
    for(i in seq_along(scheme$db_sets)) {
        set_i <- sort(unique(as.character(scheme$db_sets[[i]])))
        if(!length(set_i)) stop("Empty db_set in scheme at index ", i, call. = FALSE)
        
        # Try reuse exact-match DB if overwrite==FALSE
        reused <- FALSE
        if(!overwrite && length(existing_sets)) {
            hit <- vapply(existing_sets, function(s) identical(s, set_i), logical(1L))
            if(any(hit)) {
                list_fn <- names(existing_sets)[which(hit)[1L]]
                cds_fn <- sub("_blastdb\\.list$", ".fa", list_fn)
                # Old naming used all_cds.fa; also allow that
                if(!file.exists(cds_fn)) {
                    alt <- sub("_blastdb\\.list$", "_cds.fa", list_fn)
                    if(file.exists(alt)) cds_fn <- alt
                }
                db_prefix <- sub("\\.fa$", ".blastdb", cds_fn)
                db_map[[i]] <- list(
                    set = set_i,
                    cds_fn = cds_fn,
                    db_prefix = db_prefix,
                    blastdb_list_fn = list_fn,
                    reused = TRUE
                )
                reused <- TRUE
            }
        }
        
        if(!reused) {
            idx_offset <- idx_offset + 1L
            cds_fn <- file.path(blast_dir, paste0(idx_offset, "_group_cds.fa"))
            cds_fn_list <- input_list$cds[input_list$names %in% set_i]
            cds_fn_list <- normalizePath(cds_fn_list, mustWork = TRUE)
            
            sys <- Sys.info()[["sysname"]]
            if(sys %in% c("Linux", "Darwin")) {
                status <- system2("cat", cds_fn_list, stdout = cds_fn)
                if(!identical(status, 0L)) {
                    stop("Merge of CDS FASTA files failed for DB set ", i, " (cat exit ", status, ")", call. = FALSE)
                }
            } else {
                fa <- readDNAStringSet(cds_fn_list)
                writeXStringSet(fa, cds_fn)
            }
            
            .blastDBengine(cds_fn = cds_fn, blast_dir = blast_dir, blast_path = blast_path, overwrite = overwrite)
            blastdb_list_fn <- sub("\\.fa$", "_blastdb.list", cds_fn)
            write.table(
                x = set_i,
                file = blastdb_list_fn,
                quote = FALSE,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE
            )
            
            db_map[[i]] <- list(
                set = set_i,
                cds_fn = cds_fn,
                db_prefix = sub("\\.fa$", ".blastdb", cds_fn),
                blastdb_list_fn = blastdb_list_fn,
                reused = FALSE
            )
        }
    }
    
    names(db_map) <- paste0("db_set_", seq_along(db_map))
    db_map
}

.checkBLASTout_target <- function(blast_jobs, db_map, blast_dir, overwrite) {
    if(!length(blast_jobs)) {
        return(list(to_be_blast = data.frame(genome = character(), db = character(), stringsAsFactors = FALSE),
                    blast_out_index_offset = 0L))
    }
    
    # Build desired genome->db tasks
    rows <- lapply(blast_jobs, function(job) {
        q <- as.character(job$query)
        sid <- as.integer(job$db_set_id)
        if(!length(q)) return(NULL)
        if(is.na(sid) || sid < 1L || sid > length(db_map)) {
            stop("Invalid db_set_id in blast_jobs: ", sid, call. = FALSE)
        }
        db_prefix <- db_map[[sid]]$db_prefix
        data.frame(genome = q, db = rep(db_prefix, length(q)), stringsAsFactors = FALSE)
    })
    to_be_blast <- do.call(rbind, rows)
    if(is.null(to_be_blast) || !nrow(to_be_blast)) {
        return(list(to_be_blast = data.frame(genome = character(), db = character(), stringsAsFactors = FALSE),
                    blast_out_index_offset = 0L))
    }
    to_be_blast <- unique(to_be_blast)
    
    blast_out_list <- list.files(blast_dir, "blast\\.out\\.list$", full.names = TRUE, recursive = TRUE)
    blast_out_index_offset <- 0L
    if(length(blast_out_list) > 0L && !overwrite) {
        done <- lapply(blast_out_list, function(x) {
            x_out <- fread(file = x, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
            as.data.frame(x_out)
        })
        done <- do.call(rbind, done)
        if(ncol(done) >= 2L) {
            done <- done[, 1:2, drop = FALSE]
            names(done) <- c("genome", "db")
            to_be_id <- paste(to_be_blast$genome, to_be_blast$db, sep = "_")
            done_id <- paste(done$genome, done$db, sep = "_")
            to_be_blast <- to_be_blast[!to_be_id %in% done_id, , drop = FALSE]
        }
        idx <- suppressWarnings(as.integer(sub("_.+", "", basename(blast_out_list))))
        blast_out_index_offset <- max(idx, na.rm = TRUE)
        if(!is.finite(blast_out_index_offset)) blast_out_index_offset <- 0L
    }
    
    list(to_be_blast = to_be_blast, blast_out_index_offset = as.integer(blast_out_index_offset))
}

.blastDBengine <- function(cds_fn, blast_dir, blast_path, overwrite){
    out_prefix <- sub(".fa", "", cds_fn)
    job_finish_marker <- paste0(out_prefix, ".makeblastdb.finish")
    if(overwrite){
        exec <- TRUE
        
    } else {
        exec <- !file.exists(job_finish_marker)
    }
    unlink(job_finish_marker, force = TRUE)
    
    if(exec){
        blast_args <- paste("-in", cds_fn, 
                            "-dbtype nucl", 
                            "-out", paste0(out_prefix, ".blastdb"))
        
        out <- try({
            system2(command = file.path(blast_path, "makeblastdb"),
                    args = blast_args, 
                    stdout = FALSE)
        }, silent = TRUE)
        
        if(inherits(out, "try-error")){
            write(out, paste0(out_prefix, ".makeblastdb.error"))
            
        } else {
            write("", job_finish_marker)
        }
    }
    invisible(TRUE)
}

## Sum of on-disk byte sizes for all BLAST DB files sharing db_prefix (e.g. *.ndb, *.nin, ...).
.blastdb_disk_bytes <- function(db_prefix) {
    dir <- dirname(db_prefix)
    base <- basename(db_prefix)
    if(!nzchar(base) || !dir.exists(dir)) {
        return(0)
    }
    pat <- paste0("^", gsub("([.|()\\^{}$+\\[\\]\\\\])", "\\\\\\1", base), "\\.")
    fn <- list.files(dir, pattern = pat, full.names = TRUE)
    if(!length(fn)) {
        return(0)
    }
    sum(file.size(fn))
}

## Per-genome CDS file sizes (bytes), named by sample name.
.cds_file_sizes <- function(input_list) {
    nm <- input_list$names
    cds <- input_list$cds
    if(length(nm) != length(cds)) {
        stop("input_list$names and input_list$cds must have the same length.")
    }
    sz <- file.size(cds)
    names(sz) <- nm
    sz[is.na(sz)] <- 0
    sz
}

## Greedy packing: keep sum(S_q) <= S_q_max per chunk; order follows `genomes` order.
.chunk_genomes_by_Sq_max <- function(genomes, cds_sizes, S_q_max) {
    if(!length(genomes)) {
        return(integer(0))
    }
    if(!length(S_q_max) || !is.finite(S_q_max) || S_q_max <= 0) {
        return(seq_along(genomes))
    }
    chunk_id <- integer(length(genomes))
    cur <- 1L
    sum_s <- 0
    for(i in seq_along(genomes)) {
        g <- genomes[[i]]
        sz <- as.numeric(cds_sizes[[g]])
        if(!length(sz) || is.na(sz)) {
            sz <- 0
        }
        if(sum_s > 0 && sum_s + sz > S_q_max) {
            cur <- cur + 1L
            sum_s <- 0
        }
        chunk_id[[i]] <- cur
        sum_s <- sum_s + sz
    }
    chunk_id
}

.jobAssign <- function(check_blast_out,
                       input_list,
                       n_threads,
                       min_threads = 8L,
                       k_db = 1,
                       k_q = 1.75,
                       M_fixed_bytes = 384 * 1024^2,
                       frac_ma = 0.25,
                       max_threads_blastn = 16L,
                       verbose = TRUE) {
    ma <- .mem_available_bytes()
    margin <- max(512 * 1024^2, 0.05 * ma)
    ma_budget <- ma - margin
    if(!is.finite(ma_budget) || ma_budget <= 0) {
        warning("MemAvailable margin leaves no budget; using single-threaded BLAST.")
        ma_budget <- max(1, ma * 0.5)
    }
    M_q_max <- frac_ma * ma
    S_q_max <- M_q_max / k_q
    
    cds_sizes <- .cds_file_sizes(input_list)
    
    if(!length(check_blast_out$to_be_blast$genome)) {
        return(list(
            fasta_chunk = data.frame(
                chunk = integer(),
                genome = character(),
                db = character(),
                stringsAsFactors = FALSE
            ),
            n_parallel_jobs = 1L,
            threads_per_job = max(1L, as.integer(n_threads))
        ))
    }
    
    fasta_chunk <- tapply(seq_along(check_blast_out$to_be_blast$genome),
                          check_blast_out$to_be_blast$db,
                          function(i) {
                              x <- as.character(check_blast_out$to_be_blast$genome[i])
                              db_name <- check_blast_out$to_be_blast$db[i][1]
                              chunk_loc <- .chunk_genomes_by_Sq_max(x, cds_sizes, S_q_max)
                              out <- data.frame(chunk = as.integer(chunk_loc),
                                                genome = x,
                                                db = db_name,
                                                stringsAsFactors = FALSE)
                              out
                          })
    fasta_chunk <- do.call(rbind, fasta_chunk)
    rownames(fasta_chunk) <- NULL
    fasta_chunk_id <- paste(fasta_chunk$chunk, fasta_chunk$db, sep = "_")
    fasta_chunk_id <- factor(fasta_chunk_id)
    fasta_chunk$chunk <- as.integer(fasta_chunk_id)
    n_chunk <- max(fasta_chunk$chunk)
    
    ## M_job per global chunk, J_mem
    M_worst <- 0
    for(ch in seq_len(n_chunk)) {
        rows <- fasta_chunk$chunk == ch
        db_pref <- fasta_chunk$db[rows][1]
        S_db <- .blastdb_disk_bytes(as.character(db_pref))
        genomes_ch <- as.character(fasta_chunk$genome[rows])
        S_q <- sum(as.numeric(cds_sizes[genomes_ch]), na.rm = TRUE)
        M_job <- k_db * S_db + k_q * S_q + M_fixed_bytes
        if(M_job > M_worst) {
            M_worst <- M_job
        }
    }
    if(M_worst <= 0) {
        J_mem <- Inf
    } else {
        J_mem <- floor(ma_budget / M_worst)
    }
    J_cap_display <- if(!is.finite(J_mem)) {
        as.integer(n_chunk)
    } else {
        max(1L, min(as.integer(J_mem), as.integer(n_chunk)))
    }
    
    if(is.finite(J_mem) && J_mem < 1L) {
        warning(
            "Estimated RAM allows <1 parallel BLAST job (J_mem=", J_mem,
            "). Forcing single job; consider smaller chunks or more memory."
        )
        threads_per_job <- min(as.integer(n_threads), as.integer(max_threads_blastn))
        n_parallel_jobs <- 1L
        
    } else if(n_threads < min_threads) {
        threads_per_job <- as.integer(n_threads)
        n_parallel_jobs <- 1L
        
    } else {
        if(!is.finite(J_mem)) {
            J_cap <- as.integer(n_chunk)
        } else {
            J_cap <- max(1L, min(as.integer(J_mem), as.integer(n_chunk)))
        }
        threads_per_job <- max(
            min_threads,
            min(max_threads_blastn, floor(n_threads / J_cap))
        )
        n_parallel_jobs <- min(J_cap, max(1L, floor(n_threads / threads_per_job)))
        ## Reconcile oversubscription of logical CPUs
        while(n_parallel_jobs * threads_per_job > n_threads) {
            if(threads_per_job > min_threads) {
                threads_per_job <- threads_per_job - 1L
            } else if(n_parallel_jobs > 1L) {
                n_parallel_jobs <- n_parallel_jobs - 1L
            } else {
                break
            }
        }
    }
    
    if(isTRUE(verbose)) {
        message(sprintf(
            "[rbh] blastn schedule: MemAvailable~%.1f GiB, M_worst~%.2f MiB, J_mem=%s, J_cap=%d, n_chunk=%d, threads_per_job=%d, n_parallel_jobs=%d",
            ma / 1024^3,
            M_worst / 1024^2,
            if(M_worst <= 0) {
                "Inf"
            } else {
                as.character(J_mem)
            },
            J_cap_display,
            n_chunk,
            threads_per_job,
            n_parallel_jobs
        ))
    }
    
    out <- list(fasta_chunk = fasta_chunk,
                n_parallel_jobs = as.integer(n_parallel_jobs),
                threads_per_job = as.integer(threads_per_job))
    return(out)
}

.blastn_search <- function(index,
                           input_list, 
                           blast_dir,
                           blast_path,
                           fasta_chunk, 
                           threads_per_job,
                           index_offset){
    db_fn <- fasta_chunk$db[fasta_chunk$chunk == index][1]
    db_prefix <- sub(".ndb", "", db_fn)
    fasta_chunk_i <- fasta_chunk$genome[fasta_chunk$chunk == index]
    index <- index_offset + index
    query_cds_list <- input_list$cds[input_list$names %in% fasta_chunk_i]
    query_fn <- file.path(blast_dir, paste0(index, "_query_cds.fa"))
    
    query_cds_list <- normalizePath(query_cds_list, mustWork = TRUE)
    sys <- Sys.info()[["sysname"]]
    if(!sys %in% c("Linux", "Darwin")) {
        stop("Linux/macOS only (requires cat).")
    }
    status <- system2("cat", query_cds_list, stdout = query_fn)
    if(!identical(status, 0L)) {
        stop("Merge of query CDS FASTA files failed (cat exited with status ", status, ")")
    }
    
    n_genome <- length(input_list$names)
    blast_out <- .run_blastn(query_fn = query_fn, 
                             db_prefix = db_prefix,
                             blast_dir = blast_dir,
                             blast_path = blast_path,
                             threads_per_job = threads_per_job,
                             max_target_seqs = n_genome * 100,
                             index = index)
    if(isTRUE(blast_out)){
        blast_out_list_fn <- file.path(blast_dir, 
                                       paste0(index, "_blast.out.list"))
        blast_out_list <- data.frame(fasta_chunk_i, db_fn)
        fwrite(x = blast_out_list,
               file = blast_out_list_fn,
               quote = FALSE, 
               sep = "\t", 
               row.names = FALSE,
               col.names = FALSE)
    }
    return(blast_out)
}

.run_blastn <- function(query_fn,
                        db_prefix, 
                        blast_dir, 
                        blast_path, 
                        threads_per_job, 
                        max_target_seqs,
                        index){
    out_fn <- file.path(blast_dir, paste0(index, "_blast.out"))
    options(scipen = 9999)
    blast_args <- paste("-query", query_fn,
                        "-db", db_prefix,
                        "-task blastn", 
                        "-max_target_seqs", max_target_seqs,
                        "-evalue 1e-4",
                        "-strand plus",
                        "-outfmt '6 qseqid sseqid pident qstart qend sstart send qlen qcovs'",
                        "-num_threads", threads_per_job,
                        "-out", out_fn)
    out <- try({
        system2(command = file.path(blast_path, "blastn"),
                args = blast_args)
    }, silent = TRUE)
    
    if(inherits(out, "try-error")){
        write(out, file.path(blast_dir, paste0(index, "_blast.out.error")))
        
    } else {
        out <- TRUE
    }
    return(out)
}

rbh_extract <- function(blast_out_list,
                        out_rbh_fn,
                        genome_width = 4L,
                        n_threads = detectCores(),
                        tmpdir = tempdir(),
                        sort_tmp = tmpdir,
                        sort_mem = NULL,
                        strategy = c("auto", "global_sort", "bucket"),
                        n_buckets = 256L,
                        keep_intermediate = FALSE,
                        final_sort_by_mutual_ci = FALSE,
                        verbose = TRUE) {
    strategy <- match.arg(strategy)
    stopifnot(length(blast_out_list) >= 1L)
    stopifnot(all(file.exists(blast_out_list)))
    
    sys <- Sys.info()[["sysname"]]
    if(!sys %in% c("Linux", "Darwin")) {
        stop("Linux/macOS only (mclapply + sort/awk/xargs).")
    }
    
    dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
    dir.create(sort_tmp, showWarnings = FALSE, recursive = TRUE)
    tmpdir <- normalizePath(tmpdir, mustWork = TRUE)
    sort_tmp <- normalizePath(sort_tmp, mustWork = TRUE)
    
    # ------------ helpers ------------
    
    assert_cmd <- function(cmd) {
        if(Sys.which(cmd) == "") stop("Required command not found in PATH: ", cmd)
    }
    
    get_free_bytes <- function(path) {
        out <- suppressWarnings(system2("df", c("-B1", path), stdout = TRUE, stderr = TRUE))
        if(length(out) < 2) return(NA_real_)
        cols <- strsplit(out[2], "[[:space:]]+")[[1]]
        cols <- cols[nzchar(cols)]
        avail <- suppressWarnings(as.numeric(cols[4]))
        if(is.na(avail)) NA_real_ else avail
    }
    
    guess_sort_mem <- function(default = "4G") {
        if(identical(sys, "Linux") && file.exists("/proc/meminfo")) {
            mi <- readLines("/proc/meminfo", warn = FALSE)
            ma <- grep("^MemAvailable:", mi, value = TRUE)
            if(length(ma)) {
                kb <- suppressWarnings(as.numeric(sub("MemAvailable:\\s+([0-9]+)\\s+kB.*", "\\1", ma)))
                if(!is.na(kb)) {
                    bytes <- kb * 1024 * 0.25
                    gb <- floor(bytes / 1024^3)
                    gb <- max(1, gb)
                    return(sprintf("%dG", gb))
                }
            }
        }
        default
    }
    
    normalize_many <- function(blast_out_list, norm_dir) {
        dir.create(norm_dir, showWarnings = FALSE, recursive = TRUE)
        norm_files <- file.path(norm_dir, sprintf("norm_%05d.tsv", seq_along(blast_out_list)))
        for(i in seq_along(blast_out_list)) {
            if(verbose) message(sprintf("[norm %d/%d] %s", i, length(blast_out_list), basename(blast_out_list[i])))
            normalize_blast_tsv_cpp(
                infile = blast_out_list[i],
                outfile = norm_files[i],
                genome_width = genome_width,
                drop_self = TRUE
            )
        }
        norm_files
    }
    
    # Write NUL-separated file list for xargs -0
    write_nul_list <- function(files, list_fn) {
        con <- file(list_fn, open = "wb")
        on.exit(close(con), add = TRUE)
        for(f in files) {
            writeBin(charToRaw(f), con)
            writeBin(as.raw(0), con)
        }
        invisible(list_fn)
    }
    
    # ---- robust sort runner: cat inputs | sort > out ----
    # Use bash 'while read -d ""' to feed files to cat (avoids xargs -0 issues on some systems)
    run_sort_cmd <- function(in_files, out_file, sort_mem, sort_tmp, n_threads_sort) {
        assert_cmd("sort"); assert_cmd("cat"); assert_cmd("bash")
        key_part <- "-k1,1 -k2,2"
        in_files <- vapply(in_files, function(f) normalizePath(f, mustWork = TRUE), character(1L))
        list_fn <- tempfile("sort_inputs_", tmpdir = tmpdir, fileext = ".nul")
        write_nul_list(in_files, list_fn)
        
        sort_part <- if(identical(sys, "Linux")) {
            if(is.null(sort_mem)) sort_mem <- guess_sort_mem()
            sprintf("LC_ALL=C sort -t $'\\t' -S %s -T %s --parallel=%d %s",
                    sort_mem, shQuote(sort_tmp), as.integer(n_threads_sort), key_part)
        } else {
            sprintf("LC_ALL=C sort -t $'\\t' %s", key_part)
        }
        script <- sprintf(
            "while IFS= read -r -d '' f; do cat \"$f\"; done < %s | %s > %s",
            shQuote(list_fn), sort_part, shQuote(out_file)
        )
        on.exit(unlink(list_fn, force = TRUE), add = TRUE)
        
        rc <- system2("bash", c("-c", shQuote(script)), stdout = TRUE, stderr = TRUE)
        if(!is.null(attr(rc, "status")) && attr(rc, "status") != 0L) {
            stop("sort pipeline failed. ", paste(rc, collapse = "\n"))
        }
        
        out_file
    }
    
    final_sort_rbh <- function(rbh_fn) {
        if(!final_sort_by_mutual_ci) return(invisible(rbh_fn))
        assert_cmd("sort"); assert_cmd("bash")
        orig_size <- file.info(rbh_fn)$size
        if(is.na(orig_size) || orig_size == 0L) return(invisible(rbh_fn))
        tmp <- paste0(rbh_fn, ".tmp_sorted")
        script <- sprintf("LC_ALL=C sort -t $'\\t' -k8,8gr %s > %s",
                          shQuote(rbh_fn), shQuote(tmp))
        rc <- system2("bash", c("-c", shQuote(script)), stdout = TRUE, stderr = TRUE)
        if(!is.null(attr(rc, "status")) && attr(rc, "status") != 0L) {
            stop("final sort failed: ", paste(rc, collapse = "\n"))
        }
        tmp_size <- file.info(tmp)$size
        if(is.na(tmp_size) || tmp_size == 0L) {
            unlink(tmp, force = TRUE)
            stop("final sort produced empty output; ", rbh_fn, " left unchanged. ",
                 "Check that the file has 12 tab-separated columns (column 12 = mutual_ci).")
        }
        ok <- file.rename(tmp, rbh_fn)
        if(!ok) {
            file.copy(tmp, rbh_fn, overwrite = TRUE)
            unlink(tmp, force = TRUE)
        }
        invisible(rbh_fn)
    }
    
    bucketize_norm_files <- function(norm_files, bucket_dir, n_buckets, overwrite = TRUE) {
        assert_cmd("awk")
        dir.create(bucket_dir, showWarnings = FALSE, recursive = TRUE)
        if(overwrite) {
            old <- list.files(bucket_dir, full.names = TRUE, pattern = "\\.tsv$")
            if(length(old)) unlink(old)
        }
        
        awk_script <- sprintf(
            'BEGIN{FS="\\t"; OFS="\\t"; nb=%d; outdir="%s"}
{
  gA = substr($1,1,4) + 0;
  gB = substr($2,1,4) + 0;
  b = (gA*10000 + gB) %% nb;
  fn = sprintf("%%s/bucket_%%03d.tsv", outdir, b);
  print $0 >> fn;
}
', as.integer(n_buckets), gsub('"', '\\"', bucket_dir))
        
        awk_fn <- file.path(bucket_dir, "bucketize.awk")
        writeLines(awk_script, awk_fn)
        
        if(verbose) message(sprintf("[bucketize] distributing into %d buckets", n_buckets))
        rc <- system2("awk", c("-f", awk_fn, norm_files), stdout = TRUE, stderr = TRUE)
        if(!is.null(attr(rc, "status")) && attr(rc, "status") != 0) {
            stop("awk bucketize failed:\n", paste(rc, collapse = "\n"))
        }
        
        unlink(awk_fn)
        invisible(TRUE)
    }
    
    process_one_bucket <- function(bucket_fn, out_bucket_fn, sort_tmp, sort_mem, sort_parallel = 1L) {
        fi <- file.info(bucket_fn)
        if(is.na(fi$size) || fi$size == 0) return(NA_character_)
        sorted_fn <- sub("\\.tsv$", ".sorted.tsv", bucket_fn)
        
        run_sort_cmd(in_files = bucket_fn,
                     out_file = sorted_fn,
                     sort_mem = sort_mem,
                     sort_tmp = sort_tmp,
                     n_threads_sort = sort_parallel)
        
        rbh_from_sorted_norm_cpp(sorted_fn, out_bucket_fn)
        unlink(sorted_fn)
        out_bucket_fn
    }
    
    merge_files_xargs_cat <- function(files, out_file) {
        assert_cmd("xargs"); assert_cmd("cat")
        list_fn <- tempfile("merge_inputs_", tmpdir = tmpdir, fileext = ".nul")
        write_nul_list(files, list_fn)
        cmd <- sprintf("xargs -0 cat < %s > %s", shQuote(list_fn), shQuote(out_file))
        rc <- system(cmd, ignore.stdout = TRUE, ignore.stderr = FALSE)
        unlink(list_fn)
        if(!identical(rc, 0L)) stop("merge(cat) failed (exit code ", rc, "):\n", cmd)
        out_file
    }
    
    # ------------ decide auto strategy ------------
    
    total_bytes <- sum(file.info(blast_out_list)$size, na.rm = TRUE)
    if(is.null(sort_mem)) sort_mem <- guess_sort_mem()
    
    if(strategy == "auto") {
        free_bytes <- get_free_bytes(sort_tmp)
        need <- total_bytes * 2
        if(!is.na(free_bytes) && free_bytes < need) {
            if(verbose) message(sprintf("[auto] tmp free %.1f GB < need %.1f GB -> bucket",
                                         free_bytes/1024^3, need/1024^3))
            strategy <- "bucket"
        } else {
            if(verbose) message("[auto] using global_sort")
            strategy <- "global_sort"
        }
    }
    
    # ------------ prerequisites ------------
    
    assert_cmd("sort"); assert_cmd("awk"); assert_cmd("xargs"); assert_cmd("cat")
    
    norm_dir <- file.path(tmpdir, "norm")
    bucket_dir <- file.path(tmpdir, "buckets")
    rbh_bucket_dir <- file.path(tmpdir, "rbh_buckets")
    
    on.exit({
        if(!keep_intermediate) {
            unlink(norm_dir, recursive = TRUE, force = TRUE)
            unlink(bucket_dir, recursive = TRUE, force = TRUE)
            unlink(rbh_bucket_dir, recursive = TRUE, force = TRUE)
        }
    }, add = TRUE)
    
    # ------------ execute ------------
    
    if(strategy == "global_sort") {
        norm_files <- normalize_many(blast_out_list, norm_dir)
        norm_sizes <- file.info(norm_files)$size
        norm_nonempty <- !is.na(norm_sizes) & norm_sizes > 0L
        if(!any(norm_nonempty)) {
            stop("All norm files are empty. Normalization wrote no rows. ",
                 "Check that blast_out_list files are BLAST outfmt 6 (tab-separated, 10 columns: ",
                 "qseqid sseqid pident qcovs qstart qend sstart send qlen slen).")
        }
        norm_files <- norm_files[norm_nonempty]
        
        sorted_norm_fn <- file.path(tmpdir, "all.norm.sorted.tsv")
        if(verbose) message(sprintf("[sort] global -> %s (%d norm files)", sorted_norm_fn, length(norm_files)))
        
        run_sort_cmd(in_files = norm_files,
                     out_file = sorted_norm_fn,
                     sort_mem = sort_mem,
                     sort_tmp = sort_tmp,
                     n_threads_sort = n_threads)
        
        sz <- file.info(sorted_norm_fn)$size
        if(is.na(sz) || sz == 0L) {
            stop("Sort produced empty output file: ", sorted_norm_fn,
                 ". Norm files had total ", sum(norm_sizes, na.rm = TRUE), " bytes. ",
                 "Try running the pipeline manually: LC_ALL=C xargs -0 cat < <list_of_norm_paths> | sort -t $'\\t' ...")
        }
        if(verbose) message(sprintf("[rbh] streaming -> %s", out_rbh_fn))
        rbh_from_sorted_norm_cpp(sorted_norm_fn, out_rbh_fn)
        
        if(!keep_intermediate) unlink(sorted_norm_fn)
        final_sort_rbh(out_rbh_fn)
        return(invisible(out_rbh_fn))
    }
    
    if(strategy == "bucket") {
        norm_files <- normalize_many(blast_out_list, norm_dir)
        norm_sizes <- file.info(norm_files)$size
        norm_nonempty <- !is.na(norm_sizes) & norm_sizes > 0L
        if(!any(norm_nonempty)) {
            stop("All norm files are empty. Normalization wrote no rows. ",
                 "Check that blast_out_list files are BLAST outfmt 6 (tab-separated, 10 columns: qseqid sseqid pident qcovs qstart qend sstart send qlen slen).")
        }
        norm_files <- norm_files[norm_nonempty]
        bucketize_norm_files(norm_files, bucket_dir, n_buckets, overwrite = TRUE)
        
        bucket_files_all <- file.path(bucket_dir, sprintf("bucket_%03d.tsv", 0:(n_buckets-1L)))
        sizes <- suppressWarnings(file.info(bucket_files_all)$size)
        ok <- which(!is.na(sizes) & sizes > 0)
        if(length(ok) == 0L) {
            file.create(out_rbh_fn)
            return(invisible(out_rbh_fn))
        }
        bucket_files <- bucket_files_all[ok]
        if(verbose) message(sprintf("[bucket] non-empty buckets: %d / %d", length(bucket_files), n_buckets))
        
        dir.create(rbh_bucket_dir, showWarnings = FALSE, recursive = TRUE)
        out_bucket_files <- file.path(rbh_bucket_dir, sub("\\.tsv$", ".rbh.tsv", basename(bucket_files)))
        
        mc_cores <- max(1L, min(as.integer(n_threads), length(bucket_files)))
        
        bucket_sort_mem <- sort_mem
        if(is.null(match.call()$sort_mem) && mc_cores > 1L) bucket_sort_mem <- "2G"
        if(verbose) message(sprintf("[bucket] mclapply workers=%d, bucket_sort_mem=%s", mc_cores, bucket_sort_mem))
        
        out_done <- mclapply(seq_along(bucket_files), function(i) {
            process_one_bucket(
                bucket_fn = bucket_files[i],
                out_bucket_fn = out_bucket_files[i],
                sort_tmp = sort_tmp,
                sort_mem = bucket_sort_mem,
                sort_parallel = 1L
            )
        }, mc.cores = mc_cores)
        
        out_done <- Filter(Negate(is.na), unlist(out_done, use.names = FALSE))
        if(length(out_done) == 0L) {
            file.create(out_rbh_fn)
            return(invisible(out_rbh_fn))
        }
        
        if(verbose) message(sprintf("[merge] %d bucket outputs -> %s", length(out_done), out_rbh_fn))
        merge_files_xargs_cat(out_done, out_rbh_fn)
        
        final_sort_rbh(out_rbh_fn)
        return(invisible(out_rbh_fn))
    }
    
    stop("Unknown strategy: ", strategy)
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
