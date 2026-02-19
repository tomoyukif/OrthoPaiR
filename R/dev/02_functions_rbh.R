#' Create BLAST databases and perform reciprocal BLAST
#'
#' This function prepare BLAST databases and perform BLAST searches reciprocally
#' on given input datasets to create the list of reciprocal BLAST hits.
#'
#' @export
#'
rbh <- function(object, blast_path, target_pair, n_threads, overwrite){
    blast_dir <- file.path(object$working_dir, "blast")
    dir.create(path = blast_dir, showWarnings = FALSE, recursive = TRUE)
    
    .makeBlastDB(object = object, 
                 blast_dir = blast_dir,
                 blast_path = blast_path,
                 n_threads = n_threads,
                 overwrite = overwrite)
    
    check_blast_out <- .checkBLASTout(object, blast_dir, overwrite)
    job_assign <- .jobAssign(check_blast_out = check_blast_out,
                             n_threads = n_threads,
                             min_threads = 8)
    
    mclapply(X = unique(job_assign$fasta_chunk$chunk),
             mc.cores = job_assign$n_parallel_jobs, 
             FUN = .blastn_search, 
             object = object,
             blast_dir = blast_dir,
             blast_path = blast_path,
             fasta_chunk = job_assign$fasta_chunk,
             threads_per_job = job_assign$threads_per_job,
             index_offset = check_blast_out$blast_out_index_offset)
    
    rbh_dir <- file.path(object$working_dir, "rbh")
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
}

.splitRBHbyGenomePair <- function(rbh_fn, rbh_dir, genome_width = 4L){
    if (!file.exists(rbh_fn)) {
        warning("RBH file not found: ", rbh_fn, ". Skipping split by genome pairs.")
        return(invisible(NULL))
    }
    
    # Create rbh output directory
    dir.create(path = rbh_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Use awk for fast single-pass splitting (much faster than R line-by-line processing)
    # Columns: query_tx, subject_tx, pident, q2s_qcovs, s2q_qcovs, q2s_ci, s2q_ci, mutual_ci
    sys <- Sys.info()[["sysname"]]
    if (!sys %in% c("Linux", "Darwin")) {
        stop("Linux/macOS only (requires awk)")
    }
    
    if (Sys.which("awk") == "") {
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
  if (gA <= gB) {
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
    if (!is.null(attr(rc, "status")) && attr(rc, "status") != 0L) {
        stop("awk split failed:\n", paste(rc, collapse = "\n"))
    }
    
    invisible(TRUE)
}

.makeBlastDB <- function(object, blast_dir, blast_path, n_threads, overwrite){
    check_db <- .checkDBexist(object = object, 
                              blast_dir = blast_dir,
                              overwrite = overwrite)
    
    .makeMergedFASTA(object = object, 
                     check_db = check_db, 
                     blast_dir = blast_dir)
    
    cds_fn_list <- list.files(blast_dir, "all_cds.fa", full.names = TRUE)
    
    mclapply(X = cds_fn_list,
             mc.cores = min(n_threads, length(cds_fn_list)),
             FUN = .blastDBengine,
             blast_dir = blast_dir,
             blast_path = blast_path,
             overwrite = overwrite)
    invisible(TRUE)
}

.checkDBexist <- function(object, blast_dir, overwrite){
    to_be_db <- object$input$name
    blast_db_list <- list.files(blast_dir,
                                "_blastdb.list",
                                full.names = TRUE,
                                recursive = TRUE)
    db_index_offset <- 0
    if(length(blast_dir) > 0 & !overwrite){
        genome_in_db <- sapply(blast_db_list, function(x){
            x_out <- fread(file = x, 
                           sep = "\t",
                           header = FALSE,
                           stringsAsFactors = FALSE)
            return(as.vector(x_out))
        })
        to_be_db <- to_be_db[!to_be_db %in% genome_in_db]
        db_index_offset <- max(as.numeric(sub("_.+",
                                              "",
                                              basename(blast_db_list))))
    }
    out <- list(to_be_db = to_be_db, db_index_offset = db_index_offset)
    return(out)
}

.makeMergedFASTA <- function(object, check_db, blast_dir){
    index <- check_db$db_index_offset + 1
    fn <- file.path(blast_dir, paste(index, "all_cds.fa", sep = "_"))
    cds_fn_list <- object$input$cds[object$input$name %in% check_db$to_be_db]
    cds_fn_list <- normalizePath(cds_fn_list, mustWork = TRUE)
    
    sys <- Sys.info()[["sysname"]]
    if (sys %in% c("Linux", "Darwin")) {
        # Fast merge via cat: no R parsing, minimal memory, very fast
        status <- system2("cat", cds_fn_list, stdout = fn)
        if (status != 0L) {
            stop("Merge of CDS FASTA files failed (cat exited with status ", status, ")")
        }
    } else {
        # Fallback: read and write in R (e.g. Windows)
        fa <- readDNAStringSet(cds_fn_list)
        writeXStringSet(fa, fn)
    }
    invisible(TRUE)
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
            write(out, file.path(blast_dir, paste0(index, "_makeblastdb.error")))
            
        } else {
            write("", job_finish_marker)
        }
    }
    invisible(TRUE)
}

.checkBLASTout <- function(object, blast_dir, overwrite){
    blast_db_list <- list.files(blast_dir,
                                ".blastdb.nsq",
                                full.names = TRUE,
                                recursive = TRUE)
    to_be_blast <- expand.grid(genome = object$input$name, db = blast_db_list)
    
    blast_out_list <- list.files(blast_dir,
                                 "blast.out.list",
                                 full.names = TRUE,
                                 recursive = TRUE)
    blast_out_index_offset <- 0
    if(length(blast_out_list) > 0 & !overwrite){
        blast_out_list <- lapply(blast_out_list, function(x){
            x_out <- fread(file = x, 
                           sep = "\t",
                           header = FALSE,
                           stringsAsFactors = FALSE)
            return(x_out)
        })
        blast_out_list <- do.call(rbind, blast_out_list)
        to_be_blast_id <- apply(to_be_blast, 1, paste, collapse = "_")
        blast_out_list_id <- apply(blast_out_list, 1, paste, collapse = "_")
        to_be_blast <- to_be_blast[!to_be_blast_id %in% blast_out_list_id, ]
        blast_out_index_offset <- max(as.numeric(sub("_.+",
                                                     "",
                                                     basename(blast_out_list))))
    }
    out <- list(to_be_blast = to_be_blast, 
                blast_out_index_offset = blast_out_index_offset)
    return(out)
}

.jobAssign <- function(check_blast_out, n_threads, min_threads = 8){
    n_fasta <- length(check_blast_out$to_be_blast$genome)
    ma <- .mem_available_bytes()
    ma_per_single_fasta <- n_fasta * 2e8
    chunk_per_fasta <- ma / ma_per_single_fasta
    fasta_chunk <- tapply(seq_along(check_blast_out$to_be_blast$genome), 
                          check_blast_out$to_be_blast$db,
                          function(i){
                              x <- check_blast_out$to_be_blast$genome[i]
                              db_name <- check_blast_out$to_be_blast$db[i][1]
                              n_breaks <- ceiling(length(x) / chunk_per_fasta)
                              if(n_breaks > 1){
                                  chunk <- cut(seq_along(x), 
                                               breaks = n_breaks)
                                  chunk <- as.integer(chunk)
                              } else {
                                  chunk <- rep(1L, length(x))
                              }
                              out <- data.frame(chunk = chunk, 
                                                genome = x,
                                                db = db_name)
                              return(out)
                          })
    fasta_chunk <- do.call(rbind, fasta_chunk)
    fasta_chunk_id <- paste(fasta_chunk$chunk, fasta_chunk$db, sep = "_")
    fasta_chunk_id <- factor(fasta_chunk_id)
    fasta_chunk$chunk <- as.integer(fasta_chunk_id)
    n_chunk <- max(fasta_chunk$chunk)
    
    if(n_threads < min_threads){
        threads_per_job <- n_threads
        n_parallel_jobs <- 1
        
    } else {
        threads_per_job <- floor(n_threads / n_chunk)
        if(threads_per_job < min_threads){
            threads_per_job <- min_threads
        }
        n_parallel_jobs <- floor(n_threads / threads_per_job)
    }
    out <- list(fasta_chunk = fasta_chunk,
                n_parallel_jobs = n_parallel_jobs, 
                threads_per_job = threads_per_job)
    return(out)
}

.blastn_search <- function(index,
                           object, 
                           blast_dir,
                           blast_path,
                           fasta_chunk, 
                           threads_per_job,
                           index_offset){
    db_fn <- fasta_chunk$db[fasta_chunk$chunk == index][1]
    db_prefix <- sub(".nsq", "", db_fn)
    fasta_chunk_i <- fasta_chunk$genome[fasta_chunk$chunk == index]
    index <- index_offset + index
    query_cds_list <- object$input$cds[object$input$name %in% fasta_chunk_i]
    query_cds <- sapply(query_cds_list, function(x){
        x_out <- fread(file = x, 
                       sep = "\t",
                       header = FALSE,
                       stringsAsFactors = FALSE)
        return(x_out)
    })
    query_cds <- unlist(query_cds)
    query_fn <- file.path(blast_dir, paste0(index, "_query_cds.fa"))
    write(x = query_cds, file = query_fn, sep = "\t")
    
    blast_out <- .run_blastn(query_fn = query_fn, 
                             db_prefix = db_prefix,
                             blast_dir = blast_dir,
                             blast_path = blast_path,
                             threads_per_job = threads_per_job,
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
                        index){
    out_fn <- file.path(blast_dir, paste0(index, "_blast.out"))
    
    blast_args <- paste(paste("-query", query_fn),
                        paste("-db", db_prefix),
                        "-task blastn -max_target_seqs 100",
                        "-evalue 1e-4 -strand plus",
                        "-outfmt '6 qseqid sseqid pident qcovs qstart qend sstart send'",
                        paste("-num_threads", threads_per_job),
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
                        n_threads = parallel::detectCores(),
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
    if (!sys %in% c("Linux", "Darwin")) {
        stop("Linux/macOS only (mclapply + sort/awk/xargs).")
    }
    
    dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
    dir.create(sort_tmp, showWarnings = FALSE, recursive = TRUE)
    tmpdir <- normalizePath(tmpdir, mustWork = TRUE)
    sort_tmp <- normalizePath(sort_tmp, mustWork = TRUE)
    
    # ------------ helpers ------------
    
    assert_cmd <- function(cmd) {
        if (Sys.which(cmd) == "") stop("Required command not found in PATH: ", cmd)
    }
    
    get_free_bytes <- function(path) {
        out <- suppressWarnings(system2("df", c("-B1", path), stdout = TRUE, stderr = TRUE))
        if (length(out) < 2) return(NA_real_)
        cols <- strsplit(out[2], "[[:space:]]+")[[1]]
        cols <- cols[nzchar(cols)]
        avail <- suppressWarnings(as.numeric(cols[4]))
        if (is.na(avail)) NA_real_ else avail
    }
    
    guess_sort_mem <- function(default = "4G") {
        if (identical(sys, "Linux") && file.exists("/proc/meminfo")) {
            mi <- readLines("/proc/meminfo", warn = FALSE)
            ma <- grep("^MemAvailable:", mi, value = TRUE)
            if (length(ma)) {
                kb <- suppressWarnings(as.numeric(sub("MemAvailable:\\s+([0-9]+)\\s+kB.*", "\\1", ma)))
                if (!is.na(kb)) {
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
        for (i in seq_along(blast_out_list)) {
            if (verbose) message(sprintf("[norm %d/%d] %s", i, length(blast_out_list), basename(blast_out_list[i])))
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
        for (f in files) {
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
        
        sort_part <- if (identical(sys, "Linux")) {
            if (is.null(sort_mem)) sort_mem <- guess_sort_mem()
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
        if (!is.null(attr(rc, "status")) && attr(rc, "status") != 0L) {
            stop("sort pipeline failed. ", paste(rc, collapse = "\n"))
        }
        
        out_file
    }
    
    final_sort_rbh <- function(rbh_fn) {
        if (!final_sort_by_mutual_ci) return(invisible(rbh_fn))
        assert_cmd("sort"); assert_cmd("bash")
        orig_size <- file.info(rbh_fn)$size
        if (is.na(orig_size) || orig_size == 0L) return(invisible(rbh_fn))
        tmp <- paste0(rbh_fn, ".tmp_sorted")
        script <- sprintf("LC_ALL=C sort -t $'\\t' -k8,8gr %s > %s",
                          shQuote(rbh_fn), shQuote(tmp))
        rc <- system2("bash", c("-c", shQuote(script)), stdout = TRUE, stderr = TRUE)
        if (!is.null(attr(rc, "status")) && attr(rc, "status") != 0L) {
            stop("final sort failed: ", paste(rc, collapse = "\n"))
        }
        tmp_size <- file.info(tmp)$size
        if (is.na(tmp_size) || tmp_size == 0L) {
            unlink(tmp, force = TRUE)
            stop("final sort produced empty output; ", rbh_fn, " left unchanged. ",
                 "Check that the file has 12 tab-separated columns (column 12 = mutual_ci).")
        }
        ok <- file.rename(tmp, rbh_fn)
        if (!ok) {
            file.copy(tmp, rbh_fn, overwrite = TRUE)
            unlink(tmp, force = TRUE)
        }
        invisible(rbh_fn)
    }
    
    bucketize_norm_files <- function(norm_files, bucket_dir, n_buckets, overwrite = TRUE) {
        assert_cmd("awk")
        dir.create(bucket_dir, showWarnings = FALSE, recursive = TRUE)
        if (overwrite) {
            old <- list.files(bucket_dir, full.names = TRUE, pattern = "\\.tsv$")
            if (length(old)) unlink(old)
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
        
        if (verbose) message(sprintf("[bucketize] distributing into %d buckets", n_buckets))
        rc <- system2("awk", c("-f", awk_fn, norm_files), stdout = TRUE, stderr = TRUE)
        if (!is.null(attr(rc, "status")) && attr(rc, "status") != 0) {
            stop("awk bucketize failed:\n", paste(rc, collapse = "\n"))
        }
        
        unlink(awk_fn)
        invisible(TRUE)
    }
    
    process_one_bucket <- function(bucket_fn, out_bucket_fn, sort_tmp, sort_mem, sort_parallel = 1L) {
        fi <- file.info(bucket_fn)
        if (is.na(fi$size) || fi$size == 0) return(NA_character_)
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
        if (!identical(rc, 0L)) stop("merge(cat) failed (exit code ", rc, "):\n", cmd)
        out_file
    }
    
    # ------------ decide auto strategy ------------
    
    total_bytes <- sum(file.info(blast_out_list)$size, na.rm = TRUE)
    if (is.null(sort_mem)) sort_mem <- guess_sort_mem()
    
    if (strategy == "auto") {
        free_bytes <- get_free_bytes(sort_tmp)
        need <- total_bytes * 2
        if (!is.na(free_bytes) && free_bytes < need) {
            if (verbose) message(sprintf("[auto] tmp free %.1f GB < need %.1f GB -> bucket",
                                         free_bytes/1024^3, need/1024^3))
            strategy <- "bucket"
        } else {
            if (verbose) message("[auto] using global_sort")
            strategy <- "global_sort"
        }
    }
    
    # ------------ prerequisites ------------
    
    assert_cmd("sort"); assert_cmd("awk"); assert_cmd("xargs"); assert_cmd("cat")
    
    norm_dir <- file.path(tmpdir, "norm")
    bucket_dir <- file.path(tmpdir, "buckets")
    rbh_bucket_dir <- file.path(tmpdir, "rbh_buckets")
    
    on.exit({
        if (!keep_intermediate) {
            unlink(norm_dir, recursive = TRUE, force = TRUE)
            unlink(bucket_dir, recursive = TRUE, force = TRUE)
            unlink(rbh_bucket_dir, recursive = TRUE, force = TRUE)
        }
    }, add = TRUE)
    
    # ------------ execute ------------
    
    if (strategy == "global_sort") {
        norm_files <- normalize_many(blast_out_list, norm_dir)
        norm_sizes <- file.info(norm_files)$size
        norm_nonempty <- !is.na(norm_sizes) & norm_sizes > 0L
        if (!any(norm_nonempty)) {
            stop("All norm files are empty. Normalization wrote no rows. ",
                 "Check that blast_out_list files are BLAST outfmt 6 (tab-separated, 8 columns: ",
                 "query_tx, subject_tx, pident, qcovs, qstart, qend, sstart, send).")
        }
        norm_files <- norm_files[norm_nonempty]
        
        sorted_norm_fn <- file.path(tmpdir, "all.norm.sorted.tsv")
        if (verbose) message(sprintf("[sort] global -> %s (%d norm files)", sorted_norm_fn, length(norm_files)))
        
        run_sort_cmd(in_files = norm_files,
                     out_file = sorted_norm_fn,
                     sort_mem = sort_mem,
                     sort_tmp = sort_tmp,
                     n_threads_sort = n_threads)
        
        sz <- file.info(sorted_norm_fn)$size
        if (is.na(sz) || sz == 0L) {
            stop("Sort produced empty output file: ", sorted_norm_fn,
                 ". Norm files had total ", sum(norm_sizes, na.rm = TRUE), " bytes. ",
                 "Try running the pipeline manually: LC_ALL=C xargs -0 cat < <list_of_norm_paths> | sort -t $'\\t' ...")
        }
        if (verbose) message(sprintf("[rbh] streaming -> %s", out_rbh_fn))
        rbh_from_sorted_norm_cpp(sorted_norm_fn, out_rbh_fn)
        
        if (!keep_intermediate) unlink(sorted_norm_fn)
        final_sort_rbh(out_rbh_fn)
        return(invisible(out_rbh_fn))
    }
    
    if (strategy == "bucket") {
        norm_files <- normalize_many(blast_out_list, norm_dir)
        norm_sizes <- file.info(norm_files)$size
        norm_nonempty <- !is.na(norm_sizes) & norm_sizes > 0L
        if (!any(norm_nonempty)) {
            stop("All norm files are empty. Normalization wrote no rows. ",
                 "Check that blast_out_list files are BLAST outfmt 6 (tab-separated, 8 columns).")
        }
        norm_files <- norm_files[norm_nonempty]
        bucketize_norm_files(norm_files, bucket_dir, n_buckets, overwrite = TRUE)
        
        bucket_files_all <- file.path(bucket_dir, sprintf("bucket_%03d.tsv", 0:(n_buckets-1L)))
        sizes <- suppressWarnings(file.info(bucket_files_all)$size)
        ok <- which(!is.na(sizes) & sizes > 0)
        if (length(ok) == 0L) {
            file.create(out_rbh_fn)
            return(invisible(out_rbh_fn))
        }
        bucket_files <- bucket_files_all[ok]
        if (verbose) message(sprintf("[bucket] non-empty buckets: %d / %d", length(bucket_files), n_buckets))
        
        dir.create(rbh_bucket_dir, showWarnings = FALSE, recursive = TRUE)
        out_bucket_files <- file.path(rbh_bucket_dir, sub("\\.tsv$", ".rbh.tsv", basename(bucket_files)))
        
        mc_cores <- max(1L, min(as.integer(n_threads), length(bucket_files)))
        
        bucket_sort_mem <- sort_mem
        if (is.null(match.call()$sort_mem) && mc_cores > 1L) bucket_sort_mem <- "2G"
        if (verbose) message(sprintf("[bucket] mclapply workers=%d, bucket_sort_mem=%s", mc_cores, bucket_sort_mem))
        
        out_done <- parallel::mclapply(seq_along(bucket_files), function(i) {
            process_one_bucket(
                bucket_fn = bucket_files[i],
                out_bucket_fn = out_bucket_files[i],
                sort_tmp = sort_tmp,
                sort_mem = bucket_sort_mem,
                sort_parallel = 1L
            )
        }, mc.cores = mc_cores)
        
        out_done <- Filter(Negate(is.na), unlist(out_done, use.names = FALSE))
        if (length(out_done) == 0L) {
            file.create(out_rbh_fn)
            return(invisible(out_rbh_fn))
        }
        
        if (verbose) message(sprintf("[merge] %d bucket outputs -> %s", length(out_done), out_rbh_fn))
        merge_files_xargs_cat(out_done, out_rbh_fn)
        
        final_sort_rbh(out_rbh_fn)
        return(invisible(out_rbh_fn))
    }
    
    stop("Unknown strategy: ", strategy)
}


