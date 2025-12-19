##' @title Riparian-style plot from orthopair HDF5 data
##' @description
##' `plot_riparian` draws a riparian / braided-river style plot directly
##' from an HDF5 orthology file or an in-memory orthology object that can
##' be passed to `getOrthoPair()`.
##'
##' Internally the function calls `getOrthoPair(score = FALSE, loc = TRUE)`
##' to obtain a list of data.frames (`op`) and then builds the
##' plot from the synteny block information, similar to `devel_script.R`.
##'
##' Each element of the resulting `op` should contain at least
##' the columns defined in the `target_col` object from `devel_script.R`,
##' in particular:
##'   - `query_chr`, `query_synteny_block`
##'   - `subject_chr`, `subject_synteny_block`
##'
##' The function aggregates genes into synteny blocks on each side,
##' arranges blocks along the x–axis, and draws curved ribbons connecting
##' matched query/subject blocks between two genome rows.
##'
##' @param hdf5_fn character scalar, path to the HDF5 file used by
##'   `getOrthoPair()`. Exactly one of `hdf5_fn` or `object` must be non-NULL.
##' @param object optional in-memory object accepted by `getOrthoPair()`
##'   (see that function for details). Exactly one of `hdf5_fn` or `object`
##'   must be non-NULL.
##' @param min_block_genes integer, minimum number of unique genes that a
##'   synteny block must contain on *both* query and subject sides to be kept.
##' @param genome_names optional character vector of genome row labels.
##'   If `NULL`, the function tries to derive names from `names(op)`
##'   assuming patterns like `"GENOME1_GENOME2"`.
##' @param ribbon_alpha numeric (0–1), transparency of ribbons.
##' @param ribbon_width numeric, half–width of each block in x units.
##' @param y_gap numeric, vertical gap between the two genomes of a pair.
##' @param curve_npts integer, number of points used to generate the cosine
##'   curve in `calc_curvePolygon`.
##' @param curve_keepat integer, thinning parameter passed to
##'   `calc_curvePolygon` (controls how many points are retained).
##' @param chr_lwd numeric, line width for chromosome segments.
##' @param chr_palette function that takes an integer `n` and returns `n`
##'   colors for chromosomes.
##' @param ribbon_palette function that takes an integer `n` and returns `n`
##'   colors for block pairs.
##'
##' @return A `ggplot` object.
##'
##' @details
##' This function is meant as a light-weight alternative to the full
##' GENESPACE `plot_riparian` pipeline, for cases where orthology and
##' synteny are already summarised in an HDF5 file that can be read
##' with `getOrthoPair()`. It does **not** require GENESPACE block files
##' and instead works directly from the query/subject synteny block IDs.
##'
##' The function expects that each element of `op` corresponds
##' to a pair of genomes (query vs subject). For each pair, two horizontal
##' rows are drawn (query on top, subject below), and curved ribbons are
##' plotted between matching synteny blocks.
##'
##' @seealso `getOrthoPair`, `calc_curvePolygon`, `round_rect`
##' 
##' @import data.table
##' @import ggplot2
##' @import dplyr
##' @import tidyr
##' @import stringr
##' @import scales
##' @import ggforce
##' @import rhdf5
##' 
##' 
##' @export
plotRiparian <- function(hdf5_fn = NULL,
                          object = NULL,
                          min_block_genes = 5L,
                          genomes = NULL) {
    
    if (is.null(hdf5_fn) && is.null(object)) {
        stop("Provide either `hdf5_fn` (HDF5 path) or `object`, but not both.")
    }
    if (!is.null(hdf5_fn) && !is.null(object)) {
        stop("Provide only one of `hdf5_fn` (HDF5 path) or `object`, not both.")
    }
    
    if (!is.null(hdf5_fn)) {
        op <- getOrthoPair(hdf5_fn = hdf5_fn, score = FALSE, loc = TRUE)
        meta <- getMeta(hdf5_fn = hdf5_fn)
    } else {
        op <- getOrthoPair(object = object, score = FALSE, loc = TRUE)
        meta <- getMeta(object = object)
    }
    
    if (!is.list(op) || length(op) == 0L) {
        stop("`op` must be a non-empty list of data.frames.")
    }
    pair_blocks <- .get_block_list(op = op, 
                                   genomes = genomes,
                                   meta = meta,
                                   min_block_genes = min_block_genes)
    n_list <- length(pair_blocks)
    p <- .riparian_plot_engine(pair_blocks = pair_blocks[-n_list],
                               genomes = pair_blocks[[n_list]])
    return(p)
}

.get_block_list <- function(op, genomes, meta, min_block_genes){
    if(is.null(genomes)){
        genomes <- sort(meta$genomes)
        message("Since no genome order was specified, ",
                "plot genomes in the alphabetical order from the top.\n",
                paste(genomes, collapse = " "))
    }
    meta_pair_id <- strsplit(meta$pair_id, "_")
    target_pair_id <- NULL
    target_pair_order <- NULL
    ordered_pair_id <- NULL
    for(i in head(seq_along(genomes), -1)){
        i_pair_id <- c(genomes[i], genomes[i + 1])
        check <- lapply(meta_pair_id, setdiff, i_pair_id)
        check <- sapply(check, length) == 0
        if(sum(check) == 0){
            stop("No pairing infomation for ", i_pair_id, "!")
        }
        target_pair_id <- c(target_pair_id, meta$pair_id[check])
        i_pair_id <- paste(i_pair_id, collapse = "_")
        check <- i_pair_id == meta$pair_id[check]
        if(check){
            target_pair_order <- c(target_pair_order, 1)
        } else {
            target_pair_order <- c(target_pair_order, -1)
        }
        ordered_pair_id <- c(ordered_pair_id, i_pair_id)
    }
    
    eval_order <- match(target_pair_id, names(op))
    pair_blocks <- NULL
    for (i in seq_along(eval_order)) {
        index <- eval_order[i]
        pair_id <- ordered_pair_id[i]
        op_i <- as.data.table(op[[index]])
        needed_cols <- c("query_gene", "subject_gene",
                         "query_chr", "query_synteny_block",
                         "subject_chr", "subject_synteny_block",
                         "query_start","query_end","subject_start","subject_end")
        
        missing_cols <- setdiff(needed_cols, names(op_i))
        if (length(missing_cols) > 0L) {
            warning("Element ", index, " is missing required columns: ",
                    paste(missing_cols, collapse = ", "), "; skipping.")
            next
        }
        op_i <- subset(op_i, select = needed_cols)
        if(target_pair_order[i] == -1){
            op_i_names <- names(op_i)
            query <- grep("query_", op_i_names)
            subject <- grep("subject_", op_i_names)
            names(op_i)[query] <- sub("query_", "subject_", op_i_names[query])
            names(op_i)[subject] <- sub("subject_", "query_", op_i_names[subject])
        }
        
        # 1) sort（in-placeに近い）
        setorder(op_i, query_chr, query_start, subject_chr, subject_start)
        
        # 2) blockごとの遺伝子数（uniqueN が速い）
        q_sizes <- op_i[, .(q_n = uniqueN(query_gene)), by = query_synteny_block]
        s_sizes <- op_i[, .(s_n = uniqueN(subject_gene)), by = subject_synteny_block]
        
        # 3) joinしてフィルタ（tapply→join）
        op_i <- q_sizes[op_i, on = .(query_synteny_block)]
        op_i <- s_sizes[op_i, on = .(subject_synteny_block)]
        op_i <- op_i[q_n >= min_block_genes & s_n >= min_block_genes]
        
        # 4) q_blocks / s_blocks（tapply複数回→1回の集計）
        q_blocks <- op_i[, .(
            chr   = first(query_chr),
            start = min(query_start),
            end   = max(query_end)
        ), by = .(block = query_synteny_block)]
        q_blocks[, width := end - start]
        setorder(q_blocks, chr, start)
        q_blocks[, x_idx := .I]
        
        s_blocks <- op_i[, .(
            chr   = first(subject_chr),
            start = min(subject_start),
            end   = max(subject_end)
        ), by = .(block = subject_synteny_block)]
        s_blocks[, width := end - start]
        setorder(s_blocks, chr, start)
        s_blocks[, x_idx := .I]
        
        # 5) q,s の付与（paste+named vector→2回のjoin）
        #    blockとchrの両方でマップ（元コードの "::" と同等）
        op_i <- q_blocks[op_i,
                         on = .(block = query_synteny_block, chr = query_chr),
                         nomatch = 0L]
        setnames(op_i, "x_idx", "q")
        setnames(op_i, "block", "query_block")
        setnames(op_i, "chr", "query_block_chr")
        setnames(op_i, "start", "query_block_start")
        setnames(op_i, "end", "query_block_end")
        setnames(op_i, "width", "query_block_width")
        
        op_i <- s_blocks[op_i,
                         on = .(block = subject_synteny_block, chr = subject_chr),
                         nomatch = 0L]
        setnames(op_i, "x_idx", "s")
        setnames(op_i, "block", "subject_block")
        setnames(op_i, "chr", "subject_block_chr")
        setnames(op_i, "start", "subject_block_start")
        setnames(op_i, "end", "subject_block_end")
        setnames(op_i, "width", "subject_block_width")
        
        dt <- as.data.table(op_i)[, .(q, s)]
        # op_i の行順を保ったまま run 分割
        dt[, `:=`(dq = q - shift(q), ds = s - shift(s))]
        dt[, adjacent := !is.na(dq) & (abs(dq) <= 1L) & (abs(ds) <= 1L)]
        dt[, block_id := cumsum(!adjacent)]
        
        # block_id を op_i に戻す
        op_i$block_id <- dt$block_id
        
        # 9) block_df（tapply→unique）
        block_df <- op_i[, .(query_block_start = min(query_start), 
                             query_block_end = max(query_end), 
                             query_block_direction = tail(query_start, 1) - head(query_start, 1),
                             subject_block_start = min(subject_start), 
                             subject_block_end = max(subject_end),
                             subject_block_direction = tail(subject_start, 1) - head(subject_start, 1)),
                         by = .(block_id, query_block_chr, subject_block_chr)]
        query_reverse <- block_df$query_block_direction <= 0
        block_df$query_block_direction <- rep("+", length(block_df$query_block_direction))
        block_df$query_block_direction[query_reverse] <- "-"
        subject_reverse <- block_df$subject_block_direction <= 0
        block_df$subject_block_direction <- rep("+", length(block_df$subject_block_direction))
        block_df$subject_block_direction[subject_reverse] <- "-"
        new_blocks <- list(block_df)
        names(new_blocks) <- pair_id
        pair_blocks <- c(pair_blocks, new_blocks)
    }
    pair_blocks <- c(pair_blocks, list(genomes = genomes))
    return(pair_blocks)
}

.riparian_plot_engine <- function(pair_blocks,
                                  genomes,
                                  chr_sizes = NULL,
                                  chr_gap = 2e6,
                                  track_gap = 1.3,
                                  max_links = 50000,
                                  seed = 1,
                                  alpha_base = NULL,
                                  linewidth = 0.45,
                                  inv_darken = 0.35,
                                  inv_alpha_mult = 1.6,
                                  palette = c("#d73027","#fc8d59","#fee08b","#91bfdb","#4575b4"),
                                  min_block_width = 0,
                                  chr_bar_lw = 2.8,
                                  chr_sep_lw = 0.35,
                                  chr_label_size = 3.2) {
    
    stopifnot(length(genomes) >= 2)
    stopifnot(is.list(pair_blocks))
    
    expected_pairs <- paste0(genomes[-length(genomes)], "_", genomes[-1])
    missing <- setdiff(expected_pairs, names(pair_blocks))
    if (length(missing) > 0) stop("pair_blocks missing: ", paste(missing, collapse=", "))
    
    # --- chr_sizes 推定 or 整形 ---
    if (is.null(chr_sizes)) {
        chr_sizes <- setNames(vector("list", length(genomes)), genomes)
        for (g in genomes) chr_sizes[[g]] <- data.table(chr=character(), len=numeric())
        
        for (pnm in expected_pairs) {
            parts <- str_split(pnm, "_", simplify = TRUE)
            g1 <- parts[1]; g2 <- parts[2]
            dt <- as.data.table(pair_blocks[[pnm]])
            
            dt[, `:=`(
                query_block_chr   = norm_chr(query_block_chr),
                subject_block_chr = norm_chr(subject_block_chr),
                query_block_end   = as.numeric(query_block_end),
                subject_block_end = as.numeric(subject_block_end)
            )]
            
            cs1 <- dt[, .(len = max(query_block_end, na.rm=TRUE)), by=.(chr=query_block_chr)]
            cs2 <- dt[, .(len = max(subject_block_end, na.rm=TRUE)), by=.(chr=subject_block_chr)]
            chr_sizes[[g1]] <- rbindlist(list(chr_sizes[[g1]], cs1), use.names=TRUE, fill=TRUE)
            chr_sizes[[g2]] <- rbindlist(list(chr_sizes[[g2]], cs2), use.names=TRUE, fill=TRUE)
        }
        
        chr_sizes <- lapply(chr_sizes, function(cs) {
            cs <- as.data.table(cs)
            cs[, chr := norm_chr(chr)]
            cs[, len := as.numeric(len)]
            cs <- cs[is.finite(len)]
            cs <- cs[, .(len=max(len)), by=.(chr)]
            # 並び順固定：chr番号→文字
            key <- chr_order_key(cs$chr)
            cs <- merge(cs, key, by="chr", all.x=TRUE)
            cs <- cs[order(num, chr2)][, .(chr, len)]
            cs
        })
    } else {
        chr_sizes <- lapply(chr_sizes, function(cs) {
            cs <- as.data.table(cs)
            setnames(cs, c("chr","len"))
            cs[, chr := norm_chr(chr)]
            cs[, len := as.numeric(len)]
            cs <- cs[is.finite(len)]
            cs <- cs[, .(len=max(len)), by=.(chr)]
            key <- chr_order_key(cs$chr)
            cs <- merge(cs, key, by="chr", all.x=TRUE)
            cs <- cs[order(num, chr2)][, .(chr, len)]
            cs
        })
    }
    
    # offsets per genome
    build_offsets <- function(cs) {
        cs <- as.data.table(copy(cs))
        cs[, offset := c(0, cumsum(head(len, -1) + chr_gap))]
        cs
    }
    chr_offsets <- lapply(chr_sizes[genomes], build_offsets)
    
    map_x <- function(genome, chr, pos) {
        off <- chr_offsets[[genome]]
        chr <- norm_chr(chr)
        idx <- match(chr, off$chr)
        off$offset[idx] + as.numeric(pos)
    }
    
    # y positions
    y_map <- setNames(rev(seq_along(genomes)) * track_gap, genomes)
    
    # chromosome track segments (chr単位)
    tracks <- rbindlist(lapply(genomes, function(g) {
        off <- chr_offsets[[g]]
        off[, .(genome=g, chr=chr,
                x0=offset, x1=offset+len,
                xmid=offset+len/2,
                y=y_map[g])]
    }))
    chr_seps <- tracks[, .(x = x0, y = y, genome=genome, chr=chr)]
    
    # chr label: "chr01" -> "1"
    chr_num <- function(chr) {
        z <- suppressWarnings(as.integer(str_extract(chr, "(?<=^chr)\\d+")))
        ifelse(is.na(z), chr, as.character(z))
    }
    chr_labels <- tracks[, .(genome, chr, x=xmid, y=y, lab=chr_num(chr))]
    
    # --- links: build endpoints ---
    link_dt_list <- vector("list", length(expected_pairs))
    names(link_dt_list) <- expected_pairs
    
    for (pnm in expected_pairs) {
        parts <- str_split(pnm, "_", simplify = TRUE)
        g1 <- parts[1]; g2 <- parts[2]
        dt <- as.data.table(copy(pair_blocks[[pnm]]))
        
        needed <- c("block_id",
                    "query_block_chr","query_block_start","query_block_end","query_block_direction",
                    "subject_block_chr","subject_block_start","subject_block_end","subject_block_direction")
        miss <- setdiff(needed, names(dt))
        if (length(miss) > 0) stop(pnm, " missing columns: ", paste(miss, collapse=", "))
        
        dt[, `:=`(
            query_block_chr     = norm_chr(query_block_chr),
            subject_block_chr   = norm_chr(subject_block_chr),
            query_block_start   = as.numeric(query_block_start),
            query_block_end     = as.numeric(query_block_end),
            subject_block_start = as.numeric(subject_block_start),
            subject_block_end   = as.numeric(subject_block_end),
            qmid = (as.numeric(query_block_start) + as.numeric(query_block_end))/2,
            smid = (as.numeric(subject_block_start) + as.numeric(subject_block_end))/2,
            inv  = (query_block_direction != subject_block_direction)
        )]
        
        dt[, bw := (abs(query_block_end - query_block_start) + abs(subject_block_end - subject_block_start))/2]
        if (min_block_width > 0) dt <- dt[bw >= min_block_width]
        
        dt[, x1 := map_x(g1, query_block_chr, qmid)]
        dt[, x2 := map_x(g2, subject_block_chr, smid)]
        dt <- dt[is.finite(x1) & is.finite(x2)]
        
        dt[, `:=`(
            pair = pnm,
            group = paste0(pnm, "::", block_id),
            y1 = y_map[g1],
            y2 = y_map[g2],
            ymid = (y_map[g1] + y_map[g2]) / 2,
            colval = x1
        )]
        
        link_dt_list[[pnm]] <- dt[, .(pair, block_id, group, inv, colval, y1, y2, ymid, x1, x2, bw)]
    }
    
    links <- rbindlist(link_dt_list, use.names=TRUE, fill=TRUE)
    if (nrow(links) == 0) stop("No links after filtering. Check chr mapping / min_block_width.")
    
    # --- sampling & auto alpha ---
    n_links <- nrow(links)
    if (n_links > max_links) {
        set.seed(seed)
        links <- links[sample.int(n_links, max_links)]
        n_links <- nrow(links)
    }
    
    if (is.null(alpha_base)) {
        alpha_base <- pmin(0.45, pmax(0.06, 0.55 / (1 + log10(n_links))))
    }
    alpha_inv <- pmin(1, alpha_base * 1.6)
    
    # --- bezier control points ---
    bez <- rbindlist(list(
        links[, .(group, pt=1L, x=x1, y=y1, inv, colval)],
        links[, .(group, pt=2L, x=x1, y=ymid, inv, colval)],
        links[, .(group, pt=3L, x=x2, y=ymid, inv, colval)],
        links[, .(group, pt=4L, x=x2, y=y2, inv, colval)]
    ))
    setorder(bez, group, pt)
    
    # colval -> hex via gradient, then inv -> darken
    rng <- range(bez$colval, na.rm=TRUE)
    pal_fun <- scales::gradient_n_pal(palette)
    bez[, col_t := (colval - rng[1]) / (rng[2] - rng[1] + 1e-12)]
    bez[, col_hex := pal_fun(pmin(1, pmax(0, col_t)))]
    bez[inv == TRUE, col_hex := colorspace::darken(col_hex, amount = inv_darken)]
    
    p <- ggplot() +
        ggforce::geom_bezier(
            data = bez,
            aes(x=x, y=y, group=group, color=col_hex, alpha=inv),
            linewidth = linewidth,
            lineend = "round",
            show.legend = FALSE
        ) +
        scale_color_identity() +
        scale_alpha_manual(values = c(`FALSE` = alpha_base, `TRUE` = alpha_inv), guide = "none") +
        # chr bars
        geom_segment(
            data = tracks,
            aes(x=x0, xend=x1, y=y, yend=y),
            linewidth = chr_bar_lw,
            lineend = "round",
            color = "grey15"
        ) +
        # chr separators
        geom_segment(
            data = chr_seps,
            aes(x=x, xend=x, y=y-0.18, yend=y+0.18),
            linewidth = chr_sep_lw,
            color = "grey5"
        ) +
        # chr circle label
        geom_point(
            data = chr_labels,
            aes(x=x, y=y),
            shape = 21, fill = "white", color = "grey10", stroke = 0.3,
            size = 3.4
        ) +
        geom_text(
            data = chr_labels,
            aes(x=x, y=y, label=lab),
            size = chr_label_size/3.2,
            color = "grey10"
        ) +
        # genome labels
        geom_text(
            data = data.table(genome=genomes, y=y_map[genomes]),
            aes(x = min(tracks$x0, na.rm=TRUE) - chr_gap*0.6, y=y, label=genome),
            hjust = 1, size = 4
        ) +
        theme_void() +
        theme(legend.position = "none")
    
    attr(p, "alpha_base") <- alpha_base
    attr(p, "n_links_plotted") <- n_links
    p
}

norm_chr <- function(x, pad = 2, prefix = "chr") {
    x <- as.character(x)
    x <- trimws(x)
    
    # まず数字を抽出（最初に出てくる連続数字）
    num <- suppressWarnings(as.integer(stringr::str_extract(x, "\\d+")))
    
    out <- ifelse(
        !is.na(num),
        paste0(prefix, stringr::str_pad(num, width = pad, pad = "0")),
        # 数字がないものは "chr_" + 元文字列を簡易クリーニング
        paste0(prefix, "_", stringr::str_replace_all(tolower(x), "[^a-z0-9]+", "_"))
    )
    out
}

# chrを 1..n（数字）優先で並べるためのキー
chr_order_key <- function(chr, prefix = "chr") {
    chr <- as.character(chr)
    # "chr01" -> 1, "chr_foo" -> NA
    num <- suppressWarnings(as.integer(stringr::str_extract(chr, "(?<=^chr)\\d+")))
    # NAは大きい値にして後ろへ、同値は文字で安定化
    data.table(chr = chr,
               num = ifelse(is.na(num), 1e9, num),
               chr2 = chr)
}