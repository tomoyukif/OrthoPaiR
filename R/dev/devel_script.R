library(data.table)
library(rhdf5)
library(ggplot2)
min_block_genes = 5L
genome_names = NULL
ribbon_alpha = 0.4
ribbon_width = 0.45
y_gap = 1
curve_npts = 250L
curve_keepat = round(curve_npts / 20)
chr_lwd = 2
chr_palette = function(n) grDevices::hcl.colors(n, "Dark 3")
ribbon_palette = function(n) grDevices::hcl.colors(n, "Zissou 1")

hdf5_fn <- "~/workspace/orthology/output/benchmark/orthopair/hdf5_out/Osat_Ogla.h5"
devtools::load_all("./")
op <- getOrthoPair(hdf5_fn = hdf5_fn, score = FALSE, loc = TRUE)

# Prepare data for plotting
n_pairs <- length(op)
all_chrom <- list()
all_ribbons <- list()
y_levels <- c()
row_labels <- c()

for (i in seq_along(op)) {
    op_i <- as.data.table(op[[i]])
    needed_cols <- c("query_gene", "subject_gene",
                     "query_chr", "query_synteny_block",
                     "subject_chr", "subject_synteny_block",
                     "query_start","query_end","subject_start","subject_end")
    
    missing_cols <- setdiff(needed_cols, names(op_i))
    if (length(missing_cols) > 0L) {
        warning("Element ", i, " is missing required columns: ",
                paste(missing_cols, collapse = ", "), "; skipping.")
        next
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
    block_df <- op_i[, .(query_block_chr = unique(query_block_chr), 
                         query_block_start = min(query_start), 
                         query_block_end = max(query_end), 
                         subject_block_chr = unique(subject_block_chr), 
                         subject_block_start = min(subject_start), 
                         subject_block_end = max(subject_end)),
                         by = block_id]
    
    # Build chromosome segments for this pair
    # query side
    q_chr_op <- stats::aggregate(x_idx ~ query_chr, data = q_blocks, FUN = range)
    q_chr_out <- data.frame(
        chr   = q_chr_op$query_chr,
        xmin  = vapply(q_chr_op$x_idx, function(x) min(x) - 0.5, numeric(1)),
        xmax  = vapply(q_chr_op$x_idx, function(x) max(x) + 0.5, numeric(1)),
        y     = y_query,
        side  = "query",
        pair  = i,
        stringsAsFactors = FALSE
    )
    # subject side
    s_chr_op <- stats::aggregate(x_idx ~ subject_chr, data = s_blocks, FUN = range)
    s_chr_out <- data.frame(
        chr   = s_chr_op$subject_chr,
        xmin  = vapply(s_chr_op$x_idx, function(x) min(x) - 0.5, numeric(1)),
        xmax  = vapply(s_chr_op$x_idx, function(x) max(x) + 0.5, numeric(1)),
        y     = y_subject,
        side  = "subject",
        pair  = i,
        stringsAsFactors = FALSE
    )
    all_chrom[[length(all_chrom) + 1L]] <- rbind(q_chr_out, s_chr_out)
    
    # Build ribbons as curved polygons for each block pair
    ribbon_list_i <- vector("list", nrow(pair_summary))
    for (j in seq_len(nrow(pair_summary))) {
        ps <- pair_summary[j, ]
        
        start1 <- ps$q_x - ribbon_width
        end1   <- ps$q_x + ribbon_width
        start2 <- ps$s_x - ribbon_width
        end2   <- ps$s_x + ribbon_width
        
        poly <- calc_curvePolygon(
            start1 = start1,
            end1   = end1,
            start2 = start2,
            end2   = end2,
            y1     = y_query - 0.1,
            y2     = y_subject + 0.1,
            npts   = curve_npts,
            keepat = curve_keepat
        )
        
        poly$pair       <- i
        poly$block_pair <- paste(ps$query_synteny_block,
                                 ps$subject_synteny_block,
                                 sep = "->")
        poly$poly_id    <- paste0("pair", i, "_", poly$block_pair, "_", j)
        
        ribbon_list_i[[j]] <- poly
    }
    
    ribbon_op <- do.call(rbind, ribbon_list_i)
    all_ribbons[[length(all_ribbons) + 1L]] <- ribbon_op
    
    # Row labels
    if (is.null(genome_names)) {
        # Try to derive from list names like "NB_OL"
        nm <- if (!is.null(names(op))) names(op)[i] else ""
        if (!is.null(nm) && !is.na(nm) && nzchar(nm)) {
            sp <- strsplit(nm, "_")[[1L]]
            if (length(sp) >= 2L) {
                row_labels <- c(row_labels, sp[1L], sp[2L])
            } else {
                row_labels <- c(row_labels,
                                paste0("pair", i, "_query"),
                                paste0("pair", i, "_subject"))
            }
        } else {
            row_labels <- c(row_labels,
                            paste0("pair", i, "_query"),
                            paste0("pair", i, "_subject"))
        }
    }
    
    y_levels <- c(y_levels, y_query, y_subject)
}

if (length(all_ribbons) == 0L) {
    stop("No syntenic blocks left after filtering; try lowering 'min_block_genes'.",
         call. = FALSE)
}

chrom_op <- do.call(rbind, all_chrom)
ribbon_op <- do.call(rbind, all_ribbons)

# Final genome labels
if (!is.null(genome_names)) {
    if (length(genome_names) != length(y_levels)) {
        stop("Length of 'genome_names' (", length(genome_names),
             ") does not match number of genome rows (", length(y_levels), ").",
             call. = FALSE)
    }
    labels <- genome_names
} else {
    labels <- row_labels
}

# Build palettes
chr_levels <- sort(unique(chrom_op$chr))
chr_cols <- chr_palette(length(chr_levels))
names(chr_cols) <- chr_levels

bp_levels <- sort(unique(ribbon_op$block_pair))
bp_cols <- ribbon_palette(length(bp_levels))
names(bp_cols) <- bp_levels

# Build plot
p <- ggplot() +
    geom_segment(
        data = chrom_op,
        aes(x = xmin, xend = xmax, y = y, yend = y, colour = chr),
        linewidth = chr_lwd
    ) +
    geom_polygon(
        data = ribbon_op,
        aes(x = x, y = y, group = poly_id, fill = block_pair),
        alpha = ribbon_alpha,
        colour = NA
    ) +
    scale_colour_manual(values = chr_cols) +
    scale_fill_manual(values = bp_cols) +
    scale_y_continuous(
        breaks = sort(unique(y_levels)),
        labels = labels[order(y_levels)]
    ) +
    theme_void() +
    theme(
        legend.position = "none",
        axis.text.y = element_text(size = 10),
        plot.margin = margin(5.5, 20, 5.5, 5.5, unit = "pt")
    )

print(p)
