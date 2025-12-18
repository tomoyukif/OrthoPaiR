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
    op_i <- op[[i]]
    gene_order <- order(op_i$query_chr, 
                        op_i$query_start,
                        op_i$subject_chr, 
                        op_i$subject_start)
    op_i <- op_i[gene_order, ]
    needed_cols <- c("query_gene", "subject_gene",
                     "query_chr", "query_synteny_block",
                     "subject_chr", "subject_synteny_block")
    missing_cols <- setdiff(needed_cols, colnames(op_i))
    if (length(missing_cols) > 0L) {
        warning("Element ", i, " is missing required columns: ",
                paste(missing_cols, collapse = ", "), "; skipping.")
        next
    }
    
    # Count genes per block on each side
    q_block_sizes <- tapply(op_i$query_gene,
                            op_i$query_synteny_block,
                            function(x) length(unique(x)))
    s_block_sizes <- tapply(op_i$subject_gene,
                            op_i$subject_synteny_block,
                            function(x) length(unique(x)))
    
    # For each row, require both blocks >= threshold
    q_ok <- q_block_sizes[as.character(op_i$query_synteny_block)] >= min_block_genes
    s_ok <- s_block_sizes[as.character(op_i$subject_synteny_block)] >= min_block_genes
    keep <- q_ok & s_ok
    keep[is.na(keep)] <- FALSE
    op_i <- op_i[keep, , drop = FALSE]
    
    # Define vertical positions for this pair (top row higher y)
    # Pairs are stacked from top (first) to bottom (last)
    y_query   <- (n_pairs - i) * (y_gap + 1)
    y_subject <- y_query - y_gap
    
    # Summarize unique blocks for positioning on x-axis
    q_blocks <- data.frame(block = tapply(op_i$query_synteny_block,
                                          op_i$query_synteny_block, "[", 1),
                           chr = tapply(op_i$query_chr,
                                        op_i$query_synteny_block, "[", 1),
                           start = tapply(op_i$query_start,
                                          op_i$query_synteny_block, min),
                           end = tapply(op_i$query_end,
                                        op_i$query_synteny_block, max))
    q_blocks$width <- q_blocks$end - q_blocks$start
    q_blocks <- q_blocks[order(q_blocks$chr, q_blocks$start), ]
    q_blocks$x_idx <- seq_len(nrow(q_blocks))
    
    s_blocks <- data.frame(block = tapply(op_i$subject_synteny_block,
                                          op_i$subject_synteny_block, "[", 1),
                           chr = tapply(op_i$subject_chr,
                                        op_i$subject_synteny_block, "[", 1),
                           start = tapply(op_i$subject_start,
                                          op_i$subject_synteny_block, min),
                           end = tapply(op_i$subject_end,
                                        op_i$subject_synteny_block, max))
    s_blocks$width <- s_blocks$end - s_blocks$start
    s_blocks <- s_blocks[order(s_blocks$chr, s_blocks$start), ]
    s_blocks$x_idx <- seq_len(nrow(s_blocks))
    
    # Map block ID to x index
    q_key <- paste(q_blocks$chr, q_blocks$block, sep = "::")
    s_key <- paste(s_blocks$chr, s_blocks$block, sep = "::")
    op_i$q_key <- paste(op_i$query_chr, op_i$query_synteny_block, sep = "::")
    op_i$s_key <- paste(op_i$subject_chr, op_i$subject_synteny_block, sep = "::")
    q_map <- q_blocks$x_idx
    names(q_map) <- q_key
    s_map <- s_blocks$x_idx
    names(s_map) <- s_key
    op_i$q <- as.numeric(q_map[op_i$q_key])
    op_i$s <- as.numeric(s_map[op_i$s_key])
    lcb <- unique(subset(op_i, select = c(q, s)))
    dt <- as.data.table(lcb)
    setorder(dt, q, s)
    n <- nrow(dt)
    
    # --- Union-Find ---
    parent <- seq_len(n)
    findp <- function(x) { while (parent[x] != x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }; x }
    unionp <- function(a,b) { ra <- findp(a); rb <- findp(b); if (ra != rb) parent[rb] <<- ra }
    
    # (q,s) -> 行番号 のハッシュ
    key <- paste(dt$q, dt$s, sep=":")
    idx_map <- new.env(hash = TRUE, parent = emptyenv())
    
    # 8近傍（Δq,Δs）
    dqs <- expand.grid(dq = -1:1, ds = -1:1)
    dqs <- dqs[!(dqs$dq == 0 & dqs$ds == 0), , drop=FALSE]
    
    for (i in seq_len(n)) {
        q <- dt$q[i]; s <- dt$s[i]
        
        # 近傍が既に登録されていれば union
        for (k in seq_len(nrow(dqs))) {
            nq <- q + dqs$dq[k]
            ns <- s + dqs$ds[k]
            nk <- paste(nq, ns, sep=":")
            j <- idx_map[[nk]]
            if (!is.null(j)) unionp(i, j)
        }
        
        # 自分を登録
        idx_map[[ key[i] ]] <- i
    }
    
    # 連結成分ID
    comp <- vapply(seq_len(n), findp, integer(1))
    dt[, comp := comp]
    
    # comp（Union-Findの代表値）を 1,2,3,... に変換
    dt[, new_block_id := as.integer(factor(comp))]
    dt[, pair_id := paste(q, s, sep = "_")]
    block_summary <- dt[, .(
        q_range = sprintf("%d-%d", min(q), max(q)),
        s_range = sprintf("%d-%d", min(s), max(s)),
        n_pairs = .N
    ), by = new_block_id][order(new_block_id)]
    
    op_i$pair_id <- paste(op_i$q, op_i$s, sep = "_")
    hit <- match(op_i$pair_id, dt$pair_id)
    op_i$block_id <- dt$new_block_id[hit]
    
    block_df <- data.frame(block_id = tapply(op_i$block_id, op_i$block_id, "[", 1))
    
    # Summarize one record per block pair
    pair_summary <- unique(op_i[, c("query_chr", "subject_chr",
                                    "query_synteny_block", "subject_synteny_block",
                                    "q_x", "s_x")])
    pair_summary$pair_id <- seq_len(nrow(pair_summary))
    
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
