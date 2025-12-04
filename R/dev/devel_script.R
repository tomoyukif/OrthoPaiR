library(rhop5)
library(ggplot2)
hop5_fn <- "../orthology/Osat_Ogla.h5"
h5 <- H5Fopen(hop5_fn)
target_col <- c("original_query_gene", "original_subject_gene", 
                "query_gene", "query_tx", 
                "subject_gene", "subject_tx", 
                "query_synteny_block", "subject_synteny_block",
                "class", "SOG", 
                "query_collapse", "subject_collapse")
target_col <- c(target_col, 
                "query_chr", "query_start",
                "query_end", "query_strand",
                "subject_chr", "subject_start",
                "subject_end", "subject_strand")
op <- NULL
name <- names(h5$orthopair)
for(i in seq_along(name)){
    col_names <- names(h5$orthopair[[i]])
    out_i <- h5$orthopair[[i]][, target_col]
    op <- c(op, list(out_i))
}
names(op) <- name
H5Fclose(h5)

if(is.null(object)){
    op <- getOrthoPair(hop5_fn = hop5_fn, score = FALSE, loc = TRUE)
    
} else {
    op <- getOrthoPair(object = object, score = FALSE, loc = TRUE)
}

# Prepare data for plotting
all_chrom <- list()
all_ribbons <- list()
y_levels <- c()
row_labels <- c()


for (i in seq_along(op)) {
    # Count genes per block on each side
    q_block_sizes <- tapply(op$query_gene,
                            op$query_synteny_block,
                            function(x) length(unique(x)))
    s_block_sizes <- tapply(op$subject_gene,
                            op$subject_synteny_block,
                            function(x) length(unique(x)))
    # For each row, require both blocks >= threshold
    q_ok <- q_block_sizes[as.character(op$query_synteny_block)] >= min_block_genes
    s_ok <- s_block_sizes[as.character(op$subject_synteny_block)] >= min_block_genes
    keep <- q_ok & s_ok
    op <- op[keep, ]
    if (nrow(op) == 0L) next
    
    # Define vertical positions for this pair (top row higher y)
    y_query   <- (n_pairs - i) * 2 + 2
    y_subject <- y_query - 1
    
    # Summarize unique blocks for positioning on x-axis
    q_blocks <- unique(op[, c("query_chr", "query_synteny_block")])
    q_blocks <- q_blocks[order(q_blocks$query_chr, q_blocks$query_synteny_block), ]
    q_blocks$x_idx <- seq_len(nrow(q_blocks))
    
    s_blocks <- unique(op[, c("subject_chr", "subject_synteny_block")])
    s_blocks <- s_blocks[order(s_blocks$subject_chr, s_blocks$subject_synteny_block), ]
    s_blocks$x_idx <- seq_len(nrow(s_blocks))
    
    # Map block ID to x index
    q_key <- paste(q_blocks$query_chr, q_blocks$query_synteny_block, sep = "::")
    s_key <- paste(s_blocks$subject_chr, s_blocks$subject_synteny_block, sep = "::")
    op$q_key <- paste(op$query_chr, op$query_synteny_block, sep = "::")
    op$s_key <- paste(op$subject_chr, op$subject_synteny_block, sep = "::")
    q_map <- q_blocks$x_idx
    names(q_map) <- q_key
    s_map <- s_blocks$x_idx
    names(s_map) <- s_key
    op$q_x <- as.numeric(q_map[op$q_key])
    op$s_x <- as.numeric(s_map[op$s_key])
    
    # Summarize one record per block pair
    pair_summary <- unique(op[, c("query_chr", "subject_chr",
                                  "query_synteny_block", "subject_synteny_block",
                                  "q_x", "s_x")])
    pair_summary$pair_id <- seq_len(nrow(pair_summary))
    
    # Build chromosome segments for this pair
    # query side
    q_chr_op <- aggregate(q_x ~ query_chr, data = q_blocks, FUN = range)
    q_chr_out <- data.frame(
        chr   = q_chr_op$query_chr,
        xmin  = vapply(q_chr_op$q_x, function(x) min(x) - 0.5, numeric(1)),
        xmax  = vapply(q_chr_op$q_x, function(x) max(x) + 0.5, numeric(1)),
        y = y_query,
        side  = "query",
        pair  = i,
        stringsAsFactors = FALSE
    )
    # subject side
    s_chr_op <- aggregate(s_x ~ subject_chr, data = s_blocks, FUN = range)
    s_chr_out <- data.frame(
        chr   = s_chr_op$subject_chr,
        xmin  = vapply(s_chr_op$s_x, function(x) min(x) - 0.5, numeric(1)),
        xmax  = vapply(s_chr_op$s_x, function(x) max(x) + 0.5, numeric(1)),
        y = y_subject,
        side  = "subject",
        pair  = i,
        stringsAsFactors = FALSE
    )
    all_chrom[[length(all_chrom) + 1L]] <- rbind(q_chr_out, s_chr_out)
    
    # Build ribbons (as simple quadrilaterals) for each block pair
    rib_op <- data.frame(
        xmin = pair_summary$q_x - 0.45,
        xmax = pair_summary$q_x + 0.45,
        ymin = y_subject + 0.05,
        ymax = y_query   - 0.05,
        pair = i,
        block_pair = paste(pair_summary$query_synteny_block,
                           pair_summary$subject_synteny_block,
                           sep = "->"),
        stringsAsFactors = FALSE
    )
    all_ribbons[[length(all_ribbons) + 1L]] <- rib_op
    
    # Row labels
    if (is.null(genome_names)) {
        # Try to derive from list names like "NB_OL"
        if (!is.null(names(orthopair_list)) && names(orthopair_list)[i] != "") {
            nm <- names(orthopair_list)[i]
            sp <- strsplit(nm, "_")[[1L]]
            if (length(sp) >= 2L) {
                row_labels <- c(row_labels, sp[1L], sp[2L])
            } else {
                row_labels <- c(row_labels,
                                paste0("pair", i, "_query"),
                                paste0("pair", i, "_subject"))
            }
        } else {row_labels <- c(row_labels,
                                paste0("pair", i, "_query"),
                                paste0("pair", i, "_subject"))
        }
    }
    y_levels <- c(y_levels, y_query, y_subject)
    chrom_op <- do.call(rbind, all_chrom)
    ribbon_op <- do.call(rbind, all_ribbons)
    if (nrow(ribbon_op) == 0L) {
        stop("No syntenic blocks left after filtering; try lowering 'min_block_genes'.", call. = FALSE)}
    
    # Final genome labels
    if (!is.null(genome_names)) {
        if (length(genome_names) != length(y_levels)) {
            stop("Length of 'genome_names' (", length(genome_names), ") does not match number of genome rows (", length(y_levels), ").", call. = FALSE)
        }
        labels <- genome_names
    } else {
        labels <- row_labels
    }
}
# Build plot
suppressWarnings({
    p <- opggplot() +
        opgeom_segment(data = chrom_op, opaes(x = xmin, xend = xmax, y = y, yend = y, colour = chr),linewidth = 2) +
        opgeom_rect(data = ribbon_op, opaes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = factor(block_pair)),alpha = 0.4,colour = NA) +
        opscale_y_continuous(breaks = sort(unique(y_levels)),labels = labels[order(y_levels)]) +
        optheme_void() +optheme(legend.position = "none",axis.text.y = opelement_text(size = 10),plot.margin = opmargin(5.5, 20, 5.5, 5.5, unit = "pt"))
})
print(p)
invisible(p)
