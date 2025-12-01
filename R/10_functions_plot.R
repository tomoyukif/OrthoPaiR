#' Riparian plot of syntenic blocks
#'
#' This function draws a \"riparian\" style plot to visualize syntenic blocks
#' between genomes using the orthologous pairs table produced by OrthoPaiR
#' (e.g. `orthopair.csv` created by `graph2df()` / `orthopair()`).
#'
#' The basic use case is a single pair of genomes, where the first genome
#' (query) is shown on the top row and the second genome (subject) on the
#' second row. Syntenic blocks are defined by the integer columns
#' `query_synteny_block` and `subject_synteny_block`; for each pair of blocks
#' with the same block IDs, a ribbon is drawn connecting the corresponding
#' regions between the two genomes.
#'
#' You can optionally filter out small blocks by specifying a minimum number
#' of genes per block via the `min_block_genes` argument (default: 10).
#'
#' For more than two genomes, you can pass a named list of ortholog tables in
#' `orthopair_list`. Each element should be a data.frame for one adjacent pair
#' of genomes (e.g. NB–OL, OL–WK21, ...), having the same columns as
#' `orthopair`. The pairs are stacked vertically in the order of the list,
#' so additional genomes appear below the first two in the plot.
#'
#' @param orthopair A data.frame of ortholog pairs for a single query–subject
#'   genome pair. It must contain at least the following columns:
#'   `query_gene`, `subject_gene`, `query_synteny_block`,
#'   `subject_synteny_block`. If the chromosome columns `query_chr` and
#'   `subject_chr` are missing, they are inferred from the GFF files via the
#'   `query_gff` and `subject_gff` arguments.
#' @param orthopair_list Optional named list of ortholog tables for more than
#'   two genomes. If supplied, `orthopair` can be omitted. Each element must
#'   have the same required columns as `orthopair`. Names are used for
#'   labelling genome rows; if a name has the form \"A_B\", the top row is
#'   labelled \"A\" and the bottom row \"B\".
#' @param query_gff Path to a query genome GFF file used to obtain gene
#'   locations when `orthopair` does not already contain `query_chr`.
#' @param subject_gff Path to a subject genome GFF file used to obtain gene
#'   locations when `orthopair` does not already contain `subject_chr`.
#' @param min_block_genes Integer, minimum number of genes per block to retain
#'   (on both query and subject sides). Default is 10.
#' @param genome_names Optional character vector of labels for genome rows.
#'   For a single pair this must be length 2 (query, subject). For a list of
#'   pairs it must be length `2 * length(orthopair_list)`. If `NULL`, labels
#'   are derived from `orthopair_list` names or set to generic values.
#' @return A `ggplot` object (the plot is also printed as a side effect).
#' @export
#'
#' @importFrom ggplot2 ggplot geom_segment geom_rect aes theme_void theme
#'   scale_y_continuous scale_fill_manual guides guide_legend
plot_riparian <- function(orthopair = NULL,
                          orthopair_list = NULL,
                          query_gff = NULL,
                          subject_gff = NULL,
                          min_block_genes = 10,
                          genome_names = NULL) {
    if (is.null(orthopair) && is.null(orthopair_list)) {
        stop(\"Specify at least 'orthopair' or 'orthopair_list'.\", call. = FALSE)
    }
    if (!is.null(orthopair)) {
        if (!is.data.frame(orthopair)) {
            stop(\"'orthopair' must be a data.frame.\", call. = FALSE)
        }
        base_required <- c(\"query_gene\", \"subject_gene\",
                           \"query_synteny_block\", \"subject_synteny_block\")
        missing_cols <- base_required[!base_required %in% names(orthopair)]
        if (length(missing_cols) > 0) {
            stop(\"'orthopair' is missing required columns: \",
                 paste(missing_cols, collapse = \", \"), call. = FALSE)
        }
        # If chr columns are missing, infer them from GFF files
        need_query_chr  <- !\"query_chr\"   %in% names(orthopair)
        need_subject_chr <- !\"subject_chr\" %in% names(orthopair)
        if (need_query_chr || need_subject_chr) {
            if (is.null(query_gff) || is.null(subject_gff)) {
                stop(\"When 'orthopair' lacks 'query_chr'/'subject_chr', \",\n+                     \"'query_gff' and 'subject_gff' must be supplied.\",\n+                     call. = FALSE)\n+            }\n+            # Import GFFs and build gene -> chr map\n+            q_gff <- rtracklayer::import.gff3(query_gff)\n+            s_gff <- rtracklayer::import.gff3(subject_gff)\n+\n+            # Prefer gene-level features, otherwise use any with gene_id\n+            if (!is.null(q_gff$gene_id)) {\n+                q_gene_df <- as.data.frame(q_gff)\n+                q_gene_df <- q_gene_df[!is.na(q_gene_df$gene_id), ]\n+                q_map_chr <- as.character(q_gene_df$seqnames)\n+                names(q_map_chr) <- q_gene_df$gene_id\n+            } else {\n+                q_gene_df <- as.data.frame(q_gff)\n+                if (!\"ID\" %in% names(q_gene_df)) {\n+                    stop(\"Query GFF must have either 'gene_id' or 'ID' column.\",\n+                         call. = FALSE)\n+                }\n+                q_map_chr <- as.character(q_gene_df$seqnames)\n+                names(q_map_chr) <- q_gene_df$ID\n+            }\n+\n+            if (!is.null(s_gff$gene_id)) {\n+                s_gene_df <- as.data.frame(s_gff)\n+                s_gene_df <- s_gene_df[!is.na(s_gene_df$gene_id), ]\n+                s_map_chr <- as.character(s_gene_df$seqnames)\n+                names(s_map_chr) <- s_gene_df$gene_id\n+            } else {\n+                s_gene_df <- as.data.frame(s_gff)\n+                if (!\"ID\" %in% names(s_gene_df)) {\n+                    stop(\"Subject GFF must have either 'gene_id' or 'ID' column.\",\n+                         call. = FALSE)\n+                }\n+                s_map_chr <- as.character(s_gene_df$seqnames)\n+                names(s_map_chr) <- s_gene_df$ID\n+            }\n+\n+            if (need_query_chr) {\n+                orthopair$query_chr <- as.character(q_map_chr[orthopair$query_gene])\n+            }\n+            if (need_subject_chr) {\n+                orthopair$subject_chr <- as.character(s_map_chr[orthopair$subject_gene])\n+            }\n+        }\n+        # Now require chr columns to exist\n+        required_cols <- c(\"query_chr\", \"subject_chr\",\n+                           \"query_synteny_block\", \"subject_synteny_block\",\n+                           \"query_gene\", \"subject_gene\")\n+        missing_cols2 <- required_cols[!required_cols %in% names(orthopair)]\n+        if (length(missing_cols2) > 0) {\n+            stop(\"Unable to add chromosome information; missing columns: \",\n+                 paste(missing_cols2, collapse = \", \"), call. = FALSE)\n+        }\n         orthopair_list <- c(list(orthopair), orthopair_list)
        orthopair_list <- c(list(orthopair), orthopair_list)
    }
    if (!is.list(orthopair_list)) {
        stop(\"'orthopair_list' must be a list of data.frames.\", call. = FALSE)
    }

    n_pairs <- length(orthopair_list)
    if (n_pairs == 0L) {
        stop(\"'orthopair_list' is empty.\", call. = FALSE)
    }

    # Prepare data for plotting
    all_chrom <- list()
    all_ribbons <- list()
    y_levels <- c()
    row_labels <- c()

    for (i in seq_len(n_pairs)) {
        df <- orthopair_list[[i]]
        # Basic sanity check
        req_cols <- c(\"query_chr\", \"subject_chr\",
                      \"query_synteny_block\", \"subject_synteny_block\",
                      \"query_gene\", \"subject_gene\")
        if (!all(req_cols %in% names(df))) {
            stop(\"Element \", i,
                 \" of 'orthopair_list' is missing required columns.\", call. = FALSE)
        }
        # Remove NA blocks
        df <- df[!is.na(df$query_synteny_block) &
                     !is.na(df$subject_synteny_block), ]
        if (nrow(df) == 0L) next

        # Count genes per block on each side
        q_block_sizes <- tapply(df$query_gene,
                                df$query_synteny_block,
                                function(x) length(unique(x)))
        s_block_sizes <- tapply(df$subject_gene,
                                df$subject_synteny_block,
                                function(x) length(unique(x)))
        # For each row, require both blocks >= threshold
        q_ok <- q_block_sizes[as.character(df$query_synteny_block)] >= min_block_genes
        s_ok <- s_block_sizes[as.character(df$subject_synteny_block)] >= min_block_genes
        keep <- q_ok & s_ok
        df <- df[keep, ]
        if (nrow(df) == 0L) next

        # Define vertical positions for this pair (top row higher y)
        y_query   <- (n_pairs - i) * 2 + 2
        y_subject <- y_query - 1

        # Summarize unique blocks for positioning on x-axis
        q_blocks <- unique(df[, c(\"query_chr\", \"query_synteny_block\")])
        q_blocks <- q_blocks[order(q_blocks$query_chr, q_blocks$query_synteny_block), ]
        q_blocks$x_idx <- seq_len(nrow(q_blocks))

        s_blocks <- unique(df[, c(\"subject_chr\", \"subject_synteny_block\")])
        s_blocks <- s_blocks[order(s_blocks$subject_chr, s_blocks$subject_synteny_block), ]
        s_blocks$x_idx <- seq_len(nrow(s_blocks))

        # Map block ID to x index
        q_key <- paste(q_blocks$query_chr, q_blocks$query_synteny_block, sep = \"::\")
        s_key <- paste(s_blocks$subject_chr, s_blocks$subject_synteny_block, sep = \"::\")
        df$q_key <- paste(df$query_chr, df$query_synteny_block, sep = \"::\")
        df$s_key <- paste(df$subject_chr, df$subject_synteny_block, sep = \"::\")
        q_map <- q_blocks$x_idx
        names(q_map) <- q_key
        s_map <- s_blocks$x_idx
        names(s_map) <- s_key
        df$q_x <- as.numeric(q_map[df$q_key])
        df$s_x <- as.numeric(s_map[df$s_key])

        # Summarize one record per block pair
        pair_summary <- unique(df[, c(\"query_chr\", \"subject_chr\",
                                      \"query_synteny_block\", \"subject_synteny_block\",
                                      \"q_x\", \"s_x\")])
        pair_summary$pair_id <- seq_len(nrow(pair_summary))

        # Build chromosome segments for this pair
        # query side
        q_chr_df <- aggregate(q_x ~ query_chr, data = q_blocks, FUN = range)
        q_chr_out <- data.frame(
            chr   = q_chr_df$query_chr,
            xmin  = vapply(q_chr_df$q_x, function(x) min(x) - 0.5, numeric(1)),
            xmax  = vapply(q_chr_df$q_x, function(x) max(x) + 0.5, numeric(1)),
            y     = y_query,
            side  = \"query\",
            pair  = i,
            stringsAsFactors = FALSE
        )
        # subject side
        s_chr_df <- aggregate(s_x ~ subject_chr, data = s_blocks, FUN = range)
        s_chr_out <- data.frame(
            chr   = s_chr_df$subject_chr,
            xmin  = vapply(s_chr_df$s_x, function(x) min(x) - 0.5, numeric(1)),
            xmax  = vapply(s_chr_df$s_x, function(x) max(x) + 0.5, numeric(1)),
            y     = y_subject,
            side  = \"subject\",
            pair  = i,
            stringsAsFactors = FALSE
        )
        all_chrom[[length(all_chrom) + 1L]] <- rbind(q_chr_out, s_chr_out)

        # Build ribbons (as simple quadrilaterals) for each block pair
        rib_df <- data.frame(
            xmin = pair_summary$q_x - 0.45,
            xmax = pair_summary$q_x + 0.45,
            ymin = y_subject + 0.05,
            ymax = y_query   - 0.05,
            pair = i,
            block_pair = paste(pair_summary$query_synteny_block,
                               pair_summary$subject_synteny_block,
                               sep = \"->\"),
            stringsAsFactors = FALSE
        )
        all_ribbons[[length(all_ribbons) + 1L]] <- rib_df

        # Row labels
        if (is.null(genome_names)) {
            # Try to derive from list names like \"NB_OL\"
            if (!is.null(names(orthopair_list)) && names(orthopair_list)[i] != \"\") {
                nm <- names(orthopair_list)[i]
                sp <- strsplit(nm, \"_\")[[1L]]
                if (length(sp) >= 2L) {
                    row_labels <- c(row_labels, sp[1L], sp[2L])
                } else {
                    row_labels <- c(row_labels,
                                    paste0(\"pair\", i, \"_query\"),
                                    paste0(\"pair\", i, \"_subject\"))
                }
            } else {\n+                row_labels <- c(row_labels,
                                paste0(\"pair\", i, \"_query\"),
                                paste0(\"pair\", i, \"_subject\"))\n+            }\n+        }\n+        y_levels <- c(y_levels, y_query, y_subject)\n+    }\n+\n+    chrom_df <- do.call(rbind, all_chrom)\n+    ribbon_df <- do.call(rbind, all_ribbons)\n+\n+    if (nrow(ribbon_df) == 0L) {\n+        stop(\"No syntenic blocks left after filtering; try lowering 'min_block_genes'.\",\n+             call. = FALSE)\n+    }\n+\n+    # Final genome labels\n+    if (!is.null(genome_names)) {\n+        if (length(genome_names) != length(y_levels)) {\n+            stop(\"Length of 'genome_names' (\", length(genome_names),\n+                 \") does not match number of genome rows (\", length(y_levels), \").\",\n+                 call. = FALSE)\n+        }\n+        labels <- genome_names\n+    } else {\n+        labels <- row_labels\n+    }\n+\n+    # Build plot\n+    suppressWarnings({\n+        p <- ggplot2::ggplot() +\n+            ggplot2::geom_segment(\n+                data = chrom_df,\n+                ggplot2::aes(x = xmin, xend = xmax, y = y, yend = y,\n+                             colour = chr),\n+                linewidth = 2\n+            ) +\n+            ggplot2::geom_rect(\n+                data = ribbon_df,\n+                ggplot2::aes(xmin = xmin, xmax = xmax,\n+                             ymin = ymin, ymax = ymax,\n+                             fill = factor(block_pair)),\n+                alpha = 0.4,\n+                colour = NA\n+            ) +\n+            ggplot2::scale_y_continuous(\n+                breaks = sort(unique(y_levels)),\n+                labels = labels[order(y_levels)]\n+            ) +\n+            ggplot2::theme_void() +\n+            ggplot2::theme(\n+                legend.position = \"none\",\n+                axis.text.y = ggplot2::element_text(size = 10),\n+                plot.margin = ggplot2::margin(5.5, 20, 5.5, 5.5, unit = \"pt\")\n+            )\n+    })\n+\n+    print(p)\n+    invisible(p)\n+}\n+\n*** End Patch"}]}/>

