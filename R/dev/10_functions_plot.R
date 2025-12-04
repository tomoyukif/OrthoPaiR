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
##' @import ggplot2
##' @export
plot_riparian <- function(hdf5_fn = NULL,
                          object = NULL,
                          min_block_genes = 5L,
                          genome_names = NULL,
                          ribbon_alpha = 0.4,
                          ribbon_width = 0.45,
                          y_gap = 1,
                          curve_npts = 250L,
                          curve_keepat = round(curve_npts / 20),
                          chr_lwd = 2,
                          chr_palette = function(n) grDevices::hcl.colors(n, "Dark 3"),
                          ribbon_palette = function(n) grDevices::hcl.colors(n, "Zissou 1")) {

  if (is.null(hdf5_fn) && is.null(object)) {
    stop("Provide either `hdf5_fn` (HDF5 path) or `object`, but not both.")
  }
  if (!is.null(hdf5_fn) && !is.null(object)) {
    stop("Provide only one of `hdf5_fn` (HDF5 path) or `object`, not both.")
  }

  if (!is.null(hdf5_fn)) {
    op <- getOrthoPair(hdf5_fn = hdf5_fn, score = FALSE, loc = TRUE)
  } else {
    op <- getOrthoPair(object = object, score = FALSE, loc = TRUE)
  }

  if (!is.list(op) || length(op) == 0L) {
    stop("`op` must be a non-empty list of data.frames.")
  }

  n_pairs <- length(op)

  all_chrom <- list()
  all_ribbons <- list()
  y_levels <- c()
  row_labels <- c()

  for (i in seq_along(op)) {
    op_i <- op[[i]]

    if (!is.data.frame(op_i)) {
      warning("Element ", i, " of `op` is not a data.frame; skipping.")
      next
    }

    needed_cols <- c(
      "query_gene", "subject_gene",
      "query_chr", "query_synteny_block",
      "subject_chr", "subject_synteny_block"
    )
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
    op_i <- op_i[keep, , drop = FALSE]
    if (nrow(op_i) == 0L) {
      next
    }

    # Define vertical positions for this pair (top row higher y)
    # Pairs are stacked from top (first) to bottom (last)
    y_query   <- (n_pairs - i) * (y_gap + 1)
    y_subject <- y_query - y_gap

    # Summarize unique blocks for positioning on x-axis
    q_blocks <- unique(op_i[, c("query_chr", "query_synteny_block")])
    q_blocks <- q_blocks[order(q_blocks$query_chr, q_blocks$query_synteny_block), ]
    q_blocks$x_idx <- seq_len(nrow(q_blocks))

    s_blocks <- unique(op_i[, c("subject_chr", "subject_synteny_block")])
    s_blocks <- s_blocks[order(s_blocks$subject_chr, s_blocks$subject_synteny_block), ]
    s_blocks$x_idx <- seq_len(nrow(s_blocks))

    # Map block ID to x index
    q_key <- paste(q_blocks$query_chr, q_blocks$query_synteny_block, sep = "::")
    s_key <- paste(s_blocks$subject_chr, s_blocks$subject_synteny_block, sep = "::")
    op_i$q_key <- paste(op_i$query_chr, op_i$query_synteny_block, sep = "::")
    op_i$s_key <- paste(op_i$subject_chr, op_i$subject_synteny_block, sep = "::")
    q_map <- q_blocks$x_idx
    names(q_map) <- q_key
    s_map <- s_blocks$x_idx
    names(s_map) <- s_key
    op_i$q_x <- as.numeric(q_map[op_i$q_key])
    op_i$s_x <- as.numeric(s_map[op_i$s_key])

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
  invisible(p)
}

## -------------------------------------------------------------------------
## Local helpers (only used if not already available in the environment)
## -------------------------------------------------------------------------

if (!exists("scale_between")) {
  scale_between <- function(x, min, max) {
    rng <- range(x, na.rm = TRUE)
    if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
      return(rep((min + max) / 2, length(x)))
    }
    ( (x - rng[1]) / (rng[2] - rng[1]) ) * (max - min) + min
  }
}

if (!exists("calc_curvePolygon")) {
  calc_curvePolygon <- function(start1,
                                end1 = NULL,
                                start2,
                                end2 = NULL,
                                y1,
                                y2,
                                npts = 250,
                                keepat = round(npts / 20)) {
    cosine_points <- function(npts, keepat) {
      # initial number of points
      # grid to keep always
      grid <- seq(from = 0, to = pi, length.out = npts) # grid
      x <- (1 - cos(grid)) / max((1 - cos(grid))) # scaled cosine
      y <- grid / max(grid) # scaled grid
      # calculate slope for each point
      x1 <- x[-1];  y1 <- y[-1]
      x2 <- x[-length(x)];  y2 <- y[-length(y)]
      s <-  (y1 - y2) / (x1 - x2)
      # choose points that capture changes in slope
      ds <- cumsum(abs(diff(s))) * 5
      wh <- c(1, which(!duplicated(round(ds))), length(x))
      wh2 <- c(wh, seq(from = 0, to = length(x), by = round(keepat)))
      wh <- c(wh, wh2)[!duplicated(c(wh, wh2))]
      wh <- wh[order(wh)]
      cbind(x[wh], y[wh])
    }

    scaledCurve <- cosine_points(npts = npts, keepat = keepat)
    if (!is.null(end1) | !is.null(end2)) {
      sc1 <- scaledCurve[, 1]
      sc2 <- scaledCurve[, 2]

      tp <- rbind(
        start1 = data.frame(
          x = start1, y = y1
        ),
        poly1 = data.frame(
          x = scale_between(x = sc1, min = start1, max = start2),
          y = scale_between(x = sc2, min = y1, max = y2)
        ),
        start2 = data.frame(x = start2, y = y2),
        end2 = data.frame(
          x = end2, y = y2
        ),
        poly2 = data.frame(
          x = scale_between(x = sc1, min = end2, max = end1),
          y = scale_between(x = sc2, min = y2, max = y1)
        ),
        end1 = data.frame(
          x = end1, y = y1
        )
      )
    } else {
      tp <- data.frame(
        x = scale_between(x = scaledCurve[, 1], min = start1, max = start2),
        y = scale_between(x = scaledCurve[, 2], min = y1, max = y2)
      )
    }

    tp
  }
}
