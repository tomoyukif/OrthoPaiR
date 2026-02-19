# Ideas to Speed Up `makeOrthoGraph()` and `graph2df()`

This document outlines bottlenecks in the orthology graph pipeline and concrete ways to improve performance, including optional Rcpp usage (the package already has Rcpp in `DESCRIPTION`).

---

## 1. Profiling First

Before heavy refactoring, profile to confirm where time is spent:

```r
# Install if needed: install.packages("profvis")
library(profvis)
reorg_hdf5_fn <- "path/to/reorg_orthopair.h5"
graph <- makeOrthoGraph(hdf5_fn = reorg_hdf5_fn)
profvis({
  graph2df(hdf5_fn = reorg_hdf5_fn, graph = graph, orthopair_fn = "orthopair_list.csv")
})
```

Use `summaryRprof()` or `profvis` to see whether most time is in `ego()`, list/string operations, or I/O.

---

## 2. `makeOrthoGraph()` Optimizations

### 2.1 Avoid repeated `lapply` + `do.call("rbind")` on many data.frames

**Current pattern:**  
`orthopair_gene` is a list of data.frames; you `subset` each, then `do.call("rbind", ...)`. Building many small data.frames and rbinding is slow.

**Improvements:**

- **Single pass over the list:** Loop once, collect `query_gene`, `subject_gene`, `mutual_ci`, `class` into **vectors** (or a pre-allocated matrix), then build one data.frame at the end. This avoids many small allocations and one big rbind.
- **data.table:** If you allow `data.table` as Suggests, `rbindlist(orthopair_gene[, .(query_gene, subject_gene, mutual_ci, class)])` is much faster than `do.call("rbind", ...)`.
- **Rcpp:** One C++ function that takes the list of orthopair tables (as list of CharacterVector/ NumericVector), flattens into four vectors, and returns a DataFrame. No need for R-side rbind.

### 2.2 Single genome lookup

**Current:**  
`match(edges$query_gene, gene_list$gene)` and `match(edges$subject_gene, gene_list$gene)` (and again for vertices). Each `match` is O(n).

**Improvement:**  
Create once a named vector (or environment) `gene_to_genome <- setNames(gene_list$genome, gene_list$gene)` and use `gene_to_genome[edges$query_gene]` and `gene_to_genome[edges$subject_gene]`. Factor levels can also be used for integer indexing if you prefer.

### 2.3 Build graph from edge list in one go

**Current:**  
`make_empty_graph(n = 0)` then `add_vertices` then `add_edges`. Multiple igraph API calls.

**Improvement:**  
If igraph allows building from a full edge list and vertex attribute table in one step, use that (e.g. `graph_from_data_frame()` with a vertices data.frame and edges data.frame). That minimizes internal reallocation. Then add orphan vertices in a single `add_vertices` call if needed.

### 2.4 Typo

- Fix message: `"onyl"` → `"only"` in `makeOrthoGraph()`.

---

## 3. `graph2df()` Optimizations (Main Bottleneck)

Most cost is expected in **neighborhood enumeration** and **post-processing**. Below: algorithmic and implementation ideas.

### 3.1 Parallelize `ego()` by splitting vertices

**Current:**  
`ego(graph, order = n_genomes)` computes the `n_genomes`-neighborhood for **every** vertex in one call. Internally that is one BFS per vertex, so time is O(V × (V + E)) in the worst case and is the main bottleneck.

**Improvement:**  
Split vertex indices into chunks and call `ego(..., nodes = chunk)` per chunk in parallel:

```r
library(parallel)
n_vert <- vcount(graph)
chunk_size <- max(1L, ceiling(n_vert / n_core))  # n_core from function arg
chunks <- split(seq_len(n_vert), ceiling(seq_len(n_vert) / chunk_size))
ego_list <- mclapply(chunks, function(nodes) {
  ego(graph = graph, order = n_genomes, nodes = nodes)
}, mc.cores = n_core)
ego_list <- unlist(ego_list, recursive = FALSE)
ego_list <- lapply(ego_list, names)
# ... rest unchanged
```

- Use `n_core` (or add it to `graph2df(..., n_core = 1)`) so the pipeline can control parallelism.
- Ensure `graph` and `n_genomes` are in scope and that you don’t oversubscribe cores (e.g. respect `n_threads` from the main pipeline).

This keeps the same algorithm but parallelizes the expensive ego computations.

### 3.2 Replace `ego()` with a single batch BFS in Rcpp (largest gain)

**Idea:**  
Implement in C++ a function that, for each vertex, returns the set of vertex indices within distance `order` (e.g. `n_genomes`). Return a list of integer vectors. Then in R you only need to map vertex indices back to names and genomes.

**Why it helps:**  
- One pass over the graph; better cache locality.  
- No per-call R/igraph overhead.  
- You can optionally use a compact representation (e.g. adjacency list built once) and run BFS in a tight loop.

**Sketch:**

- **Input:**  
  - `edges`: 2-column integer matrix (0-based or 1-based vertex IDs, be consistent).  
  - `order`: integer (max distance).  
  - `n_vert`: number of vertices.
- **Output:**  
  - List of length `n_vert`; element `i` = integer vector of vertex indices in the `order`-neighborhood of vertex `i`.

**Algorithm:**  
For each vertex, run BFS (or use igraph’s C library via Rcpp if you link against it), and stop at depth `order`. Push reachable vertices into a set or vector.

**R side:**  
Build once a character vector `vnames <- V(graph)$name` and use `vnames[ego_indices[[i]]]` to get names for the i-th vertex. Then replace the current `ego_list` construction with this list of name vectors and keep the rest of `graph2df()` logic (genome extraction, full_pair, valid_group1/2, etc.) as is.

### 3.3 Avoid repeated string operations on the same data

**Current:**  
Multiple `lapply(..., sub, pattern = "&.+", replace = "")` and `sub(".+&", "", ...)` on the same or similar vectors. Also `apply(out, 2, sub, ...)` which is a slow loop over cells.

**Improvements:**

- **Do the split once:**  
  When you set `v$name <- paste(v$genome, v$name, sep = "&")`, also keep two vectors: `v_name_only` (gene) and `v_genome` (genome). When you get ego as **vertex indices**, use `v_name_only` and `v_genome` to build genome membership without repeatedly parsing `"genome&gene"` strings.
- **Use `stringi::stri_replace_*`:**  
  If you keep the current string-based workflow, replace `sub()` with `stringi::stri_replace_first_regex()` (or the appropriate variant); it’s much faster on long vectors.
- **Column-wise operations instead of `apply`:**  
  For `out_label` and `out_values`, avoid `apply(..., 2, sub, ...)`. Either:
  - Loop over columns and call `sub()` once per column (vectorized), or  
  - Keep genome and gene as separate columns from the start so you never need to parse `"genome&gene"` in the big table.

### 3.4 More efficient handling of “full” vs “rest” groups

**Current:**  
- `ego_list_genomes <- lapply(ego_list, sub, pattern = "&.+", replace = "")`  
- Then `unique`, `length`, etc. on each element.

**Improvement:**  
If you have ego as **integer vertex indices**, you can:

- Build once a vector `vertex_genome` (length = vcount(graph)).  
- For each ego set `idx`, get `vertex_genome[idx]` and then `length(unique(vertex_genome[idx]))` and `length(idx)` to decide full_pair. No string parsing in the inner loop.

### 3.5 `valid_group1`: avoid `do.call(rbind)` on many rows

**Current:**  
`valid_group1 <- lapply(valid_group1, sort)` then `do.call(rbind, valid_group1)`.

**Improvement:**  
- If each group has the same length (n_genomes), store as a matrix: one row per group, columns = genomes. Fill the matrix in a loop or with a fast Rcpp function that takes the list of sorted character vectors and writes into a pre-allocated matrix. Then `duplicated.matrix(valid_group1)` is fast.  
- Alternatively use `data.table::rbindlist` with a key and then `unique()`.

### 3.6 `rest_list` and `expand.grid`

**Current:**  
`mapply(split, ...)`, then `lapply(rest_list, expand.grid)`. For groups with many genes per genome, `expand.grid` can explode in size and time.

**Improvements:**

- **Early exit / cap:**  
  If the product of group sizes exceeds a threshold (e.g. 1e6), skip or cap (e.g. take a sample or only first N combinations) to avoid memory and CPU blow-up.
- **Rcpp for combinations:**  
  A C++ function that takes per-genome vectors of gene names and returns a matrix (or list of vectors) of all combinations, written in one pass without building a huge list in R.  
- **data.table::CJ:**  
  If you use data.table, `CJ(genome1 = x1, genome2 = x2, ..., unique = TRUE)` can be faster and more memory-efficient than `expand.grid` for building the combinations.

### 3.7 Faster CSV write

**Current:**  
`write.csv(out, file = orthopair_fn, row.names = FALSE)`.

**Improvement:**  
- Use `data.table::fwrite(out, orthopair_fn)` if you allow data.table; it is much faster for large tables.  
- Or at least use `write.table(out, orthopair_fn, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)` and avoid `write.csv`’s default quoting if you don’t need it.

### 3.8 Unused arguments

**Current:**  
`graph2df(..., n_core = 1, n_batch = 50)` — `n_core` and `n_batch` are not used.

**Improvement:**  
Use `n_core` for the parallel ego chunks (see 3.1). If you implement chunked processing (e.g. writing the table in batches), use `n_batch` there; otherwise remove it or document “reserved for future use.”

---

## 4. Optional: Rcpp Summary

| Place              | What to do in Rcpp |
|--------------------|--------------------|
| **makeOrthoGraph** | Flatten list of orthopair tables into four vectors (query_gene, subject_gene, mutual_ci, class) and return a single DataFrame. Optionally build edge list and vertex table for igraph. |
| **graph2df**       | **Batch BFS:** input = edge list + n_vert + order; output = list of integer vectors (neighborhoods). Biggest expected win. |
| **graph2df**       | Build valid_group1 matrix from list of sorted name vectors (no rbind). |
| **graph2df**       | Generate combination table (expand.grid replacement) from list of per-genome character vectors. |

---

## 5. Suggested Order of Implementation

1. **Low effort:**  
   - Fix typo in `makeOrthoGraph`; use named vector for genome lookup; single data.frame build (or rbindlist) in `makeOrthoGraph`; use `fwrite` or `write.table` in `graph2df`; use `n_core` to parallelize `ego()` by vertex chunks.

2. **Medium effort:**  
   - Replace repeated `sub()` with one-time parsing or stringi; avoid `apply(..., 2, sub)`; build valid_group1 as a matrix and use `duplicated.matrix`.

3. **High effort / max speed:**  
   - Rcpp batch BFS replacing `ego()`; optional Rcpp for flattening orthopair list and for combination generation in `graph2df`.

---

## 6. Testing

After each change:

- Compare output of `graph2df()` (and optionally `makeOrthoGraph()`) with the previous implementation on a small and a medium-sized run (same `reorg_orthopair.h5` and graph).  
- Check that row counts and column structure of `orthopair_list.csv` match, and that SOG IDs and gene memberships are equivalent (e.g. same groups, possibly reordered).  
- Use `system.time()` or `microbenchmark` on a fixed dataset to ensure improvements are real.

This should give you a clear path to faster `makeOrthoGraph()` and especially `graph2df()` while keeping results identical.
