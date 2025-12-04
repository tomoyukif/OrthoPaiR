## OrthoPaiR

[![R](https://img.shields.io/badge/R-%3E%3D4.1-blue)](https://www.r-project.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**OrthoPaiR** is an R package for **synteny‑aware ortholog detection** across two or more genomes.
It organises genome annotations, runs reciprocal sequence searches, and summarises orthologous relationships as tables and graphs.

### Features

- **Multi‑genome support**: run ortholog detection for all pairwise combinations of genomes.
- **Synteny‑aware pairing**: combine reciprocal best hits with genomic context.
- **Network‑based ortholog groups**: summarise orthologs as an orthology graph and CSV tables.
- **Genome‑wise outputs**: reorganised GFF files, ortholog lists, and orphan genes.

### Prerequisites

- **R packages** (installed via CRAN and Bioconductor):

```r
install.packages("devtools")
install.packages("BiocManager")

BiocManager::install(c(
  "GenomicRanges", "GenomeInfoDb", "GenomicFeatures",
  "rhdf5", "BiocGenerics", "Biobase", "IRanges",
  "rtracklayer", "Biostrings", "S4Vectors", "BSgenome"
))

install.packages(c("parallel", "ggplot2", "dplyr"))
```

- **External tools**
  - [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs) (`makeblastdb` in `blast_path`) – for DNA sequence similarity
  - [DIAMOND](https://github.com/bbuchfink/diamond) (`diamond` in `diamond_path`) – required when `use_prot = TRUE` for protein-based RBH
  - [Miniprot](https://github.com/lh3/miniprot) (optional; only if `run_miniprot = TRUE` to predict missing ORFs) – requires genome FASTA files in `orgInputFiles()`

### Installation

Install the development version from GitHub:

```r
install.packages("devtools")
devtools::install_github("tomoyukif/OrthoPaiR", dependencies = TRUE)
```

### Usage (simplified)

The main workflow is:

1. Prepare a list of input files for each genome with `orgInputFiles()`.
2. Run the full pipeline with `orthopair()`.
3. Use the generated HDF5, GFF, graph, and CSV files for downstream analyses.

```r
library(OrthoPaiR)

## 1. Prepare input definitions per genome
in_list <- orgInputFiles(
  name   = "Osat",
  genome = NA,                                # optional; required if run_miniprot = TRUE
  gff    = "path/to/osat.gff",
  cds    = "path/to/osat_cds.fa",
  prot   = "path/to/osat_prot.fa"            # optional; NA if unused
)

in_list <- orgInputFiles(
  object = in_list,
  name   = "Hvul",
  genome = NA,
  gff    = "path/to/hvul.gff",
  cds    = "path/to/hvul_cds.fa",
  prot   = "path/to/hvul_prot.fa"
)

## 2. Run the OrthoPaiR pipeline
wd <- "path/to/output/orthopair"

out_files <- orthopair(
  in_list      = in_list,
  working_dir  = wd,
  miniprot_path = "/path/to/miniprot/bin",
  blast_path    = "/path/to/blast/bin",
  diamond_path  = "/path/to/diamond/bin",
  n_threads     = 8,
  overwrite     = TRUE,
  verbose       = TRUE,
  use_prot      = FALSE,                      # TRUE: use DIAMOND for protein-based RBH
  target_pair   = NULL,
  orthopair     = TRUE,
  run_miniprot  = FALSE,                      # TRUE: predict missing ORFs (requires genome FASTA)
  reorg         = TRUE,
  makegraph     = TRUE,
  output_table  = TRUE
)

out_files
```

- For **two genomes**, you typically do **not** need multi‑genome summarisation, so you can leave:

  - `reorg = FALSE`, `makegraph = FALSE`, `output_table = FALSE`

  and work directly from the pairwise HDF5 output in `hdf5_out/`.

- For **more than two genomes**, setting:

  - `reorg = TRUE`, `makegraph = TRUE`, `output_table = TRUE`

  tells OrthoPaiR to:

  - reorganise gene IDs and write `<name>_orthopair.gff` files,
  - build an orthology graph (`orthopair.graphml`),
  - and export ortholog and orphan tables (`orthopair_list.csv`, `orphan_list.csv`).

**Note on parameters:**
- If `run_miniprot = TRUE` (to reciprocally predict missing ORFs), you must specify the `genome` parameter (path to genome FASTA) for each genome in `orgInputFiles()`.
- If `use_prot = TRUE`, `orthopair()` will use DIAMOND to compute reciprocal best hits based on protein sequence similarity instead of (or in addition to) DNA sequence similarity.

A complete four‑genome example is provided in `sample_data/sample_script.R`, and the corresponding output directory structure is shown in `sample_data/orthopair/`.

### Documentation

For a detailed description of input requirements, pipeline steps, and all output files, see the vignette:

```r
vignette("orthopair", package = "OrthoPaiR")
```

### License

OrthoPaiR is licensed under the **GPL‑3.0** License.  
See `LICENSE.md` for details.

