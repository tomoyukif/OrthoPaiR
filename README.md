# OrthoPaiR

[![R](https://img.shields.io/badge/R-4.1.0-blue)](https://www.r-project.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**OrthoPaiR** is an R package for conducting syntenic orthologous gene pairing between given genomes. It provides tools to identify, filter, and analyze orthologous gene pairs, facilitating comparative genomic studies.

## Table of Contents

- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [License](#license)

## Features

- Identify syntenic orthologous gene pairs between genomes.
- Filter orthologous pairs using LCB and RBH methods.
- Generate gene-wise and transcript-wise ortholog summaries.
- Map proteins using Miniprot and integrate with ortholog analysis.
- Support for handling orphan genes and split genes.

## Installation

To install the OrthoPaiR package, use the following commands in R:

```r
# Install the devtools package if not already installed
install.packages("devtools")

# Install OrthoPaiR from GitHub
devtools::install_github("tomoyukif/OrthoPaiR", dependencies = TRUE)
```

## Prerequisites

OrthoPaiR relies on several R packages and external tools. Install the prerequisite R packages using:

```r
install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "GenomeInfoDb", "GenomicFeatures", 
"rhdf5", txdbmaker, "BiocGenerics", "Biobase", "IRanges",
"rtracklayer", "Biostrings", "S4Vectors", "BSgenome"))
install.packages(c("ggplot2", "rBLAST", "parallel", "dplyr"))
```

Ensure the following external tools are installed:

- **Miniprot**: [Miniprot Installation Guide](https://github.com/lh3/miniprot)

## Usage

Here is a basic example of how to use OrthoPaiR for syntenic orthologous gene pairing.
You can go through the OrthoPaiR pipeline using the `runOrthoPaiR()` function, which is a wrapper that executes a series of functions to complete the OrthoPaiR pipeline.

```r
library(OrthoPaiR)

# Define file paths
query_genome <- "input/nb_genome.fa"
subject_genome <- "input/wk21_genome.fa"
query_gff <- "input/nb.gff"
subject_gff <- "input/wk21.gff"
query_cds <- "input/nb_cds.fa"
subject_cds <- "input/wk21_cds.fa"
query_prot <- "input/nb_prot.fa"
subject_prot <- "input/wk21_prot.fa"
hdf5_path <- "output/OrthoPaiR.h5"

# Run OrthoPaiR pipeline
object <- runOrthoPaiR(
  query_genome = query_genome,
  subject_genome = subject_genome,
  query_gff = query_gff,
  subject_gff = subject_gff,
  query_cds = query_cds,
  subject_cds = subject_cds,
  query_prot = query_prot,
  subject_prot = subject_prot,
  hdf5_path = hdf5_path,
  maf2synteny_bin = "maf2synteny",
  conda = "/home/ftom/miniforge/miniforge3/bin/conda",
  sibeliaz_condaenv = "sibeliaz",
  miniprot_bin = "miniprot",
  miniprot_condaenv = "miniprot",
  n_threads = 30
)
```

The `runOrthoPaiR()` function outputs a `OrthoPaiRDB` object containing a link to the HDF5 file that stores syntenic ortholog pairing results. You can access these results using the following functions:

```r
# Access results
genewise_summary <- summaryOrthoPaiR(object = object, gene = TRUE)
genewise_ortholog <- getOrthoPaiR(object = object, gene = TRUE)
orphan_genewise <- getOrphan(object = object, gene = TRUE)
```

## Documentation

For a detailed guide and more advanced usage, refer to the [package vignette](vignettes/orthopair.html).

## Contributing

Contributions are welcome! If you have any ideas, suggestions, or bug reports, please open an issue or submit a pull request.

## License

OrthoPaiR is licensed under the GPL-3.0 License. See the [LICENSE](LICENSE.md) file for details.
