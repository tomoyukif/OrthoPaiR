# Synog

[![R](https://img.shields.io/badge/R-4.1.0-blue)](https://www.r-project.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**Synog** is an R package for conducting syntenic orthologous gene pairing between given genomes. It provides tools to identify, filter, and analyze orthologous gene pairs, facilitating comparative genomic studies.

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

To install the Synog package, use the following commands in R:

```r
# Install the devtools package if not already installed
install.packages("devtools")

# Install Synog from GitHub
devtools::install_github("tomoyukif/Synog")
```

## Prerequisites

Synog relies on several R packages and external tools. Install the prerequisite R packages using:

```r
install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "GenomeInfoDb", "GenomicFeatures", 
"rhdf5", txdbmaker, "BiocGenerics", "Biobase", "IRanges",
"rtracklayer", "Biostrings", "S4Vectors", "BSgenome"))
install.packages(c("ggplot2", "rBLAST", "parallel", "dplyr"))
```

Ensure the following external tools are installed:

- **Miniprot**: [Miniprot Installation Guide](https://github.com/lh3/miniprot)
- **Sibeliaz**: [Sibeliaz Installation Guide](https://github.com/medvedevgroup/sibeliaz)

## Usage

Here is a basic example of how to use Synog for syntenic orthologous gene pairing.
You can go through the Synog pipeline using the `runSynog()` function, which is a wrapper that executes a series of functions to complete the Synog pipeline.

```r
library(Synog)

# Define file paths
query_genome <- "input/nb_genome.fa"
subject_genome <- "input/wk21_genome.fa"
query_gff <- "input/nb.gff"
subject_gff <- "input/wk21.gff"
query_cds <- "input/nb_cds.fa"
subject_cds <- "input/wk21_cds.fa"
query_prot <- "input/nb_prot.fa"
subject_prot <- "input/wk21_prot.fa"
hdf5_path <- "output/synog.h5"

# Run Synog pipeline
object <- runSynog(
  query_genome = query_genome,
  subject_genome = subject_genome,
  query_gff = query_gff,
  subject_gff = subject_gff,
  query_cds = query_cds,
  subject_cds = subject_cds,
  query_prot = query_prot,
  subject_prot = subject_prot,
  hdf5_path = hdf5_path,
  sibeliaz_out_dir = "./sibeliaz_out",
  sibeliaz_bin = "sibeliaz",
  maf2synteny_bin = "maf2synteny",
  conda = "/home/ftom/miniforge/miniforge3/bin/conda",
  sibeliaz_condaenv = "sibeliaz",
  miniprot_bin = "miniprot",
  miniprot_condaenv = "miniprot",
  miniprot_out_dir = "./miniprot_out",
  n_threads = 30,
  omit_chr = "chrUn|chrSy"
)
```

The `runSynog()` function outputs a `SynogDB` object containing a link to the HDF5 file that stores syntenic ortholog pairing results. You can access these results using the following functions:

```r
# Access results
genewise_summary <- summarySynog(object = object, gene = TRUE)
synog_genewise <- getSynog(object = object, gene = TRUE)
orphan_genewise <- getOrphan(object = object, gene = TRUE)
genewise_split_summary <- summarySynog(object = object, gene = TRUE, split = TRUE)
synog_genewise_split <- getSynog(object = object, gene = TRUE, split = TRUE)
orphan_genewise_split <- getOrphan(object = object, gene = TRUE, split = TRUE)
```

## Documentation

For a detailed guide and more advanced usage, refer to the [package vignette](vignettes/Synog.html).

## Contributing

Contributions are welcome! If you have any ideas, suggestions, or bug reports, please open an issue or submit a pull request.

## License

Synog is licensed under the GPL-3.0 License. See the [LICENSE](LICENSE.md) file for details.
