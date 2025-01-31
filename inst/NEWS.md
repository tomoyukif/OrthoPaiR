Changes in version 0.4.9 (2025-01-31)
+ Bug fix to accept no blast hit.

Changes in version 0.4.7 (2025-01-30)
+ Update .prepPairs() to create symlinks to input files and put all symlink files
+ into the input directory of the working directory. 

Changes in version 0.4.6 (2025-01-29)
+ Update fixInFiles() to support automatic transcript rename. 

Changes in version 0.4.3 (2025-01-07)
+ Minor bug fix in .orgOrphan(). and the main pipeline control. 

Changes in version 0.4.2 (2024-12-09)
+ Minor bug fix in syntenicOrtho() to skip .splitGene() if the first run of 
+ .splitGene() returned NULL to the 'rest' slot in the output object.

Changes in version 0.4.1 (2024-11-15)
+ Remove the LCB detection step by SibeliaZ and change the algorithm to let the anchoring step only depends on the RBBH result with filtering by the 10%-quantile of mutual-CI.

Changes in version 0.3.7 (2024-11-1)
+ Add function to fix genome and GFF files to meet the demand of OrthoPaiR.

Changes in version 0.3.5 (2024-10-30)
+ Add function to validate whether sequence levels, CDS, and protein names match
+ IDs in GFF.
+ Add function to validate entry types and gene_id in GFF.
+ Add function to resume the process.

Changes in version 0.3.2 (2024-9-30)
+ Replace wrapper functions with a single wrapper function that covers all 
+ functions that had been provided by the removed wrapper functions.
+ Minor bug fixes in the core algorithm.
+ Add a function to conduct the multiple sequence alignment.

Changes in version 0.2.2 (2024-8-2)
+ Add a wrapper function to run OrthoPaiR over multiple combinations of genomes.
+ Add a function to reorganized GFF files.
+ Add a function to create a graph object holding complex ortholog pairing information.
+ Add a function to convert the complex ortholog pairing information to a dataframe.

Changes in version 0.1.2 (2024-6-20)
+ Add a wrapper function for running the OrthoPaiR pipeline.
+ Minor bug fixes in the core algorithm.

Changes in version 0.1.1 (2024-6-6)
+ Beta release of OrthoPaiR.
+ Vignette was prepared.

Changes in version 0.1.0 (2024-2-15)
+ Initial release of OrthoPaiR.