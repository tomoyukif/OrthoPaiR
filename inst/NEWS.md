Changes in version 0.6.2 (2025-02-13)
+ Post syntenic orthologous gene pairing step was largely revised to improve speed and ram usage.

Changes in version 0.5.12 (2025-01-07)
+ Add a function to draw a riparian plot representing syntenic blocks.
+ Add columns indicating whether pairs are anchors or not in the output data

Changes in version 0.5.9 (2025-12-01)
+ Bug fix in .getRBH()

Changes in version 0.5.5 (2025-11-21)
+ Bug fix in mapProt() and graph2df()

Changes in version 0.5.5 (2025-10-21)
+ Update the handling of input files
+ Update the output directory of compareOrthoSeq().
+ Add code to remove AA sequences after the terminal codon in each AA sequence 
+ in compareOrthoSeq().

Changes in version 0.5.2 (2025-06-16)
+ Anchor genes is filtered based on if the genes are possibly TEs that have 
+ a large number of multiple hits in RBH.
+ Orthologous pairs are filtered based on mutual-CI if it is detected as
+ outliers in the distribution of the mutual-CIs observed in the detected 
+ orthologous pairs.
+ A function to identify syntenic blocks is added.

Changes in version 0.4.15 (2025-02-20)
+ Change the algorithm to set SOGs.
+ Change the arugment name "redo" to "module".

Changes in version 0.4.14 (2025-02-19)
+ Replace BLAST to DIAMOND for the protein BLAST search in rbh().

Changes in version 0.4.11 (2025-02-10)
+ Minor bug fix in the no-reorg mode of reorgOrthopair().

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