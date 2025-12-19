library(rhdf5)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(scales)
library(ggforce)


hdf5_fn <- "~/workspace/orthology/output/benchmark/orthopair/hdf5_out/Osat_Ogla.h5"
hdf5_fn <- "~/workspace/orthology/output/benchmark_plant/orthopair/reorg_out/reorg_orthopair.h5"
genomes <- c("Osat", "Ogla")
min_block_genes = 5L
genome_names = NULL
ribbon_alpha = 0.4
ribbon_width = 0.45
y_gap = 1
curve_npts = 250L
curve_keepat = round(curve_npts / 20)
chr_lwd = 2
chr_palette = function(n) hcl.colors(n, "Dark 3")
ribbon_palette = function(n) hcl.colors(n, "Zissou 1")

devtools::load_all("./")
op <- getOrthoPair(hdf5_fn = hdf5_fn, score = FALSE, loc = TRUE)
meta <- getMeta(hdf5_fn = hdf5_fn)

# Prepare data for plotting

p <- plot_riparian(hdf5_fn = hdf5_fn)

ggsave("~/workspace/orthology/output/benchmark_plant/synteny_stacked.pdf",
       p, width = 12, height = 12)
