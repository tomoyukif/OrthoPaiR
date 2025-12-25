library(rhdf5)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(scales)
library(ggforce)

source("R/dev/10_functions_plot.R")

# Prepare data for plotting
p <- plotRiparian(hdf5_fn = "~/workspace/orthology/output/benchmark/orthopair/hdf5_out/Osat_Ogla.h5",
                  genomes = c("Osat", "Ogla"), 
                  chr_sizes = NULL,
                  chr_gap = 3,
                  track_gap = 0.5,
                  max_links = 50000,
                  seed = 1,
                  alpha_base = NULL,
                  linewidth = 0.45,
                  inv_darken = 0.35,
                  inv_alpha_mult = 1.6,
                  palette = c("#d73027","#fc8d59","#fee08b","#91bfdb","#4575b4"),
                  min_block_width = 1e4,
                  chr_bar_lw = 8,
                  chr_label_size = 15,
                  normalize_genome = TRUE,
                  genome_scale = 100)

ggsave("~/workspace/orthology/output/benchmark/synteny_stacked.pdf",
       p, width = 8, height = 4)



# Prepare data for plotting
p <- plotRiparian(hdf5_fn = "~/workspace/orthology/output/benchmark_plant/orthopair/reorg_out/reorg_orthopair.h5",
                  genomes = c("Osat", "Hvul", "Atha", "Mpol"),
                  select_chr = list(Osat = c(3, 11, 12), Hvul = c(4)),
                  chr_sizes = NULL,
                  chr_gap = 3,
                  track_gap = 0.5,
                  max_links = 50000,
                  seed = 1,
                  alpha_base = NULL,
                  linewidth = 0.45,
                  inv_darken = 0.35,
                  inv_alpha_mult = 1.6,
                  palette = c("#d73027","#fc8d59","#fee08b","#91bfdb","#4575b4"),
                  min_block_width = 1e4,
                  chr_bar_lw = 8,
                  chr_label_size = 15,
                  normalize_genome = TRUE,
                  genome_scale = 200)

ggsave("~/workspace/orthology/output/benchmark_plant/synteny_stacked.pdf",
       p, width = 8, height = 8)
