library(OrthoPaiR)
hdf5_fn <- "~/hdd5/ngsData/assemble/og_pt/orthopair/pairwise_orthopair/hdf5_out/reorg_orthopair.h5"

makeOrthoGraph(hdf5_fn = hdf5_fn)

x <- grp[[4]]

# h5 <- H5Fopen("~/hdd5/ngsData/assemble/og_pt/orthopair/pairwise_orthopair/hdf5_out/WK21_C7530.h5")
# a <- h5$orthopair_gene
nb <- import.gff3("~/hdd5/ngsData/assemble/og_pt/orthopair/pairwise_orthopair/hdf5_out/NB_orthopair.gff")
nb <- nb$ID[nb$type == "gene"]
wk21 <- import.gff3("~/hdd5/ngsData/assemble/og_pt/orthopair/pairwise_orthopair/hdf5_out/WK21_orthopair.gff")
wk21 <- wk21$ID[wk21$type == "gene"]
C7530 <- import.gff3("~/hdd5/ngsData/assemble/og_pt/orthopair/pairwise_orthopair/hdf5_out/C7530_orthopair.gff")
C7530 <- C7530$ID[C7530$type == "gene"]
c8547 <- import.gff3("~/hdd5/ngsData/assemble/og_pt/orthopair/pairwise_orthopair/hdf5_out/C8547_orthopair.gff")
c8547 <- c8547$ID[c8547$type == "gene"]
gene_list <- rbind(data.frame(genome = "NB", gene = nb), 
                   data.frame(genome = "WK21", gene = wk21), 
                   data.frame(genome = "C7530", gene = C7530), 
                   data.frame(genome = "C8547", gene = c8547))
"C8547_C7530_MP000006"
#! Change the MP gene renaming rule
# The prefix for the genome + MP-IDs
146               OgC7530_01G00003600  32
147               OgC8547_01G00003300  32
148                    LOC_Os01g01320  32
149                    LOC_Os01g01330  32
150          C7530_WK21_MP000054:gene  32
151  OgWK21_01T00004300.02:split_gene  32
152  OgWK21_01T00004300.03:split_gene  32
153          C8547_WK21_MP000054:gene  32
154                OgWK21_01G00004300  32

"OgWK21_01T00004300.02:split_gene"
#################################################################################

# > mp_gff[mp_gff$ID == "NB_MP001262:gene"]
# GRanges object with 1 range and 10 metadata columns:
#     seqnames          ranges strand |   source     type     score     phase               ID
# <Rle>       <IRanges>  <Rle> | <factor> <factor> <numeric> <integer>      <character>
#     [1]    chr01 2823257-2823763      - | miniprot     gene       277      <NA> NB_MP001262:gene
# Name          gene_id          Parent                Target         oldID
# <character>      <character> <CharacterList>           <character>   <character>
#     [1] NB_MP001262:gene NB_MP001262:gene                 LOC_Os01g06910.1 1 57 MP001262:gene
# -------
#     seqinfo: 12 sequences from an unspecified genome; no seqlengths

in_list <- orgInputFiles(name = "NB",
                         genome = "~/hdd3/genomeData/rice/cultivar_sativa/os_nb_msu7/all.con",
                         gff = "~/hdd3/genomeData/rice/cultivar_sativa/nb_combined/version_2023/gff/nb_combined_all.gff",
                         cds = "~/hdd3/genomeData/rice/cultivar_sativa/nb_combined/version_2023/fasta/nb_combined_all_cds.fa",
                         prot = "~/hdd3/genomeData/rice/cultivar_sativa/nb_combined/version_2023/fasta/nb_combined_all_prot.fa")
in_list <- orgInputFiles(object = in_list,
                         name = "WK21",
                         genome = "/home/ftom/hdd2/og_gwas/genomics/01_genomeData/final_genome/wk21_genome_chr.fa",
                         gff = "/home/ftom/hdd2/og_gwas/genomics/03_reorganizeGenes/wk21_genome_reorganized.gff",
                         cds = "/home/ftom/hdd2/og_gwas/genomics/03_reorganizeGenes/wk21_reorganized_cds.fa",
                         prot = "/home/ftom/hdd2/og_gwas/genomics/03_reorganizeGenes/wk21_reorganized_prot.fa")
in_list <- orgInputFiles(object = in_list,
                         name = "C7530",
                         genome = "/home/ftom/hdd5/ngsData/assemble/og_pt/final_output/C7530_genome.fa.gz",
                         gff = "/home/ftom/hdd5/ngsData/assemble/og_pt/final_output/C7530_genome.gff",
                         cds = "/home/ftom/hdd5/ngsData/assemble/og_pt/final_output/C7530_cds.fa",
                         prot = "/home/ftom/hdd5/ngsData/assemble/og_pt/final_output/C7530_prot.fa")
in_list <- orgInputFiles(object = in_list,
                         name = "C8547",
                         genome = "/home/ftom/hdd5/ngsData/assemble/og_pt/final_output/C8547_genome.fa.gz",
                         gff = "/home/ftom/hdd5/ngsData/assemble/og_pt/final_output/C8547_genome.gff",
                         cds = "/home/ftom/hdd5/ngsData/assemble/og_pt/final_output/C8547_cds.fa",
                         prot = "/home/ftom/hdd5/ngsData/assemble/og_pt/final_output/C8547_prot.fa")

setwd("~/hdd5/ngsData/assemble/og_pt")
in_list = in_list
hdf5_out_dir = "synog/pairwise_synog/hdf5_out"
sibeliaz_out_dir = "synog/pairwise_synog/sibeliaz_out"
sibeliaz_bin = "sibeliaz"
maf2synteny_bin = "maf2synteny"
conda = "/home/ftom/miniforge/miniforge3/bin/conda"
sibeliaz_condaenv = "sibeliaz"
miniprot_bin = "miniprot"
miniprot_condaenv = "miniprot"
miniprot_out_dir = "synog/pairwise_synog/miniprot_out"
n_threads = 30
overwrite = TRUE







known_1to1$pair_id[!known_1to1$pair_id %in% synog$gene_pair_id]
# [1] "Os12g0212100_OgWK21_12G007590" "Os12g0567200_OgWK21_12G018890"

#
# seqnames          ranges strand |   source       type     score     phase                     ID
# <Rle>       <IRanges>  <Rle> | <factor>   <factor> <numeric> <integer>            <character>
#     [1]    chr08 5635831-5638392      + |     RGIR transcript      1956      <NA>  OgWK21_08T00085900.01
# [2]    chr08 5635965-5721409      + |     RGIR transcript      5548      <NA>  OgWK21_08T00085900.02
# [3]    chr08 5636066-5721409      + |     RGIR transcript     10468      <NA>  OgWK21_08T00085900.03
# [4]    chr08 5643015-5646838      + |     RGIR transcript      2147      <NA>  OgWK21_08T00085900.04
# [5]    chr08 5643015-5645015      + |     RGIR transcript       987      <NA>  OgWK21_08T00085900.05
# ...      ...             ...    ... .      ...        ...       ...       ...                    ...
# [295]    chr08 6217987-6220747      + |     RGIR transcript      2236      <NA> OgWK21_08T00085900.295
# [296]    chr08 6217987-6219904      + |     RGIR transcript      1594      <NA> OgWK21_08T00085900.296
# [297]    chr08 6230553-6230844      + |     RGIR transcript       381      <NA> OgWK21_08T00085900.297
# [298]    chr08 6237474-6239955      + |     RGIR transcript      2831      <NA> OgWK21_08T00085900.298
# [299]    chr08 6239182-6239958      + |     RGIR tran

setwd("/home/ftom/hdd2/og_gwas/genomics/02_orthology/07_synog")

"~/hdd2/og_gwas/genomics/02_orthology/07_synog/output/nb_wk21.h5"

query_genome = "input/nb_genome.fa"
subject_genome = "input/wk21_genome.fa"
query_gff = "input/nb.gff"
subject_gff = "input/wk21.gff"
query_cds = "input/nb_cds.fa"
subject_cds = "input/wk21_cds.fa"
query_prot = "input/nb_prot.fa"
subject_prot = "input/wk21_prot.fa"
hdf5_path = "output/nb_wk21.h5"
sibeliaz_out_dir = "output/sibeliaz_out"
sibeliaz_bin = "sibeliaz"
maf2synteny_bin = "maf2synteny"
conda = "/home/ftom/miniforge/miniforge3/bin/conda"
sibeliaz_condaenv = "sibeliaz"
miniprot_bin = "miniprot"
miniprot_condaenv = "miniprot"
miniprot_out_dir = "output/miniprot_out"
n_threads = 30
omit_chr = ""

runSynog(query_genome = "input/nb_genome.fa",
         subject_genome = "input/wk21_genome.fa",
         query_gff = "input/nb.gff",
         subject_gff = "input/wk21.gff",
         query_cds = "input/nb_cds.fa",
         subject_cds = "input/wk21_cds.fa",
         query_prot = "input/nb_prot.fa",
         subject_prot = "input/wk21_prot.fa",
         hdf5_path = "output/nb_wk21.h5",
         sibeliaz_out_dir = "output/sibeliaz_out",
         sibeliaz_bin = "sibeliaz",
         maf2synteny_bin = "maf2synteny",
         conda = "/home/ftom/miniforge/miniforge3/bin/conda",
         sibeliaz_condaenv = "sibeliaz",
         miniprot_bin = "miniprot",
         miniprot_condaenv = "miniprot",
         miniprot_out_dir = "output/miniprot_out",
         n_threads = 30,
         omit_chr = "")

object <- makeSynogDB(query_genome = "input/nb_genome.fa",
                      subject_genome = "input/wk21_genome.fa",
                      query_gff = "input/nb.gff",
                      subject_gff = "input/wk21.gff",
                      query_cds = "input/nb_cds.fa",
                      subject_cds = "input/wk21_cds.fa",
                      query_prot = "input/nb_prot.fa",
                      subject_prot = "input/wk21_prot.fa",
                      hdf5_path = "output/nb_wk21.h5")

subject_prot = object$subject_prot
query_prot = object$query_prot
miniprot_bin = "miniprot"
conda = "/home/ftom/miniforge/miniforge3/bin/conda"
condaenv = "miniprot"
n_threads = 30
out_dir = "output/miniprot_out"

###############################################################################
setwd("~/hdd5/ngsData/assemble/og_pt")
hdf5_fn <- c(WK21_C7530 = "~/hdd5/ngsData/assemble/og_pt/orthopair/pairwise_orthopair/hdf5_out/WK21_C7530.h5",
             WK21_C8547 = "~/hdd5/ngsData/assemble/og_pt/orthopair/pairwise_orthopair/hdf5_out/WK21_C8547.h5",
             C7530_C8547 = "~/hdd5/ngsData/assemble/og_pt/orthopair/pairwise_orthopair/hdf5_out/C7530_C8547.h5",
             NB_C7530 = "~/hdd5/ngsData/assemble/og_pt/orthopair/pairwise_orthopair/hdf5_out/NB_C7530.h5",
             NB_C8547 = "~/hdd5/ngsData/assemble/og_pt/orthopair/pairwise_orthopair/hdf5_out/NB_C8547.h5")
