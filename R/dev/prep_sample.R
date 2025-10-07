library(Biostrings)
genome <- readDNAStringSet("~/hdd2/og_gwas/genomics/02_orthology/07_synog/input/nb_genome.fa")
writeXStringSet(genome[1], "inst/extdata/query_genome.fa")

genome <- readDNAStringSet("~/hdd2/og_gwas/genomics/02_orthology/07_synog/input/wk21_genome.fa")
writeXStringSet(genome[1], "inst/extdata/subject_genome.fa")

library(rtracklayer)
gff <- import.gff3("~/hdd2/og_gwas/genomics/02_orthology/07_synog/input/nb.gff")
gff <- gff[as.character(seqnames(gff)) %in% "chr01"]
export.gff3(gff, "inst/extdata/query.gff")
cds <- readDNAStringSet("~/hdd2/og_gwas/genomics/02_orthology/07_synog/input/nb_cds.fa")
cds <- cds[names(cds) %in% gff$ID]
writeXStringSet(cds, "inst/extdata/query_cds.fa")
prot <- readAAStringSet("~/hdd2/og_gwas/genomics/02_orthology/07_synog/input/nb_prot.fa")
prot <- prot[names(prot) %in% gff$ID]
writeXStringSet(prot, "inst/extdata/query_prot.fa")

gff <- import.gff3("~/hdd2/og_gwas/genomics/02_orthology/07_synog/input/wk21.gff")
gff <- gff[as.character(seqnames(gff)) %in% "chr01"]
export.gff3(gff, "inst/extdata/subject.gff")
cds <- readDNAStringSet("~/hdd2/og_gwas/genomics/02_orthology/07_synog/input/wk21_cds.fa")
cds <- cds[names(cds) %in% gff$ID]
writeXStringSet(cds, "inst/extdata/subject_cds.fa")
prot <- readAAStringSet("~/hdd2/og_gwas/genomics/02_orthology/07_synog/input/wk21_prot.fa")
prot <- prot[names(prot) %in% gff$ID]
writeXStringSet(prot, "inst/extdata/subject_prot.fa")
