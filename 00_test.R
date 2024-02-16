object <- makeSynogDB(query_genome = "input/nb_genome.fa",
                      subject_genome = "input/wk21_genome.fa",
                      query_gff = "input/nb.gff",
                      subject_gff = "input/wk21.gff",
                      query_cds = "input/nb_cds.fa",
                      subject_cds = "input/wk21_cds.fa",
                      query_prot = "input/nb_prot.fa",
                      subject_prot = "input/wk21_prot.fa",
                      positive_list = "input/positive_list.csv",
                      negative_list = "input/negative_list.csv",
                      hdf5_path = "output/synog.h5")

################
# Run Sibeliaz #
################
runSibeliaz(object = object,
            out_dir = "sibeliaz_out",
            conda = "/home/ftom/tools/anaconda3/bin/conda",
            condaenv = "sibeliaz", run_sibeliaz = TRUE)
sibeliaLCB2DF(object = object)
lcbClassify(object = object)
getLCBpairs(object = object)

# Find ortholog pairs based on CDS RBH
rbh(object = object, n_threads = 30)

# Run LCB-based RBBH filtering
anchorOrtho(object = object)

# Synteny-based RBH filtering
syntenyOrtho(object = object, omit_chr = "chrUn|chrSy",
             pident = 90, evalue = 1e-100, qcovs = 50)

geneOrtho(object = object)

txwise_summary <- summarySynog(object = object, gene = FALSE)
genewise_summary <- summarySynog(object = object, gene = TRUE)
synog_genewise <- getSynog(object = object, gene = TRUE)
orphan_genewise <- getOrphan(object = object, gene = TRUE)
write.csv(txwise_summary, "output/synog_nonMiniprot_txwise_summary.csv")
write.csv(genewise_summary, "output/synog_nonMiniprot_genewise_summary.csv")
write.csv(synog_genewise, "output/synog_nonMiniprot_genewise_ortho.csv")
write.csv(orphan_genewise$query, "output/synog_nonMiniprot_orphan_query.csv", row.names = FALSE)
write.csv(orphan_genewise$subject, "output/synog_nonMiniprot_orphan_subject.csv", row.names = FALSE)


miniprot_bin <- "/home/ftom/tools/miniprot/miniprot"
mapProt(object = object,
        miniprot_bin = miniprot_bin,
        n_core = 20,
        out_prefix = "output/nb_wk21_")

createFASTA(object = object)

object <- updateFiles(object = object)
rbh(object = object, n_threads = 30)
anchorOrtho(object = object)
syntenyOrtho(object = object, omit_chr = "chrUn|chrSy",
             pident = 90, evalue = 1e-100, qcovs = 50)

geneOrtho(object = object)

txwise_summary <- summarySynog(object = object, gene = FALSE)
genewise_summary <- summarySynog(object = object, gene = TRUE)
synog_genewise <- getSynog(object = object, gene = TRUE)
orphan_genewise <- getOrphan(object = object, gene = TRUE)
write.csv(txwise_summary, "output/synog_txwise_summary.csv")
write.csv(genewise_summary, "output/synog_genewise_summary.csv")
write.csv(synog_genewise, "output/synog_genewise_ortho.csv")
write.csv(orphan_genewise$query, "output/synog_orphan_query.csv", row.names = FALSE)
write.csv(orphan_genewise$subject, "output/synog_orphan_subject.csv", row.names = FALSE)


splitGenes(object = object)
genewise_split_summary <- summarySynog(object = object, gene = TRUE, split = TRUE)
synog_txwise <- getSynog(object = object, gene = FALSE)
synog_genewise_split <- getSynog(object = object, gene = TRUE, split = TRUE)
orphan_genewise_split <- getOrphan(object = object, gene = TRUE, split = TRUE)
write.csv(genewise_split_summary, "output/synog_genewise_split_summary.csv")
write.csv(synog_genewise_split, "output/synog_genewise_split_ortho.csv")
write.csv(synog_txwise, "output/synog_txwise_ortho.csv")
write.csv(orphan_genewise_split$query, "output/synog_orphan_split_query.csv", row.names = FALSE)
write.csv(orphan_genewise_split$subject, "output/synog_orphan_split_subject.csv", row.names = FALSE)
