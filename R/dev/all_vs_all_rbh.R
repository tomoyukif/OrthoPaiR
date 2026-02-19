# get all CDS files
devtools::load_all("/home/ftom/01_wd/softDevel/OrthoPaiR")
working_dir <- "/home/ftom/workspace/orthology/rice/output"
blast_path <- "/home/ftom/07_tools/bin"
diamond_path <- "/home/ftom/conda/envs/diamond/bin"

all_cds_fn <- list.files(file.path(working_dir, "input"),
                         "cds.fa$",
                         recursive = TRUE,
                         full.names = TRUE)

fn <- file.path(working_dir, "input", "all_cds.fa")

fa <- readDNAStringSet(all_cds_fn)
ids <- names(fa)
ufa <- unique(fa)
gid <- Biostrings::match(fa, ufa)
gsize <- tabulate(gid, nbins = length(u))
keep <- gsize[gid] > 1L
gid2 <- gid[keep]
ids2 <- ids[keep]
o <- order(gid2)
gid2 <- gid2[o]
ids2 <- ids2[o]
r <- rle(gid2)
starts <- cumsum(c(1L, head(r$lengths, -1L)))
rep_ids <- ids2[starts]
rep_col <- rep(rep_ids, times = r$lengths)
sel <- ids2 != rep_col
dup_dt <- data.table(rep = rep_col[sel], member = ids2[sel])
saveRDS(dup_dt, file.path(working_dir, "input", "dup_cds_table.rds"))

writeXStringSet(ufa, fn)

all_taxid_map_fn <- list.files(file.path(working_dir, "input"),
                               "taxid_map.tsv$",
                               recursive = TRUE,
                               full.names = TRUE)
taxid_map <- file.path(working_dir, "input", "all_taxid_map.tsv")
out_con <- file(taxid_map, "wb")
for (f in all_taxid_map_fn) {
    in_con <- file(f, "rb")
    while (length(buf <- readBin(in_con, "raw", 1024^2)) > 0) {
        writeBin(buf, out_con)
    }
    close(in_con)
}
close(out_con)

rm(fa, ufa, dup_dt); gc(); gc()

.makeblastdb(blast_path = blast_path, fn = fn, taxid_map = taxid_map)
out <- .blast_search(fa = fn, db = fn, 
                     blast_path = blast_path, 
                     n_threads = 30, 
                     stdout = FALSE)
# All vs all BLAST took 68 hours to finish.

Rcpp::sourceCpp("inst/src/reciprocal_allhsps.cpp")

blast_out_fn <- file.path(working_dir, "input", "blast_out.txt")
rbh_out_fn <- file.path(working_dir, "input", "rbh_out.txt")
tmp_dir <- "/home/ftom/workspace/orthology/input/tmp"
dir.create(tmp_dir)
reciprocal_hits_allhsps_cpp(
    infile  = blast_out_fn,
    outfile = rbh_out_fn,
    pident_min = 0.0,     # 必要なら 60 など
    qcov_min   = 0.0,     # 必要なら 50 など
    chunk_lines = 3000000,
    tmpdir = tmp_dir,
    filter_self = TRUE
)

################################################################################
ufa <- readDNAStringSet(fn)
prot <- translate(ufa, if.fuzzy.codon = "solve")
fn <- file.path(working_dir, "input", "all_prot.fa")
writeXStringSet(prot, fn)
.makediamonddb(diamond_path = diamond_path, fn = fn)

.diamond_search(fa = fn, db = fn, n_threads = 30, diamond_path = diamond_path)
# All vs all BLAST took  hours to finish.

