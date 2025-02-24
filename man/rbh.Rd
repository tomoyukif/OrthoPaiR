% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04_functions_blast.R
\name{rbh}
\alias{rbh}
\title{Execute Reciprocal BLAST Hits (RBH) search}
\usage{
rbh(
  object,
  db1 = NULL,
  db2 = NULL,
  use_prot = FALSE,
  n_threads = 1,
  n_batch = NULL,
  makedb = TRUE,
  max_target_seqs = 1e+05,
  pident = 0,
  qcovs = 0,
  evalue = 1e-04,
  diamond_exec_path = NULL,
  diamond_out_dir = ""
)
}
\arguments{
\item{object}{A OrthoPairDB object.}

\item{db1}{Path to the BLAST database for query genome (default is NULL).}

\item{db2}{Path to the BLAST database for subject genome (default is NULL).}

\item{n_threads}{Number of threads to use for BLAST (default is 1).}

\item{n_batch}{Number of sequences to process in each batch (default is NULL).}

\item{makedb}{Logical indicating whether to create BLAST databases (default is TRUE).}

\item{max_target_seqs}{Maximum number of target sequences (default is 100000).}

\item{pident}{Minimum percentage identity for BLAST hits (default is 90).}

\item{qcovs}{Minimum query coverage for BLAST hits (default is 0).}

\item{evalue}{E-value threshold for BLAST (default is 1e-10).}
}
\description{
This function performs a Reciprocal BLAST Hits (RBH) search between query and subject genomes.
}
