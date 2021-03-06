% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BasicFunctions.R
\name{run_ADARR_partition}
\alias{run_ADARR_partition}
\title{A function for annotating the Sprint results with genomic regions and repeats (allows for partition)}
\usage{
run_ADARR_partition(
  input.file.path,
  annotation.granges.path,
  complete.annotated.path,
  partition = 10000
)
}
\arguments{
\item{input.file.path}{The input file path}

\item{annotation.granges.path}{The path to output initial annotations}

\item{complete.annotated.path}{The path to output the complete annotations}

\item{partition}{the number of RES to partition by}
}
\description{
A function for annotating the Sprint results with genomic regions and repeats (allows for partition)
}
