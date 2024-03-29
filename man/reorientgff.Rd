% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{reorientgff}
\alias{reorientgff}
\title{Reorientation of a gff file type}
\usage{
reorientgff(
  x,
  proteinID = NA,
  bplocation = bp_location,
  replicon = NA,
  Rep_size = fasta
)
}
\arguments{
\item{x}{gff file, processed with functions processProkkagff or processNCBIgff since the additional column called 'alternend' is used.}

\item{proteinID}{Supplies the start position at which file should be reoriented; Defaults to NA. If this ProteinID is not found on the biggest replicon, the replicon must be supplied as well.}

\item{bplocation}{Supplies the start position at which file should be reoriented in base pairs and needs to be identical with the start position of one of the proteins.}

\item{replicon}{Replicon to be reoriented. This option defaults to the largest replicon.}

\item{Rep_size}{Indicating the size of the Replicon to be used in base pairs. Alternatively, a genomic fasta sequence in fna format can be supplied.}
}
\value{
The reoriented gff file with three additional columns called Ostart and Oend containing the adjusted start and end base pair locations as well as a column called 'Oaltend' containing the reoriented end position based on the alternative end column, in case the function 'processNCBIgff' has been used.
}
\description{
Reorientation of the start position of a gff file type based on ProteinID or base pair position. New start and end locations are added in two additional columns called 'Ostart' and 'Oend'.
}
\examples{
gff <- reorientgff(gff, "WP_012176686.1")
gff <- reorientgff(gff, bplocation = 1866, replicon = "CP000031.2")
gff <- reorientgff(gff, proteinID = "AAV97145.1", replicon = "CP000032.1")
gff <- reorientgff(gff, bplocation = 0, Rep_size = fna_path)
}
