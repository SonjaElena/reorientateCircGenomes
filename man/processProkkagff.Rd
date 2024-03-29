% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{processProkkagff}
\alias{processProkkagff}
\title{Processing of Prokka supplied gff file type}
\usage{
processProkkagff(x)
}
\arguments{
\item{x}{Path to the unprocessed gff file generated from Prokka.}
}
\value{
The processed gff file as data.frame.
}
\description{
Process the gff file created by Prokka into a data.frame.
}
\examples{
gff_unprocessed <- path_to_gff
gff <- processProkkagff(gff_unprocessed);
}
