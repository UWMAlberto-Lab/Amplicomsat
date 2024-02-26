\name{AlleleCum}
\alias{AlleleCum}
\title{
Plots cummulative distributions of scored alleles
}
\description{
Plots cummulative distributions of scored alleles to help diagnose any possible errors
}
\usage{
AlleleCum(GenotypeA.DF="./GenotypesOut/Genotypes allele table.txt",locus,locusInfo="primers",
          ymin=NULL,ymax=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GenotypeA.DF}{
A data.frame as produced by \code{\link{alleleFinder}}
}
\item{locus}{
A string with the name of the locus to plot. The name needs to be written as in the file provided to locusInfo argument to
 other functions, e.g., \code{\link{read.size.freqs}}
}
\item{locusInfo}{
A tab delimited text file containing the information for the microsatellite loci used. Should contain four columns with the following header names:\cr
\cr
Locus: The loci names\cr
\cr
Fprimer: The sequence of the FWD primer in upper case\cr
\cr
Rprimer: The sequence of the RVS primer in upper case\cr
\cr
Motif: The sequence for the microsatellite motif (e.g., CTTG for tetranucleotide). Note that for some loci you might need to use the reverse
 complementar sequence of what your reference value is. Thus, inspect your reads first or if you are not getting data for some loci try this out.
 For example, for a motif CAGG, you might need to use CCTG.
}
 \item{ymin}{
The y lower lim for the y axis. Usefull for loci that are very polymorphic and need to be zoomed in.
}
 \item{ymax}{
The y upper lim for the y axis. Usefull for loci that are very polymorphic and need to be zoomed in.
}
}
\details{
Plots the cumulative distribution of alleles scored for the selected locus stored in the data.frame supplied to argument GenotypeA.DF. 
The sample names for specific observations can be obtained by interacting with the plot using function \code{\link{getpoints}}.
}
\value{
Plots the cumulative allele distribution for a selected locus
}
\references{
Alberto F (in prep)...
}
\author{
F. Alberto
}
\seealso{
\code{\link{getpoints}}
}
\keyword{ IO }
