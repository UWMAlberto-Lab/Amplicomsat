\name{alleleFinder}
\alias{alleleFinder}
\title{
Scores the diploid genotypes into an output table
}
\description{
Estimates the diploid genotypes based on read size (primer to primer) in bp. The results are written to two files: one with two alleles
 per column and the other with the two alleles collapsed in a single column. The files are written to a folder named GenotypesOut
created by the function. 
}
\usage{
alleleFinder(sampleDir="./samplesF",locusInfo="primers",coverageF=25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampleDir}{
The directory where the filtered reads are stored. This directory should only contain the read files and 
nothing else. By default this is a directory named samplesF written in the session working directory, 
as produced by \code{\link{read.quality.filter}}
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
\item{coverageF}{
Reads scored as alleles based on the sum of squares minimization method need to have coverage larger than this value, otherwise are written 
as missing data NA. Defaults to coverage of 25. It is possible to send a vector with different values for each locus. If you want to use the same coverage
filter for all loci you can input a single value.
  }
}
\details{
Uses the tables created by \code{\link{GenotypeModel}}, containing expected relative frequencies of read cover for different combinations (2 out of N) of alleles,
to contrast with observed data in each sample and find the diploid genotypes minimizing a sum of squares. For each observed sample, the row table that minimizes 
the sum of squares: sum((ORF(a)-ERF(a))^2), with ORF and ERF the Observed and Expected relative frequencies at each read size, is the estimated genotype.

}
\value{
A data frame with a concatenated genotype (the two alleles in the same cell) is returned. This table is also written to file, as well as a version with 
the two alleles in one locus in different columns (two columns per locus). The output files are written to a folder named GenotypesOut that is created by 
the function in the R session working directory.
}
\references{
Alberto F (in prep)...
}
\author{
F. Alberto
}
\seealso{
\code{\link{GenotypeModel}}, \code{\link{AlleleCum}}
}
\keyword{ IO }
