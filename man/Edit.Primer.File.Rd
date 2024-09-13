\name{Edit.Primer.File}
\alias{Edit.Primer.File}
\title{
Checks and corrects for mismatches in primer sequences.
}
\description{
Sequencing reads might show primer sequences that diverge from published sequences for various reasons. This function will search the reads for primer sequences and correct those that diverge from the supplied parameter file.
}
\usage{
Edit.Primer.File(sampleDir="./samplesF",locusInfo="primers",calibN=2,mismatch=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{sampleDir}{
The directory where the filtered reads are stored. This directory should only contain the read files and nothing else. By default this is a directory
named samplesF written in the session working directory, as produced by \code{\link{read.quality.filter}} 
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
 For example, for a motif CAGG, you might need to use CCTG.\cr
\cr
nmotif: Reads are retained if they have at least these many repetitions of the repeat motif.
}
\item{calibN}{
A numeric value (defaults to 2) or a character vector with sample names. The function will search primers in these many random samples or specific samples supplied in a vector of sample names. The default is two random samples. If using specific samples, we recommend using two or three samples to avoid slowing down the code too much.
}
\item{mismatch}{
A numeric value for the maximum mismatch allowed between supplied primer sequences and those searched in sequenced reads. It defaults to 1.
}
}
\details{
All forward and reverse primers supplied in locusInfo will be searched in the sequenced reads, allowing for mismatches of up to one base pair (default). If differences are found, an edited file with primerInfo is written to the R session working directory. A log file reporting the changes made to the file is also written. The function uses matchPattern from package Biostrings to perform the search with mismatches. Using this function at the start of the pipeline allows for faster read filtering in subsequent steps by searching exact sequences without mismatches.
}
\value{
The function returns a table with the edited primer information in a similar format to argument locusInfo (see above). The table is also written in the R session working directory, where it can be conveniently parsed to the other functions in the Amplicomsat pipeline. A log reporting the changes made to the primer sequences is also written into a file.
}
\references{
Alberto F (in prep)...
}
\author{
F. Alberto
}
\seealso{
\code{\link{read.size.freqs}}
}
\keyword{ IO }