\name{read.size.freqs}
\alias{read.size.freqs}
\title{
Finds reads that contain the microsatellite primers and four or more repeats of the microsatellite motif
}
\description{
Searches sequenced reads for those containing the microsatellite FWD and RVS primers and a user-defined number of microsatellite motif repeats.
 Writes files for each sample and loci combination containing the read size and the number of msat motif repeats for each read detected.
}
\usage{
read.size.freqs(locusInfo="primers",sampleDir="./samplesF",nmotif=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
\item{sampleDir}{
The directory where the filtered reads are stored. This directory should only contain the read files and nothing else. By default this is a directory
named samplesF written in the session working directory, as produced by \code{\link{read.quality.filter}} 
}
\item{nmotif}{
Reads are retained if they have at least these many repetitions of the repeat motif. It takes either a single value to be 
applied to all loci or a vector of values
}
}
\details{
Searches sequenced reads for those containing the microsatellite FWD and RVS primers and a user-defined number of microsatellite motif repeats. 
 Writes files for each sample and loci combination containing the read size and the number of msat motif repeats for each read detected.
 Running this function is the slowest step in the pipeline, but the files only need to be created once per set of new samples.
}
\value{
Files for each sample and locus combinations with the following two columns for each read detected:\cr
\cr
ReadSize: The length in bp of each read\cr
\cr
RepeatNumber: The number of microsatellite motifs for each read\cr
\cr
These files are written to a folder named ReadSizesRepeatNumbers that is created by the function in the working environment and will be used
down the pipeline by the functions \code{\link{reads.plots}}, \code{\link{GenotypeModel}}, and \code{\link{alleleFinder}}.
}
\references{
Alberto F (in prep)...
}
\author{
F. Alberto
}
\seealso{
\code{\link{reads.plots}}, \code{\link{alleleFinder}}
}
\keyword{ IO }
