\name{reads.plots}
\alias{reads.plots}
\title{
Creates frequency plots for the reads of specfic locus for each sample
}
\description{
Creates plots to visualize the read frequencies and help discerning the possible microsatellite alleles. 
}
\usage{
reads.plots(sampleDir="./samplesF",locusInfo="primers",x.range=c(90,350),
            x.label.cex=0.8,RepNumberPlot=FALSE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{sampleDir}{
A string with the path for the directory where the filtered reads are stored. This directory should only contain the read files and 
nothing else. This should be the samplesF directory written in the session working directory, 
as produced by \code{\link{read.quality.filter}}.
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
\item{x.range}{
The x-axis limits plot limits. If some sample has alleles outside this range the function will adjust the x-axis limit for that sample.
The defautl range is 90-350 bp.
}
\item{x.label.cex}{
The plot size for allele labels (sizes in bp) along the x-axis. The default size is 0.8.
}
\item{RepNumberPlot}{
A logical to decide if the second plot per sample should be plotted. It defaults to FALSE.
}
}
\details{
The function provides a visualization tool that separates samples by loci, i.e. one plot per loci, with all samples and two pages per sample. 
The first sample plot is the classical barplot with x-axis for read sizes and y-axis for their frequency. The second plot shows read size vs
number of motif repeats for reads. The plots are based on information from reads that contain the locus primers and a certain number of
repeats for the msat motif in each locus. The function reads the files created by function \code{\link{read.size.freqs}} written
to the folder ReadSizesRepeatNumbers.
}
\value{
Plots for each loci, each containing two pages per sample. The first plot is the classical barplot with x-axis for read sizes and
 y-axis for their frequency. The second plot shows read size vs number of motif repeats for reads containing both primers and a minimum number of
 motif repeats (as set by function \code{\link{read.size.freqs}}).
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
\keyword{ hplot }
